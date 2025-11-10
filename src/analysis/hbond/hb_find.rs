use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
use rayon::prelude::*;

use std::collections::HashMap;

use crate::utils::spatial::mic::mic;
use crate::utils::spatial::inverse_cell::inverse_cell;
use crate::io::reader::read_atoms_from_bin;

use crate::analysis::types::molecule::{
    MoleculeRecord as MoleculeResult,
    HydrogenBondPartner,
};

use crate::analysis::pre_process::extract_molecules_from_atoms;

fn calc_ang(v1: [f64; 3], v2: [f64; 3]) -> f64 {
    let dot = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
    dot.clamp(-1.0, 1.0).acos()
}

#[pyfunction]
pub fn find_hb_per_frame(
    frame_idx: usize,
    index_silanol: Vec<i32>,
    index_water: Vec<i32>,
    surface_normal: Vec<f64>,
    center: Option<f64>,
    z_threshold_water: Option<Vec<[f64; 2]>>,
    z_threshold_silanol: Option<f64>,
    hb_conditions: Vec<f64>,
    dir: &str
) -> PyResult<Vec<MoleculeResult>> {

    let file_path = format!("{}/timestep_{}.bin", dir, frame_idx);
    let atoms = read_atoms_from_bin(&file_path)?;

    let mut molecules: Vec<MoleculeResult> = Python::with_gil(|py| {
        py.allow_threads(|| {
            extract_molecules_from_atoms(
                &atoms,
                &index_silanol,
                &index_water,
                &surface_normal,
                center.clone(),
                &z_threshold_water,
                &z_threshold_silanol,
            )
        })
    })?;

    let nmols = molecules.len();
    if nmols < 2 { return Ok(molecules); }
    
    if hb_conditions.len() != 3 {
        return Err(PyValueError::new_err("hb_conditions must be a list of length 3: [rOO_max, rHO_max, ang_max]"));
    }
    let r_oo_max = hb_conditions[0];
    let r_ho_max = hb_conditions[1];
    let ang_max_deg = hb_conditions[2];
    
    let cell = atoms.cell.as_ref().unwrap();
    let cell_inv = inverse_cell(cell);
    // --- Neighbor List Algorithm for ANY cell geometry (including Triclinic) ---

    // 1. Determine grid dimensions based on the shortest perpendicular box heights.
    let h_inv_sq: Vec<f64> = (0..3).map(|i| cell_inv[i].iter().map(|&x| x*x).sum()).collect();
    let grid_dims = [
        (1.0 / (h_inv_sq[0].sqrt() * r_oo_max)).floor() as i32,
        (1.0 / (h_inv_sq[1].sqrt() * r_oo_max)).floor() as i32,
        (1.0 / (h_inv_sq[2].sqrt() * r_oo_max)).floor() as i32,
    ];
    let grid_dims = [grid_dims[0].max(1), grid_dims[1].max(1), grid_dims[2].max(1)];

    // 2. Create and populate the cell list using fractional coordinates
    let mut cell_map: HashMap<(i32, i32, i32), Vec<usize>> = HashMap::new();
    let frac_coords: Vec<[f64; 3]> = molecules.iter().map(|mol| {
        let p = mol.o_pos.expect("o_pos must be Some for HB search");
        let mut s = [
            p[0] * cell_inv[0][0] + p[1] * cell_inv[1][0] + p[2] * cell_inv[2][0],
            p[0] * cell_inv[0][1] + p[1] * cell_inv[1][1] + p[2] * cell_inv[2][1],
            p[0] * cell_inv[0][2] + p[1] * cell_inv[1][2] + p[2] * cell_inv[2][2],
        ];
        s[0] -= s[0].floor(); s[1] -= s[1].floor(); s[2] -= s[2].floor();
        s
    }).collect();

    let cell_indices: Vec<(i32, i32, i32)> = frac_coords.iter().map(|s| (
        (s[0] * grid_dims[0] as f64).floor() as i32,
        (s[1] * grid_dims[1] as f64).floor() as i32,
        (s[2] * grid_dims[2] as f64).floor() as i32,
    )).collect();

    for (i, &cell_idx) in cell_indices.iter().enumerate() {
        cell_map.entry(cell_idx).or_default().push(i);
    }
    
    // 3. Search for H-bonds in parallel using the fractional cell list
    let hbond_updates: Vec<(usize, HydrogenBondPartner)> = (0..nmols)
        .into_par_iter()
        .flat_map(|i| {
            let mut found_partners: Vec<(usize, HydrogenBondPartner)> = Vec::new();
            let donor: &MoleculeResult = &molecules[i];
            
            // Set variable for donor information
            let donor_o_pos = donor
                .o_pos
                .expect("donor.o_pos must exist for HB search");

            let donor_h1_pos = donor
                .h1_pos
                .expect("donor.o_pos must exist for HB search");
            
            let donor_e_oh1 = donor
                .e_oh1
                .expect("donor.e_oh1 must exist for HB angle");

            let (has_h2, donor_h2_pos, donor_e_oh2) = if donor.h2_idx.is_some() {
                (true,
                Some(donor
                    .h2_pos
                    .expect("donor.h2_pos must exist when h2_idx is Some")),
                Some(donor
                    .e_oh2
                    .expect("donor.e_oh2 must exist when h2_idx is Some")),
                )
            } else {
                (false, None, None)
            };

            let donor_cell_idx = cell_indices[i];

            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        let neighbor_cell_idx = (
                            (donor_cell_idx.0 + dx + grid_dims[0]) % grid_dims[0],
                            (donor_cell_idx.1 + dy + grid_dims[1]) % grid_dims[1],
                            (donor_cell_idx.2 + dz + grid_dims[2]) % grid_dims[2],
                        );

                        if let Some(acceptor_indices) = cell_map.get(&neighbor_cell_idx) {
                            for &j in acceptor_indices {
                                if i == j { continue; }

                                let acceptor: &MoleculeResult = &molecules[j];

                                let acceptor_o_pos = acceptor
                                    .o_pos
                                    .expect("acceptor.o_pos must exist for HB search");

                                let oo_vec = mic([acceptor_o_pos[0] - donor_o_pos[0],
                                                                    acceptor_o_pos[1] - donor_o_pos[1],
                                                                    acceptor_o_pos[2] - donor_o_pos[2]],
                                                                    cell, &cell_inv);
                                if oo_vec.iter().map(|&x| x*x).sum::<f64>() < r_oo_max.powi(2) {
                                    let d_oo = oo_vec.iter().map(|&x| x*x).sum::<f64>().sqrt();
                                    
                                    
                                    let e_oo = [oo_vec[0]/d_oo, oo_vec[1]/d_oo, oo_vec[2]/d_oo]; // Vector from donor O to acceptor O 

                                    // Check donor's H1
                                    let ho_vec = mic([acceptor_o_pos[0] - donor_h1_pos[0],
                                                                        acceptor_o_pos[1] - donor_h1_pos[1], 
                                                                        acceptor_o_pos[2] - donor_h1_pos[2]],
                                                                        cell,
                                                                        &cell_inv);
                                    let angle1 = calc_ang(donor_e_oh1, e_oo).to_degrees();

                                    // if ho_vec.iter().map(|&x| x*x).sum::<f64>().sqrt() < r_ho_max && calc_ang(donor.e_oh1, e_oo).to_degrees() < ang_max_deg {
                                    if angle1 < ang_max_deg {
                                        found_partners.push((i, HydrogenBondPartner { partner_o_idx: acceptor.o_idx, partner_moltype: acceptor.mol_type.clone(), h_bond_type: Some("donor".to_string()) }));
                                        found_partners.push((j, HydrogenBondPartner { partner_o_idx: donor.o_idx, partner_moltype: donor.mol_type.clone(), h_bond_type: Some("acceptor".to_string()) }));
                                    }

                                    // Check donor's H2
                                    if has_h2 {
                                        let donor_h2_pos = donor_h2_pos.expect("checked above");
                                        let donor_e_oh2  = donor_e_oh2.expect("checked above");
                                        let ho_vec2 = mic([acceptor_o_pos[0] - donor_h2_pos[0], 
                                                                            acceptor_o_pos[1] - donor_h2_pos[1], 
                                                                            acceptor_o_pos[2] - donor_h2_pos[2]], 
                                                                            cell, &cell_inv);
                                        let angle2 = calc_ang(donor_e_oh2, e_oo).to_degrees();

                                        // if ho_vec2.iter().map(|&x| x*x).sum::<f64>().sqrt() < r_ho_max && calc_ang(donor.e_oh2, e_oo).to_degrees() < ang_max_deg {
                                        if angle2 < ang_max_deg {
                                            found_partners.push((i, HydrogenBondPartner { partner_o_idx: acceptor.o_idx, partner_moltype: acceptor.mol_type.clone(), h_bond_type: Some("donor".to_string()) }));
                                            found_partners.push((j, HydrogenBondPartner { partner_o_idx: donor.o_idx, partner_moltype: donor.mol_type.clone(), h_bond_type: Some("acceptor".to_string()) }));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            found_partners
        })
        .collect();

    // 4. Apply all updates sequentially
    for (mol_idx, partner_info) in hbond_updates {
        molecules[mol_idx]
            .h_bond_partners
            .get_or_insert_with(Vec::new)  // Option<Vec<...>>
            .push(partner_info);
    }
    
    Ok(molecules)
}