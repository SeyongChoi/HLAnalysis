use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
use rayon::prelude::*;
use std::collections::HashMap;

use crate::atoms::Atoms;
use crate::utils::spatial::mic::mic;
use crate::utils::spatial::pycell3::PyCell3;
use crate::utils::spatial::inverse_cell::inverse_cell;
use crate::io::reader::read_atoms_from_bin;

use crate::analysis::types::oh::OHRecord as OHResult;
use crate::analysis::types::molecule::MoleculeRecord as MoleculeResult;


#[inline]
fn effect_normal(surface_normal: &[f64], o_z: f64, center: Option<f64>) -> [f64; 3] {
    if let Some(zc) = center {
        if o_z <= zc {
            [surface_normal[0], surface_normal[1], surface_normal[2]]
        } else {
            [-surface_normal[0], -surface_normal[1], -surface_normal[2]]
        }
    } else {
        [surface_normal[0], surface_normal[1], surface_normal[2]]
    }
}

#[inline]
fn unit(v: [f64; 3]) -> [f64; 3] {
    let n = (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]).sqrt();
    if n > 1e-10 { [v[0]/n, v[1]/n, v[2]/n] } else { [0.0, 0.0, 0.0] }
}

#[inline]
fn rz_from_unit(e: [f64; 3], n: [f64; 3]) -> [f64; 3] {
    // p_z = dot(e, n), p_r = sqrt(1 - p_z^2), psi is unused -> 0.0
    let p_z = (e[0]*n[0] + e[1]*n[1] + e[2]*n[2]).clamp(-1.0, 1.0);
    let p_r = (1.0 - p_z*p_z).sqrt();
    [p_r, 0.0, p_z]
}

pub fn precompute_indices(
    atoms: &Atoms,
    index_silanol: &[i32],
    index_water: &[i32],
) -> (Vec<i32>, Vec<i32>, Vec<i32>, Vec<i32>) {
    let symbols = atoms.symbols.as_ref().unwrap();
    
    let oxygen_indices: Vec<_> = symbols.iter().enumerate()
        .filter_map(|(i, s)| if s == "O" { Some(i as i32) } else { None })
        .collect();

    let (silanol_oxygen, silanol_hydrogen): (Vec<_>, Vec<_>) = symbols.iter().enumerate()
        .filter(|(i, s)| s == &"O" && index_silanol.contains(&(*i as i32)) || s == &"H" && index_silanol.contains(&(*i as i32)))
        .partition(|(_, s)| s == &"O");

    let silanol_oxygen_indices = silanol_oxygen.into_iter().map(|(i, _)| i as i32).collect();
    let silanol_hydrogen_indices = silanol_hydrogen.into_iter().map(|(i, _)| i as i32).collect();

    let water_oxygen_indices: Vec<_> = symbols.iter().enumerate()
        .filter_map(|(i, s)| if s == "O" && index_water.contains(&(i as i32)) { Some(i as i32) } else { None })
        .collect();

    (oxygen_indices, silanol_oxygen_indices, silanol_hydrogen_indices, water_oxygen_indices)
}

pub fn process_silanol(
    o_idx: i32,
    o_pos: [f64; 3],
    silanol_map: &HashMap<i32, i32>,
    cell: &[[f64; 3]; 3],
    cell_inv: &[[f64; 3]; 3],
    atoms: &Atoms,
    surface_normal: &[f64], 
    center: Option<f64>,               // z-center to flip normal
    z_threshold_silanol: &Option<f64>, // Option<z_center of substrate>
) -> Option<MoleculeResult> {
    let o_vel = atoms.velocities.as_ref().unwrap().get(o_idx as usize)?;

    let h_idx = silanol_map.get(&o_idx)?;
    let h_pos_raw = atoms.positions.as_ref()?.get(*h_idx as usize)?;
    let h_vel = atoms.velocities.as_ref().unwrap().get(*h_idx as usize)?;

    // OH vector with MIC; make unit
    let oh = mic(
        [h_pos_raw[0] - o_pos[0], h_pos_raw[1] - o_pos[1], h_pos_raw[2] - o_pos[2]],
        cell, cell_inv,
    );
    let e_oh = unit(oh);
    let h_pos = [o_pos[0] + oh[0], o_pos[1] + oh[1], o_pos[2] + oh[2]];

    // effective normal & r–z projection
    let effective_normal = effect_normal(surface_normal, o_pos[2], center);
    let rz_oh1 = rz_from_unit(e_oh, effective_normal);

    // angle from p_z
    let oh_theta = rz_oh1[2].acos().to_degrees();

    // label
    let mut mol_type = "S".to_string();
    if let Some(zc) = z_threshold_silanol {
        mol_type = if o_pos[2] > *zc { "S-U".into() } else { "S-L".into() };
        let state = if (5.0..=45.0).contains(&oh_theta) { "D" }
                    else if (75.0..=125.0).contains(&oh_theta) { "A" }
                    else { "NA" };
        mol_type = format!("{}-{}", mol_type, state);
    }
    
    // build MoleculeResult: only fields needed here are set
    let rec = MoleculeResult::new(o_idx)
        .with_h_indices(Some(*h_idx), None)
        .with_positions(Some(o_pos), Some(h_pos), Some(h_pos))
        .with_velocities(Some(*o_vel), Some(*h_vel), Some(*h_vel))
        .with_orient(Some(e_oh), None, None, None)
        .with_labels(Some(mol_type));

    Some(rec)
}

pub fn process_water(
    o_idx: i32,
    o_pos: [f64; 3],
    cell: &[[f64; 3]; 3],
    cell_inv: &[[f64; 3]; 3],
    atoms: &Atoms, 
    surface_normal: &[f64], 
    center: Option<f64>,  // z_center coordinate   
    z_threshold_water: &Option<Vec<[f64; 2]>>, // Option<Vec<[z_min, z_max]>
) -> Option<MoleculeResult> {
    
    let o_vel = atoms.velocities.as_ref().unwrap().get(o_idx as usize)?;

    let h1_idx = o_idx + 1;
    let h2_idx = o_idx + 2;
    
    let h1_pos_raw = atoms.positions.as_ref()?.get(h1_idx as usize)?;
    let h2_pos_raw = atoms.positions.as_ref()?.get(h2_idx as usize)?;

    let h1_vel = atoms.velocities.as_ref().unwrap().get(h1_idx as usize)?;
    let h2_vel = atoms.velocities.as_ref().unwrap().get(h2_idx as usize)?;

    // OH vector with MIC; make unit
    let oh1 = mic(
        [h1_pos_raw[0] - o_pos[0], h1_pos_raw[1] - o_pos[1], h1_pos_raw[2] - o_pos[2]],
        cell, cell_inv,
    );
    let oh2 = mic(
        [h2_pos_raw[0] - o_pos[0], h2_pos_raw[1] - o_pos[1], h2_pos_raw[2] - o_pos[2]],
        cell, cell_inv,
    );

    let e_oh1 = unit(oh1);
    let e_oh2 = unit(oh2);

    let h1_pos = [o_pos[0] + oh1[0], o_pos[1] + oh1[1], o_pos[2] + oh1[2]];
    let h2_pos = [o_pos[0] + oh2[0], o_pos[1] + oh2[1], o_pos[2] + oh2[2]];

    // unit HH and HOH(bisector/dipole-like)
    let hh = unit([h1_pos[0] - h2_pos[0], h1_pos[1] - h2_pos[1], h1_pos[2] - h2_pos[2]]);
    let hoh = unit([h1_pos[0] + h2_pos[0] - 2.0 * o_pos[0],
                              h1_pos[1] + h2_pos[1] - 2.0 * o_pos[1],
                              h1_pos[2] + h2_pos[2] - 2.0 * o_pos[2]]);

    // effective normal & r–z projection
    let effective_normal = effect_normal(surface_normal, o_pos[2], center);
    let rz_hh = rz_from_unit(hh, effective_normal);
    let rz_hoh = rz_from_unit(hoh, effective_normal);

    // angles in degrees from p_z
    let mut hh_angle    = rz_hh[2].acos().to_degrees();
    if hh_angle < 90.0 { hh_angle = 180.0 - hh_angle; } // symmetry
    let dipole_angle = rz_hoh[2].acos().to_degrees();

    // layer/state labeling
    let mut mol_type = "W".to_string();
    if let Some(ths) = z_threshold_water {
        for (i, &[z_min, z_max]) in ths.iter().enumerate() {
            if o_pos[2] >= z_min && o_pos[2] < z_max {
                mol_type = match i {
                    0 => if o_pos[2] < z_min + 1.6 { "W-1a-L" } else { "W-1b-L" }.into(),
                    1 => "W-2-L".into(),
                    2 => "W-3-L".into(),
                    3 => "W-B".into(),
                    4 => "W-3-U".into(),
                    5 => "W-2-U".into(),
                    6 => if o_pos[2] > z_max - 1.6 { "W-1a-U" } else { "W-1b-U" }.into(),
                    _ => "W".into(),
                };
                let state = if (90.0..=110.0).contains(&hh_angle) && (155.0..=180.0).contains(&dipole_angle) { "S2" }
                            else if (110.0..=155.0).contains(&hh_angle) && (110.0..=155.0).contains(&dipole_angle) { "S1" }
                            else if (155.0..=180.0).contains(&hh_angle) && (70.0..=110.0).contains(&dipole_angle) { "SB" }
                            else if (90.0..=110.0).contains(&hh_angle) && (70.0..=110.0).contains(&dipole_angle) { "P" }
                            else if (110.0..=155.0).contains(&hh_angle) && (25.0..=70.0).contains(&dipole_angle) { "B1" }
                            else if (90.0..=110.0).contains(&hh_angle) && (0.0..=25.0).contains(&dipole_angle) { "B2" }
                            else { "NA" };
                mol_type = format!("{}-{}", mol_type, state);
                break;
            }
        }
    }
    
    // build MoleculeResult: only fields needed here are set
    let rec = MoleculeResult::new(o_idx)
        .with_h_indices(Some(h1_idx), Some(h2_idx))
        .with_positions(Some(o_pos), Some(h1_pos), Some(h2_pos))
        .with_velocities(Some(*o_vel), Some(*h1_vel), Some(*h2_vel))
        .with_orient(Some(e_oh1), Some(e_oh2), Some(hh), Some(hoh))
        .with_labels(Some(mol_type))
        .with_hbonds(None);

    Some(rec)
}


#[derive(Debug, Clone, Copy)]
pub enum ReadMode {
    Spectrum,
    CrossCorr,
}

#[derive(Debug)]
pub enum ReadFrameResult {
    Spectrum {
        n_ohs: i32,
        cell: PyCell3,
        mol_atom_map: HashMap<i32, Vec<i32>>,
        oh_infos: Vec<OHResult>,
        mol_infos: Vec<MoleculeResult>,
    },
    CrossCorr {
        results: Vec<Option<MoleculeResult>>,
        e_vecs: Vec<[f64; 3]>,
        mol_types: Vec<String>,
    },
}


// #[pyfunction]
pub fn read_frame(
    mode: String,
    frame_idx: usize,
    index_silanol: Vec<i32>,
    index_water: Vec<i32>,
    surface_normal: Vec<f64>,
    center: Option<f64>, // z_center coordinate
    z_threshold_water: Option<Vec<[f64; 2]>>,
    z_threshold_silanol: Option<f64>,
    dir: &str, // Optional directory path
) -> PyResult<ReadFrameResult> {    
    
    // -----------------------
    // 1) Load atoms & cell
    // -----------------------
    let directory = dir.to_string();
    let file_path = format!("{}/timestep_{}.bin", directory, frame_idx);
    let atoms = read_atoms_from_bin(&file_path)
        .or_else(|_| Err(PyErr::new::<PyValueError, _>("Failed to read atoms")))?;
    let cell = *atoms.cell.as_ref().unwrap();
    let cell_inv = inverse_cell(&cell);

    // -----------------------
    // 2) Precompute indices
    // -----------------------
    let (oxygen_indices, silanol_oxygen_indices, silanol_hydrogen_indices, water_oxygen_indices) = 
        precompute_indices(&atoms, &index_silanol, &index_water);
    let silanol_group_map: HashMap<i32, i32> = silanol_oxygen_indices
        .iter()
        .zip(silanol_hydrogen_indices.iter())
        .map(|(&o, &h)| (o, h))
        .collect();

    // -----------------------
    // 3) Per-oxygen processing (parallel)
    // -----------------------
    let results: Vec<MoleculeResult> = oxygen_indices.par_iter()
        .filter_map(|&o_idx| {
            let o_pos = atoms.positions.as_ref().unwrap()[o_idx as usize];
            let is_water = water_oxygen_indices.contains(&o_idx);
            let is_silanol = silanol_oxygen_indices.contains(&o_idx);

            if is_water {
                return process_water(
                    o_idx, 
                    o_pos, 
                    &cell, 
                    &cell_inv, 
                    &atoms,
                    &surface_normal,
                    center.clone(),
                    &z_threshold_water

                );
            }

            if is_silanol{
                return process_silanol(
                    o_idx, 
                    o_pos,
                    &silanol_group_map, 
                    &cell, 
                    &cell_inv, 
                    &atoms,
                    &surface_normal,
                    center.clone(),
                    &z_threshold_silanol

                );
            }
            None
        })
        .collect();
    // -----------------------
    // 4) Mode-specific aggregation
    // -----------------------
    match mode.as_str() {
        "spectrum" =>{
            let mut n_ohs: i32 = 0;
            let mut oh_infos: Vec<OHResult> = Vec::new();
            let mut mol_atom_map: HashMap<i32, Vec<i32>> = HashMap::new();

            for r in &results {
                let mt_ok = r.mol_type.as_ref().map(|s| s.contains('W')).unwrap_or(false);
                if mt_ok {
                    let o_idx = r.o_idx;
                    let h1 = r.h1_idx.expect("water h1_idx must be Some");
                    let h2 = r.h2_idx.expect("water h2_idx must be Some");

                    let o_pos = r.o_pos.expect("water o_pos must be Some");
                    let h1_pos = r.h1_pos.expect("water h1_pos must be Some");
                    let h2_pos = r.h2_pos.expect("water h2_pos must be Some");

                    let o_vel = r.o_vel.expect("water o_vel must be Some");
                    let h1_vel = r.h1_vel.expect("water h1_vel must be Some");
                    let h2_vel = r.h2_vel.expect("water h2_vel must be Some");

                    let mt = r.mol_type.as_ref().unwrap().clone();

                    // OH1
                    mol_atom_map.insert(n_ohs, vec![o_idx, h1]);
                    oh_infos.push(OHResult {
                        o_idx,
                        h_idx: h1,
                        o_pos: Some(o_pos),
                        h_pos: Some(h1_pos),
                        o_vel: Some(o_vel),
                        h_vel: Some(h1_vel),
                        mol_type: Some(mt.clone()),
                    });
                    n_ohs += 1;

                    // OH2
                    mol_atom_map.insert(n_ohs, vec![o_idx, h2]);
                    oh_infos.push(OHResult {
                        o_idx,
                        h_idx: h2,
                        o_pos: Some(o_pos),
                        h_pos: Some(h2_pos),
                        o_vel: Some(o_vel),
                        h_vel: Some(h2_vel),
                        mol_type: Some(mt.clone())
                    });
                    n_ohs += 1;
                } else {
                    let mt_ok_s = r.mol_type.as_ref().map(|s| s.contains('S')).unwrap_or(false);
                    if mt_ok_s {
                        let o_idx = r.o_idx;
                        let h1 = r.h1_idx.expect("silanol h1_idx must be Some");

                        let o_pos = r.o_pos.expect("silanol o_pos must be Some");
                        let h1_pos = r.h1_pos.expect("silanol h1_pos must be Some");

                        let o_vel = r.o_vel.expect("silanol o_vel must be Some");
                        let h1_vel = r.h1_vel.expect("silanol h1_vel must be Some");

                        let mt = r.mol_type.as_ref().unwrap().clone();

                        mol_atom_map.insert(n_ohs, vec![o_idx, h1]);
                        oh_infos.push(OHResult {
                            o_idx,
                            h_idx: h1,
                            o_pos: Some(o_pos),
                            h_pos: Some(h1_pos),
                            o_vel: Some(o_vel),
                            h_vel: Some(h1_vel),
                            mol_type: Some(mt),
                        });
                        n_ohs += 1;
                    }
                }
            }

            Ok(ReadFrameResult::Spectrum {
                n_ohs,
                cell: PyCell3(cell),
                mol_atom_map,
                oh_infos,
                mol_infos: results,
            })
        }
        "crosscorr" =>{
            Err(PyErr::new::<pyo3::exceptions::PyNotImplementedError, _>(
                    "CrossCorr mode is not implemented yet",
                ))        
            }
        other => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            format!("unknown mode: {other}")
        )),
    }
 }