use pyo3::prelude::*;
use rayon::prelude::*;
use std::collections::HashMap;

use crate::atoms::Atoms;
use crate::utils::spatial::mic::mic;
use crate::io::reader::read_atoms_from_bin;

use crate::analysis::types::OHRecord as OHResult;
use crate::analysis::types::MoleculeRecord as MoleculeResult;
use crate::analysis::types::DipolePolarRecord as DipolePolarResult;

#[inline]
fn effect_normal(surface_normal: &[f64; 3], o_z: f64, center: Option<f64>) -> [f64; 3] {
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
    surface_normal: &[f64;3], 
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
    let effective_normal = effective_normal(surface_normal, o_pos[2], center);
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
        .with_h_indices(*h_idx, None)
        .with_positions(o_pos, h1_pos, h1_pos)
        .with_velocities(*o_vel, *h_vel, *h_vel)
        .with_labels(mol_type, None);

    Some(rec)
}


/* To do
   - process_water
   - read_frame
 */
// pub fn process_water

// #[pyfunction]
// pub fn read_frame