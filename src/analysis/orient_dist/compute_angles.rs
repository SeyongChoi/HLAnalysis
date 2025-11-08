use rayon::prelude::*;
use std::collections::HashMap;

use crate::atoms::Atoms;
use crate::utils::spatial::mic::mic;

/// Compute O–H angles (in degrees) w.r.t. a surface normal for a given frame.
/// - Applies MIC to bond vectors.
/// - If `center` is provided, flips the normal above `z_center`.
pub fn compute_oh_angles(
    atoms: &Atoms,
    oh_bond_map: &HashMap<i32, (i32, i32)>,
    surface_normal: &[f64],
    cell: &[[f64; 3];3],
    inv_cell: &[[f64; 3];3],
    center: Option<f64>,  // z-center coordinate
) -> Vec<(i32, i32, f64)> {
    let mut local_angles = Vec::new();
    let surface_normal_inv = [-surface_normal[0], -surface_normal[1], -surface_normal[2]];

    if let Some(positions) = &atoms.positions {
        // Parallel over O–H pairs
        local_angles = oh_bond_map
            .par_iter()
            .map(|(&_idx, &(o_idx, h_idx))| {
                let o_pos = positions[o_idx as usize];
                let h_pos = positions[h_idx as usize];

                // Raw O->H vector
                let mut oh_vec = [h_pos[0] - o_pos[0], h_pos[1] - o_pos[1], h_pos[2] - o_pos[2]];

                // Apply MIC to bring H near O
                let mic_corr = mic(oh_vec, cell, inv_cell);
                let h_pos_corr = [mic_corr[0] + o_pos[0], mic_corr[1] + o_pos[1], mic_corr[2] + o_pos[2]];

                // Recompute O->H after MIC
                oh_vec = [h_pos_corr[0] - o_pos[0], h_pos_corr[1] - o_pos[1], h_pos_corr[2] - o_pos[2]];

                // Magnitude
                let dist = (oh_vec[0].powi(2) + oh_vec[1].powi(2) + oh_vec[2].powi(2)).sqrt();

                // Dot with proper normal (flip if above z_center)
                let dot = if let Some(zc) = center {
                    if o_pos[2] <= zc {
                        oh_vec[0]*surface_normal[0] + oh_vec[1]*surface_normal[1] + oh_vec[2]*surface_normal[2]
                    } else {
                        oh_vec[0]*surface_normal_inv[0] + oh_vec[1]*surface_normal_inv[1] + oh_vec[2]*surface_normal_inv[2]
                    }
                } else {
                    oh_vec[0]*surface_normal[0] + oh_vec[1]*surface_normal[1] + oh_vec[2]*surface_normal[2]
                };

                // Safe acos
                let cosv = (dot / dist).clamp(-1.0, 1.0);
                let angle_deg = cosv.acos().to_degrees();

                (o_idx, h_idx, angle_deg)
            })
            .collect();
    }

    local_angles
}

/// Compute dipole–normal angles (in degrees) for water molecules in one frame.
/// - Dipole ≈ (r_H1 + r_H2 - 2 r_O).
/// - Applies MIC to both O–H vectors before building the dipole.
/// - Flips the reference normal above `z_center` if provided.
pub fn compute_dipole_angles(
    atoms: &Atoms,
    water_map: &HashMap<i32, (i32, i32, i32)>,
    surface_normal: &[f64],
    cell: &[[f64; 3];3],
    inv_cell: &[[f64; 3];3],
    center: Option<f64>,  // z_center coordinate
) -> Vec<(i32, i32, i32, f64)> {
    let mut local_angles = Vec::new();
    let surface_normal_inv = [-surface_normal[0], -surface_normal[1], -surface_normal[2]];

    if let Some(positions) = &atoms.positions {
        // Parallel over water triplets
        local_angles = water_map
            .par_iter()
            .map(|(&_idx, &(o_idx, h1_idx, h2_idx))| {
                let o_pos  = positions[o_idx  as usize];
                let h1_pos = positions[h1_idx as usize];
                let h2_pos = positions[h2_idx as usize];

                // Raw O->H vectors
                let oh1 = [h1_pos[0] - o_pos[0], h1_pos[1] - o_pos[1], h1_pos[2] - o_pos[2]];
                let oh2 = [h2_pos[0] - o_pos[0], h2_pos[1] - o_pos[1], h2_pos[2] - o_pos[2]];

                // MIC-correct H positions around O
                let h1_corr = {
                    let c = mic(oh1, cell, inv_cell);
                    [c[0] + o_pos[0], c[1] + o_pos[1], c[2] + o_pos[2]]
                };
                let h2_corr = {
                    let c = mic(oh2, cell, inv_cell);
                    [c[0] + o_pos[0], c[1] + o_pos[1], c[2] + o_pos[2]]
                };

                // Approximate dipole vector
                let dip = [
                    h1_corr[0] + h2_corr[0] - 2.0 * o_pos[0],
                    h1_corr[1] + h2_corr[1] - 2.0 * o_pos[1],
                    h1_corr[2] + h2_corr[2] - 2.0 * o_pos[2],
                ];

                // Magnitude
                let dist = (dip[0].powi(2) + dip[1].powi(2) + dip[2].powi(2)).sqrt();

                // Dot with proper normal
                let dot = if let Some(zc) = center {
                    if o_pos[2] <= zc {
                        dip[0]*surface_normal[0] + dip[1]*surface_normal[1] + dip[2]*surface_normal[2]
                    } else {
                        dip[0]*surface_normal_inv[0] + dip[1]*surface_normal_inv[1] + dip[2]*surface_normal_inv[2]
                    }
                } else {
                    dip[0]*surface_normal[0] + dip[1]*surface_normal[1] + dip[2]*surface_normal[2]
                };

                // Safe acos
                let cosv = (dot / dist).clamp(-1.0, 1.0);
                let angle_deg = cosv.acos().to_degrees();

                (o_idx, h1_idx, h2_idx, angle_deg)
            })
            .collect();
    }

    local_angles
}

/// Compute (θ_dipole, θ_HH) **or** (cos θ_dipole, cos θ_HH) per water molecule.
/// - θ_HH is adjusted to be in [90°, 180°] by reflecting values < 90°.
/// - If `mode == "cos"`, returns cosines instead of angles (with HH cosine made non-positive).
/// - MIC is applied to O–H vectors prior to building dipole/HH vectors.
/// - Normal is flipped above `z_center` if `center` is provided.
pub fn compute_hh_hoh_angles(
    atoms: &Atoms,
    water_map: &HashMap<i32, (i32, i32, i32)>,
    surface_normal: &[f64],
    cell: &[[f64; 3];3],
    inv_cell: &[[f64; 3];3],
    center: Option<f64>,  // z_center coordinate
    mode: &str,           // "angle" or "cos"
) -> Vec<(i32, i32, i32, f64, f64)> {
    let mut local_angles = Vec::new();
    let surface_normal_inv = [-surface_normal[0], -surface_normal[1], -surface_normal[2]];

    if let Some(positions) = &atoms.positions {
        // Parallel over water triplets
        local_angles = water_map
            .par_iter()
            .map(|(&_idx, &(o_idx, h1_idx, h2_idx))| {
                let o_pos  = positions[o_idx  as usize];
                let h1_pos = positions[h1_idx as usize];
                let h2_pos = positions[h2_idx as usize];

                // Raw O->H vectors
                let oh1 = [h1_pos[0] - o_pos[0], h1_pos[1] - o_pos[1], h1_pos[2] - o_pos[2]];
                let oh2 = [h2_pos[0] - o_pos[0], h2_pos[1] - o_pos[1], h2_pos[2] - o_pos[2]];

                // MIC-correct H positions
                let h1_corr = {
                    let c = mic(oh1, cell, inv_cell);
                    [c[0] + o_pos[0], c[1] + o_pos[1], c[2] + o_pos[2]]
                };
                let h2_corr = {
                    let c = mic(oh2, cell, inv_cell);
                    [c[0] + o_pos[0], c[1] + o_pos[1], c[2] + o_pos[2]]
                };

                // Dipole and H–H vectors
                let dip = [
                    h1_corr[0] + h2_corr[0] - 2.0 * o_pos[0],
                    h1_corr[1] + h2_corr[1] - 2.0 * o_pos[1],
                    h1_corr[2] + h2_corr[2] - 2.0 * o_pos[2],
                ];
                let hh = [
                    h1_corr[0] - h2_corr[0],
                    h1_corr[1] - h2_corr[1],
                    h1_corr[2] - h2_corr[2],
                ];

                // Magnitudes
                let dip_norm = (dip[0].powi(2) + dip[1].powi(2) + dip[2].powi(2)).sqrt();
                let hh_norm  = (hh[0].powi(2)  + hh[1].powi(2)  + hh[2].powi(2)).sqrt();

                // Dot with flipped-or-not normal
                let dip_dot = if let Some(zc) = center {
                    if o_pos[2] <= zc {
                        dip[0]*surface_normal[0] + dip[1]*surface_normal[1] + dip[2]*surface_normal[2]
                    } else {
                        dip[0]*surface_normal_inv[0] + dip[1]*surface_normal_inv[1] + dip[2]*surface_normal_inv[2]
                    }
                } else {
                    dip[0]*surface_normal[0] + dip[1]*surface_normal[1] + dip[2]*surface_normal[2]
                };

                let hh_dot = if let Some(zc) = center {
                    if o_pos[2] <= zc {
                        hh[0]*surface_normal[0] + hh[1]*surface_normal[1] + hh[2]*surface_normal[2]
                    } else {
                        hh[0]*surface_normal_inv[0] + hh[1]*surface_normal_inv[1] + hh[2]*surface_normal_inv[2]
                    }
                } else {
                    hh[0]*surface_normal[0] + hh[1]*surface_normal[1] + hh[2]*surface_normal[2]
                };

                // Safe cosines
                let dip_cos = (dip_dot / dip_norm).clamp(-1.0, 1.0);
                let mut hh_cos = (hh_dot / hh_norm).clamp(-1.0, 1.0);

                // Convert to degrees if needed; ensure HH angle in [90°, 180°]
                let dip_ang = dip_cos.acos().to_degrees();
                let mut hh_ang = hh_cos.acos().to_degrees();
                if hh_ang < 90.0 {
                    hh_ang = 180.0 - hh_ang;
                }

                // For "cos" mode, make HH cosine non-positive (align with ≥ 90°)
                let hh_cos_corrected = if hh_cos > 0.0 { -hh_cos } else { hh_cos };

                if mode == "cos" {
                    (o_idx, h1_idx, h2_idx, dip_cos, hh_cos_corrected)
                } else {
                    (o_idx, h1_idx, h2_idx, dip_ang, hh_ang)
                }
            })
            .collect();
    }

    local_angles
}
