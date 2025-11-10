use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
use rayon::prelude::*;
use std::sync::Arc;
use crate::io::reader::read_dipole_from_bin;

use crate::utils::spatial::mic::mic;
use crate::utils::spatial::pycell3::PyCell3;
use crate::utils::spatial::inverse_cell::inverse_cell;

use crate::analysis::types::dipole_polar::DipolePolarRecord as DipolePolarResult;

use crate::analysis::spectroscopy::tools::norm;

/// Compute VVACF by streaming dipole–polar frames from bincode files on disk.
///
/// Reads `{bin_file}/dipole_polar_{t}.bin` for `t = t0 .. t0+tcfl-1` at each time origin `t0`,
/// applies optional OH type filtering, and accumulates auto / intra / inter terms.
/// Supports IR vs SFG by toggling the second factor (`rel_vel_z` vs `v_proj`).
#[pyfunction]
pub fn compute_vvacf_from_bin(
    bin_file: &str,
    n_ohs: i32,
    start: Option<usize>,
    end: Option<usize>,
    frame_interval: Option<usize>,
    tcfl: usize,
    cell: PyCell3,
    auto: bool,
    intra: bool,
    inter: bool,
    rcut_inter: f64,
    oh_type: Option<Vec<String>>, // e.g., Some(vec!["W".into(), "1".into()])
    spectrum: Option<String>,     // "IR" or "SFG"
) -> PyResult<(f64, Vec<f64>)> {
    // --- Setup & initialization ---
    let cell = &cell.0;
    let cell_inv = inverse_cell(cell);

    let start_frame = start.unwrap_or(0);
    let end_frame = end.unwrap_or(usize::MAX);
    let frame_interval = frame_interval.unwrap_or(1);

    // Generate list of time origins (t0 values)
    let t0_values: Vec<usize> = (start_frame..(end_frame.saturating_sub(tcfl + 1)))
        .step_by(frame_interval)
        .collect();

    let oh_type_filter = oh_type.map(|v| Arc::new(v));

    // Spectrum type: IR or SFG
    let spec_type = spectrum.unwrap_or_else(|| "IR".to_string());
    let is_sfg = spec_type.eq_ignore_ascii_case("SFG");
    match spec_type.as_str() {
        "SFG" | "IR" => {}
        _ => {
            return Err(PyErr::new::<PyValueError, _>(
                "spectrum must be 'IR' or 'SFG'",
            ))
        }
    }

    // --- Parallelized loop over t0 using Rayon global thread pool ---
    let (vvacf, total_mean, total_count) = t0_values
        .into_par_iter()
        .map(|t0| {
            // Local accumulators per time origin
            let mut local_vvacf = vec![0.0; tcfl];
            let mut local_mean = 0.0f64;
            let mut local_count = 0usize;

            // Read a contiguous block of tcfl frames
            let start_time = std::time::Instant::now();
            let mut dipole_frames: Vec<Vec<DipolePolarResult>> = Vec::with_capacity(tcfl);
            for t in 0..tcfl {
                let file_path = format!("{}/dipole_polar_{}.bin", bin_file, t0 + t);
                let dipole_polars = match read_dipole_from_bin(&file_path) {
                    Ok(data) => data,
                    Err(_) => return (vec![0.0; tcfl], 0.0, 0),
                };
                dipole_frames.push(dipole_polars);
            }
            let elapsed = start_time.elapsed();
            let start_time_2 = std::time::Instant::now();

            // --- Loop over OH bonds ---
            for i in 0..n_ohs as usize {
                // Extract type string as &str ("" if None)
                let oh0: &str = dipole_frames[0][i].oh_type.as_deref().unwrap_or("");

                // Apply optional OH-type filter
                let should_process = match &oh_type_filter {
                    Some(filter) => filter.iter().all(|tok| oh0.contains(tok)),
                    None => true,
                };
                if !should_process {
                    continue;
                }

                local_mean += 1.0;

                // --- Time correlation accumulation ---
                for t_ in 0..tcfl {
                    let ref_term = &dipole_frames[0][i];
                    let next_term = &dipole_frames[t_][i];
                    let sign = if oh0.contains('L') { 1.0 } else { -1.0 };
                    // Safe numeric unwrapping (missing -> 0.0)
                    let rvz_ref    = ref_term.rel_vel_z.unwrap_or(0.0);
                    let rvz_next   = next_term.rel_vel_z.unwrap_or(0.0);
                    let vproj_next = next_term.v_proj.unwrap_or(0.0);

                    // Auto-correlation term
                    if auto {
                        let contrib = if is_sfg {
                            sign * rvz_ref * vproj_next
                        } else {
                            sign * rvz_ref * rvz_next
                        };
                        local_vvacf[t_] += sign * contrib;
                    }

                    // Intra-molecular term (same O index, different H)
                    if intra {
                        for j in 0..n_ohs as usize {
                            if i != j && ref_term.o_idx == dipole_frames[0][j].o_idx {
                                let contrib = if is_sfg {
                                    sign * rvz_ref * dipole_frames[t_][j].v_proj.unwrap_or(0.0)
                                } else {
                                    sign * rvz_ref * dipole_frames[t_][j].rel_vel_z.unwrap_or(0.0)
                                };
                                local_vvacf[t_] += sign * contrib;
                            }
                        }
                    }

                    // Inter-molecular term (different O index within cutoff)
                    if inter {
                        // --- i OH-center (safe unwrap) ---
                        let (o_i, h_i) = match (ref_term.o_pos, ref_term.h_pos) {
                            (Some(o), Some(h)) => (o, h),
                            _ => {
                                // missing positions → skip this i
                                continue;
                            }
                        };
                        let oh_i_pos = [
                            0.5 * (o_i[0] + h_i[0]),
                            0.5 * (o_i[1] + h_i[1]),
                            0.5 * (o_i[2] + h_i[2]),
                        ];

                        for j in 0..n_ohs as usize {
                            if i == j || ref_term.o_idx == dipole_frames[0][j].o_idx {
                                continue;
                            }
                            // --- j OH-center (safe unwrap) ---
                            let (o_j_opt, h_j_opt) = (dipole_frames[0][j].o_pos, dipole_frames[0][j].h_pos);
                            let (o_j, h_j) = match (o_j_opt, h_j_opt) {
                                (Some(o), Some(h)) => (o, h),
                                _ => {
                                    // missing positions → skip this j
                                    continue;
                                }
                            };

                            let oh_j_pos = [
                                0.5 * (o_j[0] + h_j[0]),
                                0.5 * (o_j[1] + h_j[1]),
                                0.5 * (o_j[2] + h_j[2]),
                            ];

                            let rij_vec = mic(
                                [
                                    oh_j_pos[0] - oh_i_pos[0],
                                    oh_j_pos[1] - oh_i_pos[1],
                                    oh_j_pos[2] - oh_i_pos[2],
                                ],
                                cell,
                                &cell_inv,
                            );

                            if norm(rij_vec) <= rcut_inter {
                                let contrib = if is_sfg {
                                    sign * rvz_ref * dipole_frames[t_][j].v_proj.unwrap_or(0.0)
                                } else {
                                    sign * rvz_ref * dipole_frames[t_][j].rel_vel_z.unwrap_or(0.0)
                                };
                                local_vvacf[t_] += sign * contrib;
                            }
                        }
                    }
                }
            }
            let elapsed_2 = start_time_2.elapsed();
            println!("Time taken to whole process vvacf for t0 = {}: {:?}", t0, elapsed+elapsed_2);

            local_count += 1;
            (local_vvacf, local_mean, local_count)
        })
        .reduce(
            || (vec![0.0; tcfl], 0.0, 0usize),
            |(mut acc_vvacf, acc_mean, acc_count),
             (local_vvacf, local_mean, local_count)| {
                for t in 0..tcfl {
                    acc_vvacf[t] += local_vvacf[t];
                }
                (acc_vvacf, acc_mean + local_mean, acc_count + local_count)
            },
        );

    // --- Final normalization ---
    let mean = if total_count > 0 {
        total_mean / total_count as f64
    } else {
        0.0
    };

    Ok((mean, vvacf))
}