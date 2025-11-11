use pyo3::prelude::*;
use pyo3::exceptions::{PyIOError, PyValueError};
use rayon::prelude::*;

use std::fs::{self, File};
use std::sync::{Arc, Mutex};
use std::io::BufWriter;
use std::time::Instant;
use crate::analysis::pre_process::{read_frame, ReadFrameResult};
use crate::analysis::types::molecule::MoleculeRecord as MoleculeResult;
use crate::analysis::types::time_corr::TimeCorrRecord as TimeCorrResult;



fn calc_cn(
    e_o: [f64; 3],
    e_t: [f64; 3],
) -> (f64, f64, f64){

    let edot = e_o.iter().zip(e_t.iter()).map(|(a, b)| a * b).sum::<f64>();
    let c1 = edot;
    let c2 = 0.5 * (3.0 * edot.powi(2) - 1.0);
    let c3 = 0.5 * (5.0 * edot.powi(3) - 3.0 * edot);
    (c1, c2, c3)
}

#[pyfunction]
pub fn init_e_oh(
    frame_idx: usize,
    index_silanol: Vec<i32>,
    index_water: Vec<i32>,
    z_threshold_water: Option<Vec<[f64; 2]>>,
    z_threshold_silanol: Option<f64>,
    dir: &str,
) -> PyResult<(Vec<MoleculeResult>, Vec<Vec<[f64; 3]>>, Vec<String>)> {
    // Call read_frame in "timecorr" mode
    let res = read_frame(
        "timecorr".to_string(),           // mode: String
        frame_idx,
        index_silanol,
        index_water,
        vec![0.0, 0.0, 1.0],              // surface_normal: Vec<f64>
        None,                              // center
        z_threshold_water,
        z_threshold_silanol,
        dir,
    )?;

    // Extract fields from the TimeCorr variant
    let (results, e_ohs, mol_types) = match res {
        ReadFrameResult::TimeCorr { results, e_vecs, mol_types } => (results, e_vecs, mol_types),
        _ => {
            return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                "init_e_oh: read_frame did not return TimeCorr variant",
            ));
        }
    };

    // Duplicate e_ohs → [ e_ohs, e_ohs ]
    let eOHs = vec![e_ohs.clone(), e_ohs];

    Ok((results, eOHs, mol_types))
}

#[pyfunction]
pub fn time_corr(
    start_frame: i32,
    end_frame: i32,
    frame_interval: i32,
    corr_times: i32,
    corr_interval: i32,
    index_silanol: Vec<i32>,
    index_water: Vec<i32>,
    z_threshold_water: Option<Vec<[f64; 2]>>,
    z_threshold_silanol: Option<f64>,
    dir: Option<&str>,
    // Only include OHs whose moltype string contains ALL tokens
    moltype_filter: Option<Vec<String>>,
) -> PyResult<()> {
    // Shared inputs
    let directory = dir.unwrap_or("./tmp_atoms/");
    let index_silanol = Arc::new(index_silanol);
    let index_water = Arc::new(index_water);
    let z_threshold_water = z_threshold_water.map(Arc::new);
    let z_threshold_silanol = z_threshold_silanol.map(Arc::new);

    // Output directory
    fs::create_dir_all("./time_corr_results").map_err(|e| {
        PyIOError::new_err(format!("Failed to create output directory: {}", e))
    })?;

    // t0 list
    let t0_values: Vec<usize> = (start_frame as usize .. (end_frame - (corr_times + 1)) as usize)
        .step_by(frame_interval as usize)
        .collect();

    // time offsets per t0
    let t_values: Vec<usize> = (0 .. corr_times as usize)
        .step_by(corr_interval as usize)
        .collect();

    let overall_timer = Instant::now();

    // Parallel over t0 blocks
    t0_values.into_par_iter().try_for_each(|t0| -> PyResult<()> {
        let frame_timer = Instant::now();
        println!("time_corr: starting frame {}", t0);

        // Initialize at t0
        let (_mol_result_orig, mut e_OH, mol_types_orig) = init_e_oh(
            t0,
            (*index_silanol).clone(),
            (*index_water).clone(),
            z_threshold_water.as_ref().map(|arc| (**arc).clone()),
            z_threshold_silanol.as_ref().map(|arc| (**arc).clone()),
            directory,
        )?;

        let num_OHs = e_OH[0].len();
        let num_steps = t_values.len();

        // Build i-mask from moltype_filter at t0
        let i_mask: Vec<bool> = (0..num_OHs)
            .map(|i| {
                if let Some(ref filt) = moltype_filter {
                    filt.iter().all(|tok| mol_types_orig[i].contains(tok))
                } else {
                    true
                }
            })
            .collect();

        if !i_mask.iter().any(|&x| x) {
            return Err(PyValueError::new_err(format!(
                "time_corr(t0={}): no OH indices matched the provided moltype_filter",
                t0
            )));
        }

        // Result buffer (shared with single lock per step)
        let result = Arc::new(Mutex::new(TimeCorrResult::new(num_steps)));

        // Iterate over time offsets (sequential to preserve order)
        for (step_idx, t_offset) in t_values.iter().enumerate() {
            let t = t0 + *t_offset;

            // Read e_OH at time t
            let res_new = read_frame(
                "timecorr".to_string(),
                t,
                (*index_silanol).clone(),
                (*index_water).clone(),
                vec![0.0, 0.0, 1.0],
                None,
                z_threshold_water.as_ref().map(|arc| (**arc).clone()),
                z_threshold_silanol.as_ref().map(|arc| (**arc).clone()),
                directory,
            )?;

            let e_oh_new = match res_new {
                ReadFrameResult::TimeCorr { e_vecs, .. } => e_vecs,
                _ => {
                    return Err(PyValueError::new_err(format!(
                        "time_corr: read_frame returned non-TimeCorr variant at t={}",
                        t
                    )));
                }
            };

            // Update slot-1
            e_OH[1] = e_oh_new;

            // Parallel reduction over filtered i's (auto-correlation: i with itself)
            let (total_c1, total_c2, total_c3, total_norm) = (0..num_OHs)
                .into_par_iter()
                .map(|i| {
                    if !i_mask[i] {
                        return (0.0, 0.0, 0.0, 0.0);
                    }
                    let (c1t, c2t, c3t) = calc_cn(e_OH[0][i], e_OH[1][i]);
                    (c1t, c2t, c3t, 1.0)
                })
                .reduce(
                    || (0.0, 0.0, 0.0, 0.0),
                    |a, b| (a.0 + b.0, a.1 + b.1, a.2 + b.2, a.3 + b.3),
                );

            // Single lock update
            let mut guard = result.lock().unwrap();
            guard.c1[step_idx] += total_c1;
            guard.c2[step_idx] += total_c2;
            guard.c3[step_idx] += total_c3;
            guard.norm_t[step_idx] += total_norm;
        }

        // Finalize averages and serialize
        let mut result = Arc::try_unwrap(result)
            .map_err(|_| PyValueError::new_err("Mutex still has multiple owners"))?
            .into_inner()
            .unwrap();

        result.normalize_in_place();

        let tag = if let Some(ref f) = moltype_filter {
            let joined = f.join("+"); // TODO: sanitize if necessary
            format!("time_corr_{}_{}.bin", t0, joined)
        } else {
            format!("time_corr_{}.bin", t0)
        };

        let file_path = format!("./time_corr_results/{}", tag);
        let file = File::create(&file_path)
            .map_err(|e| PyIOError::new_err(format!("Failed to create {}: {}", file_path, e)))?;
        let writer = BufWriter::new(file);
        bincode::serialize_into(writer, &result)
            .map_err(|e| PyValueError::new_err(format!("Serialization failed (frame {}): {}", t0, e)))?;

        println!(
            "time_corr: frame {} saved → {} (elapsed: {:.2?})",
            t0, file_path, frame_timer.elapsed()
        );
        Ok(())
    })?;

    println!(
        "time_corr: completed all frames (total elapsed: {:.2?})",
        overall_timer.elapsed()
    );
    Ok(())
}
