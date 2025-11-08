use pyo3::prelude::*;
use rayon::prelude::*;

use std::sync::Mutex;
use std::io::BufReader; 
use std::path::PathBuf;
use std::fs::{self, File};
use std::collections::HashMap;

use bincode::deserialize_from;

use crate::atoms::Atoms;
use crate::utils::spatial::inverse_cell::inverse_cell;
use crate::utils::histogram::histogram_1d::histogram_1d;
use crate::utils::histogram::generate_bins::generate_angle_bins_along_normal;

use crate::analysis::orient_dist::compute_angles::compute_oh_angles;

/// Compute OH-bond angles relative to a given surface normal for each timestep.
///
/// Notes:
/// - Input `.bin` files are detected by the pattern `timestep_<N>.bin` and filtered by [start, end] & interval.
/// - Files are stably **sorted by timestep index** to keep a deterministic order.
#[pyfunction]
pub fn cal_oh_bond_angles(
    oh_bond_map: HashMap<i32, (i32, i32)>, // {idx: (o_idx, h_idx)}
    surface_normal: Vec<f64>,
    parallel: bool,                         // kept for API compatibility; ignored
    dir: Option<&str>,                      // input directory
    start_ts: Option<usize>,                // inclusive
    end_ts: Option<usize>,                  // inclusive
    interval_ts: Option<usize>,             // sampling interval
    center: Option<f64>,                    // optional z-center
) -> PyResult<Vec<Vec<(i32, i32, f64)>>> {
    // Basic validation
    if oh_bond_map.is_empty() || surface_normal.is_empty() {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "Input data must not be empty (oh_bond_map, surface_normal).",
        ));
    }

    let directory = dir.unwrap_or("./tmp_atoms/");
    // Read all the bin files in the directory
    let start_t = start_ts.unwrap_or(1);
    let end_t = end_ts.unwrap_or(usize::MAX);
    let interval = interval_ts.unwrap_or(1);

    // Collect candidate files, keep only timestep_* bins in the requested range, and sort by timestep index.
    let mut files: Vec<(usize, PathBuf)> = fs::read_dir(directory)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Failed to read dir: {e}")))?
        .filter_map(|entry| entry.ok().map(|e| e.path()))
        .filter_map(|path| {
            let is_bin = path.extension().map(|e| e == "bin").unwrap_or(false);
            if !is_bin { return None; }
            let stem = path.file_stem()?.to_str()?;
            let n = stem.strip_prefix("timestep_")?.parse::<usize>().ok()?;
            Some((n, path))
        })
        .filter(|(timestep, _)| *timestep >= start_t && *timestep <= end_t && (*timestep - start_t) % interval == 0)
        .collect();

    // Stable sort by timestep number to enforce deterministic processing order.
    files.sort_by_key(|(n, _)| *n);

    // Parallel processing of all selected files.
    let all_angles: Vec<Vec<(i32, i32, f64)>> = files
        .par_iter()
        .filter_map(|(_, file_path)| {
            let file   = File::open(file_path).ok()?;
            let reader = BufReader::new(file);
            let atoms: Atoms = deserialize_from(reader).ok()?;

            let cell = atoms.cell.as_ref()?;           // Option<[[f64;3];3]> → &[[f64;3];3]
            let inv_cell = inverse_cell(cell);

            // Compute angles for this frame
            Some(compute_oh_angles(
                &atoms,
                &oh_bond_map,
                &surface_normal,
                cell,
                &inv_cell,
                center,
            ))
        })
        .collect();

    Ok(all_angles)
}

/// Bin the previously computed angles along a surface-normal coordinate for each frame.
///
/// Notes:
/// - The `all_angles` vector is zipped with the sorted file list by index order.
#[pyfunction]
pub fn angle_bins_along_normal(
    all_angles: Vec<Vec<(i32, i32, f64)>>,
    bins_normal: Vec<f64>,
    centers_normal: Vec<f64>,
    parallel: bool,                  // kept for API compatibility; ignored
    dir: Option<&str>,               // input directory
    start_ts: Option<usize>,         // inclusive
    end_ts: Option<usize>,           // inclusive
    interval_ts: Option<usize>,      // sampling interval
) -> PyResult<Vec<Vec<Vec<f64>>>> {
    if all_angles.is_empty() || bins_normal.is_empty() || centers_normal.is_empty() {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "Input vectors must not be empty (all_angles, bins_normal, centers_normal).",
        ));
    }

    let directory = dir.unwrap_or("./tmp_atoms/");
    let start_t  = start_ts.unwrap_or(1);
    let end_t    = end_ts.unwrap_or(usize::MAX);
    let interval = interval_ts.unwrap_or(1);

    // Collect, filter, and sort files by timestep index (same as above).
    let mut files: Vec<(usize, PathBuf)> = fs::read_dir(directory)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Failed to read dir: {e}")))?
        .filter_map(|entry| entry.ok().map(|e| e.path()))
        .filter_map(|path| {
            let is_bin = path.extension().map(|e| e == "bin").unwrap_or(false);
            if !is_bin { return None; }
            let stem = path.file_stem()?.to_str()?;
            let n = stem.strip_prefix("timestep_")?.parse::<usize>().ok()?;
            Some((n, path))
        })
        .filter(|(timestep, _)| *timestep >= start_t && *timestep <= end_t && (*timestep - start_t) % interval == 0)
        .collect();

    files.sort_by_key(|(n, _)| *n);

    // Zip files with all_angles by index; parallel map into per-frame 2D bins.
    let all_angle_bins: Vec<Vec<Vec<f64>>> = files
        .par_iter()
        .zip(all_angles.par_iter())
        .filter_map(|((_, file_path), angles)| {
            let file   = File::open(file_path).ok()?;
            let reader = BufReader::new(file);
            let atoms: Atoms = deserialize_from(reader).ok()?;

            Some(generate_angle_bins_along_normal(
                &atoms,
                angles,
                &bins_normal,
                &centers_normal,
            ))
        })
        .collect();

    Ok(all_angle_bins)
}

/// Accumulate angle histograms across frames, returning a 2D array:
/// `[normal_bin_index][angle_bin_index] -> accumulated counts`
///
/// Notes:
/// - The histogram over angles is computed per normal bin using `histogram_1d`.
#[pyfunction]
pub fn accumulate_angle_bins_along_normal(
    all_angles_bins: Vec<Vec<Vec<f64>>>,
    centers_normal: Vec<f64>,
    bins_angles: Vec<f64>,
    centers_angles: Vec<f64>,
) -> Vec<Vec<f64>> {
    let num_normal = centers_normal.len();
    let num_angles = centers_angles.len();

    // Global accumulator across frames
    let angle_avg_bins = Mutex::new(vec![vec![0.0; num_angles]; num_normal]);

    // Parallel over frames
    all_angles_bins.into_par_iter().for_each(|angle_bins| {
        // Per-frame temporary accumulator
        let mut angle_bins_temp = vec![vec![0.0; num_angles]; num_normal];

        for (bin_idx, angle_bin) in angle_bins.iter().enumerate() {
            if !angle_bin.is_empty() {
                let (counts, _) = histogram_1d(angle_bin.to_vec(), bins_angles.clone());
                for (j, &count) in counts.iter().enumerate() {
                    angle_bins_temp[bin_idx][j] = count as f64;
                }
            }
        }

        // Reduce into the global accumulator
        let mut lock = angle_avg_bins.lock().expect("accumulator poisoned");
        for i in 0..num_normal {
            for j in 0..num_angles {
                lock[i][j] += angle_bins_temp[i][j];
            }
        }
    });

    angle_avg_bins.into_inner().expect("accumulator poisoned")
}

