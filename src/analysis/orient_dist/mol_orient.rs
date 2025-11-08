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
use crate::utils::histogram::generate_bins::generate_angle_pair_bins_for_normal_range;

use crate::analysis::orient_dist::compute_angles::compute_hh_hoh_angles;


/// Compute molecular orientation (θ_HH and θ_HOH) for each timestep.
///
/// Notes:
/// - Input `.bin` files are detected by the pattern `timestep_<N>.bin` and filtered by [start, end] & interval.
/// - Files are stably **sorted by timestep index** to keep a deterministic order.
#[pyfunction]
pub fn cal_mol_orient(
    water_map: HashMap<i32, (i32, i32, i32)>, // {idx: (o_idx, h1_idx, h2_idx)}
    surface_normal: Vec<f64>,
    parallel: bool,                           // kept for API compatibility; ignored
    dir: Option<&str>,                        // input directory
    start_ts: Option<usize>,                  // inclusive
    end_ts: Option<usize>,                    // inclusive
    interval_ts: Option<usize>,               // sampling interval
    center: Option<f64>,                      // optional z-center
    mode: Option<&str>,                       // mode -> 'cos' or 'angle'
) -> PyResult<Vec<Vec<(i32, i32, i32, f64, f64)>>> {
    // Basic validation
    if water_map.is_empty() || surface_normal.is_empty() {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "Input data must not be empty (water_map, surface_normal).",
        ));
    }
    
    let directory = dir.unwrap_or("./tmp_atoms/");
    let mode = mode.unwrap_or("angle");
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

    // Parallel processing of all selected files
    let all_angles: Vec<Vec<(i32, i32, i32, f64, f64)>> = files
        .par_iter()
        .filter_map(|(_, file_path)| {
            let file = File::open(file_path).ok()?;
            let reader = BufReader::new(file);
            let atoms: Atoms = deserialize_from(reader).ok()?;

            let cell = atoms.cell.as_ref()?;           // Option<[[f64;3];3]> → &[[f64;3];3]
            let inv_cell = inverse_cell(cell); 

            // Compute angles for this frame
            Some(compute_hh_hoh_angles(
                &atoms,
                &water_map,
                &surface_normal, 
                cell, 
                &inv_cell, 
                center, 
                mode))
        })
        .collect();
    
    Ok(all_angles)

}

/// Generate 2D angle–angle histograms (θ_dipole, θ_HH) within z-range slices.
#[pyfunction]
pub fn angle_pair_bins_for_normal_range(
    all_angles: Vec<Vec<(i32, i32, i32, f64, f64)>>,
    bins_dipole: Vec<f64>,
    bins_hh: Vec<f64>,
    z_range: Vec<[f64; 2]>,       // [[z_min, z_max], ...]
    parallel: bool,               // kept for API compatibility; ignored
    dir: Option<&str>,            // input directory
    start_ts: Option<usize>,      // inclusive
    end_ts: Option<usize>,        // inclusive
    interval_ts: Option<usize>,   // sampling interval
) -> PyResult<Vec<Vec<Vec<usize>>>> {

    if all_angles.is_empty() || bins_dipole.is_empty() || bins_hh.is_empty() {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "Input vectors must not be empty (all_angles, bins_dipole, bins_hh).",
        ));
    }

    let directory = dir.unwrap_or("./tmp_atoms/");
    let start_t = start_ts.unwrap_or(1);
    let end_t = end_ts.unwrap_or(usize::MAX);
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
    let all_angle_bins: Vec<Vec<Vec<usize>>> = files
        .par_iter()
        .zip(all_angles.par_iter())
        .filter_map(|((_, file_path), angles)| {
            let file = File::open(file_path).ok()?;
            let reader = BufReader::new(file);
            let atoms: Atoms = deserialize_from(reader).ok()?;

            // Atoms 데이터가 유효하다면 각도 계산
            Some(generate_angle_pair_bins_for_normal_range(
                &atoms,
                angles,
                &z_range,
                &bins_dipole,
                &bins_hh,
            ))
        })
        .collect();
    
    Ok(all_angle_bins)
}