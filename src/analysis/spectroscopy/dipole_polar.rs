use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
use rayon::prelude::*;

use std::fs::{self, File};
use std::io::BufWriter;
use std::path::Path;
use std::sync::Arc;
use std::time::Instant;

use crate::utils::spatial::mic::mic;
use crate::utils::spatial::pycell3::PyCell3;
use crate::utils::spatial::inverse_cell::inverse_cell;

use crate::analysis::types::oh::OHRecord as OHResult;
use crate::analysis::types::dipole_polar::DipolePolarRecord as DipolePolarResult;

use crate::analysis::pre_process::{read_frame, ReadFrameResult};
use crate::analysis::spectroscopy::tools::{dp, norm};


/// Compute dipole–polarizability terms for a set of OH bonds from one frame.
#[pyfunction]
pub fn cal_dipole_polar_term(
    oh_infos: Vec<OHResult>,
    cell: PyCell3,
    cell_inv: PyCell3,
) -> PyResult<Vec<DipolePolarResult>> {
    let cell = &cell.0;
    let cell_inv = &cell_inv.0;

    let dipole_polar_terms: Vec<DipolePolarResult> = 
        oh_infos
        .par_iter()
        .filter_map(|oh_info| {
            let o_pos = oh_info.o_pos?;
            let h_pos = oh_info.h_pos?;
            let o_vel = oh_info.o_vel?;
            let h_vel = oh_info.h_vel?;
            let h_index = oh_info.h_idx;
            let o_index = oh_info.o_idx;
            let oh_type = oh_info.mol_type.clone();
            
            // OH bond distance / norm
            let oh = [h_pos[0]-o_pos[0], h_pos[1]-o_pos[1], h_pos[2]-o_pos[2]];
            let oh_mic = mic(oh, cell, cell_inv);
            let oh_norm = norm(oh_mic);

            // OH bond relative velocity --> O-H stretching is main property for SFG
            let rel_vel = [h_vel[0] - o_vel[0], h_vel[1] - o_vel[1], h_vel[2] - o_vel[2]];

            // polarizability terms for the "2nd order - Susceptibility"
            let v_proj = dp([oh_mic[0] / oh_norm, oh_mic[1] / oh_norm, oh_mic[2] / oh_norm], rel_vel);

            Some(DipolePolarResult{
                rel_vel_z: Some(rel_vel[2]),
                v_proj: Some(v_proj),
                h_pos: Some(h_pos),
                o_pos: Some(o_pos),
                h_idx: h_index,
                o_idx: o_index,
                oh_type: oh_type,
        })
        })
        .collect();

    Ok(dipole_polar_terms)
}

/// Batch compute dipole–polarizability terms across frames and save per-frame results as bincode.
#[pyfunction]
pub fn compute_dipole_polar_terms_saving(
    start: Option<usize>,
    end: Option<usize>,
    frame_interval: Option<usize>,
    index_silanol: Vec<i32>,
    index_water: Vec<i32>,
    surface_normal: Vec<f64>,
    center: Option<f64>, // z_center coordinate
    z_threshold_water: Option<Vec<[f64; 2]>>,
    z_threshold_silanol: Option<f64>,
    dir: Option<&str>,
    output_dir: Option<&str>,
) -> PyResult<()> {

    let index_silanol = Arc::new(index_silanol);
    let index_water = Arc::new(index_water);
    let surface_normal = Arc::new(surface_normal);

    let z_threshold_water = z_threshold_water.map(Arc::new);
    let z_threshold_silanol = z_threshold_silanol.map(Arc::new);

    let start_frame = start.unwrap_or(0);
    let end_frame = end.unwrap_or(usize::MAX);
    let frame_interval = frame_interval.unwrap_or(1);
    let directory = dir.unwrap_or("./tmp_atoms/");
    let output_directory = output_dir.unwrap_or("./dipole_polar_terms/");
    let frames: Vec<usize> = (start_frame..end_frame).step_by(frame_interval).collect();

    // Ensure output directory exists
    fs::create_dir_all(output_directory).map_err(|e| {
        PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Failed to create output directory: {}", e))
    })?;

    frames.into_par_iter().try_for_each(|frame_idx|-> Result<(), PyErr> {
        let timer = Instant::now(); 

        // Read frame in "spectrum" mode (enum result)
        let res = read_frame(
            "spectrum".to_string(),
            frame_idx,
            (*index_silanol).clone(),
            (*index_water).clone(),
            (*surface_normal).clone(),
            center.clone(),
            z_threshold_water.as_ref().map(|arc| (**arc).clone()),
            z_threshold_silanol.as_ref().map(|arc| (**arc).clone()),
            directory,
        ).map_err(|e| {
            PyErr::new::<PyValueError, _>(format!("Frame read error: {:?}", e))
        })?;

        // Unwrap the Spectrum variant and grab fields we need
        let (_n_ohs, cell, oh_infos): (i32, PyCell3, Vec<OHResult>) = match res {
            ReadFrameResult::Spectrum { n_ohs, cell, oh_infos, .. } => (n_ohs, cell, oh_infos),
            other => {
                return Err(PyErr::new::<PyValueError, _>(format!(
                    "Expected Spectrum variant, but got {:?}",
                    other
                )));
            }
        };

        // Compute dipole–polar terms for the frame
        let cell_inv = PyCell3(inverse_cell(&cell.0));
        let dipole_polar_terms = cal_dipole_polar_term(oh_infos, cell, cell_inv).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Dipole term error: {:?}", e))
        })?;

        // Serialize to file
        let file_path = format!("{}/dipole_polar_{}.bin", output_directory, frame_idx);
        let path = Path::new(&file_path);
        let file = File::create(path).map_err(|e| {
            eprintln!("Failed to create file for frame {}: {}", frame_idx, e); // 여기 추가
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("File create error: {}", e))
        })?;
        
        let writer = BufWriter::new(file);

        bincode::serialize_into(writer, &dipole_polar_terms).map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Serialize error: {}", e))
        })?;

        let duration = timer.elapsed(); 
        println!(
            "Frame {} done – elapsed: {}.{:03}s",
            frame_idx,
            duration.as_secs(),
            duration.subsec_millis()
        );

        Ok(())
    })?;

    Ok(())
}