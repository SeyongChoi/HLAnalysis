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

use crate::analysis::hbond::hb_find::find_hb_per_frame;

use crate::analysis::types::molecule::{
    MoleculeRecord as MoleculeResult,
    HydrogenBondPartner,
};


#[inline]
fn matches_target_type(mol_type: &str, filter: &Option<Vec<String>>) -> bool {
    match filter {
        None => true,
        Some(v) if v.is_empty() => true,
        Some(v) => v.iter().all(|t| mol_type.contains(t)),
    }
}

#[pyfunction]
pub fn hbond_statistics(
    start_frame: i32,
    end_frame: i32,
    frame_interval: i32,
    index_silanol: Vec<i32>,
    index_water: Vec<i32>,
    surface_normal: Vec<f64>,
    center: Option<f64>,
    target_type: Option<Vec<String>>,
    z_threshold_water: Option<Vec<[f64; 2]>>,
    z_threshold_silanol: Option<f64>,
    hb_conditions: Vec<f64>,
    dir: Option<&str>,
) -> PyResult<()> {
    let overall_timer = Instant::now();

    // Share immutable inputs across threads
    let index_silanol = Arc::new(index_silanol);
    let index_water = Arc::new(index_water);
    let surface_normal = Arc::new(surface_normal);
    let z_threshold_water = Arc::new(z_threshold_water);
    let hb_conditions = Arc::new(hb_conditions);
    let target_type_filter = Arc::new(target_type);
    let directory = Arc::new(dir.unwrap_or("./tmp_atoms/").to_string());

    let frame_indices: Vec<usize> = (start_frame as usize..=end_frame as usize)
        .step_by(frame_interval as usize)
        .collect();

    let num_frames = frame_indices.len();
    if num_frames == 0 {
        println!("No frames to process.");
        return Ok(());
    }

    println!(
        "Processing {} frames from {} to {}...",
        num_frames, start_frame, end_frame
    );

    // Per-frame: (4 HB counts, # of target molecules)
    let results: PyResult<Vec<(usize, usize, usize, usize, usize)>> = Python::with_gil(|py| {
        py.allow_threads(|| {
            frame_indices
                .par_iter()
                .map(|&t0| -> PyResult<(usize, usize, usize, usize, usize)> {
                    Python::with_gil(|_py| {
                        let molecules = find_hb_per_frame(
                            t0,
                            (*index_silanol).clone(),
                            (*index_water).clone(),
                            (*surface_normal).clone(),
                            center,
                            (*z_threshold_water).clone(),
                            z_threshold_silanol,
                            (*hb_conditions).clone(),
                            &directory,
                        )?;

                        let mut donated_to_silanol = 0;
                        let mut accepted_from_silanol = 0;
                        let mut donated_to_water = 0;
                        let mut accepted_from_water = 0;
                        let mut target_molecule_count = 0; // Count of molecules matching the filter

                        for molecule in &molecules {
                            let should_process =
                                matches_target_type(molecule.mol_type.as_deref().unwrap_or(""), &*target_type_filter);

                            if should_process {
                                target_molecule_count += 1;

                                if let Some(partners) = molecule.h_bond_partners.as_ref() {
                                    for partner in partners {
                                        let hb_ty = partner.h_bond_type.as_deref().unwrap_or("");
                                        let is_s = partner
                                            .partner_moltype
                                            .as_deref()
                                            .map(|s| s.starts_with('S'))
                                            .unwrap_or(false);
                                        let is_w = partner
                                            .partner_moltype
                                            .as_deref()
                                            .map(|s| s.starts_with('W'))
                                            .unwrap_or(false);

                                        match hb_ty {
                                            "donor" => {
                                                if is_s { donated_to_silanol += 1; }
                                                else if is_w { donated_to_water += 1; }
                                            }
                                            "acceptor" => {
                                                if is_s { accepted_from_silanol += 1; }
                                                else if is_w { accepted_from_water += 1; }
                                            }
                                            _ => {}
                                        }
                                    }
                                }
                            }
                        }

                        Ok((
                            donated_to_silanol,
                            accepted_from_silanol,
                            donated_to_water,
                            accepted_from_water,
                            target_molecule_count,
                        ))
                    })
                })
                .collect()
        })
    });

    let collected_counts = results?;

    // Sum over all frames: total HB counts and total # of target molecules
    let totals = collected_counts.iter().fold(
        (0usize, 0usize, 0usize, 0usize, 0usize),
        |acc, val| (acc.0 + val.0, acc.1 + val.1, acc.2 + val.2, acc.3 + val.3, acc.4 + val.4),
    );

    let total_target_molecules = totals.4;

    println!(
        "\n--- Hydrogen Bond Statistics (Target filter = {:?}) ---",
        &*target_type_filter
    );
    println!("Processed Frames: {}", num_frames);
    println!(
        "Total Target Molecules (summed over frames): {}",
        total_target_molecules
    );
    println!(
        "Average Target Molecules / Frame: {:.2}",
        total_target_molecules as f64 / num_frames as f64
    );
    println!("---------------------------------------------------");

    // [A] Molecule-weighted mean (denominator = total # of target molecules across frames)
    if total_target_molecules > 0 {
        let avg_d_s = totals.0 as f64 / total_target_molecules as f64;
        let avg_a_s = totals.1 as f64 / total_target_molecules as f64;
        let avg_d_w = totals.2 as f64 / total_target_molecules as f64;
        let avg_a_w = totals.3 as f64 / total_target_molecules as f64;

        println!("[Molecule-weighted] Avg D(->S) / molecule: {:.3}", avg_d_s);
        println!("[Molecule-weighted] Avg A(<-S) / molecule: {:.3}", avg_a_s);
        println!("[Molecule-weighted] Avg D(->W) / molecule: {:.3}", avg_d_w);
        println!("[Molecule-weighted] Avg A(<-W) / molecule: {:.3}", avg_a_w);
    } else {
        println!("No molecules of the target type were found.");
    }

    // [B] Frame-averaged mean: average of per-frame means (each frame weighted equally)
    let mut framewise_means: Vec<(f64, f64, f64, f64)> = Vec::with_capacity(num_frames);
    for (d_s, a_s, d_w, a_w, n) in &collected_counts {
        let n = *n as f64;
        if n > 0.0 {
            framewise_means.push((
                *d_s as f64 / n,
                *a_s as f64 / n,
                *d_w as f64 / n,
                *a_w as f64 / n,
            ));
        }
    }
    if !framewise_means.is_empty() {
        let (sum_ds, sum_as, sum_dw, sum_aw) = framewise_means.iter().fold(
            (0.0, 0.0, 0.0, 0.0),
            |acc, v| (acc.0 + v.0, acc.1 + v.1, acc.2 + v.2, acc.3 + v.3),
        );
        let denom = framewise_means.len() as f64;
        println!(
            "[Frame-averaged]   Avg D(->S) / molecule: {:.3}",
            sum_ds / denom
        );
        println!(
            "[Frame-averaged]   Avg A(<-S) / molecule: {:.3}",
            sum_as / denom
        );
        println!(
            "[Frame-averaged]   Avg D(->W) / molecule: {:.3}",
            sum_dw / denom
        );
        println!(
            "[Frame-averaged]   Avg A(<-W) / molecule: {:.3}",
            sum_aw / denom
        );
    }

    println!("---------------------------------------------------");
    println!("\nTotal elapsed time: {:.2?}", overall_timer.elapsed());
    Ok(())
}