use pyo3::prelude::*;
use rayon::prelude::*;

use std::path::PathBuf;
use std::fs::{self, File};
use std::io::BufReader; 
use bincode::deserialize_from;

use crate::atoms::Atoms;
use crate::utils::histogram::histogram_1d::histogram_1d;


#[pyfunction]
pub fn accumulate_counts(
    target_symbol: &str,   // Atom symbol to filter (e.g., "O")
    target_index: Vec<i32>,// Optional index whitelist (0-based). Empty = all indices.
    bins: Vec<f64>,        // Bin edges for histogram (right-open intervals)
    parallel: bool,        // Kept for API compatibility; CURRENTLY IGNORED (always parallel)
    dir: Option<&str>,     // Directory containing `timestep_*.bin` files
    start: Option<usize>,  // Optional lower bound of timestep number to include
    end: Option<usize>,    // Optional upper bound of timestep number to include
) -> (Vec<f64>, Vec<usize>) {
    // -------------------------------------------------------------------------
    // 1) Prepare z-bin centers from the provided bin edges.
    //    Centers have length (bins.len() - 1).
    // -------------------------------------------------------------------------
    let z_centers: Vec<f64> = bins.windows(2).map(|w| (w[0] + w[1]) / 2.0).collect();

    // Default directory where serialized Atoms are stored.
    let directory = dir.unwrap_or("./tmp_atoms/");

    // Timestep filtering bounds (inclusive).
    let start_t = start.unwrap_or(1);
    let end_t = end.unwrap_or(usize::MAX);

    // -------------------------------------------------------------------------
    // 2) Discover all `timestep_*.bin` files within [start_t, end_t] (inclusive).
    // -------------------------------------------------------------------------
    let files: Vec<PathBuf> = fs::read_dir(directory)
        .unwrap()
        .filter_map(|entry| entry.ok().map(|e| e.path()))
        .filter(|path| {
            // Keep only files with `.bin` extension AND matching `timestep_{N}.bin`.
            path.extension().map(|e| e == "bin").unwrap_or(false)
                && path
                    .file_stem()
                    .and_then(|stem| stem.to_str())
                    .and_then(|s| {
                        if let Some(num_str) = s.strip_prefix("timestep_") {
                            num_str.parse::<usize>().ok()
                        } else {
                            None
                        }
                    })
                    .map(|timestep| (start_t..=end_t).contains(&timestep))
                    .unwrap_or(false)
        })
        .collect();
    // -------------------------------------------------------------------------
    // 3) Accumulate z-positions of the target atoms **in parallel**.
    //
    //    NOTE:
    //    - The `parallel` flag is kept for API compatibility, but not used.
    //    - Each file is deserialized to `Atoms`, then we filter and collect
    //      z-coordinates (position[2]) of atoms that match `target_symbol`
    //      and (optionally) `target_index`.
    // -------------------------------------------------------------------------
    // Accumulate z-positions of the target symbol atoms from the files
    let z_positions: Vec<f64> = files
            .into_par_iter()
            .filter_map(|file_path| {
                // Open and deserialize each file; skip on error.
                let file = File::open(file_path.clone()).ok()?;
                let reader = BufReader::new(file);
                let atoms: Atoms = deserialize_from(reader).ok()?;

                println!("Processing file: {:?}", file_path.clone());

                // We require both positions and symbols to be present.
                if let (Some(positions), Some(symbols)) = (atoms.positions, atoms.symbols) {
                    // Determine the index list:
                    // - If `atoms.indices` is Some but empty, generate 1..=N.
                    // - If `atoms.indices` is None, generate 1..=N.
                    // - Otherwise, use the provided indices.
                    let indices = if let Some(indices) = &atoms.indices {
                        if indices.is_empty() {
                            (1..=symbols.len()).map(|i| i as i32).collect::<Vec<i32>>()
                        } else {
                            indices.clone()
                        }
                    } else {
                        (1..=symbols.len()).map(|i| i as i32).collect::<Vec<i32>>()
                    };
                    
                    // Collect z for atoms that match:
                    // - symbol equals `target_symbol`
                    // - index is in `target_index` (0-based), if a whitelist is provided
                    let mut result = Vec::new();

                    for (index, (symbol, position)) in indices.iter().zip(symbols.iter().zip(positions.iter())) {
                        if symbol == target_symbol && (target_index.is_empty() || target_index.contains(&(index - 1))){
                            result.push(position[2]); // Get the z-coordinate
                        }
                    }
                    if !result.is_empty() {
                        Some(result)
                    } else {
                        None
                    }
                } else {
                    // Skip if atoms data is missing
                    None 
                }
            })
            .flatten() // Flatten Vec<Vec<f64>> -> Vec<f64>
            .collect();

    // -------------------------------------------------------------------------
    // 4) Build a 1D histogram over the collected z-positions.
    //    We return (z_centers, counts) for convenient plotting/use downstream.
    // -------------------------------------------------------------------------
    let (counts, _) = histogram_1d(z_positions, bins.clone());

    (z_centers, counts)
}
