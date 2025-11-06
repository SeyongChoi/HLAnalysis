use pyo3::prelude::*;

use std::fs::{create_dir_all, File};
use std::io::{BufRead, BufReader, Write};

use crate::atoms::Atoms;
use bincode;



/// Reads a LAMMPS dump file and saves specified timesteps as binary files
/// - Reads the file line-by-line with a large buffer (64 MB)
/// - Splits content into timestep blocks
/// - For blocks within `[start, end)` matching the `interval`, serialize to `{output_dir}/timestep_<idx>.bin`
///
/// # Arguments
/// * `file_path`  - Path to the LAMMPS dump file
/// * `start`      - Starting block index (default: 0)
/// * `end`        - End block index, non-inclusive (default: usize::MAX)
/// * `interval`   - Keep every `interval`-th block after `start` (default: 1)
/// * `output_dir` - Output directory (default: "./tmp_atoms")

#[pyfunction]
pub fn read_lammps_dump(
    file_path: &str,
    start: Option<usize>, 
    end: Option<usize>, 
    interval: Option<usize>, 
    output_dir: Option<&str>
) -> PyResult<()> {
    // Open the LAMMPS dump file
    let file = File::open(file_path)
        .map_err(|e| {PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("File open error: {}", e))})?;
    
    let reader = BufReader::with_capacity(64*1024*1024,file);

    // Ensure the output directory exists, if not create it
    let output_dir = output_dir.unwrap_or("./tmp_atoms");
    if let Err(e) = create_dir_all(output_dir) {
        eprintln!("Error creating directory {}: {}", output_dir, e);
        return Err(PyErr::new::<pyo3::exceptions::PyIOError, _>(format!(
            "Failed to create directory: {}",
            e
        )));
    }
    // Defaults for range and sampling
    let start_idx = start.unwrap_or(0);
    let end_idx = end.unwrap_or(usize::MAX);
    let interval_size = interval.unwrap_or(1);

    // Accumulate lines for the current timestep block
    let mut current_block: Vec<String> = Vec::new();
    let mut block_idx: usize = 0;

    // Scan file line-by-line and detect timestep boundaries
    for line_result in reader.lines() {
        let line = line_result.map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("File read error: {}", e))
        })?;

        // New timestep header → process the previous block if any
        if line.starts_with("ITEM: TIMESTEP") {
            if !current_block.is_empty() {
                if block_idx >= start_idx 
                   && block_idx < end_idx
                   && (block_idx - start_idx) % interval_size == 0
                {
                    // Parse and serialize the selected timestep block
                    if let Ok(timestep) = parse_timestep(&current_block) {
                        let output_file_name = format!("{}/timestep_{}.bin", output_dir, block_idx);
                        match bincode::serialize(&timestep) {
                            Ok(binary_data) => {
                                if let Err(e) = File::create(&output_file_name)
                                    .and_then(|mut file| file.write_all(&binary_data))
                                {
                                    eprintln!("Error writing to file {}: {}", output_file_name, e);
                                }
                            }
                            Err(e) => eprintln!("Error serializing timestep: {}", e),
                        }
                        // Free memory explicitly (not strictly necessary)
                        std::mem::drop(timestep);
                    } else {
                        eprintln!("Error parsing timestep at block {}", block_idx);
                    }
                }

                // Advance to next block; early-exit if we've reached the end
                block_idx += 1;
                if block_idx >= end_idx {
                    break;
                }
                current_block.clear();
            }
        }
        // Keep accumulating lines for the current block
        current_block.push(line);
    }

    // Handle the final block after EOF
    if !current_block.is_empty()
       && block_idx >= start_idx
       && block_idx < end_idx 
       && (block_idx - start_idx) % interval_size == 0 
    {
        if let Ok(timestep) = parse_timestep(&current_block) {
            let output_file_name = format!("{}/timestep_{}.bin", output_dir, block_idx);
            match bincode::serialize(&timestep) {
                Ok(binary_data) => {
                    if let Err(e) = File::create(&output_file_name)
                        .and_then(|mut file| file.write_all(&binary_data))
                    {
                        eprintln!("Error writing to file {}: {}", output_file_name, e);
                    }
                }
                Err(e) => eprintln!("Error serializing timestep: {}", e),
            }
        } else {
            eprintln!("Error parsing timestep at block {}", block_idx);
        }
    }

    Ok(())

}


// Parse a single LAMMPS timestep block (lines) into an 'Atoms' instance
pub fn parse_timestep(lines: &[String]) -> PyResult<Atoms> {
    // Per-timestep accumulators
    let mut indices = Vec::new();
    let mut symbols = Vec::new();
    let mut masses = Vec::new();
    let mut charges = Vec::new();
    let mut positions = Vec::new();
    let mut positions_unwrap = Vec::new();
    let mut velocities = Vec::new();
    let mut forces = Vec::new();
    let mut cell = None;
    let mut pbc = None;

    let mut n_atoms: i32 = 0;
    let mut lines_iter = lines.iter(); 

    while let Some(line) = lines_iter.next() {
        // Get the number of atoms
        if line.starts_with("ITEM: NUMBER OF ATOMS") {
            if let Some(n_atoms_line) = lines_iter.next() {
                n_atoms = n_atoms_line.trim().parse()?;
            }
        }
        // Periodic boundary conditions and cell dimensions
        else if line.starts_with("ITEM: BOX BOUNDS") {
            // Detect periodic boundary conditions
            pbc = if line.contains("pp pp pp") {
                Some([true, true, true])
            } else if line.contains("pp pp ff") {
                Some([true, true, false])
            } else if line.contains("pp ff pp") {
                Some([true, false, true])
            } else if line.contains("ff pp pp") {
                Some([false, true, true])
            } else if line.contains("pp ff ff") {
                Some([true, false, false])
            } else if line.contains("ff pp ff") {
                Some([false, true, false])
            } else if line.contains("ff ff pp") {
                Some([false, false, true])
            } else {
                Some([false, false, false])
            };

            // Variables to hold bounds and tilt factors
            let mut xlo_bound = 0.0;
            let mut xhi_bound = 0.0;
            let mut ylo_bound = 0.0;
            let mut yhi_bound = 0.0;
            let mut zlo_bound = 0.0;
            let mut zhi_bound = 0.0;
            let mut xy = 0.0;
            let mut xz = 0.0;
            let mut yz = 0.0;

            // Read the next three lines for the bounds and tilt factors
            for i in 0..3 {
                if let Some(bounds_line) = lines_iter.next() {
                    let bounds: Vec<f64> = bounds_line
                        .split_whitespace()
                        .map(|v| v.parse().unwrap_or(0.0))
                        .collect();

                    if bounds.len() == 2 {
                        // Orthogonal case: only min and max are present
                        match i {
                            0 => {
                                xlo_bound = bounds[0];
                                xhi_bound = bounds[1];
                            }
                            1 => {
                                ylo_bound = bounds[0];
                                yhi_bound = bounds[1];
                            }
                            2 => {
                                zlo_bound = bounds[0];
                                zhi_bound = bounds[1];
                            }
                            _ => {}
                        }
                    } else if bounds.len() == 3 {
                        // Triclinic case: min, max, and tilt factor
                        match i {
                            0 => {
                                xlo_bound = bounds[0];
                                xhi_bound = bounds[1];
                                xy = bounds[2];
                            }
                            1 => {
                                ylo_bound = bounds[0];
                                yhi_bound = bounds[1];
                                xz = bounds[2];
                            }
                            2 => {
                                zlo_bound = bounds[0];
                                zhi_bound = bounds[1];
                                yz = bounds[2];
                            }
                            _ => {}
                        }
                    }
                }
            }

            // Calculate adjusted bounds for triclinic cells
            let xlo = xlo_bound - f64::min(0.0, f64::min(xy, f64::min(xz, xy + xz)));
            let xhi = xhi_bound - f64::max(0.0, f64::max(xy, f64::max(xz, xy + xz)));
            let ylo = ylo_bound - f64::min(0.0, yz);
            let yhi = yhi_bound - f64::max(0.0, yz);
            let zlo = zlo_bound;
            let zhi = zhi_bound;

            // Construct the cell matrix for the triclinic cell
            let cell_matrix = [
                [xhi - xlo, 0.0, 0.0], // a vector
                [xy, yhi - ylo, 0.0],  // b vector
                [xz, yz, zhi - zlo],   // c vector
            ];

            cell = Some(cell_matrix);
        }
        // Atom data section
        else if line.starts_with("ITEM: ATOMS") {
            // Set attributes to parse atom data
            let attributes: Vec<String> = line
                .split_whitespace()
                .skip(2)
                .map(|s| s.to_string())
                .collect();

            for _ in 0..n_atoms {
                if let Some(atom_line) = lines_iter.next() {
                    let values: Vec<&str> = atom_line.split_whitespace().collect();
                    let mut pos = [0.0; 3];
                    let mut pos_unwrap = [0.0; 3];
                    let mut vel = [0.0; 3];
                    let mut f = [0.0; 3];
                    let mut charge = 0.0;

                    // Parse each atom's data based on the attributes
                    for (i, attribute) in attributes.iter().enumerate() {
                        match attribute.as_str() {
                            "id" => indices.push(values[i].parse().unwrap_or_default()),
                            "element" => symbols.push(values[i].to_string()),
                            "mass" => masses.push(values[i].parse().unwrap_or_default()),
                            "charge" => charge = values[i].parse().unwrap_or_default(),
                            "x" => pos[0] = values[i].parse().unwrap_or_default(),
                            "y" => pos[1] = values[i].parse().unwrap_or_default(),
                            "z" => pos[2] = values[i].parse().unwrap_or_default(),
                            "xu" => pos_unwrap[0] = values[i].parse().unwrap_or_default(),
                            "yu" => pos_unwrap[1] = values[i].parse().unwrap_or_default(),
                            "zu" => pos_unwrap[2] = values[i].parse().unwrap_or_default(),
                            "vx" => vel[0] = values[i].parse().unwrap_or_default(),
                            "vy" => vel[1] = values[i].parse().unwrap_or_default(),
                            "vz" => vel[2] = values[i].parse().unwrap_or_default(),
                            "fx" => f[0] = values[i].parse().unwrap_or_default(),
                            "fy" => f[1] = values[i].parse().unwrap_or_default(),
                            "fz" => f[2] = values[i].parse().unwrap_or_default(),
                            _ => {}
                        }
                    }

                    charges.push(charge);
                    forces.push(f);
                    velocities.push(vel);
                    // Use `pos_unwrap` as fallback if `pos` is `[0.0, 0.0, 0.0]`
                    if pos == [0.0; 3] && pos_unwrap != [0.0; 3] {
                        positions.push(pos_unwrap);
                    } else {
                        positions.push(pos);
                    }
                    positions_unwrap.push(pos_unwrap);
                }
            }
        }
    }

    Ok(Atoms::new(
        Some(indices),
        Some(symbols),
        Some(masses),
        Some(charges),
        Some(positions),
        Some(positions_unwrap),
        Some(forces),
        Some(velocities),
        cell,
        pbc,
    ))
}
