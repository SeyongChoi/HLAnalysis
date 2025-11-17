use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;

use std::fs::{create_dir_all, File};
use std::io::{BufRead, BufReader, Write};

use crate::atoms::Atoms;
use bincode;

#[pyfunction]
pub fn read_xyz(
    file_path: &str,
    start: Option<usize>, 
    end: Option<usize>, 
    interval: Option<usize>, 
    output_dir: Option<&str>
) -> PyResult<()> {
    // Open the XYZ file
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
    for line_result in reader.lines(){
        let line = line_result.map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("File read error: {}", e))
        })?;

        
        let trimmed = line.trim();
        let is_header = !trimmed.is_empty() && trimmed.parse::<usize>().is_ok();
        
        // New timestep header → process the previous block if any
        if is_header{
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

pub fn parse_timestep(lines: &[String]) -> PyResult<Atoms> {
    // 1. Parse the number of atoms from the first line
    let n_atoms: usize = lines[0]
        .trim()
        .parse()
        .map_err(|e| PyErr::new::<PyValueError, _>(format!("Failed to parse atom count: {}", e)))?;

    // 2. Read the comment line (not used here, but can be stored if needed)
    let _comment = lines[1].trim().to_string();

    // 3. Validate that the block contains enough atom lines
    let expected_len = 2 + n_atoms;
    if lines.len() < expected_len {
        return Err(PyErr::new::<PyValueError, _>(format!(
            "Block too short: expected {} atom lines, but got {}",
            n_atoms,
            lines.len() - 2
        )));
    }

    // 4. Parse atomic symbols and positions
    let mut symbols: Vec<String> = Vec::with_capacity(n_atoms);
    let mut positions: Vec<[f64; 3]> = Vec::with_capacity(n_atoms);

    for (i, line) in lines[2..(2 + n_atoms)].iter().enumerate() {
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 4 {
            return Err(PyErr::new::<PyValueError, _>(format!(
                "Invalid atom line at index {}: '{}'",
                i,
                line
            )));
        }

        let symbol = parts[0].to_string();
        let x: f64 = parts[1].parse().map_err(|e| {
            PyErr::new::<PyValueError, _>(format!("Failed to parse x in '{}': {}", line, e))
        })?;
        let y: f64 = parts[2].parse().map_err(|e| {
            PyErr::new::<PyValueError, _>(format!("Failed to parse y in '{}': {}", line, e))
        })?;
        let z: f64 = parts[3].parse().map_err(|e| {
            PyErr::new::<PyValueError, _>(format!("Failed to parse z in '{}': {}", line, e))
        })?;

        symbols.push(symbol);
        positions.push([x, y, z]);
    }

    // 5. Generate consecutive atom indices (0..n_atoms-1)
    let indices: Vec<i32> = (0..n_atoms as i32).collect();

    // 6. Construct the Atoms struct
    Ok(Atoms::new(
        Some(indices),   // indices
        Some(symbols),   // symbols
        None,            // masses
        None,            // charges
        Some(positions), // positions
        None,            // positions_unwrap
        None,            // forces
        None,            // velocities
        None,            // cell
        None,            // pbc
    ))
}