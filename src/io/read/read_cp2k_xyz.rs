use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;

use std::fs::{create_dir_all, File};
use std::io::{BufRead, BufReader, Write};

use crate::atoms::Atoms;
use bincode;

#[pyfunction]
pub fn read_cp2k_xyz(
    file_path: &str,
    start: Option<usize>, 
    end: Option<usize>, 
    interval: Option<usize>, 
    output_dir: Option<&str>,
    velocity: bool,
    vel_file_path: Option<&str>,
    cell: Option<[[f64;3];3]>,
) -> PyResult<()> {
    // Open the position XYZ file
    let file = File::open(file_path)
        .map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(
                format!("File open error (positions): {}", e)
            )
        })?;
    let reader = BufReader::with_capacity(64 * 1024 * 1024, file);

    // Optionally open the velocity XYZ file
    let mut vel_reader_opt: Option<BufReader<File>> = None;
    if velocity {
        let vel_path = vel_file_path.ok_or_else(|| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(
                "Velocity flag is true, but vel_file_path is None",
            )
        })?;

        let vel_file = File::open(vel_path)
            .map_err(|e| {
                PyErr::new::<pyo3::exceptions::PyIOError, _>(
                    format!("File open error (velocities): {}", e)
                )
            })?;
        vel_reader_opt = Some(BufReader::with_capacity(64 * 1024 * 1024, vel_file));
    }

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

    // Accumulate lines for the current timestep block (positions / velocities)
    let mut current_block_pos: Vec<String> = Vec::new();
    let mut current_block_vel: Vec<String> = Vec::new();
    let mut block_idx: usize = 0;

    // Helper closure to process one timestep block
    let mut process_block = |block_idx: usize,
                             block_pos: &Vec<String>,
                             block_vel: Option<&Vec<String>>|
        -> PyResult<()> {
        if block_pos.is_empty() {
            return Ok(());
        }

        // Range and sampling filter
        if block_idx < start_idx || block_idx >= end_idx {
            return Ok(());
        }
        if ((block_idx - start_idx) % interval_size) != 0 {
            return Ok(());
        }

        // Parse positions into Atoms
        let mut atoms = parse_timestep(block_pos)?;

        // If velocities are provided, parse and attach them
        if let Some(v_block) = block_vel {
            // Determine atom count from symbols or positions
            let natoms = if let Some(ref syms) = atoms.symbols {
                syms.len()
            } else if let Some(ref pos) = atoms.positions {
                pos.len()
            } else {
                0
            };

            if natoms == 0 {
                return Err(PyErr::new::<pyo3::exceptions::PyIOError, _>(
                    "Atoms object has no symbols or positions to infer atom count",
                ));
            }

            let velocities = parse_velocity_block(v_block, natoms)?;
            atoms.velocities = Some(velocities);
        }

        // Set cell and PBC if provided
        if cell.is_some() {
            atoms.cell = cell;
            atoms.pbc = Some([true, true, true]);
        }

        // Serialize to .bin
        let output_file_name = format!("{}/timestep_{}.bin", output_dir, block_idx);
        match bincode::serialize(&atoms) {
            Ok(binary_data) => {
                if let Err(e) = File::create(&output_file_name)
                    .and_then(|mut file| file.write_all(&binary_data))
                {
                    eprintln!("Error writing to file {}: {}", output_file_name, e);
                }
            }
            Err(e) => eprintln!("Error serializing timestep {}: {}", block_idx, e),
        }

        // Explicit drop (optional)
        std::mem::drop(atoms);
        Ok(())
    };

    // Scan position file line-by-line and detect timestep boundaries
    // Velocity file is read line-by-line in sync with the position file
    let mut vel_reader_opt = vel_reader_opt; // make mutable
    for line_result in reader.lines() {
        let line_pos = line_result.map_err(|e| {
            PyErr::new::<pyo3::exceptions::PyIOError, _>(format!(
                "File read error (positions): {}",
                e
            ))
        })?;

        // Read the corresponding line from the velocity file, if enabled
        let line_vel_opt = if let Some(ref mut vel_reader) = vel_reader_opt {
            use std::io::Read;
            use std::io::BufRead;

            let mut buf = String::new();
            let n = vel_reader
                .read_line(&mut buf)
                .map_err(|e| {
                    PyErr::new::<pyo3::exceptions::PyIOError, _>(format!(
                        "File read error (velocities): {}",
                        e
                    ))
                })?;
            if n == 0 {
                return Err(PyErr::new::<pyo3::exceptions::PyIOError, _>(
                    "Velocity file ended before position file (line count mismatch)",
                ));
            }
            Some(buf)
        } else {
            None
        };

        let trimmed = line_pos.trim();
        let is_header = !trimmed.is_empty() && trimmed.parse::<usize>().is_ok();

        // New timestep header → process the previous block if any
        if is_header {
            if !current_block_pos.is_empty() {
                let vel_block_ref = if velocity {
                    Some(&current_block_vel)
                } else {
                    None
                };

                if let Err(e) = process_block(block_idx, &current_block_pos, vel_block_ref) {
                    eprintln!("Error processing timestep {}: {}", block_idx, e);
                }

                // Advance to next block; early-exit if we've reached the end
                block_idx += 1;
                if block_idx >= end_idx {
                    break;
                }
                current_block_pos.clear();
                current_block_vel.clear();
            }
        }

        // Accumulate position line
        current_block_pos.push(line_pos);

        // Accumulate velocity line (if enabled)
        if let Some(line_vel) = line_vel_opt {
            current_block_vel.push(line_vel);
        }
    }

    // Handle the final block after EOF
    if !current_block_pos.is_empty()
        && block_idx >= start_idx
        && block_idx < end_idx
    {
        let vel_block_ref = if velocity {
            Some(&current_block_vel)
        } else {
            None
        };

        if let Err(e) = process_block(block_idx, &current_block_pos, vel_block_ref) {
            eprintln!("Error processing final timestep {}: {}", block_idx, e);
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


/// Parse a single CP2K-style XYZ block as velocities (vx, vy, vz).
/// The format is assumed to be:
///   line 0: number of atoms (N)
///   line 1: comment
///   line 2..(2+N): "Symbol vx vy vz" (or "vx vy vz")
fn parse_velocity_block(lines: &[String], expected_n_atoms: usize) -> PyResult<Vec<[f64; 3]>> {
    if lines.len() < 3 {
        return Err(PyErr::new::<PyValueError, _>(format!(
            "Velocity block too short: got {} lines",
            lines.len()
        )));
    }

    // Parse the number of atoms from the first line
    let n_atoms: usize = lines[0]
        .trim()
        .parse()
        .map_err(|e| PyErr::new::<PyValueError, _>(format!(
            "Failed to parse atom count in velocity block: {}",
            e
        )))?;

    if n_atoms != expected_n_atoms {
        return Err(PyErr::new::<PyValueError, _>(format!(
            "Mismatch in atom count between positions ({}) and velocities ({})",
            expected_n_atoms, n_atoms
        )));
    }

    let expected_len = 2 + n_atoms;
    if lines.len() < expected_len {
        return Err(PyErr::new::<PyValueError, _>(format!(
            "Velocity block too short: expected {} atom lines, but got {}",
            n_atoms,
            lines.len() - 2
        )));
    }

    let mut velocities: Vec<[f64; 3]> = Vec::with_capacity(n_atoms);

    for (i, line) in lines[2..(2 + n_atoms)].iter().enumerate() {
        let parts: Vec<&str> = line.split_whitespace().collect();

        // Two possible formats:
        // 1) "Sym vx vy vz"  -> len == 4, offset = 1
        // 2) "vx vy vz"      -> len == 3, offset = 0
        let (offset, min_len) = (1usize, 4usize);
        let (offset, min_len) = if parts.len() == 4 {
            (1usize, 4usize)
        } else if parts.len() == 3 {
            (0usize, 3usize)
        } else {
            return Err(PyErr::new::<PyValueError, _>(format!(
                "Invalid velocity line at index {}: '{}'",
                i,
                line
            )));
        };

        if parts.len() < min_len {
            return Err(PyErr::new::<PyValueError, _>(format!(
                "Invalid velocity line at index {}: '{}'",
                i,
                line
            )));
        }

        let vx: f64 = parts[offset]
            .parse()
            .map_err(|e| PyErr::new::<PyValueError, _>(format!(
                "Failed to parse vx in '{}': {}",
                line, e
            )))?;
        let vy: f64 = parts[offset + 1]
            .parse()
            .map_err(|e| PyErr::new::<PyValueError, _>(format!(
                "Failed to parse vy in '{}': {}",
                line, e
            )))?;
        let vz: f64 = parts[offset + 2]
            .parse()
            .map_err(|e| PyErr::new::<PyValueError, _>(format!(
                "Failed to parse vz in '{}': {}",
                line, e
            )))?;

        velocities.push([vx, vy, vz]);
    }

    Ok(velocities)
}