use pyo3::prelude::*;
use std::fs::File;
use std::io::BufReader;
use bincode::deserialize_from;

use crate::atoms::Atoms;
use crate::io::read::read_lammps_dump::read_lammps_dump;

/// Main entry point for reading atomic structure files.
/// Selects the correct reader depending on the specified format.
#[pyfunction]
pub fn reader(
    file_path: &str, 
    format: String, 
    start:Option<usize>, 
    end:Option<usize>, 
    interval:Option<usize>, 
    output_dir:Option<&str>
) -> PyResult<()> {
    match format.as_str() {
        // LAMMPS-dump format
        "lammps-dump" => read_lammps_dump(file_path, start, end, interval, output_dir),

        // Placeholder for other formats (future support)
        // "xyz" => ...
        // "extxyz" => ...
        // "cp2k-xyz" => ...
        // "pdb" => ...

        // Unsupported format
        _ => Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
            "Unsupported format: {}",
            format
        ))),
    }
}

/// Reads a binary file (.bin) and deserializes it into an `Atoms` struct.
#[pyfunction]
pub fn read_atoms_from_bin(bin_file: &str) -> PyResult<Atoms> {
    // Open binary file for reading
    let file = File::open(bin_file)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(format!("Error opening file: {}", e)))?;
    let reader = BufReader::new(file);

    // Deserialize binary data into Atoms struct
    let atoms: Atoms = deserialize_from(reader)
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Error deserializing: {}", e)))?;
    
    Ok(atoms)
}