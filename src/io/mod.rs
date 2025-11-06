pub mod read;
pub mod reader;

use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

use crate::io::reader::__pyo3_get_function_reader;
use crate::io::reader::__pyo3_get_function_read_atoms_from_bin;

use crate::io::read::read_lammps_dump::__pyo3_get_function_read_lammps_dump;

pub fn register(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    // ---- functions defined in src/io/reader.rs ----
    m.add_function(wrap_pyfunction!(reader, m)?)?;
    m.add_function(wrap_pyfunction!(read_atoms_from_bin, m)?)?;
    
    // ---- function defined in src/io/read/read_lammps_dump.rs ----
    m.add_function(wrap_pyfunction!(read_lammps_dump, m)?)?;

    Ok(())
}
