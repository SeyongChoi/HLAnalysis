pub mod read;
pub mod reader;

use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

pub fn register(py: Python<'_>, m: &PyModule) -> PyResult<()> {
    // ---- functions defined in src/io/reader.rs ----
    m.add_function(wrap_pyfunction!(reader, m)?)?;
    m.add_function(wrap_pyfunction!(read_atoms_from_bin, m)?)?;
    
    // ---- function defined in src/io/read/read_lammps_dump.rs ----
    m.add_function(wrap_pyfunction!(read_lammps_dump, m)?)?;

    Ok(())
}
