mod atom;
mod atoms;
use pyo3::prelude::*;

/*------------------------------------------------------------------------------------ */
#[pymodule]
fn _hlanalysis(_py: Python, m: &PyModule) -> PyResult<()> {
    /* ------------------------------------------------------------------------ */
    m.add_class::<atom::Atom>()?;
    m.add_class::<atoms::Atoms>()?;
    Ok(())
}
