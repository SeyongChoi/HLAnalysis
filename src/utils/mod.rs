pub mod histogram;
pub mod physchem;
pub mod spatial;

use pyo3::prelude::*;

pub fn register(py: Python<'_>, m: &PyModule) -> PyResult<()> {
    // Create submodule: utils.histogram
    let submod = PyModule::new(py, "histogram")?;
    histogram::histogram(py, submod)?; // attach its functions
    m.add_submodule(submod)?;

    // Create submodule: utils.physchem
    let submod = PyModule::new(py, "physchem")?;
    physchem::physchem(py, submod)?; // attach its functions
    m.add_submodule(submod)?;
    
    Ok(())
}
