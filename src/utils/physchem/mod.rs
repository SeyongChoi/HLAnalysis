pub mod atomic_data;
pub mod constants;
pub mod converter;

use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

use crate::utils::physchem::atomic_data::__pyo3_get_function_atomic_masses;
use crate::utils::physchem::atomic_data::__pyo3_get_function_atomic_numbers;


#[pymodule]
pub fn physchem(py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(atomic_masses, m)?)?;
    m.add_function(wrap_pyfunction!(atomic_numbers, m)?)?;

    attach(py, m, "constants", constants::register)?;
    attach(py, m, "converter", converter::register)?;
    Ok(())
}

// generic attach helper (reuse from lib.rs)
fn attach<F>(py: Python<'_>, parent: &PyModule, name: &str, reg: F) -> PyResult<()>
where
    F: Fn(Python<'_>, &PyModule) -> PyResult<()>,
{
    let sub = PyModule::new(py, name)?;
    reg(py, &sub)?;
    parent.add_submodule(&sub)?;
    py.import("sys")?
        .getattr("modules")?
        .set_item(format!("hlanalysis.utils.physchem.{name}"), &sub)?;
    Ok(())
}