pub mod oh_reorient;

use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

use crate::analysis::reorient_dynamics::oh_reorient::{
    __pyo3_get_function_init_e_oh,
    __pyo3_get_function_time_corr,
};

#[pymodule]
pub fn reorient_dynamics(py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(init_e_oh, m)?)?;
    m.add_function(wrap_pyfunction!(time_corr, m)?)?;
    Ok(())

}