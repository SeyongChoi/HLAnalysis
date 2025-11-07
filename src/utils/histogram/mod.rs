pub mod generate_bins;
pub mod histogram_1d;

use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

use crate::utils::histogram::histogram_1d::__pyo3_get_function_histogram_1d;
use crate::utils::histogram::generate_bins::__pyo3_get_function_generate_bins;

#[pymodule]
pub fn histogram(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(histogram_1d, m)?)?;
    m.add_function(wrap_pyfunction!(generate_bins, m)?)?;
    Ok(())
}