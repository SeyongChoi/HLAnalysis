pub mod accumulate_counts;

use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

use crate::analysis::density_profile::accumulate_counts::__pyo3_get_function_accumulate_counts;

#[pymodule]
pub fn density_profile(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(accumulate_counts, m)?)?;
    Ok(())
}