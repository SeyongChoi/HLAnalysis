use pyo3::prelude::*;
use serde::{Serialize, Deserialize};

#[pyclass]
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct OHRecord {
    #[pyo3(get, set)]
    pub o_idx: i32,
    #[pyo3(get, set)]
    pub h_idx: i32,

    #[pyo3(get, set)]
    pub o_pos: Option<[f64;3]>,
    #[pyo3(get, set)]
    pub h_pos: Option<[f64;3]>,

    #[pyo3(get, set)]
    pub o_vel: Option<[f64;3]>,
    #[pyo3(get, set)]
    pub h_vel: Option<[f64;3]>,

    #[pyo3(get, set)]
    pub mol_type: Option<String>,
}

#[pymethods]
impl OHRecord {
    #[new]
    pub fn new(o_idx: i32, h_idx: i32) -> Self {
        Self { o_idx, h_idx, ..Default::default() }
    }
    pub fn __repr__(&self) -> PyResult<String> {
        Ok(serde_json::to_string_pretty(self).unwrap_or_else(|_| "OHRecord serialization failed".into()))
    }
}
