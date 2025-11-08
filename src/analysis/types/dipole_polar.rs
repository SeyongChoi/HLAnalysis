use pyo3::prelude::*;
use serde::{Serialize, Deserialize};

#[pyclass]
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct DipolePolarRecord {
    #[pyo3(get, set)]
    pub rel_vel_z: Option<f64>,
    #[pyo3(get, set)]
    pub v_proj: Option<f64>,
    #[pyo3(get, set)]
    pub h_pos: Option<[f64; 3]>,
    #[pyo3(get, set)]
    pub o_pos: Option<[f64; 3]>,
    #[pyo3(get, set)]
    pub h_idx: i32,
    #[pyo3(get, set)]
    pub o_idx: i32,
    #[pyo3(get, set)]
    pub oh_type: Option<String>,
}

#[pymethods]
impl DipolePolarRecord {
    pub fn __repr__(&self) -> PyResult<String> {
        Ok(serde_json::to_string_pretty(self).unwrap_or_else(|_| "DipolePolarRecord serialization failed".into()))
    }
}

#[pyclass]
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct DipolePolarRecordList { inner: Vec<DipolePolarRecord> }

#[pymethods]
impl DipolePolarRecordList {
    #[new] pub fn new(data: Vec<DipolePolarRecord>) -> Self { Self { inner: data } }
    pub fn __len__(&self) -> usize { self.inner.len() }
    pub fn __getitem__(&self, idx: usize) -> Option<DipolePolarRecord> { self.inner.get(idx).cloned() }
    pub fn __repr__(&self) -> PyResult<String> { Ok(format!("DipolePolarRecordList with {} items", self.inner.len())) }
}