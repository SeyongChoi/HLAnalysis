use pyo3::prelude::*;
use serde::{Serialize, Deserialize};


#[pyclass]
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct TimeCorrRecord{
    #[pyo3(get, set)]
    pub c1: Vec<f64>,
    #[pyo3(get, set)]
    pub c2: Vec<f64>,
    #[pyo3(get, set)]
    pub c3: Vec<f64>,
    #[pyo3(get, set)]
    pub norm_t: Vec<f64>,
}


#[pymethods]
impl TimeCorrRecord {
    #[new]
    pub fn py_new() -> Self {
        Self::default()
    }

    pub fn __repr__(&self) -> PyResult<String> {
        Ok(serde_json::to_string_pretty(self)
            .unwrap_or_else(|_| "TimeCorrRecord serialization failed".into()))
    }
}

// Rust-side helpers (crate 내부에서 사용)
impl TimeCorrRecord {
    pub fn new(num_steps: usize) -> Self {
        Self {
            c1: vec![0.0; num_steps],
            c2: vec![0.0; num_steps],
            c3: vec![0.0; num_steps],
            norm_t: vec![0.0; num_steps],
        }
    }

    pub fn normalize_in_place(&mut self) {
        for i in 0..self.c1.len() {
            if self.norm_t[i] > 0.0 {
                self.c1[i] /= self.norm_t[i];
                self.c2[i] /= self.norm_t[i];
                self.c3[i] /= self.norm_t[i];
            }
        }
    }
}