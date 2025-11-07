use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

// =====================
// Physical Constants
// =====================
pub const NUM_AVO: f64 = 6.0221408e23;       // Avogadro's number [ea/mol]
pub const C: f64 = 2.99792458e-5;            // Speed of light [cm/fs]
pub const K_B: f64 = 1.380649e-23;           // Boltzmann constant [J/K]
pub const H: f64 = 6.62607015e-34;           // Planck constant [J·s]
pub const H_BAR: f64 = H / (2.0 * std::f64::consts::PI); // Reduced Planck constant [J·s]

// =====================
// Python-callable getters
// =====================

#[pyfunction]
pub fn num_avo() -> f64 {
    NUM_AVO
}

#[pyfunction]
pub fn c() -> f64 {
    C
}

#[pyfunction]
pub fn k_b() -> f64 {
    K_B
}

#[pyfunction]
pub fn h() -> f64 {
    H
}

#[pyfunction]
pub fn h_bar() -> f64 {
    H_BAR
}

// =====================
// Register this module
// =====================

pub fn register(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(num_avo, m)?)?;
    m.add_function(wrap_pyfunction!(c, m)?)?;
    m.add_function(wrap_pyfunction!(k_b, m)?)?;
    m.add_function(wrap_pyfunction!(h, m)?)?;
    m.add_function(wrap_pyfunction!(h_bar, m)?)?;
    Ok(())
}
