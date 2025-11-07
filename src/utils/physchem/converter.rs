use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

/* -------------------------------------------
   Conversion constants
------------------------------------------- */
// Time
pub const NS2FS: f64 = 1.0e6;               // ns → fs
pub const NS2PS: f64 = 1.0e3;               // ns → ps
// Length
pub const A2CM: f64 = 1.0e-8;               // Å → cm
pub const A2M: f64 = 1.0e-10;               // Å → m
pub const BOHR2A: f64 = 5.29177208590000E-01; // Bohr → Å
// Energy
pub const HARTREE2EV: f64 = 2.72113838565563E+01;  // Hartree → eV
pub const HARTREE2KJ: f64 = 2.62549961709828E+03;  // Hartree → kJ/mol
pub const HARTREE2KCAL: f64 = 6.27509468713739E+02; // Hartree → kcal/mol

/* -------------------------------------------
   Time conversions
------------------------------------------- */
#[pyfunction] pub fn fs_to_ps() -> f64 { 1.0 / NS2PS }
#[pyfunction] pub fn fs_to_ns() -> f64 { 1.0 / NS2FS }
#[pyfunction] pub fn ps_to_ns() -> f64 { 1.0 / NS2PS }
#[pyfunction] pub fn ps_to_fs() -> f64 { NS2PS }
#[pyfunction] pub fn ns_to_fs() -> f64 { NS2FS }
#[pyfunction] pub fn ns_to_ps() -> f64 { NS2PS }

/* -------------------------------------------
   Length conversions
------------------------------------------- */
#[pyfunction] pub fn ang_to_cm() -> f64 { A2CM }
#[pyfunction] pub fn cm_to_ang() -> f64 { 1.0 / A2CM }
#[pyfunction] pub fn ang_to_m() -> f64 { A2M }
#[pyfunction] pub fn m_to_ang() -> f64 { 1.0 / A2M }
#[pyfunction] pub fn bohr_to_ang() -> f64 { BOHR2A }
#[pyfunction] pub fn ang_to_bohr() -> f64 { 1.0 / BOHR2A }

/* -------------------------------------------
   Energy conversions
------------------------------------------- */
#[pyfunction] pub fn hartree_to_ev() -> f64 { HARTREE2EV }
#[pyfunction] pub fn ev_to_hartree() -> f64 { 1.0 / HARTREE2EV }
#[pyfunction] pub fn hartree_to_kj() -> f64 { HARTREE2KJ }
#[pyfunction] pub fn kj_to_hartree() -> f64 { 1.0 / HARTREE2KJ }
#[pyfunction] pub fn hartree_to_kcal() -> f64 { HARTREE2KCAL }
#[pyfunction] pub fn kcal_to_hartree() -> f64 { 1.0 / HARTREE2KCAL }

/* -------------------------------------------
   PyO3 registration
------------------------------------------- */
pub fn register(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    // Time
    m.add_function(wrap_pyfunction!(fs_to_ps, m)?)?;
    m.add_function(wrap_pyfunction!(fs_to_ns, m)?)?;
    m.add_function(wrap_pyfunction!(ps_to_ns, m)?)?;
    m.add_function(wrap_pyfunction!(ps_to_fs, m)?)?;
    m.add_function(wrap_pyfunction!(ns_to_fs, m)?)?;
    m.add_function(wrap_pyfunction!(ns_to_ps, m)?)?;
    // Length
    m.add_function(wrap_pyfunction!(ang_to_cm, m)?)?;
    m.add_function(wrap_pyfunction!(cm_to_ang, m)?)?;
    m.add_function(wrap_pyfunction!(ang_to_m, m)?)?;
    m.add_function(wrap_pyfunction!(m_to_ang, m)?)?;
    m.add_function(wrap_pyfunction!(bohr_to_ang, m)?)?;
    m.add_function(wrap_pyfunction!(ang_to_bohr, m)?)?;
    // Energy
    m.add_function(wrap_pyfunction!(hartree_to_ev, m)?)?;
    m.add_function(wrap_pyfunction!(ev_to_hartree, m)?)?;
    m.add_function(wrap_pyfunction!(hartree_to_kj, m)?)?;
    m.add_function(wrap_pyfunction!(kj_to_hartree, m)?)?;
    m.add_function(wrap_pyfunction!(hartree_to_kcal, m)?)?;
    m.add_function(wrap_pyfunction!(kcal_to_hartree, m)?)?;
    Ok(())
}
