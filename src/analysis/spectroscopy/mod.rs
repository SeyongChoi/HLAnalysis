pub mod tools;
pub mod dipole_polar;
pub mod vvacf;
pub mod spectrum;

use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

use crate::analysis::spectroscopy::dipole_polar::{
    __pyo3_get_function_cal_dipole_polar_term,
    __pyo3_get_function_compute_dipole_polar_terms_saving,
};

use crate::analysis::spectroscopy::vvacf::__pyo3_get_function_compute_vvacf_from_bin;
use crate::analysis::spectroscopy::spectrum::__pyo3_get_function_compute_spectra;

 
/// Python module: `hlanalysis.analysis.spectroscopy`
///
/// Submodules attached:
/// - `hlanalysis.analysis.spectroscopy.dipole_polar`
/// - `hlanalysis.analysis.spectroscopy.vvacf`
/// - `hlanalysis.analysis.spectroscopy.spectrum`
/// 
#[pymodule]
pub fn spectroscopy(py: Python<'_>, m: &PyModule) -> PyResult<()> {
     /* ---------------------------------------------------------------------
       1. Create submodule `dipole_polar` and register dipole_polar calculation functions
       --------------------------------------------------------------------- */
    let dipole_polar = PyModule::new(py, "dipole_polar")?;
    dipole_polar.add_function(wrap_pyfunction!(cal_dipole_polar_term, dipole_polar)?)?;
    dipole_polar.add_function(wrap_pyfunction!(compute_dipole_polar_terms_saving, dipole_polar)?)?;
    m.add_submodule(&dipole_polar)?;
    py.import("sys")?
        .getattr("modules")?
        .set_item("hlanalysis.analysis.spectroscopy.dipole_polar", &dipole_polar)?;

    /* ---------------------------------------------------------------------
       2. Create submodule `vvacf` and register vvacf functions
       --------------------------------------------------------------------- */
    let vvacf = PyModule::new(py, "vvacf")?;
    vvacf.add_function(wrap_pyfunction!(compute_vvacf_from_bin, vvacf)?)?;
    m.add_submodule(&vvacf)?;
    py.import("sys")?
        .getattr("modules")?
        .set_item("hlanalysis.analysis.spectroscopy.vvacf", &vvacf)?;

    /* ---------------------------------------------------------------------
       3. Create submodule `spectrum` and register spectrum functions
       --------------------------------------------------------------------- */
    let spectrum = PyModule::new(py, "spectrum")?;
    spectrum.add_function(wrap_pyfunction!(compute_spectra, spectrum)?)?;
    m.add_submodule(&spectrum)?;
    py.import("sys")?
        .getattr("modules")?
        .set_item("hlanalysis.analysis.spectroscopy.spectrum", &spectrum)?;

    Ok(())
}