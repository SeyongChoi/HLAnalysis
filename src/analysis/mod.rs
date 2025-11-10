pub mod density_profile;
pub mod orient_dist;
pub mod types;
pub mod spectroscopy;

use pyo3::prelude::*;

/// Register `hlanalysis.analysis` and attach its submodules:
/// - hlanalysis.analysis.density_profile
/// - hlanalysis.analysis.orient_dist
/// - hlanalysis.analysis.spectroscopy
pub fn register(py: Python<'_>, m: &PyModule) -> PyResult<()> {
    /* -------------------------------------------------------
       1. Create & attach `density_profile` submodule
       ------------------------------------------------------- */
    let density_profile_mod = PyModule::new(py, "density_profile")?;
    density_profile::density_profile(py, &density_profile_mod)?; // call its #[pymodule] fn
    m.add_submodule(&density_profile_mod)?;

    // register full Python path in sys.modules
    py.import("sys")?
        .getattr("modules")?
        .set_item("hlanalysis.analysis.density_profile", &density_profile_mod)?;

    /* -------------------------------------------------------
       2. Create & attach `orient_dist` submodule
       ------------------------------------------------------- */
    let orient_dist_mod = PyModule::new(py, "orient_dist")?;
    orient_dist::orient_dist(py, &orient_dist_mod)?; // call its #[pymodule] fn
    m.add_submodule(&orient_dist_mod)?;

    py.import("sys")?
        .getattr("modules")?
        .set_item("hlanalysis.analysis.orient_dist", &orient_dist_mod)?;

    /* -------------------------------------------------------
       3. Create & attach `spectroscopy` submodule
       ------------------------------------------------------- */
    let spectroscopy_mod = PyModule::new(py, "spectroscopy")?;
    spectroscopy::spectroscopy(py, &spectroscopy_mod)?; // call its #[pymodule] fn
    m.add_submodule(&spectroscopy_mod)?;

    py.import("sys")?
        .getattr("modules")?
        .set_item("hlanalysis.analysis.spectroscopy", &spectroscopy_mod)?;

    Ok(())
}