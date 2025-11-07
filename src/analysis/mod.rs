pub mod density_profile;

use pyo3::prelude::*;

pub fn register(py: Python<'_>, m: &PyModule) -> PyResult<()> {
    // Create submodule for density_profile
    let density_profile_mod = PyModule::new(py, "density_profile")?;
    density_profile::density_profile(py, density_profile_mod)?; // attach functions
    m.add_submodule(density_profile_mod)?;

    
    Ok(())
}
