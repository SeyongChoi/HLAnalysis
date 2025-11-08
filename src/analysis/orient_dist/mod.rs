pub mod compute_angles;
pub mod oh_bond_orient;
pub mod dipole_orient;
pub mod mol_orient;

use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

use crate::analysis::orient_dist::oh_bond_orient::{
    __pyo3_get_function_cal_oh_bond_angles,
    __pyo3_get_function_angle_bins_along_normal,
    __pyo3_get_function_accumulate_angle_bins_along_normal
};

use crate::analysis::orient_dist::dipole_orient::{
    __pyo3_get_function_cal_dipole_angles,
    __pyo3_get_function_dipole_angle_bins_along_normal,
    __pyo3_get_function_accumulate_dipole_angle_bins_along_normal
};

use crate::analysis::orient_dist::mol_orient::{
    __pyo3_get_function_cal_mol_orient,
    __pyo3_get_function_angle_pair_bins_for_normal_range,
};


/// Python module: `hlanalysis.analysis.orient_dist`
///
/// Submodules attached:
/// - `hlanalysis.analysis.orient_dist.oh_bond`
/// - `hlanalysis.analysis.orient_dist.dipole`
/// - `hlanalysis.analysis.orient_dist.mol_orient`
#[pymodule]
pub fn orient_dist(py: Python<'_>, m: &PyModule) -> PyResult<()> {
     /* ---------------------------------------------------------------------
       1. Create submodule `oh_bond` and register OH-bond orientation functions
       --------------------------------------------------------------------- */
    let oh_bond = PyModule::new(py, "oh_bond")?;
    oh_bond.add_function(wrap_pyfunction!(cal_oh_bond_angles, oh_bond)?)?;
    oh_bond.add_function(wrap_pyfunction!(angle_bins_along_normal, oh_bond)?)?;
    oh_bond.add_function(wrap_pyfunction!(accumulate_angle_bins_along_normal, oh_bond)?)?;
    m.add_submodule(&oh_bond)?;
    py.import("sys")?
        .getattr("modules")?
        .set_item("hlanalysis.analysis.orient_dist.oh_bond", &oh_bond)?;

    /* ---------------------------------------------------------------------
       2. Create submodule `dipole` and register dipole orientation functions
       --------------------------------------------------------------------- */
    let dipole = PyModule::new(py, "dipole")?;
    dipole.add_function(wrap_pyfunction!(cal_dipole_angles, dipole)?)?;
    dipole.add_function(wrap_pyfunction!(dipole_angle_bins_along_normal, dipole)?)?;
    dipole.add_function(wrap_pyfunction!(accumulate_dipole_angle_bins_along_normal, dipole)?)?;
    m.add_submodule(&dipole)?;
    py.import("sys")?
        .getattr("modules")?
        .set_item("hlanalysis.analysis.orient_dist.dipole", &dipole)?;

    /* ---------------------------------------------------------------------
       3. Create submodule `mol_orient` and register molecular orientation functions
       --------------------------------------------------------------------- */
    let mol_orient = PyModule::new(py, "mol_orient")?;
    mol_orient.add_function(wrap_pyfunction!(cal_mol_orient, mol_orient)?)?;
    mol_orient.add_function(wrap_pyfunction!(angle_pair_bins_for_normal_range, mol_orient)?)?;
    m.add_submodule(&mol_orient)?;
    py.import("sys")?
        .getattr("modules")?
        .set_item("hlanalysis.analysis.orient_dist.mol_orient", &mol_orient)?;

    Ok(())
}