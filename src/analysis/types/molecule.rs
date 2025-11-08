use pyo3::prelude::*;
use serde::{Serialize, Deserialize};

use super::enums::{MolKind, HBondKind};

#[pyclass]
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct HydrogenBondPartner {
    #[pyo3(get, set)]
    pub partner_o_idx: i32,             // Oxygen index of the partner molecule
    #[pyo3(get, set)]
    pub partner_moltype: Option<String>, // Molecule type of the partner
    #[pyo3(get, set)]
    pub h_bond_type: HBondKind,          // "donor" or "acceptor"
}

#[pyclass]
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct MoleculeRecord {
    // --- Always valid minimal fields ---
    #[pyo3(get, set)]
    pub o_idx: i32,
    #[pyo3(get, set)]
    pub h1_idx: Option<i32>,
    #[pyo3(get, set)]
    pub h2_idx: Option<i32>,

    // --- Geometry (if needed) ---
    #[pyo3(get, set)]
    pub o_pos: Option<[f64; 3]>,
    #[pyo3(get, set)]
    pub h1_pos: Option<[f64; 3]>,
    #[pyo3(get, set)]
    pub h2_pos: Option<[f64; 3]>,

    // --- Velocity (used in correlation/spectroscopy analysis) ---
    #[pyo3(get, set)]
    pub o_vel: Option<[f64; 3]>,
    #[pyo3(get, set)]
    pub h1_vel: Option<[f64; 3]>,
    #[pyo3(get, set)]
    pub h2_vel: Option<[f64; 3]>,

    // --- Direction vectors (used in orientation/spectroscopy analysis) ---
    #[pyo3(get, set)]
    pub e_oh1: Option<[f64; 3]>,
    #[pyo3(get, set)]
    pub e_oh2: Option<[f64; 3]>,
    #[pyo3(get, set)]
    pub e_hh:  Option<[f64; 3]>,
    #[pyo3(get, set)]
    pub e_hoh: Option<[f64; 3]>,

    // --- r–z projection (if needed) ---
    #[pyo3(get, set)]
    pub rz_oh1: Option<[f64; 3]>, // [p_r, psi, p_z]
    #[pyo3(get, set)]
    pub rz_hh:  Option<[f64; 3]>,
    #[pyo3(get, set)]
    pub rz_hoh: Option<[f64; 3]>,

    // --- Classification / labels ---
    #[pyo3(get, set)]
    pub mol_type: Option<String>,  // e.g. "W-1a-L-S1"
    #[pyo3(get, set)]
    pub mol_kind: Option<MolKind>, // Enum: W, S, etc.

    // --- Hydrogen bond partners (only for HB analysis) ---
    #[pyo3(get, set)]
    pub h_bond_partners: Option<Vec<HydrogenBondPartner>>,
}

#[pymethods]
impl MoleculeRecord {
    #[new]
    pub fn new(o_idx: i32) -> Self {
        Self { o_idx, ..Default::default() }
    }

    pub fn __repr__(&self) -> PyResult<String> {
        Ok(serde_json::to_string_pretty(self).unwrap_or_else(|_| "MoleculeRecord serialization failed".into()))
    }
}

// Builder methods for setting only the necessary fields (chainable)
impl MoleculeRecord {
    pub fn with_h_indices(mut self, h1: Option<i32>, h2: Option<i32>) -> Self { 
        self.h1_idx = h1; 
        self.h2_idx = h2; 
        self 
    }
    pub fn with_positions(mut self, o: Option<[f64;3]>, h1: Option<[f64;3]>, h2: Option<[f64;3]>) -> Self { 
        self.o_pos = o; 
        self.h1_pos = h1; 
        self.h2_pos = h2; 
        self 
    }
    pub fn with_velocities(mut self, o: Option<[f64;3]>, h1: Option<[f64;3]>, h2: Option<[f64;3]>) -> Self { 
        self.o_vel = o; 
        self.h1_vel = h1; 
        self.h2_vel = h2; 
        self 
    }
    pub fn with_orient(mut self, oh1: Option<[f64;3]>, oh2: Option<[f64;3]>, hh: Option<[f64;3]>, hoh: Option<[f64;3]>) -> Self {
        self.e_oh1 = oh1; 
        self.e_oh2 = oh2; 
        self.e_hh = hh; 
        self.e_hoh = hoh; 
        self
    }
    pub fn with_rz(mut self, oh1: Option<[f64;3]>, hh: Option<[f64;3]>, hoh: Option<[f64;3]>) -> Self {
        self.rz_oh1 = oh1; 
        self.rz_hh = hh; 
        self.rz_hoh = hoh; 
        self
    }
    pub fn with_labels(mut self, mol_type: Option<String>, mol_kind: Option<MolKind>) -> Self {
        self.mol_type = mol_type; 
        self.mol_kind = mol_kind; 
        self
    }
    pub fn with_hbonds(mut self, partners: Option<Vec<HydrogenBondPartner>>) -> Self {
        self.h_bond_partners = partners; 
        self
    }
}
