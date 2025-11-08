use pyo3::prelude::*;
use serde::{Serialize, Deserialize};

#[pyclass]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum MolKind { W, S, Other(String) }

#[pyclass]
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum HBondKind { Donor, Acceptor }
