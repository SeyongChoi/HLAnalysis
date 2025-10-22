// use bincode;
use pyo3::prelude::*;
use serde::{Deserialize, Serialize};

#[pyclass]
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct Atoms {
    // Index
    #[pyo3(get, set)]
    pub indices: Option<Vec<i32>>,
    // Chemical formula, a list of chemical symbols
    #[pyo3(get, set)]
    pub symbols: Option<Vec<String>>,
    // Atomic positions in Cartesian coordinates
    #[pyo3(get, set)]
    pub positions: Option<Vec<[f64; 3]>>,
    // Atomic unwrapped positions in Cartesian coordinates
    #[pyo3(get, set)]
    pub positions_unwrap: Option<Vec<[f64; 3]>>,
    // Forces for all atoms in Cartesian coordinates
    #[pyo3(get, set)]
    pub forces: Option<Vec<[f64; 3]>>,
    // Velocities for all atoms in Cartesian coordinates
    #[pyo3(get, set)]
    pub velocities: Option<Vec<[f64; 3]>>,
    // Atomic masses in atomic units
    #[pyo3(get, set)]
    pub masses: Option<Vec<f64>>,
    // Initial atomic charges
    #[pyo3(get, set)]
    pub charges: Option<Vec<f64>>,
    // Unit cell vectors.
    pub cell: Option<[[f64; 3]; 3]>,
    // Periodic boundary conditions flags
    #[pyo3(get, set)]
    pub pbc: Option<[bool; 3]>,
}

#[pymethods]
impl Atoms {
    #[new]
    pub fn new(
        indices: Option<Vec<i32>>,
        symbols: Option<Vec<String>>,
        masses: Option<Vec<f64>>,
        charges: Option<Vec<f64>>,
        positions: Option<Vec<[f64; 3]>>,
        positions_unwrap: Option<Vec<[f64; 3]>>,
        forces: Option<Vec<[f64; 3]>>,
        velocities: Option<Vec<[f64; 3]>>,
        cell: Option<[[f64; 3]; 3]>,
        pbc: Option<[bool; 3]>,
    ) -> Self {
        Atoms {
            indices,
            symbols,
            masses,
            charges,
            positions,
            positions_unwrap,
            forces,
            velocities,
            cell,
            pbc,
        }
    }

    // Method to get a reference to indices (returning a clone of the Option<Vec>)
    #[getter]
    pub fn indices(&self) -> Option<Vec<i32>> {
        self.indices.clone()
    }

    // Method to get a reference to masses (returning a clone of the Option<Vec>)
    #[getter]
    pub fn masses(&self) -> Option<Vec<f64>> {
        self.masses.clone()
    }

    // Method to get a reference to charges (returning a clone of the Option<Vec>)
    #[getter]
    pub fn charges(&self) -> Option<Vec<f64>> {
        self.charges.clone()
    }

    // Method to get a reference to symbols (returning a clone of the Option<Vec>)
    #[getter]
    pub fn symbols(&self) -> Option<Vec<String>> {
        self.symbols.clone()
    }

    // Method to get a reference to positions (returning a clone of the Option<Vec>)
    #[getter]
    pub fn positions(&self) -> Option<Vec<[f64; 3]>> {
        self.positions.clone()
    }

    // Method to get a reference to positions (returning a clone of the Option<Vec>)
    #[getter]
    pub fn positions_unwrap(&self) -> Option<Vec<[f64; 3]>> {
        self.positions_unwrap.clone()
    }

    // Method to get a reference to velocities (returning a clone of the Option<Vec>)
    #[getter]
    pub fn velocities(&self) -> Option<Vec<[f64; 3]>> {
        self.velocities.clone()
    }

    // Method to get a reference to forces (returning a clone of the Option<Vec>)
    #[getter]
    pub fn forces(&self) -> Option<Vec<[f64; 3]>> {
        self.forces.clone()
    }

    // Method to get a referenc to cell (retruning a clone of the Option<Vec<Vec<f64>>>)
    #[getter]
    pub fn cell(&self) -> Option<Vec<Vec<f64>>> {
        // Convert `Option<[[f64; 3]; 3]>` into `Option<Vec<Vec<f64>>>` for Python compatibility
        self.cell
            .as_ref()
            .map(|cell_array| cell_array.iter().map(|&row| row.to_vec()).collect())
    }

    // Method to get a referenc to cell (retruning a clone of the Option<[[f64; 3]; 3]>)
    #[getter]
    pub fn pbc(&self) -> Option<[bool; 3]> {
        self.pbc.clone()
    }

    // Method to get positions (Option<Vec<[f64; 3]>>)
    #[getter]
    pub fn get_positions(&self) -> Option<Vec<[f64; 3]>> {
        self.positions.clone()
    }

    // Method to get symbols (Option<Vec<String>>)
    #[getter]
    pub fn get_symbols(&self) -> Option<Vec<String>> {
        self.symbols.clone()
    }
}
