use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;

#[pyclass]
#[derive(Debug, Clone, Default)]
pub struct Atom {
    #[pyo3(get, set)]
    pub symbol: String,
    #[pyo3(get, set)]
    pub position: [f64; 3],
    #[pyo3(get, set)]
    pub tag: Option<i32>,
    #[pyo3(get, set)]
    pub momentum: Option<[f64; 3]>,
    #[pyo3(get, set)]
    pub mass: Option<f64>,
    #[pyo3(get, set)]
    pub magmom: Option<f64>,
    #[pyo3(get, set)]
    pub charge: Option<f64>,
}

#[pymethods]
impl Atom {
    /// Python-friendly constructor with keyword args & defaults.
    #[new]
    pub fn new(
        mut symbol: String,
        position: [f64; 3],
        tag: Option<i32>,
        momentum: Option<[f64; 3]>,
        mass: Option<f64>,
        magmom: Option<f64>,
        charge: Option<f64>,
    ) -> PyResult<Self> {
        // (선택) 간단 검증/정규화
        symbol.make_ascii_uppercase();
        if position.iter().any(|v| !v.is_finite()) {
            return Err(PyValueError::new_err("position must be finite numbers"));
        }
        if let Some(p) = &momentum {
            if p.iter().any(|v| !v.is_finite()) {
                return Err(PyValueError::new_err("momentum must be finite numbers"));
            }
        }
        Ok(Self { symbol, position, tag, momentum, mass, magmom, charge })
    }

    /// Nice __repr__ for debugging.
    fn __repr__(&self) -> String {
        format!(
            "Atom(symbol='{}', position=[{:.6}, {:.6}, {:.6}], tag={:?}, momentum={:?}, mass={:?}, magmom={:?}, charge={:?})",
            self.symbol, self.position[0], self.position[1], self.position[2],
            self.tag, self.momentum, self.mass, self.magmom, self.charge
        )
    }

    /// (Optional) Pickle support: return a tuple state
    fn __getstate__(&self) -> (String, [f64;3], Option<i32>, Option<[f64;3]>, Option<f64>, Option<f64>, Option<f64>) {
        (
            self.symbol.clone(),
            self.position,
            self.tag,
            self.momentum,
            self.mass,
            self.magmom,
            self.charge,
        )
    }

    /// (Optional) Unpickle from tuple
    fn __setstate__(&mut self, state: (String, [f64;3], Option<i32>, Option<[f64;3]>, Option<f64>, Option<f64>, Option<f64>)) {
        let (symbol, position, tag, momentum, mass, magmom, charge) = state;
        self.symbol = symbol;
        self.position = position;
        self.tag = tag;
        self.momentum = momentum;
        self.mass = mass;
        self.magmom = magmom;
        self.charge = charge;
    }

    // (선택) 커스텀 setter로 검증 넣고 싶으면 이렇게:
    // #[setter]
    // fn set_symbol(&mut self, mut s: String) {
    //     s.make_ascii_uppercase();
    //     self.symbol = s;
    // }
}
