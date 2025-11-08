
use pyo3::{
    FromPyObject, PyAny, PyResult, 
    exceptions::PyTypeError,
    types::PyList, IntoPy, PyObject,Python  // PyList 추가
};
use serde::{Serialize, Deserialize};



// 1. 3차원 배열을 위한 래퍼 타입 정의
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct PyCell3(pub [[f64; 3]; 3]);

// 2. FromPyObject 구현 (3x3 행렬 처리)
impl<'source> FromPyObject<'source> for PyCell3 {
    fn extract(ob: &'source PyAny) -> PyResult<Self> {
        let py_list = ob.downcast::<PyList>()?;
        if py_list.len() != 3 {
            return Err(PyTypeError::new_err("3x3 행렬이 필요합니다"));
        }

        let mut cell = [[0.0; 3]; 3];
        for (i, item) in py_list.iter().enumerate() {
            let inner_list = item.downcast::<PyList>()?;
            if inner_list.len() != 3 {
                return Err(PyTypeError::new_err("각 행은 3개의 요소를 가져야 합니다"));
            }
            for (j, val) in inner_list.iter().enumerate() {
                cell[i][j] = val.extract()?;
            }
        }
        Ok(PyCell3(cell))
    }
}

impl IntoPy<PyObject> for PyCell3 {
    fn into_py(self, py: Python) -> PyObject {
        // 각 행을 PyList로 변환한 후, 전체를 다시 PyList로 만듭니다.
        let rows: Vec<PyObject> = self.0.iter()
            .map(|row| PyList::new(py, row).into_py(py))
            .collect();
        PyList::new(py, rows).into_py(py)
    }
}