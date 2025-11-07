mod atom;
mod atoms;
mod io;
mod utils;
mod analysis;

use pyo3::prelude::*;

const PKG: &str = "hlanalysis"; // sys.modules 등록에 사용할 패키지 접두사

/*------------------------------------------------------------------------------------ */
#[pymodule]
fn _hlanalysis(py: Python, root: &PyModule) -> PyResult<()> {
    /* ------------------------------------------------------------------------ */
    // top-level class
    root.add_class::<atom::Atom>()?;
    root.add_class::<atoms::Atoms>()?;

    // sub modulue attach
    attach(py, root, "io", io::register)?;
    attach(py, root, "utils", utils::register)?;
    attach(py, root, "analysis", analysis::register)?;
    Ok(())
}

// 공통 attach 유틸
fn attach<F>(py: Python<'_>, parent: &PyModule, name: &str, reg: F) -> PyResult<()>
where
    F: Fn(Python<'_>, &PyModule) -> PyResult<()>,
{
    let sub = PyModule::new(py, name)?;
    reg(py, &sub)?;                   // 서브모듈 내부에서 함수/하위모듈 등록
    parent.add_submodule(&sub)?;
    // sys.modules 에 hlanalysis.<name> 경로로 심기
    py.import("sys")?
        .getattr("modules")?
        .set_item(format!("{PKG}.{name}"), &sub)?;
    Ok(())
}