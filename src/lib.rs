mod atom;
mod atoms;
mod io;
mod utils;
mod analysis;

use pyo3::prelude::*;

const PKG: &str = "hlanalysis"; // Package prefix for sys.modules registration

/*------------------------------------------------------------------------------------ */
#[pymodule]
fn _hlanalysis(py: Python, root: &PyModule) -> PyResult<()> {
    /* ------------------------------------------------------------------------ */
    // Register top-level classes
    root.add_class::<atom::Atom>()?;
    root.add_class::<atoms::Atoms>()?;

    // Attach submodules
    attach(py, root, "io", io::register)?;
    attach(py, root, "utils", utils::register)?;
    attach(py, root, "analysis", analysis::register)?;
    Ok(())
}

/*------------------------------------------------------------------------------------ */
// Common utility to attach submodules dynamically.
//
// This function:
// 1. Creates a new PyModule instance for the given submodule name.
// 2. Calls its corresponding `register()` function to populate Python functions/classes.
// 3. Adds it as a submodule to the parent PyModule.
// 4. Inserts it into `sys.modules` as `hlanalysis.<submodule>`,
//    so it can be imported directly from Python.
fn attach<F>(py: Python<'_>, parent: &PyModule, name: &str, reg: F) -> PyResult<()>
where
    F: Fn(Python<'_>, &PyModule) -> PyResult<()>,
{
    let sub = PyModule::new(py, name)?;
    reg(py, &sub)?; // Register all functions/classes within this submodule
    parent.add_submodule(&sub)?;
    // Add the submodule to sys.modules under the path "hlanalysis.<name>"
    py.import("sys")?
        .getattr("modules")?
        .set_item(format!("{PKG}.{name}"), &sub)?;
    Ok(())
}
