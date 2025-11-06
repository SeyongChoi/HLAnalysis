from . import io
from .atoms import Atoms
from typing import Optional

def reader(
    file_path: str,
    format: str = "lammps-dump",
    start: Optional[int] = None,
    end: Optional[int] = None,
    interval: Optional[int] = None,
    output_dir: Optional[str] = None,
) -> None:
    """
    Wrapper for hlanalysis.io.reader(...)

    Parameters
    ----------
    file_path : str
        Path to the trajectory / dump file.
    format : str, optional
        Format selector. Default: "lammps-dump".
    start, end, interval : Optional[int]
        Timestep block range and sampling interval.
    output_dir : Optional[str]
        Where to store per-timestep binary outputs (if the reader writes them).

    Raises
    ------
    ValueError, IOError propagated from native layer.
    """
    return io.reader(file_path, format, start, end, interval, output_dir)

def read_atoms_from_bin(path: str) -> Atoms:
    """
    Wrapper for hlanalysis.io.read_atoms_from_bin(...)

    Parameters
    ----------
    path : str
        Path to the binary file generated from LAMMPS dump reading.
    """
    return io.read_atoms_from_bin(path)