"""HLAnalysis: Hydration Layer Analysis Toolkit (Rust + Python)

High-level Python entrypoints and lazy access to Rust-registered submodules:
- hlanalysis.io        : I/O helpers (readers, serializers)
- hlanalysis.utils     : unit converters, constants, atomic data, MIC tools
- hlanalysis.analysis  : SFG, Hbond, Jump, ReorientDynamics, DensityProfile, BondOrientDist
"""

from . import _hlanalysis as _ext 
from ._version import __version__

Atom = _ext.Atom
Atoms = _ext.Atoms

io = _ext.io
utils = _ext.utils
analysis = _ext.analysis

from .io_wrapper import reader, read_atoms_from_bin


__all__ = [
    'Atom', 'Atoms',
    'io', 'utils', 'analysis',
    'reader', 'read_atoms_from_bin',
    '__version__'
    ]
