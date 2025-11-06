"""HLAnalysis: Hydration Layer Analysis Toolkit (Rust + Python)

High-level Python entrypoints and lazy access to Rust-registered submodules:
- hlanalysis.io        : I/O helpers (readers, serializers)
- hlanalysis.utils     : unit converters, constants, atomic data, MIC tools
- hlanalysis.analysis  : SFG, Hbond, Jump, ReorientDynamics, DensityProfile, BondOrientDist
"""

from . import _hlanalysis as _ext 
from ._version import __version__

from .io_wrapper import reader, read_atoms_from_bin

Atom = _ext.Atom
Atoms = _ext.Atoms

__all__ = [
    'Atom', 'Atoms',
    'reader', 'read_atoms_from_bin',
    '__version__'
    ]
