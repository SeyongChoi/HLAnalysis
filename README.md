# HLAnalysis

**Hydration Layer Analysis Toolkit**  
A hybrid **Rust + Python** package for analyzing structural and dynamical properties of interfacial water at solid/water interfaces.

---

## Overview
HLAnalysis provides a high-performance backend written in **Rust** (with [PyO3](https://github.com/PyO3/pyo3) bindings)  
and a user-friendly Python interface.  
It is designed for **hydration layer–resolved analysis** of atomistic simulations (e.g., LAMMPS, CP2K, AIMD), including:

- Density profiles along surface normal  
- OH and molecular orientation distributions  
- Reorientation and cross-correlation dynamics  
- Hydrogen-bond statistics and residence times  
- Vibrational spectroscopy (IR, SFG) based on time-correlation functions  
---
## Development Status
**Status:** *In progress*  
- [x] Atom / Atoms Rust module  
- [x] I/O Rust module  
- [ ] Density Profile Rust module  
- [ ] Orientation Distribution Rust module  
- [ ] Reorientation Dynamics Rust module  
- [ ] Hbond Analysis Rust module  
- [ ] Spectroscopy Rust module  

---
## Requirements

HLAnalysis uses **PyO3** for Rust–Python binding and **maturin** as a build backend.  
Before installation, please install **Rust** from the official website:  
👉 [https://www.rust-lang.org/tools/install](https://www.rust-lang.org/tools/install)

## Usage
### 1. Installation (dev)
```bash
# Create & activate environment (Python ≥ 3.8)
conda create -n py38-hla python==3.8.*
conda activate py38-hla

# Install maturin
pip install maturin 
# Clone repo
git clone https://github.com/SeyongChoi/HLAnalysis.git
cd HLAnalysis

# Build & install
maturin develop -r

# (Optional) python dependencies
pip install -r requirement.txt
```

### 2. Verify Installation
```bash
python -c "import hlanalysis; print('HLAnalysis successfully loaded!')"
```