import os
import copy
import time
import numpy as np
import hlanalysis as hla

# utils namespace: physchem (contains converter/atomic_data submodules), histogram
from hlanalysis.utils import physchem, histogram
# from hlanalysis.analysis import density_profile   # ← to be added later


def main():
    # -----------------------------
    # User parameters
    # -----------------------------
    dz = 0.1
    z_lower, z_upper = 0.0, 84.0
    start_ts, end_ts = 0, 999
    target_symbol = "O"

    # -----------------------------
    # Absolute paths
    # -----------------------------
    base_dir = os.path.abspath(os.path.dirname(__file__))
    output_dir = os.path.join(base_dir, "example_input", "tmp_atoms")
    bin_file = os.path.join(output_dir, "timestep_0.bin")

    if not os.path.isdir(output_dir):
        raise FileNotFoundError(f"Output directory not found: {output_dir}")
    if not os.path.isfile(bin_file):
        raise FileNotFoundError(f"First bin file not found: {bin_file}")

    # -----------------------------
    # Load one frame and extract cell information
    # Atoms.cell() -> list[list[float]] (getter method)
    # -----------------------------
    atoms = hla.io.read_atoms_from_bin(bin_file)

    # NOTE: Atoms.cell is a getter method → must be called with parentheses
    cell_mat = atoms.cell
    if cell_mat is None:
        raise RuntimeError("No cell information found (cell == None)")

    cell = np.array(cell_mat, dtype=float)  # 3×3 cell matrix
    unit_cell = copy.deepcopy(cell)
    unit_cell[2, 2] = dz  # set slab thickness

    # Convert Å³ → cm³
    ang_to_cm = physchem.converter.ang_to_cm()
    unit_cell_vol_cm3 = float(np.linalg.det(unit_cell)) * (ang_to_cm ** 3)

    # -----------------------------
    # Count number of frames
    # -----------------------------
    num_frames = sum(
        1
        for f in os.listdir(output_dir)
        if f.startswith("timestep_") and f.endswith(".bin")
        and f[9:-4].isdigit()
        and start_ts <= int(f[9:-4]) <= end_ts
    )

    # -----------------------------
    # Compute H₂O molecular mass (g/mol)
    # -----------------------------
    mass_water_g_per_mol = (
        physchem.atomic_masses("O") + 2.0 * physchem.atomic_masses("H")
    )

    # -----------------------------
    # Generate z-bin edges and centers
    # -----------------------------
    edges, centers = histogram.generate_bins(z_lower, z_upper + dz, dz)

    # -----------------------------
    # Print diagnostic summary
    # -----------------------------
    print("Base directory:", base_dir)
    print("Dump/bin directory:", output_dir)
    print(f"\nFrames in range [{start_ts}, {end_ts}]: {num_frames}")
    print(f"Cell (Å):\n{cell}")
    print(f"Slab volume with dz={dz:.3f} Å: {unit_cell_vol_cm3:.3e} cm³")
    print(f"Mass of H₂O: {mass_water_g_per_mol:.3f} g/mol")
    print(f"z-bins: {len(centers)} (from {edges[0]:.2f} to {edges[-1]:.2f} Å, dz={dz})")

    # ------------------------------------------------------------------
    # TODO (future work):
    # Integrate the density-profile computation
    # Example:
    # from hlanalysis.analysis import density_profile
    # z_centers, counts = density_profile.accumulate_counts(
    #     target_symbol, [], list(edges), True, output_dir, start_ts, end_ts
    # )
    # Then convert counts → number density → mass density (g/cm³)
    # ------------------------------------------------------------------


if __name__ == "__main__":
    main()
