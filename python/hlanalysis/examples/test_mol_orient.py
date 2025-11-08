import os
import copy
import time
import numpy as np
import hlanalysis as hla
import matplotlib.pyplot as plt

# utils namespace: physchem (contains converter/atomic_data submodules), histogram
from hlanalysis.utils import histogram
from hlanalysis.analysis import orient_dist


def main():
    # -----------------------------
    # User parameters
    # -----------------------------
    dz = 0.1
    d_angle = 0.02
    start_ts, end_ts = 10000, 100001
    ts_interval = 1000
    mode='cos'   # or 'angle'
    z_range=[[12.9, 14.4], [69.6, 71.1]]
    # Index selection example (silanol + water ranges)
    index_silanol = []  # fill if needed
    index_water = [(42, 173), (216, 347), (390, 521), (564, 695), (738, 869),
                   (912, 1043), (1086, 1217), (1260, 1391), (1434, 1565)]
    indices_water = [i for start, end in index_water for i in range(start, end + 1)] if len(index_water) > 0  else []
    surface_normal=[0.0, 0.0, 1.0]
    # -----------------------------
    # Absolute paths
    # -----------------------------
    base_dir = os.path.abspath(os.path.dirname(__file__))             # where examples dir
    output_dir = os.path.join(base_dir, "example_input", "tmp_atoms") # where timestep_*.bin exist
    ori_result_dir = os.path.join(base_dir, "example_ori_result")     # where to return orientation distribution results
    bin_file = os.path.join(output_dir, "timestep_1.bin")

    if not os.path.isdir(output_dir):
        raise FileNotFoundError(f"Output directory not found: {output_dir}")
    if not os.path.isdir(ori_result_dir):
        os.makedirs(ori_result_dir, exist_ok=True)
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
    z_center = cell[2][2] / 2


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
    # Identify Atomic Groups (Water / Silanol)
    # -----------------------------
    oxygen_indices = []
    hydrogen_indices = []
    water_oxygen_indices = []
    water_hydrogen_indices = []
    silanol_oxygen_indices = []
    silanol_hydrogen_indices = []

    for i, symbol in enumerate(atoms.symbols):
        if symbol == 'O':
            oxygen_indices.append(i)
            if i in indices_water:
                water_oxygen_indices.append(i)
            elif i in index_silanol:
                silanol_oxygen_indices.append(i)

        elif symbol == 'H':
            hydrogen_indices.append(i)
            if i in indices_water:
                water_hydrogen_indices.append(i)
            elif i in index_silanol:
                silanol_hydrogen_indices.append(i)

    # Mapping of O–H pairs
    silanol_group_map = {o_idx:h_idx for o_idx, h_idx in zip(silanol_oxygen_indices, silanol_hydrogen_indices)}
    water_mol_map = {idx: (ox, h1x, h2x) for idx, (ox, h1x, h2x) in enumerate((ox_i, ox_i+1, ox_i+2) for ox_i in water_oxygen_indices)}
    
    oh_pair_silanol_map = {idx:(o_idx, h_idx) for idx, (o_idx, h_idx) in enumerate((ox_i, silanol_group_map[ox_i]) for ox_i in silanol_oxygen_indices)}
    oh_pair_water_map = {idx:(o_idx, h_idx) for idx, (o_idx, h_idx) in enumerate((ox_i, ox_i + k + 1) for ox_i in water_oxygen_indices for k in range(2))}
    
    # Merge maps into a single OH-pair map
    silanol_offset = len(oh_pair_silanol_map)
    oh_pair_water_map = {idx + silanol_offset:val for idx,val in oh_pair_water_map.items()}
    oh_pair_map = {**oh_pair_silanol_map, **oh_pair_water_map}

    
    # -----------------------------
    # Print diagnostic summary
    # -----------------------------
    print("Base directory:", base_dir)
    print("Dump/bin directory:", output_dir)
    print(f"\nFrames in range [{start_ts}, {end_ts}]: {num_frames}")
    print(f"Cell (Å):\n{cell}")

    # -----------------------------
    # Orientation Distribution from Rust Backend
    # -----------------------------
    angles = orient_dist.mol_orient.cal_mol_orient(water_map=water_mol_map,
                                                  surface_normal=surface_normal,
                                                  parallel=True,
                                                  dir=output_dir,
                                                  start_ts=start_ts,
                                                  end_ts=end_ts,
                                                  interval_ts=ts_interval,
                                                  center=z_center,
                                                  mode=mode)
    
    # Generate bins for z and angle
    z_edges, z_centers = histogram.generate_bins(0, cell[2][2] + dz, dz)
    if mode == "angle":
        angle_edges_hh, angle_centers_hh = histogram.generate_bins(90, 180 + d_angle, d_angle)
        angle_edges_dipole, angle_centers_dipole = histogram.generate_bins(0, 180 + d_angle, d_angle)
    elif mode == "cos":
        angle_edges_hh, angle_centers_hh = histogram.generate_bins(-1.0, 0.0 + d_angle, d_angle)
        angle_edges_dipole, angle_centers_dipole = histogram.generate_bins(-1.0, 1.0 + d_angle, d_angle)
    
    # Bin orientation data along z
    angle_pair_bins_along_z = orient_dist.mol_orient.angle_pair_bins_for_normal_range(
                                                     all_angles=angles,
                                                     bins_dipole=angle_edges_dipole,
                                                     bins_hh=angle_edges_hh,
                                                     z_range=z_range,
                                                     parallel=True,
                                                     dir=output_dir,
                                                     start_ts=start_ts,
                                                     end_ts=end_ts,
                                                     interval_ts=ts_interval)
    
    avg_angle_pair_bins_for_normal_range = np.sum(np.array(angle_pair_bins_along_z), axis=0, dtype=np.float64)

    # ---------------------------------------------------------------------------------------------------------------- #
    # Save Orientation Distribution
    # ---------------------------------------------------------------------------------------------------------------- #
    output_data = np.zeros((len(angle_centers_hh), len(angle_centers_dipole) + 1))  
    output_data[:, 0] = angle_centers_hh  
    output_data[:, 1:] = avg_angle_pair_bins_for_normal_range  
    np.savetxt(
        os.path.join(ori_result_dir, "mol_orient_dist.txt"),
        output_data,
        fmt="%.6f",
        header="angle_center_hh " + " ".join(f"angle_center_dipole_{i}" for i in angle_centers_dipole)
    )
    
    # ---------------------------------------------------------------------------------------------------------------- #
    # Visualization
    # ---------------------------------------------------------------------------------------------------------------- #
    avg_angle_pair_bins_for_normal_range /= num_frames
    avg_den_angle_pair_bins_for_normal_range = avg_angle_pair_bins_for_normal_range / np.sum(avg_angle_pair_bins_for_normal_range)

    def plot_angle_distribution(distribution, angle_centers, z_centers, save=False):
        """Visualize 2D orientation probability distribution."""
        angle_grid, z_grid = np.meshgrid(angle_centers, z_centers)

        plt.figure(figsize=(12, 8), dpi=300)
        plt.contourf(angle_grid, z_grid, distribution, levels=25, cmap=plt.cm.turbo)
        cbar = plt.colorbar()
        cbar.ax.tick_params(labelsize=18)

        plt.xlabel(r'$\psi_{D}$ (°)', fontsize=20)
        plt.ylabel('z (Å)', fontsize=20)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.tight_layout()

        if save:
            plt.savefig(os.path.join(ori_result_dir, "dipole_ori_dist.png"))
        else:
            plt.show()

    plt.figure(figsize=(12, 6), dpi=300)
    plt.pcolormesh(angle_edges_dipole, angle_edges_hh, avg_den_angle_pair_bins_for_normal_range, cmap='coolwarm', shading='auto')
    plt.gca().invert_yaxis()
    plt.colorbar(label='Counts or Average Value')
    plt.xlabel('Dipole angle (cos θ)')
    plt.ylabel('HH angle (cos θ)')
    plt.tight_layout()
    plt.savefig(os.path.join(ori_result_dir, "mol_ori_dist.png"))
    # plt.show()

if __name__ == "__main__":
    main()
