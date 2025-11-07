import hlanalysis as hla
import os


def main():
    """
    Test the Rust-backed I/O interface from HLAnalysis:
    - hlanalysis.io.reader(): Read and parse LAMMPS dump files
    - hlanalysis.io.read_atoms_from_bin(): Load serialized timestep data
    """

    # Get absolute path of this script
    base_dir = os.path.abspath(os.path.dirname(__file__))

    # Construct absolute paths for input and output
    dump_file = os.path.join(base_dir, "example_input", "test.lammpstrj")
    output_dir = os.path.join(base_dir, "example_input", "tmp_atoms")

    print(f"Base directory: {base_dir}")
    print(f"Dump file path: {dump_file}")
    print(f"Output directory: {output_dir}")

    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    # Check that required submodules exist
    print("\n[Module Checks]")
    print(f"Has io module: {hasattr(hla, 'io')}")
    print(f"Has reader: {hasattr(hla.io, 'reader')}")
    print(f"Has read_atoms_from_bin: {hasattr(hla.io, 'read_atoms_from_bin')}")

    # Run the LAMMPS dump reader if file exists
    if os.path.exists(dump_file):
        print("\n[Running reader()]")
        try:
            hla.io.reader(
                file_path=dump_file,
                format="lammps-dump",
                start=0,
                end=1000,           # test range (reduce for speed)
                interval=1,
                output_dir=output_dir,
            )
            print("Successfully executed reader()")
        except Exception as e:
            print("Error in reader():", e)
    else:
        print(f"Dump file not found: {dump_file}")

    # Test binary loading
    bin_file = os.path.join(output_dir, "timestep_0.bin")
    if os.path.exists(bin_file):
        print("\n[Running read_atoms_from_bin()]")
        try:
            atoms = hla.io.read_atoms_from_bin(bin_file)
            n_atoms = len(atoms.positions or [])
            print(f"Loaded Atoms object with {n_atoms} atoms")
        except Exception as e:
            print("Error loading .bin file:", e)
    else:
        print(f"Binary file not found: {bin_file}")


if __name__ == "__main__":
    main()
