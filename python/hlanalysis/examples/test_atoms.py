from hlanalysis import Atom, Atoms

def main():
    '''
    Create "Atom" instances for each atom in a water molecule
    &
    Create an "Atoms" instance representing the full water molecule.
    '''
    # 'Atom' instances for each atom in a water molecule
    o_atom = Atom(symbol="O", position=[0.0, 0.0, 0.0], tag=1,
                  momentum=None, mass=15.999, magmom=None, charge=-0.834)

    h1_atom = Atom(symbol="H", position=[0.9572, 0.0, 0.0], tag=2,
                   momentum=None, mass=1.008, magmom=None, charge=0.417)

    h2_atom = Atom(symbol="H", position=[-0.2399872, 0.927297, 0.0], tag=3,
                   momentum=None, mass=1.008, magmom=None, charge=0.417)
    
    print("Created individual Atom instances:")
    for atom in [o_atom, h1_atom, h2_atom]:
        print(f"  {atom.symbol} at {atom.position} (charge={atom.charge})")

    # 'Atoms' instance for a water molecule
    symbols, positions, masses, charges = zip(*[(a.symbol, a.position, a.mass, a.charge)
                                              for a in (o_atom, h1_atom, h2_atom)])
    
    atoms = Atoms(
        indices=None,
        symbols=list(symbols),
        masses=list(masses),
        charges=list(charges),
        positions=list(positions),
        positions_unwrap=None,
        forces=None,
        velocities=None,
        cell=None,
        pbc=None,
    )

    print("Created Atoms instance:")
    print("Symbols:", atoms.symbols)
    print("Positions:", atoms.positions)


if __name__ == "__main__":
    main()
