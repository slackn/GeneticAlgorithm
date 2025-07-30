from ase import Atoms
from ase.db import connect
from ase.ga.startgenerator import StartGenerator
from ase.ga.utilities import closest_distances_generator
from ase.ga.data import PrepareDB
import numpy as np

# Parameters
n_atoms = 8
charge = 0
n_to_generate = 20
db_name = f'Na{n_atoms}_q{charge}.db'

# Compute bounding cube side length
side = 2 * (0.5 + ((3 * n_atoms) / (4 * np.pi * np.sqrt(2))) ** (1/3))

# 1. Define a fake slab (empty unit cell with no atoms)
# This is the simulation box
buffer=5    #extra vacuum space
slab_side= side+buffer
slab = Atoms(cell=[slab_side, slab_side, slab_side], pbc=False)

# 2. Define what you're placing: Ge atoms
blocks = [('Na', n_atoms)]

# 3. Set minimum interatomic distances
# Sodium atomic number is 11
blmin = closest_distances_generator(atom_numbers=[11], ratio_of_covalent_radii=0.7)


# Define the box in which to place atoms (inside the slab)
#This s the placement box 
origin = [(slab_side - side) / 2] * 3  # Center the box inside the slab
box = [origin,
       [[side, 0, 0],
        [0, side, 0],
        [0, 0, side]]]

# 5. Prepare the database (overwrite if exists)
d = PrepareDB(db_name, simulation_cell=slab, stoichiometry=[11] * n_atoms)

# 6. Create the StartGenerator
sg = StartGenerator(
    slab=slab,
    blocks=blocks,
    blmin=blmin,
    box_to_place_in=box,
)

# 7. Generate and write the clusters
for _ in range(n_to_generate):
    atoms = sg.get_new_candidate()
    atoms.charge = charge
    atoms.set_initial_charges([charge / n_atoms] * n_atoms)
    d.add_unrelaxed_candidate(atoms)
