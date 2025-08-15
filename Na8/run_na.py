from random import random

from ase.calculators.emt import EMT  # Replace with GPAW or Turbomole later
from ase.ga.cutandsplicepairing import CutAndSplicePairing
from ase.ga.data import DataConnection
from ase.ga.offspring_creator import OperationSelector
from ase.ga.population import Population
from ase.ga.standard_comparators import InteratomicDistanceComparator
from ase.ga.standardmutations import MirrorMutation, RattleMutation, PermutationMutation, StrainMutation
from ase.ga.utilities import closest_distances_generator, get_all_atom_types
from ase.io import write
from ase.optimize import BFGS, LBFGS
from ase.calculators.turbomole import Turbomole
from ase.calculators.singlepoint import SinglePointCalculator
import numpy as np
import shutil
from pathlib import Path
import os


# === Parameters ===
population_size = 20
mutation_probability = 0.5
n_to_test = 30
db_file = 'Na8_q0.db'  # use the sodium DB created by StartGenerator

# === Connect to database and retrieve info ===
da = DataConnection(db_file)
atom_numbers = [11] * 8  
n_atoms = len(atom_numbers)

# === Define chemical and geometric constraints ===
slab = da.get_slab()
all_atom_types = get_all_atom_types(slab, atom_numbers)
blmin = closest_distances_generator(all_atom_types, ratio_of_covalent_radii=0.7)

# === Comparator: how to decide if two candidates are similar ===
comparator = InteratomicDistanceComparator(
    n_top=n_atoms,
    pair_cor_cum_diff=0.015,
    pair_cor_max=0.7,
    dE=0.02,
    mic=False
)

# === Operators ===
pairing = CutAndSplicePairing(slab=slab, n_top=n_atoms, blmin=blmin)
mutations = OperationSelector(
    [1.0, 1.0 ],
    [
        MirrorMutation(blmin, n_atoms),
        RattleMutation(blmin, n_atoms)
        #StrainMutation(blmin, cellbounds)

        #PermutationMutation(n_atoms)   #permutation with one atom type is not valid
    ]
)

# === Turbomole parameters ===
params = {
    'total charge': 0,
    'multiplicity': 1,
    'scf iterations': 1000,
    #'basis set name': 'def2-TZVP',
    #'density functional': 'pbe'
    #'basis set name':'dhf-TZVP',
    #'density functional':'tpss'
}
#to mimic the paper I use their parameters

def get_dft_calculator():
    #calc_dir= "calc_dir"
    #restart_flag=os.path.exists(calc_dir)
    return Turbomole( **params)



# === Relax all unrelaxed initial structures ===

while da.get_number_of_unrelaxed_candidates() > 0:
    atoms = da.get_an_unrelaxed_candidate()
    confid = atoms.info.get('confid', 'N/A')

    
    atoms.calc = get_dft_calculator()
    dyn = BFGS(atoms, logfile=None)

    try:
        print(f"[GA] Relaxing offspring confid={confid}")
        dyn.run(fmax=0.05, steps=500)
        print(f"[GA] Offspring {confid} converged in {dyn.nsteps} steps.")
    except RuntimeError as e:
        print(f"[Relax] Skipping structure due to SCF failure: {e}")
        print(f"[GA] SCF failure for confid={confid}: {e}")
         # Mark as bad
        atoms.info['key_value_pairs']['raw_score'] = 1  # very low score
        atoms.calc = SinglePointCalculator(atoms, energy=1e6)
        da.add_relaxed_step(atoms)  # moves it out of the "unrelaxed" pool
        #input(f"\n[Pause] Candidate {confid} cannot be relaxed. Press Enter to continue...")  
        continue

    energy = atoms.get_potential_energy()
    print(f"[Relax] Energy from calculator = {energy}")
    atoms.info['key_value_pairs']['raw_score'] = -energy

    # === Pause here ===
    #input(f"\n[Pause] Candidate {confid} relaxed. Press Enter to continue...")  

    da.add_relaxed_step(atoms)

# Print relaxed candidates
for row in da.get_all_relaxed_candidates():
    print(f"confid={row.info['confid']}, raw_score={row.info['key_value_pairs']['raw_score']}")

# === Build population from relaxed structures ===
population = Population(
    data_connection=da,
    population_size=population_size,
    comparator=comparator
)

print(f"\n[Population] Current population size: {len(population.pop)}")

for i, candidate in enumerate(population.pop):
    confid = candidate.info.get('confid', 'N/A')
    energy = candidate.info['key_value_pairs'].get('raw_score', 'N/A')
    print(f"  {i+1}: confid={confid}, raw_score={energy}")

'''
p1, p2=population.get_two_candidates()
confid = p1.info.get('confid', 'N/A')
score = p2.info['key_value_pairs'].get('raw_score', 'N/A')
print(f"[Population] confid={confid}, raw_score={score}")
'''

# === Main GA loop ===
for i in range(n_to_test):
    print(f"\n[GA] Starting candidate {i + 1}/{n_to_test}")

    parent1, parent2 = population.get_two_candidates()
    # --- Print parent info ---
    confid1 = parent1.info.get('confid', 'N/A')
    confid2 = parent2.info.get('confid', 'N/A')
    print(f"[GA] Selected parents:")
    print(f"   Parent 1: confid={confid1}")
    print(f"   Parent 2: confid={confid2}")

    child, description = pairing.get_new_individual([parent1, parent2])
    confid = child.info.get('confid', 'N/A')

    if child is None:
        print("[GA] Pairing failed. Skipping.")
        continue
    da.add_unrelaxed_candidate(child, description)

    if random() < mutation_probability:
        mutated_child, desc = mutations.get_new_individual([child])
        if mutated_child is not None:
            da.add_unrelaxed_step(mutated_child, desc)
            child = mutated_child
            print("[GA] Mutation applied.")

    child.calc = get_dft_calculator()
    dyn = BFGS(child, logfile=None)
    try:
        print('Offspring relaxation starts.')
        dyn.run(fmax=0.005, steps=500)
        child.info['key_value_pairs']['raw_score'] = -child.get_potential_energy()
        da.add_relaxed_step(child)
    except RuntimeError as e:
        print(f"[GA] Skipping candidate due to SCF failure: {e}")
        print(f"[GA] SCF failure for confid={confid}: {e}")
         # Mark as bad
        child.info['key_value_pairs']['raw_score'] = 1  # very high energy
        child.calc = SinglePointCalculator(child, energy=1e6)
        da.add_relaxed_step(child)  # moves it out of the "unrelaxed" pool
        continue

    population.update()

# === Save results ===
#write('all_candidates.traj', da.get_all_relaxed_candidates())



# Get all relaxed candidates from the database
candidates = list(da.get_all_relaxed_candidates())
# Sort them by confid
candidates.sort(key=lambda atoms: atoms.info.get('confid', 0))
# Write to a trajectory file
write('all_candidates_sorted.traj', candidates)

print("\n GA run completed. Relaxed structures saved to all_candidates_sorted.traj")

# Print confid and energy for each candidate
for atoms in candidates:
    confid = atoms.info.get('confid', 'N/A')
    raw_score = atoms.info['key_value_pairs'].get('raw_score', None)
    energy = -raw_score if raw_score is not None else 'N/A'
    
    print(f"confid={confid}, energy={energy:.6f} eV")


# === Remove failed structures (raw_score == 1) ===
real_candidates = [
    atoms for atoms in candidates
    if atoms.info['key_value_pairs'].get('raw_score', 1) != 1
]

print(f"[Filter] {len(real_candidates)} candidates remaining after removing failed ones.")

# === Deduplicate using the comparator manually ===
distinct_candidates = []
for atoms in real_candidates:
    is_duplicate = False
    for unique_atoms in distinct_candidates:
        if comparator.looks_like(atoms, unique_atoms):
            # If they are similar, keep the lower energy one
            energy_atoms = -atoms.info['key_value_pairs']['raw_score']
            energy_unique = -unique_atoms.info['key_value_pairs']['raw_score']
            if energy_atoms < energy_unique:
                # Replace the higher-energy duplicate
                distinct_candidates[distinct_candidates.index(unique_atoms)] = atoms
            is_duplicate = True
            break

    if not is_duplicate:
        distinct_candidates.append(atoms)

# === Sort by energy (lowest first) ===
distinct_candidates.sort(
    key=lambda atoms: -atoms.info['key_value_pairs']['raw_score']  # raw_score = -energy
)

# === Save final candidates ===
write('distinct_candidates_sorted.traj', distinct_candidates)

print(f"\n[Done] Saved {len(distinct_candidates)} unique, valid candidates to distinct_candidates_sorted.traj")
for atoms in distinct_candidates:
    confid = atoms.info.get('confid', 'N/A')
    energy = -atoms.info['key_value_pairs']['raw_score']
    print(f"confid={confid}, energy={energy:.6f} eV")




'''
# === Select distinct structures ===
# Remove failed structures 
real_candidates=[
    atoms for atoms in candidates
    if atoms.info['key_value_pairs'].get('raw_score',1)!=1
]


final_population= Population(data_connection=None,
                            population_size=len(real_candidates),
                            comparator=comparator)


for atoms in real_candidates:
    final_population.add_candidate(atoms)

distinct_candidates=final_population.pop

#Save
write('distinct_candidates.traj', distinct_candidates)
'''

'''
# === Select distinct structures ===
# Remove failed structures (raw_score == 1)
real_candidates = [
    atoms for atoms in candidates
    if atoms.info['key_value_pairs'].get('raw_score', 1) != 1
]

print(f"[Filter] {len(real_candidates)} candidates remaining after removing failed ones.")

# === Create a dummy in-memory database for Population ===
dummy_db = DataConnection(':memory:')

# Add all real candidates to this temporary database
for atoms in real_candidates:
    dummy_db.add_relaxed_step(atoms)

# === Build a population that uses the comparator to filter duplicates ===
final_population = Population(
    data_connection=dummy_db,
    population_size=len(real_candidates),  # ensure enough space for all unique ones
    comparator=comparator
)

# Extract distinct candidates
distinct_candidates = final_population.pop

# === Sort by energy (lowest first) ===
distinct_candidates.sort(
    key=lambda atoms: -atoms.info['key_value_pairs']['raw_score']  # raw_score = -energy
)

# === Save the final distinct candidates ===
write('distinct_candidates_sorted.traj', distinct_candidates)

print(f"\n[Done] Saved {len(distinct_candidates)} unique, valid candidates to distinct_candidates_sorted.traj")
for atoms in distinct_candidates:
    confid = atoms.info.get('confid', 'N/A')
    energy = -atoms.info['key_value_pairs']['raw_score']
    print(f"confid={confid}, energy={energy:.6f} eV")

'''