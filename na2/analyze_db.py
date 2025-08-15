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





db_file = 'Na8_q0.db'  # use the oxygen DB created by StartGenerator


# === Connect to database and retrieve info ===
da = DataConnection(db_file)


# === Comparator: how to decide if two candidates are similar ===
comparator = InteratomicDistanceComparator(
    n_top=8,
    pair_cor_cum_diff=0.015,
    pair_cor_max=0.7,
    dE=0.02,
    mic=False
)


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

