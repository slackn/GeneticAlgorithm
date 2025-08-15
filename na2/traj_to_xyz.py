from ase.io import read, write

# Load all structures from trajectory
structures = read('distinct_candidates_sorted.traj', index=':')

# Save all frames to a multi-frame XYZ file
write('distinct_candidates.xyz', structures)
