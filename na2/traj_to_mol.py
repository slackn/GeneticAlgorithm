from ase.io import read, write

atoms = read('distinct_candidates_sorted.traj', index='0')
write('na8_single.cml', atoms)
