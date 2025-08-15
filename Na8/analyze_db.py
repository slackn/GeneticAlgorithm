from ase.db import connect
from ase.ga.standard_comparators import InteratomicDistanceComparator

db_file = 'Na8_q0.db'
db = connect(db_file)

print(f"\n=== All structures with raw_score ===\n")
for row in db.select():
    confid = row.get('confid', 'N/A')
    raw_score = row.key_value_pairs.get('raw_score', 'N/A')
    print(f"confid={confid}, raw_score={raw_score}")

# --- Filter valid candidates (raw_score != 1) ---
candidates = [row.toatoms() for row in db.select() 
              if row.key_value_pairs.get('raw_score', 1) != 1]

print(f"\nTotal valid candidates (raw_score != 1): {len(candidates)}")

# --- Deduplicate ---
if len(candidates) > 0:
    comparator = InteratomicDistanceComparator(
        n_top=len(candidates[0]),
        pair_cor_cum_diff=0.015,
        pair_cor_max=0.7,
        dE=0.02,
        mic=False
    )

    distinct = []
    for atoms in candidates:
        if not any(comparator.looks_like(atoms, unique) for unique in distinct):
            distinct.append(atoms)

    print(f"Distinct candidates: {len(distinct)}")
else:
    print("No valid candidates found.")
