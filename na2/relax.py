from ase.db import connect
db = connect('Na8_q0.db')

unrelaxed = [row for row in db.select(relaxed=False)]
relaxed = [row for row in db.select(relaxed=True)]

print(f"Unrelaxed: {len(unrelaxed)}")
print(f"Relaxed: {len(relaxed)}")
print(f"Total rows: {len(list(db.select()))}")
