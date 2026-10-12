"""Aim 3 planted-signal benchmark on RPE1-K562 under plan v1.3 (500 reps/level, B = 999).

Decoy similarities, k and every test use the v1.3 f(C) (arm effects). Writes
results/aim3/planted_RPE1_vs_K562_v13{_replicates.parquet,_info.json,.csv}.
Run: PYTHONPATH=src python scripts/run_planted_dev_v13.py
"""
import time

from cdtest import planted

t0 = time.time()
df = planted.run("RPE1", "K562", reps=500, B=999, n_jobs=10, arm_effects=True)
s = planted.summarize(df)
s.to_csv(planted.ROOT / "results/aim3/planted_RPE1_vs_K562_v13.csv", index=False)
print(s.to_string(), flush=True)
print("detection limit (block_x4, nominal 5%):", planted.detection_limit(s, "block_x4"))
print("detection limit (block_v1.0, nominal 5%):", planted.detection_limit(s, "block_v1.0"))
print(f"done {time.time() - t0:.0f}s")
