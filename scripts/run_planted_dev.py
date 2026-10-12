"""Aim 3 planted-signal benchmark on the development pair RPE1-K562 (500 reps/level, B = 999).

Run: PYTHONPATH=src python scripts/run_planted_dev.py
"""
import time

from cdtest import planted

t0 = time.time()
df = planted.run("RPE1", "K562", reps=500, B=999, n_jobs=10)
s = planted.summarize(df)
s.to_csv(planted.ROOT / "results/aim3/planted_RPE1_vs_K562.csv", index=False)
print(s.to_string(), flush=True)
print("detection limit (block_x4, nominal 5%):", planted.detection_limit(s, "block_x4"))
print("detection limit (block_v1.0, nominal 5%):", planted.detection_limit(s, "block_v1.0"))
print(f"done {time.time() - t0:.0f}s")
