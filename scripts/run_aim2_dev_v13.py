"""Aim 2 test on the development pair RPE1-K562 under plan v1.3 (arm fixed effects in C).

Writes results/aim2/RPE1_vs_K562_test_v13.json and _null_v13.npz (v1.2 files untouched).
Run: PYTHONPATH=src python scripts/run_aim2_dev_v13.py
"""
import time

from cdtest import aim2

t0 = time.time()
aim2.run("RPE1", "K562", B=9999, n_boot=500, arm_effects=True)
print(f"done {time.time() - t0:.0f}s", flush=True)
