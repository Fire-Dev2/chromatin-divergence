"""Development pair RPE1-K562: naive analysis (Aim 1 step 7) and the Aim 2 test.

Run: PYTHONPATH=src python scripts/run_aim2_dev.py
"""
import json
import time

from cdtest import aim2

t0 = time.time()
P = aim2.load_pair("RPE1", "K562")
naive = aim2.naive_analysis(P)
(aim2.AIM1 / "naive_RPE1_vs_K562.json").write_text(json.dumps(naive, indent=1, default=aim2._json_default))
print(json.dumps(naive, indent=1, default=aim2._json_default), flush=True)
res = aim2.run("RPE1", "K562", B=9999, n_boot=500)
print(f"done {time.time() - t0:.0f}s", flush=True)
