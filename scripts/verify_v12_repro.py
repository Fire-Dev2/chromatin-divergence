"""Check that arm_effects=False reproduces the frozen plan-v1.2 results exactly.

Reruns the Aim 2 development-pair test (B = 9,999) and the planted benchmark (500 reps/level,
B = 999) with arm_effects=False into a temporary folder and compares them to
results/aim2/RPE1_vs_K562_test.json and results/aim3/planted_RPE1_vs_K562_replicates.parquet.
Writes results/aim2/v12_reproduction.json.

Run: PYTHONPATH=src python scripts/verify_v12_repro.py
"""
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

from cdtest import aim2, planted

ROOT = Path(__file__).resolve().parents[1]
TMP = ROOT / "data/tmp_v12_repro"
TMP.mkdir(parents=True, exist_ok=True)


def diff(a, b, path="", out=None):
    """Max abs numeric difference and list of non-numeric mismatches between two JSON trees."""
    out = out if out is not None else {"max_abs": 0.0, "mismatch": []}
    if isinstance(a, dict) and isinstance(b, dict):
        if set(a) != set(b):
            out["mismatch"].append(f"{path}: keys {sorted(set(a) ^ set(b))}")
        for k in set(a) & set(b):
            diff(a[k], b[k], f"{path}/{k}", out)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            out["mismatch"].append(f"{path}: len")
        for i, (u, v) in enumerate(zip(a, b)):
            diff(u, v, f"{path}[{i}]", out)
    elif isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool):
        if not (np.isnan(a) and np.isnan(b)):
            out["max_abs"] = max(out["max_abs"], abs(a - b))
    elif a != b:
        out["mismatch"].append(f"{path}: {a!r} != {b!r}")
    return out


aim2.run("RPE1", "K562", B=9999, n_boot=500, out_dir=TMP, arm_effects=False)
old = json.loads((ROOT / "results/aim2/RPE1_vs_K562_test.json").read_text())
new = json.loads((TMP / "RPE1_vs_K562_test.json").read_text())
d2 = diff(old, new)
on, nn = np.load(ROOT / "results/aim2/RPE1_vs_K562_null.npz"), np.load(TMP / "RPE1_vs_K562_null.npz")
null_max = max(float(np.max(np.abs(on[k] - nn[k]))) for k in on.files)

df = planted.run("RPE1", "K562", reps=500, B=999, n_jobs=10, out_dir=TMP, arm_effects=False)
od = pd.read_parquet(ROOT / "results/aim3/planted_RPE1_vs_K562_replicates.parquet")
num = od.select_dtypes("number").columns
pl_max = float(np.nanmax(np.abs(od[num].to_numpy(float) - df[num].to_numpy(float))))
same_meta = bool((od.select_dtypes(exclude="number").astype(str).values
                  == df.select_dtypes(exclude="number").astype(str).values).all())
oc = pd.read_csv(ROOT / "results/aim3/planted_RPE1_vs_K562.csv")
nc = planted.summarize(df)
csv_max = float(np.nanmax(np.abs(oc.select_dtypes("number").to_numpy(float) - nc[oc.select_dtypes("number").columns].to_numpy(float))))
rep = {"aim2_json_max_abs_diff": d2["max_abs"], "aim2_json_mismatches": d2["mismatch"],
       "aim2_null_max_abs_diff": null_max, "planted_replicates_max_abs_diff": pl_max,
       "planted_replicates_nonnumeric_identical": same_meta, "planted_summary_max_abs_diff": csv_max}
rep["reproduced"] = bool(d2["max_abs"] < 1e-9 and not d2["mismatch"] and null_max < 1e-12 and pl_max < 1e-9
                         and same_meta and csv_max < 1e-9)
(ROOT / "results/aim2/v12_reproduction.json").write_text(json.dumps(rep, indent=1))
print(json.dumps(rep, indent=1))
shutil.rmtree(TMP)
