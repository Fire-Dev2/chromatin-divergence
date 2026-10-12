"""Same-line baseline (plan §0): K562 essential (A) vs K562 genome-wide (B).

Run only after data/fetch_gwps_parallel.log prints GWPS_OK.
    PYTHONPATH=src python scripts/run_same_line.py
Steps: verify the genome-wide h5ad has the essential screen's obs/var fields and layout;
streaming pseudobulk pass restricted to the essential screen's knockdowns
(screens.RESTRICT_TO); aim1.run_pair('K562', 'K562gw'); write
results/aim1/same_line_baseline.json next to the cross-line RPE1-K562 values.
"""
import json
import time

import h5py
import numpy as np
import pandas as pd

from cdtest import aim1, screens

ROOT = screens.ROOT
log = (ROOT / "data/fetch_gwps_parallel.log").read_text() if (ROOT / "data/fetch_gwps_parallel.log").exists() else ""
if "GWPS_OK" not in log:
    raise SystemExit("GWPS_OK not in data/fetch_gwps_parallel.log; not touching the file")

t0 = time.time()
check = {}
for ln in ("K562", "K562gw"):
    with h5py.File(screens.SCREEN_FILES[ln], "r") as h:
        X = h["X"]
        check[ln] = {"obs": sorted(h["obs"].keys()), "var": sorted(h["var"].keys()),
                     "X_type": type(X).__name__, "X_shape": list(X.shape) if hasattr(X, "shape") else None,
                     "X_dtype": str(X.dtype) if hasattr(X, "dtype") else None,
                     "X_chunks": list(X.chunks) if getattr(X, "chunks", None) else None,
                     "X_compression": getattr(X, "compression", None)}
need = ["gene", "gene_id", "gem_group", "UMI_count", "mitopercent"]
missing = [c for c in need if c not in check["K562gw"]["obs"]]
print(json.dumps(check, indent=1), flush=True)
if missing or check["K562gw"]["X_type"] != "Dataset":
    raise SystemExit(f"K562gw layout differs: missing obs {missing}, X {check['K562gw']['X_type']}")
obs_gw, _ = screens.read_obs_var(screens.SCREEN_FILES["K562gw"])
n_ctrl = int((obs_gw.gene == screens.CONTROL).sum())
print(f"K562gw: {len(obs_gw):,} cells, {n_ctrl:,} '{screens.CONTROL}' cells, "
      f"{obs_gw.gene.nunique():,} gene labels", flush=True)
del obs_gw

screens.run("K562gw")
print(f"pseudobulk done {time.time() - t0:.0f}s", flush=True)
same = aim1.run_pair("K562", "K562gw")
cross = json.loads((aim1.OUT / "RPE1_vs_K562" / "summary.json").read_text())

gs_same = pd.read_csv(aim1.OUT / "K562_vs_K562gw" / "gene_outcome.tsv", sep="\t")
gs_cross = pd.read_csv(aim1.OUT / "RPE1_vs_K562" / "gene_outcome.tsv", sep="\t")
common = gs_same.merge(gs_cross, on="gene_id", suffixes=("_same", "_cross"))
pb = json.loads((screens.OUT / "K562gw" / "summary.json").read_text())


def pick(s):
    return {k: s[k] for k in ("n_genes", "n_kd_total", "n_kd_kept", "split_r_median", "reliability",
                              "generic_share_median_A", "generic_share_median_B", "M_mean", "M_sd")}


out = {
    "description": "Same-line baseline: A = K562 essential screen, B = K562 genome-wide screen (Replogle 2022). "
                   "Knockdowns restricted to those in both screens; genes to those passing the filter in both.",
    "layout_check": check, "n_ctrl_K562gw_raw": n_ctrl, "pseudobulk_K562gw": pb,
    "same_line": pick(same), "cross_line_RPE1_vs_K562": pick(cross),
    "ratio_M_sd_same_over_cross": same["M_sd"] / cross["M_sd"],
    "n_genes_common": len(common),
    "corr_M_same_vs_cross_common_genes": float(np.corrcoef(common.M_same, common.M_cross)[0, 1]),
    "M_sd_cross_on_common_genes": float(common.M_cross.std(ddof=0)),
    "M_sd_same_on_common_genes": float(common.M_same.std(ddof=0)),
    "seconds": round(time.time() - t0, 1),
}
(aim1.OUT / "same_line_baseline.json").write_text(json.dumps(out, indent=1, default=float))
print(json.dumps({k: v for k, v in out.items() if k not in ("layout_check", "pseudobulk_K562gw")}, indent=1,
                 default=float))
