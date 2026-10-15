"""Reproduce the paper's confirmatory numbers (Table 2) and Figures 3-4 from the deposited
intermediates -- no raw screens (~95 GB), no ATAC-seq pipeline. Under a minute on a laptop.

    python scripts/reproduce_from_intermediates.py           # everything
    python scripts/reproduce_from_intermediates.py --quick   # B = 999 instead of 9,999 (p-values then
                                                             # differ in the last digits; theta-hat, k,
                                                             # partial R2 are exact)

Steps
1. Install the deposited intermediates (intermediates/ -> data/, results/aim1/; missing files only).
2. Main test (frozen code: aim2.load_pair, aim2.test_feature; block permutation, x4 rule, arm terms,
   B = 9,999, frozen seeds) for the combined CPM feature in all six pairs, plus the FC feature in the
   two matched pairs; Holm across the two matched pairs (one-sided theta > 0).
3. Pooled matched and pooled mixed contrasts (aim4.pooled_shared_common, frozen seeds).
4. Same-line control: each pair's accessibility difference against M = K562 essential - K562
   genome-wide (as scripts/run_same_line_control.py).
5. Compare every recomputed number with the committed result files; write
   results_repro/COMPARISON.md and results_repro/table2.{tsv,md}; exit 1 on any mismatch.
6. Draw Figures 3 and 4 from the recomputed files (results_repro/figures/).

Not recomputed (taken from the committed result files, marked in COMPARISON.md): the bootstrap
interval of partial R2/reliability (needs per-knockdown terms from the raw screens), the
planted-signal detection limits (a separate 1-2 h benchmark: scripts/run_aim4.py planted A B), and
the exploratory RPE1 diagnostic in Figure 4C and the 4DN fold-change rows of Figure 4B (need the
RPE1 insertion files and the 4DN bigWig from the ATAC-seq pipeline).
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from cdtest import aim2, aim4, crt, features, intermediates  # noqa: E402
from cdtest import covariates as cov  # noqa: E402

OUT = ROOT / "results_repro"
R4 = ROOT / "results/aim4"
SAME_TAG = "K562_vs_K562gw"
#: tolerances for the comparison: theta-hat, CI ends, partial R2 are deterministic up to floating point;
#: p-values use the frozen seeds and must match exactly when B = 9,999
TOL_REL = 1e-6


def jl(p: Path) -> dict:
    return json.loads(Path(p).read_text())


def jw(p: Path, d: dict):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(d, indent=1, default=aim2._json_default))


def holm(ps):
    ps = np.asarray(ps, float)
    o = np.argsort(ps)
    adj, run = np.empty_like(ps), 0.0
    for i, j in enumerate(o):
        run = max(run, min(1.0, (len(ps) - i) * ps[j]))
        adj[j] = run
    return adj


def main_test(P: dict, label: str, B: int) -> dict:
    """theta-hat, k, CI and p of the frozen main test only (the comparators are not needed here)."""
    M, x, Q, gen = P["M"], P["X"][label], P["Q"], P["gen"]
    prep = crt.prepare(M, x, Q=Q)
    chrom, pos, arm = gen.chrom.to_numpy(), gen.pos.to_numpy(), gen.arm.to_numpy()
    order = crt.genome_order(chrom, pos)
    k = crt.choose_k(prep.r_x[order], arm[order], rule="x4")
    grid = crt.theta_grid(prep.theta_hat, float(np.std(x)))
    rng = aim2.job_seed("aim2", P["tag"], label)
    res = crt.block_permutation_test(M, x, chrom=chrom, pos=pos, arm=arm, k=k, B=B, rng=rng, prep=prep, order=order)
    pr2 = crt.partial_r2(M, x, prep=prep)
    rel = P["summary"]["reliability"]
    return {"feature": label, "G": prep.G, "theta_hat": prep.theta_hat, "k_x4": int(k), "sd_x": float(np.std(x)),
            "main_block_x4": aim2._summ(res, grid), "partial_r2_raw": pr2, "reliability": rel,
            "partial_r2_share_reproducible": pr2 / rel}


def same_line_pair() -> dict:
    gs = pd.read_csv(ROOT / f"results/aim1/{SAME_TAG}/gene_outcome.tsv", sep="\t")
    genes = gs.gene_id.tolist()
    C = cov.build_covariates("K562", "K562gw", genes, SAME_TAG)
    Z = cov.design(C)
    P = dict(a="K562", b="K562gw", tag=SAME_TAG, genes=genes, M=gs.M.to_numpy(float), gs=gs, C=C, Z=Z,
             gen=cov.genome_info(genes), Q=crt.orthonormal_basis(Z),
             summary=jl(ROOT / f"results/aim1/{SAME_TAG}/summary.json"), X={})
    for a, b, _ in aim4.PAIRS:
        x = features.pair_feature(a, b, "combined", signal="cpm").reindex(genes).to_numpy(float)
        P["X"][aim4.tag(a, b)] = np.where(np.isfinite(x), x, np.nanmean(x))
    return P


class Compare:
    def __init__(self):
        self.rows, self.bad = [], 0

    def add(self, what: str, new, old, exact: bool = False, note: str = ""):
        if old is None:
            ok, d = True, float("nan")
        elif exact:
            ok, d = new == old, abs(new - old)
        else:
            d = abs(new - old)
            ok = d <= TOL_REL * max(1.0, abs(old)) or d <= 1e-12
        self.bad += not ok
        self.rows.append((what, new, old, d, "ok" if ok else "MISMATCH", note))

    def write(self, path: Path, header: list[str]):
        L = header + ["", "| quantity | recomputed | committed | abs. difference | status |",
                      "|---|---|---|---|---|"]
        for w, n, o, d, s, note in self.rows:
            L.append(f"| {w} | {n:.10g} | {'' if o is None else f'{o:.10g}'} | {d:.2g} | {s}{' ' + note if note else ''} |")
        L += ["", f"**{len(self.rows) - self.bad} of {len(self.rows)} numbers agree**"
              + ("" if not self.bad else f"; {self.bad} mismatch(es)") + "."]
        path.write_text("\n".join(L) + "\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--quick", action="store_true", help="B = 999 (fast check; p-values not comparable)")
    ap.add_argument("--no-figures", action="store_true")
    a = ap.parse_args(argv)
    B = 999 if a.quick else 9999
    t0 = time.time()
    intermediates.install()
    OUT.mkdir(exist_ok=True)
    cmp = Compare()
    exact_p = not a.quick

    # ---- 2. main test, six pairs (+ FC in the matched pairs)
    Ps, tests = {}, {}
    for A, Bn, role in aim4.PAIRS:
        t = aim4.tag(A, Bn)
        P = aim2.load_pair(A, Bn, arm_effects=True)
        Ps[t] = P
        feats = {"combined_cpm": main_test(P, "combined_cpm", B)}
        if role == "matched":
            feats["combined_fc"] = main_test(P, "combined_fc", B)
        tests[t] = feats
        old_p = R4 / f"tests/{t}_test_v13.json"
        old = jl(old_p if old_p.exists() else ROOT / f"results/aim2/{t}_test_v13.json")["features"]
        for lab, r in feats.items():
            o = old[lab]
            cmp.add(f"{t} {lab} θ̂", r["theta_hat"], o["theta_hat"])
            cmp.add(f"{t} {lab} k (x4)", r["k_x4"], o["k_x4"], exact=True)
            cmp.add(f"{t} {lab} partial R²/reliability", r["partial_r2_share_reproducible"], o["partial_r2_share_reproducible"])
            for i, end in enumerate(("lower", "upper")):
                cmp.add(f"{t} {lab} 95% CI {end}", r["main_block_x4"]["ci95"][i], o["main_block_x4"]["ci95"][i])
            for side in ("p_one_sided", "p_two_sided"):
                cmp.add(f"{t} {lab} {side}", r["main_block_x4"][side], o["main_block_x4"][side] if exact_p else None,
                        exact=True, note="" if exact_p else "(B = 999, not compared)")
        jw(OUT / f"aim4/tests/{t}_test_v13.json", {"pair": t, "A": A, "B": Bn, "B_perm": B, "features": feats,
                                                    "source": "reproduce_from_intermediates.py"})
        print(f"[{time.time() - t0:5.0f}s] {t}: θ̂ = {feats['combined_cpm']['theta_hat']:.5g}, "
              f"p₂ = {feats['combined_cpm']['main_block_x4']['p_two_sided']}", flush=True)

    matched = [aim4.tag(x, y) for x, y, r in aim4.PAIRS if r == "matched"]
    mixed = [aim4.tag(x, y) for x, y, r in aim4.PAIRS if r == "mixed"]
    p1 = [tests[t]["combined_cpm"]["main_block_x4"]["p_one_sided"] for t in matched]
    adj = holm(p1)
    h1_old = jl(R4 / "h1_primary.json")
    h1 = json.loads(json.dumps(h1_old))
    h1["B"] = B
    h1["source"] = ("reproduce_from_intermediates.py; recomputed: theta_hat, k_x4, p, Holm p, CI, partial R2/"
                    "reliability; carried over from the committed file: bootstrap interval, detection limits")
    for t, pa in zip(matched, adj):
        f = tests[t]["combined_cpm"]
        m = f["main_block_x4"]
        h1["pairs"][t].update(theta_hat=f["theta_hat"], k_x4=f["k_x4"], p_one_sided=m["p_one_sided"],
                              p_two_sided=m["p_two_sided"], p_holm=float(pa), ci95_refined=m["ci95"],
                              partial_r2_raw=f["partial_r2_raw"],
                              partial_r2_share_reproducible=f["partial_r2_share_reproducible"])
        cmp.add(f"{t} Holm p (H1)", float(pa), h1_old["pairs"][t]["p_holm"] if exact_p else None, exact=True)
    jw(OUT / "aim4/h1_primary.json", h1)

    # ---- 3. pooled contrasts
    pool_old = jl(R4 / "mixed_pooled.json")
    pool = {}
    for key, tags, seed in (("matched", matched, "pooled_matched_common"), ("mixed", mixed, "pooled_mixed_common")):
        r = aim4.pooled_shared_common([Ps[t] for t in tags], B=B, seed_key=seed)
        r.pop("_Tstar_pairs", None)
        r.pop("_cstar_pairs", None)
        pool[key] = r
        o = pool_old[key]
        cmp.add(f"pooled {key} θ̂", r["T_pooled"], o["T_pooled"])
        cmp.add(f"pooled {key} genes", r["n_common_genes"], o["n_common_genes"], exact=True)
        for i, end in enumerate(("lower", "upper")):
            cmp.add(f"pooled {key} 95% CI {end}", r["ci95"][i], o["ci95"][i])
        cmp.add(f"pooled {key} p_two_sided", r["p_two_sided"], o["p_two_sided"] if exact_p else None, exact=True)
    jw(OUT / "aim4/mixed_pooled.json", pool)
    print(f"[{time.time() - t0:5.0f}s] pooled contrasts done", flush=True)

    # ---- 4. same-line control
    SP = same_line_pair()
    sl_old = jl(R4 / "same_line_control.json")
    sl = {}
    for lab in SP["X"]:
        r = main_test(SP, lab, B)
        m = r["main_block_x4"]
        sl[lab] = dict(theta_hat=r["theta_hat"], k=r["k_x4"], p_two_sided=m["p_two_sided"],
                       p_one_sided=m["p_one_sided"], partial_r2_share=r["partial_r2_share_reproducible"], ci95=m["ci95"])
        cmp.add(f"same-line M ~ {lab} θ̂", r["theta_hat"], sl_old[lab]["theta_hat"])
        cmp.add(f"same-line M ~ {lab} p_two_sided", m["p_two_sided"], sl_old[lab]["p_two_sided"] if exact_p else None,
                exact=True)
    jw(OUT / "aim4/same_line_control.json", sl)
    print(f"[{time.time() - t0:5.0f}s] same-line control done", flush=True)

    # sensitivity file for Fig. 4B: committed file with the recomputed CPM and FC rows
    sens = jl(R4 / "sensitivity_h1.json")
    for t in matched:
        for row, lab in (("main (v1.3, x4, CPM)", "combined_cpm"), ("FC features (combined_fc)", "combined_fc")):
            f = tests[t][lab]
            sens[t][row] = {"theta_hat": f["theta_hat"], "p_one_sided": f["main_block_x4"]["p_one_sided"],
                            "p_two_sided": f["main_block_x4"]["p_two_sided"], "ci95": f["main_block_x4"]["ci95"]}
    jw(OUT / "aim4/sensitivity_h1.json", sens)

    # ---- 5. Table 2 and comparison report
    rows = []
    for t in matched:
        f, d = tests[t]["combined_cpm"], h1["pairs"][t]
        lo, hi = f["main_block_x4"]["ci95"]
        b = d["partial_r2_share_boot_ci95"]
        rows.append({"pair": t.replace("_vs_", " − "), "theta_hat": f"{f['theta_hat']:.4f}",
                     "ci95": f"{lo:.4f} to {hi:.4f}", "p_two_sided": f"{f['main_block_x4']['p_two_sided']:.4f}",
                     "p_holm": f"{d['p_holm']:.4g}", "k_x4": f["k_x4"],
                     "partial_r2_over_reliability": f"{100 * f['partial_r2_share_reproducible']:.3f}% "
                                                    f"[{100 * b[0]:.3f}%, {100 * b[1]:.3f}%]*"})
    tab = pd.DataFrame(rows)
    tab.to_csv(OUT / "table2.tsv", sep="\t", index=False)
    md = ["| " + " | ".join(tab.columns) + " |", "|" + "---|" * len(tab.columns)]
    md += ["| " + " | ".join(map(str, r)) + " |" for r in tab.itertuples(index=False)]
    md += ["", "* bootstrap interval carried over from results/aim4/h1_primary.json (needs the raw screens)."]
    (OUT / "table2.md").write_text("\n".join(md) + "\n")
    print("\n".join(md), flush=True)
    cmp.write(OUT / "COMPARISON.md", [
        "# Reproduction from the deposited intermediates", "",
        f"B = {B:,}; frozen seeds; recomputed with scripts/reproduce_from_intermediates.py in "
        f"{time.time() - t0:.0f} s. Committed values: results/aim4/ and results/aim2/.",
        "Carried over, not recomputed: bootstrap interval of partial R²/reliability, planted-signal detection limits, "
        "exploratory RPE1 diagnostic (Fig. 4C) and 4DN fold-change rows (Fig. 4B)."])
    print(f"comparison: {len(cmp.rows) - cmp.bad}/{len(cmp.rows)} agree -> results_repro/COMPARISON.md", flush=True)

    # ---- 6. figures from the recomputed files
    if not a.no_figures:
        env = dict(os.environ, CDTEST_RESULTS_OVERLAY=str(OUT), CDTEST_FIG_OUT=str(OUT / "figures"),
                   PYTHONPATH=str(ROOT / "src"))
        subprocess.run([sys.executable, str(ROOT / "scripts/make_figures.py"), "Fig3", "Fig4"], env=env, check=True)
    print(f"done in {time.time() - t0:.0f} s", flush=True)
    return 0 if cmp.bad == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
