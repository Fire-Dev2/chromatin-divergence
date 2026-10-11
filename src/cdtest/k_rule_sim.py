"""Plan v1.1 block-size rule selection (protocol: results/aim3/k_rule_selection_protocol.md).

Targeted simulation of the shared-structure settings that drove the v1.0 excess,
across four module designs, comparing the k rules in :data:`cdtest.crt.K_RULES`.
Generator: :func:`cdtest.simulate.simulate_replicate` (no decoys).

CLI
---
python -m cdtest.k_rule_sim run --n-jobs 8          # resumable (results/aim3/_parts_krule)
python -m cdtest.k_rule_sim report                  # csv, md tables, figure
"""
from __future__ import annotations

import argparse
import itertools
import os
import time
from dataclasses import asdict

import numpy as np
import pandas as pd

from . import crt
from .simulate import MASTER_SEED, SPLINE_COLS, Setting, simulate_replicate

STUDY_TAG = 11
RULES = crt.K_RULES
SIMPLICITY = ("v1.0", "acf01", "x2", "x3", "x4", "tau1", "tau2", "pw")  # fixed in the protocol
REFS = ("oracle",)
DESIGNS = {"geom5": ("geom", 5.0), "geom10": ("geom", 10.0), "geom25": ("geom", 25.0),
           "mixed": ("mixed", 10.0)}
RHO_M, RHO_X, CONFOUND, NOISE = (0.3, 0.5), (0.3, 0.6), (0.0, 0.3), ("homo", "hetero")
EFFECTS = {"zero": 1000, "small": 300, "moderate": 300}
KEYS = ["design", "rho_M", "rho_x", "confound", "noise"]
FPR_MAX05, FPR_MAX01, TIE = 0.065, 0.02, 0.02


def cells():
    """(cell index, design name, Setting) for all 64 x 3 effect cells."""
    out = []
    for i, (d, rm, rx, cf, nz, ef) in enumerate(itertools.product(
            DESIGNS, RHO_M, RHO_X, CONFOUND, NOISE, EFFECTS)):
        md, mm = DESIGNS[d]
        out.append((i, d, Setting(rm, rx, cf, nz, ef, 2000, mm, md)))
    return out


def analyze(sim: dict, B: int, rng) -> dict:
    M, x, C, gen = sim["M"], sim["x"], sim["C"], sim["genome"]
    theta = sim["theta"]
    prep = crt.prepare(M, x, crt.build_design(C, SPLINE_COLS, df=4))
    chrom, arm, pos = gen["chrom"].to_numpy(), gen["arm"].to_numpy(), gen["pos"].to_numpy()
    order = crt.genome_order(chrom, pos)
    r_sorted, grp = prep.r_x[order], arm[order]
    out = {"theta_hat": prep.theta_hat}
    cache = {}
    for rule in RULES:
        k = crt.choose_k(r_sorted, grp, rule=rule)
        if k not in cache:  # same k -> same procedure; reuse the draw
            cache[k] = crt.block_permutation_test(M, x, chrom=chrom, pos=pos, arm=arm, k=k, B=B,
                                                  rng=rng, prep=prep, order=order)
        out[f"k_{rule}"] = k
        res = cache[k]
        out[f"p_{rule}"], out[f"p2_{rule}"] = res.p, res.p_two_sided
        out[f"cover_{rule}"] = res.covers(theta)
    res = crt.block_permutation_test(M, x, B=B, rng=rng, prep=prep, order=order,
                                     blocks=gen["module"].to_numpy())
    out["p_oracle"], out["p2_oracle"], out["cover_oracle"] = res.p, res.p_two_sided, res.covers(theta)
    sz = np.bincount(gen["module"].to_numpy())
    out["module_mean_realized"] = float(sz.mean())
    out["module_wmean_realized"] = float((sz ** 2).sum() / sz.sum())  # gene-weighted
    out["module_max_realized"] = int(sz.max())
    return out


def _batch(ci: int, design: str, s: Setting, reps, B: int, seed: int, part_dir: str):
    part = os.path.join(part_dir, f"c{ci:03d}_B{B}_r{reps[0]:05d}-{reps[-1]:05d}.parquet")
    if os.path.exists(part):
        return pd.read_parquet(part)
    rows = []
    for r in reps:
        rng = np.random.default_rng(np.random.SeedSequence([seed, STUDY_TAG, ci, r]))
        sim = simulate_replicate(s, rng, n_decoys=0)
        row = analyze(sim, B, rng)
        row.update(asdict(s))
        row.update({"design": design, "cell": ci, "rep": r, "theta": s.theta, "B": B})
        rows.append(row)
    df = pd.DataFrame(rows)
    df.to_parquet(part + ".tmp", index=False)
    os.replace(part + ".tmp", part)
    return df


def run(B=999, n_jobs=8, seed=MASTER_SEED, batch=50, part_dir="results/aim3/_parts_krule",
        out="results/aim3/k_rule_grid.parquet", scale=1.0):
    from joblib import Parallel, delayed
    os.makedirs(part_dir, exist_ok=True)
    tasks = []
    for ci, d, s in cells():
        n = max(1, int(EFFECTS[s.effect] * scale))
        tasks += [(ci, d, s, list(range(i, min(i + batch, n)))) for i in range(0, n, batch)]
    # null cells first so the bulk of the criterion's evidence arrives early
    tasks.sort(key=lambda t: t[2].effect != "zero")
    t0 = time.time()
    dfs = Parallel(n_jobs=n_jobs, verbose=5)(
        delayed(_batch)(ci, d, s, ids, B, seed, part_dir) for ci, d, s, ids in tasks)
    df = pd.concat(dfs, ignore_index=True)
    df.to_parquet(out, index=False)
    print(f"{len(df)} replicates in {time.time() - t0:.0f} s -> {out}")
    return df


# --------------------------------------------------------------------------
# Scoring and selection
# --------------------------------------------------------------------------

def _rate(v):
    v = np.asarray(v, float)
    p = v.mean()
    return p, np.sqrt(p * (1 - p) / v.size)


def per_cell(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for kv, g in df.groupby(KEYS):
        null = g[g.effect == "zero"]
        for rule in RULES + REFS:
            pn = null[f"p_{rule}"].to_numpy()
            thr = np.quantile(pn, 0.05, method="lower")
            row = dict(zip(KEYS, kv), rule=rule, n_null=len(pn), matched_thr=thr)
            row["fpr05"], row["fpr05_se"] = _rate(pn <= 0.05)
            row["fpr01"], row["fpr01_se"] = _rate(pn <= 0.01)
            row["cover_null"], row["cover_null_se"] = _rate(null[f"cover_{rule}"])
            for ef in ("small", "moderate"):
                ge = g[g.effect == ef]
                p = ge[f"p_{rule}"].to_numpy()
                row[f"n_{ef}"] = len(p)
                row[f"power05_{ef}"], row[f"power05_{ef}_se"] = _rate(p <= 0.05)
                row[f"power_matched_{ef}"], row[f"power_matched_{ef}_se"] = _rate(p <= thr)
                row[f"cover_{ef}"], _ = _rate(ge[f"cover_{rule}"])
            if rule in RULES:
                k = g[f"k_{rule}"]
                row["k_median"], row["k_q25"], row["k_q75"] = k.median(), k.quantile(.25), k.quantile(.75)
            rows.append(row)
    return pd.DataFrame(rows)


def per_rule(pc: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for rule, g in pc.groupby("rule", sort=False):
        n = len(g)
        iw = g.fpr05.idxmax()
        r = {"rule": rule, "n_cells": n,
             "fpr05_mean": g.fpr05.mean(), "fpr05_mean_se": np.sqrt((g.fpr05_se ** 2).sum()) / n,
             "fpr05_worst": g.fpr05.max(), "fpr05_worst_se": g.fpr05_se[iw],
             "fpr05_worst_cell": "/".join(str(g.loc[iw, k]) for k in KEYS),
             "fpr01_mean": g.fpr01.mean(), "fpr01_mean_se": np.sqrt((g.fpr01_se ** 2).sum()) / n,
             "fpr01_worst": g.fpr01.max(), "fpr01_worst_se": g.fpr01_se[g.fpr01.idxmax()],
             "coverage_mean": g[["cover_null", "cover_small", "cover_moderate"]].to_numpy().mean(),
             "coverage_min": g[["cover_null", "cover_small", "cover_moderate"]].to_numpy().min()}
        for ef in ("small", "moderate"):
            r[f"power_matched_{ef}"] = g[f"power_matched_{ef}"].mean()
            r[f"power_matched_{ef}_se"] = np.sqrt((g[f"power_matched_{ef}_se"] ** 2).sum()) / n
            r[f"power05_{ef}"] = g[f"power05_{ef}"].mean()
        for d in DESIGNS:
            gd = g[g.design == d]
            r[f"fpr05_worst_{d}"] = gd.fpr05.max()
            if rule in RULES:
                r[f"k_median_{d}"] = gd.k_median.median()
        rows.append(r)
    return pd.DataFrame(rows)


def select_rule(pr: pd.DataFrame) -> dict:
    """Apply the pre-set criterion (protocol section 'Pre-set criterion')."""
    cand = pr[pr.rule.isin(RULES)].set_index("rule")
    elig = cand[(cand.fpr05_worst <= FPR_MAX05) & (cand.fpr01_worst <= FPR_MAX01)]
    if elig.empty:
        best = cand.fpr05_worst.idxmin()
        return {"criterion_met": False, "chosen": None, "eligible": [], "best_fpr_rule": best}
    top = elig.power_matched_small.max()
    tied = [r for r in SIMPLICITY if r in elig.index and elig.power_matched_small[r] >= top - TIE]
    return {"criterion_met": True, "chosen": tied[0], "eligible": list(elig.index),
            "tied": tied, "best_power_rule": elig.power_matched_small.idxmax(), "best_power": top}


def plot(pc: pd.DataFrame, pr: pd.DataFrame, path: str, chosen=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    order = list(RULES) + list(REFS)
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.8))
    cols = plt.cm.tab10(np.arange(len(order)))
    xs = np.arange(len(DESIGNS))
    w = 0.8 / len(order)
    for i, rule in enumerate(order):
        g = pc[pc.rule == rule].groupby("design").fpr05
        worst = g.max().reindex(DESIGNS).to_numpy()
        mean = g.mean().reindex(DESIGNS).to_numpy()
        ax[0].bar(xs + i * w, worst, w, color=cols[i], label=rule + (" (v1.1)" if rule == chosen else ""))
        ax[0].scatter(xs + i * w, mean, color="k", s=6, zorder=3)
    ax[0].axhline(0.05, color="k", lw=0.8)
    ax[0].axhline(FPR_MAX05, color="r", ls="--", lw=0.8)
    ax[0].set_xticks(xs + 0.4 - w / 2, list(DESIGNS))
    ax[0].set_ylabel("FPR at alpha 0.05 (bar: worst cell; dot: mean)")
    ax[0].set_title("Size by module design (16 null cells each, n=1000)")
    ax[0].set_ylim(0, 0.21)
    ax[0].legend(fontsize=7, ncol=3, loc="upper center")
    p = pr.set_index("rule").reindex(order)
    for i, rule in enumerate(order):
        ax[1].errorbar(p.fpr05_worst[rule], p.power_matched_small[rule],
                       xerr=p.fpr05_worst_se[rule], yerr=1.96 * p.power_matched_small_se[rule],
                       fmt="o", color=cols[i], ms=7 if rule == chosen else 5)
        ax[1].annotate(rule, (p.fpr05_worst[rule], p.power_matched_small[rule]), fontsize=8,
                       xytext=(4, 3), textcoords="offset points")
    ax[1].axvline(FPR_MAX05, color="r", ls="--", lw=0.8)
    ax[1].set_xlabel("worst-case FPR at 0.05 (64 null cells)")
    ax[1].set_ylabel("mean power at matched FPR, partial R2 0.5%")
    ax[1].set_title("Selection: eligible = left of red line")
    for i, rule in enumerate(RULES):
        km = [pr.set_index("rule").loc[rule, f"k_median_{d}"] for d in DESIGNS]
        ax[2].plot(xs, km, marker="o", color=cols[i], label=rule)
    ax[2].set_xticks(xs, list(DESIGNS))
    ax[2].set_ylabel("median k chosen")
    ax[2].set_title("Block size chosen")
    ax[2].legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(path, dpi=150)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["run", "report"])
    ap.add_argument("--n-jobs", type=int, default=8)
    ap.add_argument("--B", type=int, default=999)
    ap.add_argument("--scale", type=float, default=1.0, help="fraction of replicates (smoke runs)")
    ap.add_argument("--grid", default="results/aim3/k_rule_grid.parquet")
    ap.add_argument("--part-dir", default="results/aim3/_parts_krule")
    a = ap.parse_args(argv)
    if a.cmd == "run":
        run(a.B, a.n_jobs, part_dir=a.part_dir, out=a.grid, scale=a.scale)
    else:
        df = pd.read_parquet(a.grid)
        pc = per_cell(df)
        pr = per_rule(pc)
        sel = select_rule(pr)
        pc.to_csv("results/aim3/k_rule_selection_by_cell.csv", index=False)
        pr.to_csv("results/aim3/k_rule_selection.csv", index=False)
        # criterion met -> chosen rule; otherwise the protocol's fallback (lowest worst-case FPR)
        plot(pc, pr, "figures/aim3_k_rule_selection.png", sel.get("chosen") or sel.get("best_fpr_rule"))
        print(sel)


if __name__ == "__main__":
    main()
