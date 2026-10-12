"""Figures for the development pair (Aim 1 M distribution, Aim 2 null, Aim 3 planted power)."""
from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import screens

ROOT = screens.ROOT
FIG = ROOT / "figures"


def fig_m_distribution(tag: str = "RPE1_vs_K562"):
    gs = pd.read_csv(ROOT / f"results/aim1/{tag}/gene_outcome.tsv", sep="\t")
    annot = pd.read_csv(ROOT / "data/annot/genes_gencode_v50.tsv", sep="\t", low_memory=False,
                        usecols=["gene_id", "gene_name", "p53_fischer2017_core"])
    p53 = set(annot.gene_id[annot.p53_fischer2017_core.astype(str) == "True"])
    names = dict(zip(annot.gene_id, annot.gene_name))
    gate = json.loads((ROOT / f"results/aim1/{tag}/p53_gate.json").read_text())
    is53 = gs.gene_id.isin(p53).to_numpy()
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    lo, hi = np.quantile(gs.M, [0.001, 0.999])
    bins = np.linspace(lo, hi, 120)
    ax[0].hist(gs.M.clip(lo, hi), bins=bins, color="0.6", label=f"all genes (n={len(gs):,})")
    for m in gs.M[is53]:
        ax[0].axvline(min(max(m, lo), hi), color="crimson", lw=0.8, alpha=0.7)
    ax[0].axvline(0, color="k", lw=0.8, ls=":")
    ax[0].set_xlabel("M (RPE1 minus K562)")
    ax[0].set_ylabel("genes")
    ax[0].set_title(f"Aim 1 outcome M, {tag.replace('_vs_', ' vs ')} (red: p53 core targets)")
    ax[0].legend(frameon=False)
    top = gs[is53].nlargest(6, "M")
    for i, (_, r) in enumerate(top.iterrows()):
        ax[0].annotate(names.get(r.gene_id, r.gene_id), (min(r.M, hi), ax[0].get_ylim()[1] * (0.75 - 0.2 * (i % 3))),
                       rotation=90, fontsize=7, color="crimson", ha="right")
    # p53 gate
    o = np.sort(gs.M[is53].to_numpy())
    ax[1].hist(gs.M[~is53].clip(lo, hi), bins=bins, density=True, color="0.7", label="other genes")
    ax[1].hist(o.clip(lo, hi), bins=30, density=True, color="crimson", alpha=0.6,
               label=f"p53 targets (n={len(o)})")
    ax[1].set_title(f"p53 gate: mean M {gate['mean_M_p53']:.3f} vs matched null "
                    f"{gate['null_mean']:.3f}±{gate['null_sd']:.3f}, p={gate['p_one_sided']:.3g}", fontsize=9)
    ax[1].set_xlabel("M")
    ax[1].legend(frameon=False)
    fig.tight_layout()
    FIG.mkdir(exist_ok=True)
    fig.savefig(FIG / "aim1_M_distribution.png", dpi=150)
    plt.close(fig)


def fig_dev_null(tag: str = "RPE1_vs_K562", suffix: str = ""):
    res = json.loads((ROOT / f"results/aim2/{tag}_test{suffix}.json").read_text())
    nulls = np.load(ROOT / f"results/aim2/{tag}_null{suffix}.npz")
    labels = list(res["features"])
    fig, axes = plt.subplots(1, len(labels), figsize=(4 * len(labels), 3.6))
    for ax, lab in zip(axes, labels):
        d = res["features"][lab]
        Ts = nulls[lab]
        T = d["theta_hat"]
        m = d["main_block_x4"]
        lim = max(np.abs(Ts).max(), abs(T)) * 1.1
        ax.hist(Ts, bins=80, range=(-lim, lim), color="steelblue", alpha=0.8,
                label=f"block null (k={d['k_x4']}, B={m['B']:,})")
        ax.axvline(T, color="crimson", lw=2, label=f"observed T={T:.4f}")
        ax.axvline(0, color="k", lw=0.6, ls=":")
        ex = d["exchangeability_x4"]
        ax.text(0.02, 0.97, f"null mean {m['null_mean']:.4f}\n{ex['frac_genes_movable']:.0%} of genes in\nmovable blocks",
                transform=ax.transAxes, fontsize=7, va="top")
        ax.set_title(f"{lab}\none-sided p={m['p_one_sided']:.4f}; two-sided p={m['p_two_sided']:.4f}", fontsize=9)
        ax.set_xlabel("T* (OLS coefficient on x | f(C))")
        ax.legend(frameon=False, fontsize=7, loc="upper right")
    plan = " plan v1.3 (arm effects in C)" if suffix == "_v13" else ""
    fig.suptitle(f"Aim 2, development pair {tag.replace('_vs_', ' - ')} (mixed pair, descriptive){plan}", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / f"aim2_dev_pair_null{suffix}.png", dpi=150)
    plt.close(fig)


def fig_planted(tag: str = "RPE1_vs_K562", suffix: str = ""):
    s = pd.read_csv(ROOT / f"results/aim3/planted_{tag}{suffix}.csv")
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    styles = {"block_x4": ("k", "-", 2.2), "block_v1.0": ("0.4", "--", 1.6)}
    for t, g in s.groupby("test", sort=False):
        g = g.sort_values("level_partial_r2")
        c, ls, lw = styles.get(t, (None, "-", 1.0))
        x = g.level_partial_r2 * 100
        ax[0].errorbar(x, g.rej05, yerr=1.96 * g.rej05_se_cluster.fillna(g.rej05_se), color=c, ls=ls, lw=lw,
                       marker="o", ms=3, capsize=2, label=t)
        ax[1].plot(x, g.power_matched, color=c, ls=ls, lw=lw, marker="o", ms=3, label=t)
    for a in ax:
        a.axhline(0.8, color="0.7", lw=0.6, ls=":")
        a.set_xlabel("planted partial R² (%)")
        a.set_ylim(-0.02, 1.02)
    ax[0].axhline(0.05, color="crimson", lw=0.6, ls=":")
    ax[0].set_ylabel("rejection rate at nominal 5% (level 0 = FPR)")
    ax[1].set_ylabel("power at matched 5% FPR")
    ax[0].legend(frameon=False, fontsize=7)
    plan = ", plan v1.3" if suffix == "_v13" else ""
    ax[0].set_title(f"Planted signal on real M ({tag.replace('_vs_', ' - ')}{plan}); bars = 95% cluster MC", fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG / f"aim3_planted_power{suffix}.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    import sys
    which = sys.argv[1:] or ["m", "null", "planted"]
    if "m" in which:
        fig_m_distribution()
    if "null" in which:
        fig_dev_null()
    if "planted" in which:
        fig_planted()
    if "v13" in which:
        fig_dev_null(suffix="_v13")
        fig_planted(suffix="_v13")
