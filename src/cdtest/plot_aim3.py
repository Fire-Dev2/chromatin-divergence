"""Central Aim 3 figure: false-positive rate vs strength of shared structure, one line per test.

python -m cdtest.plot_aim3 --grid results/aim3/synthetic_grid.parquet --out figures/aim3_fpr_vs_structure
"""
from __future__ import annotations

import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

STYLE = {  # test: (label, colour, linestyle, marker)
    "naive_ols": ("OLS, genes independent", "#d62728", "-", "o"),
    "unrestricted": ("Unrestricted permutation", "#ff7f0e", "-", "s"),
    "expr_bin": ("Expression-bin permutation", "#bcbd22", "-", "^"),
    "stratified": ("Stratified (expr x CN)", "#8c564b", "-", "v"),
    "camera": ("CAMERA-style VIF", "#9467bd", "-", "D"),
    "decoy": ("Decoy features", "#17becf", "-", "P"),
    "block_main": ("Block permutation (main, k rule)", "#000000", "-", "o"),
    "block_oracle": ("Block = true modules (oracle)", "#2ca02c", "-", "*"),
    "block_k_div4": ("Wrong blocks: k/4", "#555555", "--", "x"),
    "block_k_x4": ("Wrong blocks: 4k", "#555555", ":", "+"),
    "block_random_order": ("Wrong blocks: random order", "#aaaaaa", "--", "1"),
}


def fpr_table(df: pd.DataFrame, level: float = 0.05) -> pd.DataFrame:
    null = df[df.effect.isin(["zero", "null"])]
    rows = []
    for (rm, rx, cf), g in null.groupby(["rho_M", "rho_x", "confound"]):
        for t in STYLE:
            rej = (g[f"p_{t}"] <= level).to_numpy()
            p = rej.mean()
            rows.append({"rho_M": rm, "rho_x": rx, "confound": cf, "test": t, "fpr": p,
                         "se": np.sqrt(p * (1 - p) / rej.size), "n": rej.size})
    return pd.DataFrame(rows)


def make_figure(df: pd.DataFrame, out: str, level: float = 0.05):
    tab = fpr_table(df, level)
    rhos_x = sorted(tab.rho_x.unique())
    confs = sorted(tab.confound.unique())
    fig, axes = plt.subplots(len(confs), len(rhos_x), figsize=(4.2 * len(rhos_x), 3.6 * len(confs)),
                             sharex=True, sharey=True, squeeze=False)
    ymax = min(1.0, max(0.12, tab.fpr.max() * 1.1))
    for i, cf in enumerate(confs):
        for j, rx in enumerate(rhos_x):
            ax = axes[i, j]
            sub = tab[(tab.confound == cf) & (tab.rho_x == rx)]
            ax.axhspan(level - 1.96 * np.sqrt(level * (1 - level) / sub.n.max()),
                       level + 1.96 * np.sqrt(level * (1 - level) / sub.n.max()),
                       color="0.9", zorder=0)
            ax.axhline(level, color="0.5", lw=0.8, zorder=1)
            for t, (lab, col, ls, mk) in STYLE.items():
                s = sub[sub.test == t].sort_values("rho_M")
                lw = 2.2 if t == "block_main" else 1.2
                ax.errorbar(s.rho_M, s.fpr, yerr=1.96 * s.se, color=col, ls=ls, marker=mk, lw=lw,
                            ms=4, capsize=2, label=lab)
            ax.set_title(f"feature shares modules: rho_x = {rx:g}\nconfounding = {cf:g}", fontsize=9)
            ax.set_ylim(0, ymax)
            if i == len(confs) - 1:
                ax.set_xlabel("response module correlation rho_M")
            if j == 0:
                ax.set_ylabel(f"false-positive rate at {level:g}")
    n = int(tab.n.max())
    fig.suptitle(f"Aim 3 synthetic grid: false-positive rate vs shared structure "
                 f"(theta = 0; {n} replicates per point, both noise models pooled; bars = 95% MC interval)",
                 fontsize=10)
    h, l = axes[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=4, fontsize=8, frameon=False)
    fig.tight_layout(rect=(0, 0.1, 1, 0.95))
    for ext in ("png", "pdf"):
        fig.savefig(f"{out}.{ext}", dpi=200)
    plt.close(fig)
    return tab


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--grid", default="results/aim3/synthetic_grid.parquet")
    ap.add_argument("--out", default="figures/aim3_fpr_vs_structure")
    ap.add_argument("--table", default="results/aim3/fpr_vs_structure.csv")
    a = ap.parse_args(argv)
    df = pd.read_parquet(a.grid) if a.grid.endswith(".parquet") else pd.read_csv(a.grid)
    tab = make_figure(df, a.out)
    tab.to_csv(a.table, index=False)
    print(f"wrote {a.out}.png/.pdf and {a.table}")


if __name__ == "__main__":
    main()
