#!/usr/bin/env python
"""Final manuscript figures: main Figs 1-4 and supplementary Figs S1-S12.

Every number drawn is read from a result file under results/ (nothing is typed by hand); the only
computation done here is the per-gene table behind Fig. S2 (local vs ENCODE K562 promoter signal),
which is written once to results/atac_local/validation_scatter.tsv and re-used afterwards.

    unset PYTHONPATH; source env/venv/bin/activate; export PYTHONPATH=src
    python scripts/make_figures.py            # all figures
    python scripts/make_figures.py Fig1 FigS2 # a subset

Outputs: figures/final/<name>.pdf (vector) and figures/final/<name>.png (300 dpi).
Style: one sans-serif font (Arial; mathtext mapped to Arial), Okabe-Ito colour-blind-safe palette,
matched pairs blue, mixed pairs orange, same-line pair grey, RPE1-containing pairs marked with a
dagger (†), panel letters A, B, ...
Fig. S10 (Jurkat GSE302221 sensitivity) is skipped while results/atac_local/jurkat_GSE302221.qc.json
is absent.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Patch  # noqa: E402
from scipy import stats  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from cdtest.aim4 import PAIRS as _PAIRS  # noqa: E402

RES = ROOT / "results"
OUT = ROOT / "figures/final"

# ----------------------------------------------------------------------------------------- style
FONT = "Arial"
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": [FONT], "font.size": 7,
    "mathtext.fontset": "custom", "mathtext.rm": FONT, "mathtext.it": f"{FONT}:italic",
    "mathtext.bf": f"{FONT}:bold", "mathtext.sf": FONT, "mathtext.cal": FONT, "mathtext.default": "regular",
    "axes.titlesize": 7, "axes.labelsize": 7, "xtick.labelsize": 6, "ytick.labelsize": 6,
    "legend.fontsize": 6, "legend.frameon": False, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 1.0,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none", "savefig.dpi": 300,
})
# Okabe-Ito
OI = dict(black="#000000", orange="#E69F00", sky="#56B4E9", green="#009E73", yellow="#F0E442",
          blue="#0072B2", vermillion="#D55E00", purple="#CC79A7", grey="#999999")
C_MATCHED, C_MIXED, C_SAME = OI["blue"], OI["orange"], OI["grey"]
C_PROM, C_DIST = OI["green"], OI["purple"]
TH = r"$\hat{\theta}$"
MINUS = "\u2212"
W1, W2 = 3.46, 7.2  # one- and two-column widths (inches)

PAIRS = [f"{a}_vs_{b}" for a, b, _ in _PAIRS]          # frozen order: matched first
ROLE = {f"{a}_vs_{b}": r for a, b, r in _PAIRS}
MATCHED = [t for t in PAIRS if ROLE[t] == "matched"]
MIXED = [t for t in PAIRS if ROLE[t] == "mixed"]
RPE1_PAIRS = ("RPE1_vs_K562", "Jurkat_vs_RPE1", "RPE1_vs_HepG2")


def jl(p):
    return json.loads((ROOT / p).read_text())


def pl(t, dagger=True):
    """Display label for a pair tag, e.g. 'RPE1 − HepG2 †'."""
    s = t.replace("_vs_", f" {MINUS} ")
    return s + (" †" if dagger and "RPE1" in t else "")


def pcol(t):
    return C_MATCHED if ROLE.get(t) == "matched" else C_MIXED


def fmt_p(p):
    return f"p = {p:.2g}" if p >= 0.001 else f"p = {p:.1g}"


def letter(ax, s, dx=-22, dy=6):
    ax.annotate(s, xy=(0, 1), xycoords="axes fraction", xytext=(dx, dy), textcoords="offset points",
                fontsize=9, fontweight="bold", ha="left", va="bottom")


def ptitle(ax, s, dx=-112, dy=8, color="k", weight="normal", size=6.3):
    """Left-aligned panel title placed at a fixed offset from the axes' top-left corner."""
    ax.annotate(s, xy=(0, 1), xycoords="axes fraction", xytext=(dx, dy), textcoords="offset points", fontsize=size,
                color=color, fontweight=weight, ha="left", va="bottom", annotation_clip=False)


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{name}.pdf", metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(OUT / f"{name}.png", dpi=300)
    plt.close(fig)
    print(f"wrote figures/final/{name}.pdf/.png", flush=True)


def role_legend(ax, loc="lower right", extra=()):
    h = [Line2D([], [], color=C_MATCHED, marker="o", ls="", ms=4, label="matched pair"),
         Line2D([], [], color=C_MIXED, marker="o", ls="", ms=4, label="mixed pair"), *extra]
    ax.legend(handles=h, loc=loc, handletextpad=0.3)


def forest(ax, rows, xlabel=TH, zero=True):
    """rows: list of dict(label, est, ci (or None), color, marker, filled). Top row first."""
    n = len(rows)
    for i, r in enumerate(rows):
        y = n - 1 - i
        if r.get("ci") is not None:
            ax.plot(r["ci"], [y, y], color=r["color"], lw=1.2, solid_capstyle="butt")
        ax.plot(r["est"], y, r.get("marker", "o"), color=r["color"], ms=r.get("ms", 4.5),
                mfc=r["color"] if r.get("filled", True) else "white", mew=0.9)
    if zero:
        ax.axvline(0, color="k", lw=0.6, ls="--", zorder=0)
    ax.set_yticks(range(n)[::-1], [r["label"] for r in rows])
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_xlabel(xlabel)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)


# ===================================================================================== Figure 1
def fig1_schematic(ax, B):
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 4.6)
    ax.axis("off")
    pos = {"K562": (0.75, 3.05), "Jurkat": (3.15, 3.05), "RPE1": (0.75, 0.95), "HepG2": (3.15, 0.95)}
    grp = {"K562": "hem", "Jurkat": "hem", "RPE1": "epi", "HepG2": "epi"}
    fill = {"hem": "#E8E8E8", "epi": "#FFFFFF"}

    def box(xy, txt, w=1.0, h=0.5, fc="white", ec="k", lw=0.8, fs=7, bold=True, ls="-"):
        x, y = xy
        ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.03,rounding_size=0.08",
                                    fc=fc, ec=ec, lw=lw, ls=ls, zorder=3))
        ax.text(x, y, txt, ha="center", va="center", fontsize=fs, fontweight="bold" if bold else "normal", zorder=4)

    def link(a, b, color, ls="-", lw=1.6, rad=0.0, label=None, lpos=0.5, loff=(0, 0.13), lc=None):
        (x1, y1), (x2, y2) = pos[a], pos[b]
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-", color=color, lw=lw, ls=ls,
                                     connectionstyle=f"arc3,rad={rad}", shrinkA=13, shrinkB=13, zorder=2))
        if label:
            ax.text(x1 + (x2 - x1) * lpos + loff[0], y1 + (y2 - y1) * lpos + loff[1], label, ha="center",
                    va="center", fontsize=5.8, color=lc or color, zorder=4,
                    bbox=dict(fc="white", ec="none", pad=0.4))

    link("K562", "Jurkat", C_MATCHED, label="matched", loff=(0, 0.2))
    link("RPE1", "HepG2", C_MATCHED, label="matched", loff=(0, -0.2))
    link("K562", "RPE1", C_MIXED, ls="--", label="mixed\n(development pair)", loff=(-0.02, 0), lc=C_MIXED)
    link("Jurkat", "HepG2", C_MIXED, ls="--")
    link("K562", "HepG2", C_MIXED, ls="--")
    link("Jurkat", "RPE1", C_MIXED, ls="--", label="mixed", lpos=0.5, loff=(0, 0))
    for k, xy in pos.items():
        box(xy, k, fc=fill[grp[k]])
    ax.text(1.95, 3.62, "haematopoietic · p53-deficient · suspension", ha="center", fontsize=5.8)
    ax.text(1.95, 0.38, "epithelial · p53 wild-type · adherent", ha="center", fontsize=5.8)
    # same-line pair
    box((0.75, 4.25), "K562 genome-wide", w=1.2, h=0.36, fs=5.6, bold=False, ec=C_SAME, ls="--")
    ax.add_patch(FancyArrowPatch((0.75, 4.07), (0.75, 3.32), arrowstyle="-", color=C_SAME, lw=1.4, ls=":", zorder=2))
    ax.text(1.55, 4.25, "same-line pair\n(K562 essential " + MINUS + " K562 genome-wide)", fontsize=5.4,
            color="#555555", va="center", ha="left", linespacing=1.2)
    ax.text(1.95, 0.02, "Pair A " + MINUS + " B; Perturb-seq knockdowns shared by A and B", ha="center", fontsize=5.6,
            color="#444444")

    # outcome / feature / test boxes
    def tbox(x, y, w, h, title, lines, ec="k"):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.04,rounding_size=0.1", fc="white",
                                    ec=ec, lw=0.8, zorder=3))
        ax.text(x + 0.1, y + h - 0.12, title, fontsize=6.6, fontweight="bold", va="top", zorder=4)
        ax.text(x + 0.1, y + h - 0.42, "\n".join(lines), fontsize=5.7, va="top", zorder=4, linespacing=1.35)

    tbox(4.35, 2.55, 3.05, 1.85, "Outcome M (per gene g)",
         ["M = mean over knockdowns of the noise-corrected",
          "squared knockdown-specific effect on gene g",
          "(squared effect minus its sampling variance;",
          "generic response removed), line A " + MINUS + " line B"])
    tbox(4.35, 0.25, 3.05, 1.95, "Feature x (per gene g)",
         ["baseline ATAC accessibility, A " + MINUS + " B",
          "(Tn5 insertion CPM; fold change = sensitivity)",
          "three features:",
          "  promoter (TSS ± 1 kb)",
          "  distal (2.5" + "\u2013" + "100 kb)",
          "  combined (promoter and distal z, equal weight)"])
    tbox(7.75, 1.0, 2.2, 2.65, "Test",
         [r"$M = f(C) + \theta\, x + \varepsilon$",
          "C: expression, copy number,",
          "chromosome-arm terms",
          "feature-side block permutation",
          "of C-residualised x within",
          f"chromosome arms (B = {B:,})",
          r"H1: $\theta > 0$ in matched pairs"])
    for y0 in (3.45, 1.25):
        ax.add_patch(FancyArrowPatch((7.42, y0), (7.73, 2.33 + (0.3 if y0 > 2 else -0.3)), arrowstyle="-|>",
                                     mutation_scale=7, color="k", lw=0.8, zorder=2))
    for y1 in (3.3, 1.3):
        ax.add_patch(FancyArrowPatch((3.8, 2.0), (4.3, y1), arrowstyle="-|>", mutation_scale=7, color="k", lw=0.8))


def fig1():
    h1 = jl("results/aim4/h1_primary.json")
    gs_ = pd.read_csv(RES / "aim1/RPE1_vs_K562/gene_outcome.tsv", sep="\t", usecols=["gene_id", "M"])
    annot = pd.read_csv(ROOT / "data/annot/genes_gencode_v50.tsv", sep="\t", low_memory=False,
                        usecols=["gene_id", "p53_fischer2017_core"])
    p53 = set(annot.gene_id[annot.p53_fischer2017_core.astype(str) == "True"])
    gate = jl("results/aim1/RPE1_vs_K562/p53_gate.json")
    a1 = pd.DataFrame(jl("results/aim4/aim1_pairs.json")["pairs"]).set_index("pair")
    sl = jl("results/aim1/same_line_baseline.json")["same_line"]
    base = sl["M_sd"] ** 2 * sl["reliability"]

    fig = plt.figure(figsize=(W2, 5.3))
    g = fig.add_gridspec(2, 3, height_ratios=[1.25, 1], width_ratios=[1.35, 0.55, 1.1], hspace=0.38,
                         wspace=0.45, left=0.07, right=0.985, top=0.97, bottom=0.12)
    axA = fig.add_subplot(g[0, :])
    fig1_schematic(axA, h1["B"])
    letter(axA, "A", dx=-18, dy=-8)

    # B: M distribution + p53 gate
    axB = fig.add_subplot(g[1, 0])
    is53 = gs_.gene_id.isin(p53).to_numpy()
    lo, hi = np.quantile(gs_.M, [0.001, 0.999])
    axB.hist(gs_.M.clip(lo, hi), bins=np.linspace(lo, hi, 100), color="#BBBBBB", lw=0)
    ymax = axB.get_ylim()[1]
    for m in gs_.M[is53]:
        axB.plot([min(max(m, lo), hi)] * 2, [0, ymax * 0.12], color=OI["vermillion"], lw=0.7)
    axB.axvline(0, color="k", lw=0.5, ls=":")
    axB.set_xlabel(f"M (RPE1 {MINUS} K562)")
    axB.set_ylabel("genes")
    axB.text(0.98, 0.95, f"{len(gs_):,} genes\n{int(is53.sum())} Fischer core p53 targets (red ticks)",
             transform=axB.transAxes, ha="right", va="top", fontsize=6)
    letter(axB, "B")
    axG = fig.add_subplot(g[1, 1])
    axG.errorbar([0], [gate["null_mean"]], yerr=[[2 * gate["null_sd"]]], fmt="o", color="#777777", ms=4,
                 capsize=2.5, lw=1)
    axG.plot([1], [gate["mean_M_p53"]], "D", color=OI["vermillion"], ms=5)
    pc = jl("results/aim4/p53_control.json")["RPE1_vs_K562"]  # same gate, stored with the number of sets
    assert np.isclose(pc["null_mean_oriented"], gate["null_mean"]) and np.isclose(pc["mean_M_p53_raw"], gate["mean_M_p53"])
    axG.set_xticks([0, 1], [f"random sets\n(n = {pc['n_sets']:,};\nmean ± 2 SD)",
                            f"p53 targets\n(n = {gate['n_p53_targets']})"], fontsize=5.2)
    axG.set_xlim(-0.75, 1.75)
    axG.set_ylim(0, None)
    axG.set_ylabel("mean M")
    axG.set_title(f"p53 gate: z = {gate['z']:.1f}, p = {gate['p_one_sided']:.3f}", fontsize=6.3)

    # C: reproducible variance relative to same-line baseline
    axC = fig.add_subplot(g[1, 2])
    rv = (a1.M_sd ** 2 * a1.reliability)
    # cross-check against the stored value (field holds SD^2 x reliability)
    assert np.allclose(rv, a1.reproducible_var_vs_same_line)
    ratio = (rv / base).reindex(PAIRS)
    xs = np.arange(len(PAIRS) + 1)
    vals = [1.0] + list(ratio.values)
    cols = [C_SAME] + [pcol(t) for t in PAIRS]
    axC.bar(xs, vals, color=cols, width=0.7)
    axC.set_yscale("log")
    axC.set_ylim(0.5, 6000)
    for x, v in zip(xs, vals):
        axC.text(x, v * 1.12, f"{v:.1f}×", ha="center", va="bottom", fontsize=5.6)
    rel = [sl["reliability"]] + list(a1.reliability.reindex(PAIRS))
    for x, r in zip(xs, rel):
        axC.text(x, 0.58, f"{r:.2f}", ha="center", va="bottom", fontsize=5.2, color="white" if x else "k")
    axC.set_xticks(xs, ["same line"] + [pl(t) for t in PAIRS], rotation=40, ha="right", fontsize=5.8)
    axC.set_ylabel("reproducible variance of M\n(SD² × reliability) / same-line")
    axC.set_title("white numbers: split-half reliability of M", fontsize=5.8, color="#444444")
    letter(axC, "C", dx=-30)
    role_legend(axC, loc="upper left", extra=[Patch(color=C_SAME, label="same-line baseline")])
    axC.legend_.set_ncols(3)
    save(fig, "Fig1")


# ===================================================================================== Figure 2
TEST_STYLE = {  # label, colour, linestyle, marker, lw
    "naive_ols": ("OLS, genes independent", OI["vermillion"], "-", "o", 1.0),
    "unrestricted": ("unrestricted permutation", OI["orange"], "-", "s", 1.0),
    "expr_bin": ("expression-bin permutation", OI["sky"], "-", "^", 1.0),
    "stratified": ("stratified permutation", OI["blue"], "-", "v", 1.0),
    "camera": ("CAMERA-style VIF", OI["purple"], "-", "D", 1.0),
    "decoy": ("decoy features", "#999933", "-", "P", 1.0),
    "block_main": ("block permutation (v1.0 rule)", "#555555", "-", "o", 1.4),
    "block_oracle": ("block = true modules (oracle)", OI["green"], "-", "*", 1.0),
    "block_k_div4": ("block, k/4", "#999999", "--", "x", 0.8),
    "block_k_x4": ("block, 4k", OI["black"], ":", "+", 1.2),
    "block_random_order": ("block, random order", "#BBBBBB", "--", "1", 0.8),
}
PLANT_STYLE = {
    "block_x4": ("main test (block, x4 rule)", OI["black"], "-", "o", 1.6),
    "block_v1.0": ("block, v1.0 rule", "#888888", "--", "o", 1.0),
    "naive_ols": TEST_STYLE["naive_ols"], "unrestricted": TEST_STYLE["unrestricted"],
    "expr_bin": TEST_STYLE["expr_bin"], "stratified": TEST_STYLE["stratified"], "camera": TEST_STYLE["camera"],
}


def fig2a(fig, gs):
    tab = pd.read_csv(RES / "aim3/fpr_vs_structure.csv")
    rx_, cf_ = sorted(tab.rho_x.unique()), sorted(tab.confound.unique())
    sub = gs.subgridspec(len(cf_), len(rx_), hspace=0.35, wspace=0.12)
    axes = []
    ymax = min(1.0, tab.fpr.max() * 1.08)
    for i, cf in enumerate(cf_):
        for j, rx in enumerate(rx_):
            ax = fig.add_subplot(sub[i, j])
            axes.append(ax)
            s0 = tab[(tab.confound == cf) & (tab.rho_x == rx)]
            n = s0.n.max()
            ax.axhspan(0.05 - 1.96 * np.sqrt(0.0475 / n), 0.05 + 1.96 * np.sqrt(0.0475 / n), color="0.92", zorder=0)
            ax.axhline(0.05, color="0.5", lw=0.6)
            for t, (lab, col, ls, mk, lw) in TEST_STYLE.items():
                s = s0[s0.test == t].sort_values("rho_M")
                ax.plot(s.rho_M, s.fpr, color=col, ls=ls, marker=mk, lw=lw, ms=2.6, label=lab, mew=0.7)
            ax.set_ylim(0, ymax)
            ax.set_xticks(sorted(tab.rho_M.unique()))
            ax.set_title(rf"$\rho_x$ = {rx:g}, confounding = {cf:g}", fontsize=6)
            if i == len(cf_) - 1:
                ax.set_xlabel(r"response module correlation $\rho_M$")
            else:
                ax.set_xticklabels([])
            if j == 0:
                ax.set_ylabel("false-positive rate (α = 0.05)")
            else:
                ax.set_yticklabels([])
    return axes


LEVELS = np.array([0.0, 0.0025, 0.005, 0.01, 0.02])


def planted_axes(ax_r, ax_p, tables, title):
    """tables: dict test -> DataFrame(level, rej05, se, power_matched); key '__excl' = main test, clean decoys."""
    for t, d in tables.items():
        if t == "__excl":
            lab, col, ls, mk, lw = (r"main test, decoys with significant own $\theta_0$ excluded", OI["black"], ":", "s", 1.2)
        else:
            lab, col, ls, mk, lw = PLANT_STYLE[t]
        x = np.searchsorted(LEVELS, np.round(d.level.values, 6))
        ax_r.errorbar(x, d.rej05, yerr=1.96 * d.se, color=col, ls=ls, marker=mk, lw=lw, ms=2.6, capsize=1.5,
                      elinewidth=0.6, label=lab, zorder=3 if t in ("block_x4", "__excl") else 2)
        if t != "__excl":
            ax_p.plot(x, d.power_matched, color=col, ls=ls, marker=mk, lw=lw, ms=2.6, label=lab,
                      zorder=3 if t == "block_x4" else 2)
    for a in (ax_r, ax_p):
        a.set_ylim(-0.03, 1.03)
        a.set_xlabel("planted partial R² (%)")
        a.set_xticks(range(len(LEVELS)), [f"{v * 100:g}" for v in LEVELS], fontsize=5.3)
        a.set_xlim(-0.4, len(LEVELS) - 0.6)
    ax_r.axhline(0.05, color=OI["vermillion"], lw=0.5, ls=":")
    ax_r.set_ylabel("rejection rate, nominal 5%")
    ax_p.set_ylabel("power at matched 5% FPR")
    ax_r.set_title(title, fontsize=6.5, loc="left")


def planted_from_csv(path, excl_path=None):
    s = pd.read_csv(path)
    out = {}
    for t in PLANT_STYLE:
        d = s[s.test == t].sort_values("level_partial_r2")
        out[t] = pd.DataFrame({"level": d.level_partial_r2.values, "rej05": d.rej05.values,
                               "se": d.rej05_se_cluster.fillna(d.rej05_se).values,
                               "power_matched": d.power_matched.values})
    if excl_path is not None:
        e = pd.read_csv(excl_path)
        d = e[(e.test == "block_x4") & e.excluded.str.startswith("theta0_block_x4")].sort_values("level_partial_r2")
        out["__excl"] = pd.DataFrame({"level": d.level_partial_r2.values, "rej05": d.rej05.values,
                                      "se": d.rej05_se_cluster.values, "power_matched": d.power_matched.values})
    return out


def planted_from_json(path):
    j = jl(path)

    def tab(t, key):
        d = j[key][t]
        lv = sorted(d, key=float)
        return pd.DataFrame({"level": [float(v) for v in lv], "rej05": [d[v]["rej05"] for v in lv],
                             "se": [d[v]["se_cluster"] for v in lv], "power_matched": [d[v]["power_matched"] for v in lv]})

    out = {t: tab(t, "table_all") for t in PLANT_STYLE}
    out["__excl"] = tab("block_x4", "table_excl_main")
    return out, j


def fig2():
    fig = plt.figure(figsize=(W2, 6.9))
    g = fig.add_gridspec(3, 1, height_ratios=[1.35, 0.3, 0.85], hspace=0.08, left=0.08, right=0.985,
                         top=0.965, bottom=0.03)
    axes = fig2a(fig, g[0])
    letter(axes[0], "A", dx=-28)
    h, l = axes[0].get_legend_handles_labels()
    lax = fig.add_subplot(g[1]); lax.axis("off")
    lax.legend(h, l, loc="center", ncol=4, fontsize=5.8, columnspacing=1.2, handlelength=2.2)

    low = g[2].subgridspec(2, 6, height_ratios=[1, 0.22], wspace=0.6, hspace=0.62)
    dev = planted_from_csv(RES / "aim3/planted_RPE1_vs_K562_v13.csv", RES / "aim3/planted_RPE1_vs_K562_v13_excl.csv")
    kj, jkj = planted_from_json("results/aim4/planted_K562_vs_Jurkat.json")
    rh, jrh = planted_from_json("results/aim4/planted_RPE1_vs_HepG2.json")
    info = jl("results/aim3/planted_RPE1_vs_K562_v13_info.json") if (RES / "aim3/planted_RPE1_vs_K562_v13_info.json").exists() else {}
    sets = [("B", dev, f"RPE1 {MINUS} K562 (development pair, mixed)"),
            ("C", kj, f"K562 {MINUS} Jurkat (matched); {len(jkj['decoys_used'])} decoys"),
            ("D", rh, f"RPE1 {MINUS} HepG2 (matched); {len(jrh['decoys_used'])} decoys")]
    first = None
    for k, (L, tabs, title) in enumerate(sets):
        ar, ap = fig.add_subplot(low[0, 2 * k]), fig.add_subplot(low[0, 2 * k + 1])
        planted_axes(ar, ap, tabs, "")
        ar.set_title(title, fontsize=6.2, loc="left", x=-0.25)
        letter(ar, L, dx=-28)
        if first is None:
            first = ar
    h, l = first.get_legend_handles_labels()
    lax2 = fig.add_subplot(low[1, :]); lax2.axis("off")
    lax2.legend(h, l, loc="center", ncol=4, fontsize=5.8, columnspacing=1.2)
    _ = info
    save(fig, "Fig2")


# ===================================================================================== Figure 3
def fig3():
    T = {t: jl(f"results/aim4/tests/{t}_test_v13.json")["features"]["combined_cpm"] for t in PAIRS
         if (RES / f"aim4/tests/{t}_test_v13.json").exists()}
    if "RPE1_vs_K562" not in T:  # development pair re-used from Aim 2 (v1.3)
        T["RPE1_vs_K562"] = jl("results/aim2/RPE1_vs_K562_test_v13.json")["features"]["combined_cpm"]
    pool = jl("results/aim4/mixed_pooled.json")
    H1 = jl("results/aim4/h1_primary.json")
    h1 = H1["pairs"]
    rows = []
    for t in MATCHED + MIXED:
        m = T[t]["main_block_x4"]
        rows.append(dict(label=pl(t), est=T[t]["theta_hat"], ci=m["ci95"], color=pcol(t), p2=m["p_two_sided"], tag=t))
    rows.append(dict(label="pooled matched", est=pool["matched"]["T_pooled"], ci=pool["matched"]["ci95"],
                     color=C_MATCHED, marker="s", p2=pool["matched"]["p_two_sided"], tag="pm"))
    rows.append(dict(label="pooled mixed", est=pool["mixed"]["T_pooled"], ci=pool["mixed"]["ci95"],
                     color=C_MIXED, marker="s", p2=pool["mixed"]["p_two_sided"], tag="px"))
    fig, ax = plt.subplots(figsize=(W2, 3.0))
    fig.subplots_adjust(left=0.13, right=0.53, bottom=0.2, top=0.86)
    forest(ax, rows, xlabel=f"{TH} (combined CPM accessibility difference), 95% CI by test inversion")
    n = len(rows)
    ax.axhline(n - len(MATCHED) - 0.5, color="0.8", lw=0.5)
    ax.axhline(1.5, color="0.8", lw=0.5)
    cols_x = (1.04, 1.25, 1.45, 1.68)
    heads = ("two-sided p", "Holm p (H1)", "detection limit*", "pR² / reliability")
    for x, h in zip(cols_x, heads):
        ax.annotate(h, xy=(x, n - 0.45), xycoords=("axes fraction", "data"), fontsize=6, fontweight="bold",
                    annotation_clip=False)
    for i, r in enumerate(rows):
        y = n - 1 - i
        vals = [f"{r['p2']:.4f}"]
        if r["tag"] in h1:
            d = h1[r["tag"]]
            vals += [f"{d['p_holm']:.2g}",
                     f"{d['detection_limit_all'] * 100:g}% / {d['detection_limit_excl_sig_theta0'] * 100:g}%",
                     f"{d['partial_r2_share_reproducible'] * 100:.3f}%"]
        for x, v in zip(cols_x, vals):
            ax.annotate(v, xy=(x, y), xycoords=("axes fraction", "data"), va="center", fontsize=6,
                        annotation_clip=False)
    fig.text(0.13, 0.95, f"H1, main test: block permutation of C-residualised x (x4 rule, chromosome-arm terms), "
             f"B = {H1['B']:,}; pooled contrasts on the {pool['mixed']['n_common_genes']:,} genes common to all pairs",
             fontsize=6.3)
    fig.text(0.13, 0.025, "† pair contains RPE1.   * planted partial R² detected with power ≥ 0.8: all decoys / decoys "
             "without significant own association with M (Fig. 2C–D)", fontsize=5.8, color="#444444")
    role_legend(ax, loc="lower right", extra=[Line2D([], [], color="0.4", marker="s", ls="", ms=4, label="pooled contrast")])
    save(fig, "Fig3")


# ===================================================================================== Figure 4
DIAG_ROWS = [  # (key, label) -- rows of manuscript Table 6
    ("primary", "primary (reproduced)"),
    ("cn_rpe1_depmap", "RPE1 CN = DepMap clone consensus"),
    ("rpe1_ds_jurkat_depth_116M", "RPE1 thinned to Jurkat depth (3 seeds)"),
    ("excl_chr10_and_chrX", "exclude chr10 and chrX"),
    ("excl_p53_union", "exclude both p53 gene sets"),
    ("promoter", "promoter feature"),
    ("distal", "distal feature"),
    ("side_rpe1_alone", "RPE1 feature alone"),
    ("side_partner_alone", "partner feature alone"),
    ("rpe1_4dn_fc_hybrid", "RPE1 = 4DN FC (others CPM)"),
    ("rpe1_local_fc_hybrid", "RPE1 = local FC (others CPM)"),
    ("partner_fc_hybrid", "partner = FC, RPE1 CPM"),
    ("rpe1_cpm_shared_with_fc", "RPE1 CPM part shared with FC"),
    ("rpe1_cpm_specific_part", "RPE1 CPM part not shared with FC"),
    ("rpe1_expr_top_tercile", "top tercile of RPE1 expression"),
    ("rpe1_expr_bottom_tercile", "bottom tercile of RPE1 expression"),
]
DIAG_FULL = [  # (key, label, group) -- all variants (Fig. S11), as in scripts/aim4x_rpe1_report.py
    ("primary", "primary (frozen main test; reproduced)", "primary"),
    ("cn_rpe1_depmap", "CN: RPE1 = DepMap clone consensus (features + cn_diff)", "copy number"),
    ("cn_rpe1_depmap__features_only", "CN: DepMap, features only", "copy number"),
    ("cn_rpe1_depmap__covariate_only", "CN: DepMap, cn_diff covariate only", "copy number"),
    ("cn_rpe1_none", "CN: no CN correction for RPE1", "copy number"),
    ("rpe1_rep1_only", "RPE1 CPM from rep1 only", "accessibility source"),
    ("rpe1_rep2_only", "RPE1 CPM from rep2 only", "accessibility source"),
    ("rpe1_4dn_fc_hybrid", "RPE1 = 4DN-processed FC (others CPM)", "accessibility source"),
    ("rpe1_local_fc_hybrid", "RPE1 = local-pipeline FC (others CPM)", "accessibility source"),
    ("rpe1_ds_jurkat_depth_116M", "RPE1 thinned to Jurkat depth (mean of 3 seeds)", "accessibility source"),
    ("rpe1_ds_58M_insertions", "RPE1 thinned to 58M insertions (mean of 3 seeds)", "accessibility source"),
    ("fc_all_lines_primary", "all lines FC (frozen FC sensitivity)", "accessibility source"),
    ("fc_all_lines_rpe1_4dn", "all lines FC, RPE1 = 4DN FC + peaks", "accessibility source"),
    ("partner_fc_hybrid", "partner = FC, RPE1 kept CPM", "accessibility source"),
    ("both_fc_hybrid", "partner and RPE1 = FC (others CPM)", "accessibility source"),
    ("rpe1_cpm_shared_with_fc", "RPE1 CPM part shared with RPE1 FC", "accessibility source"),
    ("rpe1_cpm_specific_part", "RPE1 CPM part not shared with RPE1 FC", "accessibility source"),
    ("promoter", "promoter feature", "location"),
    ("distal", "distal feature", "location"),
    ("excl_chr10", "exclude chr10", "chromosomes"),
    ("excl_chrX", "exclude chrX", "chromosomes"),
    ("excl_chr10_and_chrX", "exclude chr10 and chrX", "chromosomes"),
    ("excl_fischer_p53", "exclude Fischer p53 core targets", "p53 / expression"),
    ("excl_hallmark_p53", "exclude HALLMARK_P53_PATHWAY", "p53 / expression"),
    ("excl_p53_union", "exclude both p53 sets", "p53 / expression"),
    ("rpe1_expr_bottom_tercile", "bottom tercile of RPE1 expression", "p53 / expression"),
    ("rpe1_expr_middle_tercile", "middle tercile of RPE1 expression", "p53 / expression"),
    ("rpe1_expr_top_tercile", "top tercile of RPE1 expression", "p53 / expression"),
    ("side_rpe1_alone", "RPE1 feature alone (oriented, resid. on C)", "side"),
    ("side_partner_alone", "partner feature alone (oriented, resid. on C)", "side"),
]
DIAG_GCOL = {"primary": OI["black"], "copy number": OI["blue"], "accessibility source": OI["vermillion"],
             "location": OI["green"], "chromosomes": OI["purple"], "p53 / expression": OI["orange"],
             "side": OI["sky"]}


def diag_panel(fig, gs, rows, colour_of, label_size=6.0, title_prefix=""):
    D = jl("results/aim4_exploratory/rpe1_diagnostic.json")["pairs"]
    sub = gs.subgridspec(1, 3, wspace=0.08)
    axes = []
    for k, p in enumerate(RPE1_PAIRS):
        ax = fig.add_subplot(sub[0, k], sharey=axes[0] if axes else None)
        axes.append(ax)
        rr = []
        for key, lab, *grp in rows:
            r = D[p][key]
            rr.append(dict(label=lab, est=r["theta_hat"], ci=r["ci95"], color=colour_of(key, grp),
                           filled=r["p_two_sided"] <= 0.05, ms=3.6))
        forest(ax, rr)
        ax.axvline(D[p]["primary"]["theta_hat"], color="k", lw=0.6, ls=":")
        ax.set_title(title_prefix + pl(p, dagger=False) + f" ({ROLE[p]})", fontsize=6.5)
        ax.set_xlabel(f"{TH} (95% CI)")
        if k:
            plt.setp(ax.get_yticklabels(), visible=False)
        else:
            ax.tick_params(axis="y", labelsize=label_size)
    return axes


def fig4():
    slc = jl("results/aim4/same_line_control.json")
    diag = jl("results/aim4_exploratory/rpe1_diagnostic.json")
    link = diag["same_line_link"]["tests"]
    link_ci = {"RPE1_vs_K562": "x_RPE1_vs_K562 (CPM; reproduces same_line_control)",
               "RPE1_vs_HepG2": "x_RPE1_vs_HepG2 (CPM; reproduces same_line_control)"}
    sens = jl("results/aim4/sensitivity_h1.json")
    D = diag["pairs"]
    L, R = 0.25, 0.9

    fig = plt.figure(figsize=(W2, 9.0))
    g = fig.add_gridspec(3, 1, height_ratios=[0.62, 0.78, 1.55], hspace=0.42, left=L, right=R, top=0.955,
                         bottom=0.075)

    # A: same-line control (+ exploratory single-line partial r)
    top = g[0].subgridspec(1, 2, width_ratios=[1.0, 0.45], wspace=0.55)
    axA = fig.add_subplot(top[0, 0])
    rows = []
    for t in PAIRS:
        r = slc[t]
        ci = None
        if t in link_ci:
            lk = link[link_ci[t]]
            assert np.isclose(lk["theta_hat"], r["theta_hat"])
            ci = lk["ci95"]
        rows.append(dict(label=pl(t), est=r["theta_hat"], ci=ci, color=pcol(t), filled=r["p_two_sided"] <= 0.05,
                         p=r["p_two_sided"]))
    forest(axA, rows, xlabel=f"{TH} against same-line M")
    for i, r in enumerate(rows):
        axA.annotate(f"p = {r['p']:.3g}", xy=(1.02, len(rows) - 1 - i), xycoords=("axes fraction", "data"),
                     va="center", fontsize=5.8, annotation_clip=False)
    ptitle(axA, f"Same-line control: M = K562 essential {MINUS} K562 genome-wide ({diag['same_line_link']['G']:,} "
           "genes); filled: two-sided p ≤ 0.05")
    letter(axA, "A", dx=-128, dy=8)
    axAi = fig.add_subplot(top[0, 1])
    single = [("K562 CPM alone", "K562"), ("HepG2 CPM alone", "HepG2"), ("Jurkat CPM alone", "Jurkat"),
              ("RPE1 CPM alone", "RPE1")]
    vals = [link[k]["partial_corr"] for k, _ in single]
    axAi.barh(range(len(single))[::-1], vals, color="#777777", height=0.6)
    axAi.set_yticks(range(len(single))[::-1], [f"{s} CPM" for _, s in single])
    axAi.set_xlabel("partial r with same-line M")
    axAi.axvline(0, color="k", lw=0.5)
    axAi.set_title("single-line features\nEXPLORATORY (post hoc)", fontsize=5.8, color=OI["vermillion"])
    axAi.tick_params(axis="y", length=0)
    axAi.spines["left"].set_visible(False)

    # B: signal-track sensitivity
    axB = fig.add_subplot(g[1])
    mk = {"cpm": ("o", "CPM (primary)"), "fc": ("s", "FC, all lines"), "hyb": ("^", "only RPE1 = 4DN FC")}
    rowsB, groups = [], []

    def add(t, r, kind):
        rowsB.append(dict(label=f"{pl(t)}: {mk[kind][1]}", est=r["theta_hat"], ci=r["ci95"], color=pcol(t),
                          marker=mk[kind][0], p=r["p_two_sided"], filled=r["p_two_sided"] <= 0.05))
    for t in ("K562_vs_Jurkat", "RPE1_vs_HepG2"):
        add(t, sens[t]["main (v1.3, x4, CPM)"], "cpm")
        add(t, sens[t]["FC features (combined_fc)"], "fc")
        if "RPE1" in t:
            add(t, D[t]["rpe1_4dn_fc_hybrid"], "hyb")
        groups.append(len(rowsB))
    for t in ("RPE1_vs_K562", "Jurkat_vs_RPE1"):
        add(t, D[t]["primary"], "cpm")
        add(t, D[t]["rpe1_4dn_fc_hybrid"], "hyb")
        groups.append(len(rowsB))
    forest(axB, rowsB, xlabel=f"{TH} (95% CI)")
    for i, r in enumerate(rowsB):
        axB.annotate(f"p = {r['p']:.3g}", xy=(1.02, len(rowsB) - 1 - i), xycoords=("axes fraction", "data"),
                     va="center", fontsize=5.8, annotation_clip=False)
    for gi in groups[:-1]:
        axB.axhline(len(rowsB) - gi - 0.5, color="0.85", lw=0.5)
    ptitle(axB, "Signal track: Tn5 insertion CPM (primary) vs MACS2 fold change (FC); filled: two-sided p ≤ 0.05")
    letter(axB, "B", dx=-128, dy=8)

    # C: condensed RPE1 diagnostic
    def colour_of(key, _g):
        if key == "primary":
            return OI["black"]
        if key.startswith("side"):
            return OI["sky"]
        return OI["vermillion"] if ("fc" in key or "cpm" in key) else "#555555"
    axes = diag_panel(fig, g[2], DIAG_ROWS, colour_of)
    letter(axes[0], "C", dx=-128, dy=22)
    ptitle(axes[0], "EXPLORATORY (post hoc): RPE1 diagnostic, selected variants (manuscript Table 6)", dy=24,
           color=OI["vermillion"], weight="bold")
    fig.legend(handles=[Line2D([], [], color="k", marker="o", ls="", ms=4, label="two-sided p ≤ 0.05"),
                        Line2D([], [], color="k", marker="o", ls="", ms=4, mfc="white", label="p > 0.05"),
                        Line2D([], [], color="k", ls=":", lw=0.8, label=f"primary {TH}"),
                        Line2D([], [], color="#555555", lw=1.2, label="CN / depth / chromosomes / location / expression"),
                        Line2D([], [], color=OI["vermillion"], lw=1.2, label="signal-track variants"),
                        Line2D([], [], color=OI["sky"], lw=1.2, label="one side only")],
               loc="lower center", ncol=6, fontsize=5.6, bbox_to_anchor=(0.55, 0.0))
    save(fig, "Fig4")


# ====================================================================================== Fig S1
def md_table(text, heading):
    """Parse the first markdown table after the line that starts with `heading`."""
    lines = text.splitlines()
    i = next(k for k, l in enumerate(lines) if l.startswith(heading))
    while not lines[i].startswith("|"):
        i += 1
    rows = []
    while i < len(lines) and lines[i].startswith("|"):
        rows.append([c.strip() for c in lines[i].strip().strip("|").split("|")])
        i += 1
    return pd.DataFrame(rows[2:], columns=rows[0])


def figS1():
    txt = (RES / "features/sanity.md").read_text()
    fig = plt.figure(figsize=(W2, 6.6))
    g = fig.add_gridspec(2, 2, height_ratios=[1.15, 1], hspace=0.42, wspace=0.3, left=0.08, right=0.97,
                         top=0.95, bottom=0.08)
    for k, sig in enumerate(("cpm", "fc")):
        t = md_table(txt, f"## Promoter feature correlation across lines ({sig};")
        lines = list(t.columns[1:])
        M = t.set_index(t.columns[0])[lines].astype(float).loc[lines]
        ax = fig.add_subplot(g[0, k])
        im = ax.imshow(M.values, vmin=0.5, vmax=1, cmap="viridis")
        ax.set_xticks(range(len(lines)), lines, rotation=60, ha="right", fontsize=5.6)
        ax.set_yticks(range(len(lines)), lines, fontsize=5.6)
        for i in range(len(lines)):
            for j in range(len(lines)):
                ax.text(j, i, f"{M.values[i, j]:.2f}", ha="center", va="center", fontsize=4.2,
                        color="white" if M.values[i, j] < 0.8 else "black")
        off = M.values[~np.eye(len(lines), dtype=bool)]
        ax.set_title(("Insertion CPM (primary)" if sig == "cpm" else "MACS2 fold change (sensitivity)")
                     + f": promoter r {off.min():.2f}–{off.max():.2f}", fontsize=6.5)
        for s in ax.spines.values():
            s.set_visible(False)
        letter(ax, "AB"[k], dx=-30)
    cb = fig.colorbar(im, ax=fig.axes[:2], shrink=0.75, pad=0.02)
    cb.set_label("Pearson r, promoter feature (TSS ± 1 kb)")

    # C: promoter accessibility vs baseline expression
    axC = fig.add_subplot(g[1, 0])
    x = np.arange(2)
    for k, (sig, col) in enumerate((("cpm", "#444444"), ("fc", OI["vermillion"]))):
        t = md_table(txt, f"## Promoter feature vs baseline expression ({sig};").set_index("line")
        v = t.loc[["RPE1", "K562"], "rho promoter"].astype(float).values
        axC.bar(x + (k - 0.5) * 0.36, v, 0.36, color=col, label="CPM (primary)" if sig == "cpm" else "FC (sensitivity)")
        for xi, vi in zip(x + (k - 0.5) * 0.36, v):
            axC.text(xi, vi + (0.006 if vi > 0 else -0.006), f"{vi:.3f}".replace("-", MINUS), ha="center", va="bottom" if vi > 0 else "top",
                     fontsize=5.6)
        ngen = t.loc["RPE1", "n genes"]
    axC.axhline(0, color="k", lw=0.6)
    axC.set_xticks(x, ["RPE1", "K562"])
    axC.set_ylabel("Spearman ρ, promoter feature vs\nbaseline expression")
    axC.set_title(f"Promoter accessibility vs expression ({ngen} Aim 1 genes)", fontsize=6.5)
    axC.legend(loc="upper right")
    axC.set_ylim(-0.2, 0.2)
    letter(axC, "C", dx=-30)

    # D: housekeeping percentiles
    axD = fig.add_subplot(g[1, 1])
    rng = np.random.default_rng(0)
    for k, (sig, col) in enumerate((("cpm", "#444444"), ("fc", OI["vermillion"]))):
        t = md_table(txt, f"## Housekeeping promoters ({sig};").set_index("gene")
        lines = [c for c in t.columns if c not in ("lines >= 90th", "pass (all)")]
        for i, gene in enumerate(t.index):
            v = t.loc[gene, lines].astype(float).values
            xx = i + (k - 0.5) * 0.36 + rng.uniform(-0.08, 0.08, len(v))
            axD.scatter(xx, v, s=7, color=col, lw=0, alpha=0.85,
                        label=("CPM (primary)" if sig == "cpm" else "FC (sensitivity)") if i == 0 else None)
    axD.axhline(90, color="k", lw=0.6, ls="--")
    axD.set_xticks(range(len(t.index)), list(t.index))
    axD.set_ylabel(f"within-line percentile of promoter feature\n(one dot per line, n = {len(lines)})")
    axD.set_ylim(40, 101)
    axD.set_title("Housekeeping promoters (expected ≥ 90th percentile)", fontsize=6.5)
    axD.legend(loc="lower left")
    letter(axD, "D", dx=-30)
    save(fig, "FigS1")


# ====================================================================================== Fig S2
VAL_TSV = RES / "atac_local/validation_scatter.tsv"


def compute_validation_table():
    """Per-gene K562 promoter signal (TSS +/- 1 kb), same windows and gene filter as
    scripts/atac_validate_encode.py: ENCODE FC bigWig (ENCFF102ARJ), local-pipeline FC bigWig, and the
    pooled insertion CPM from data/atac_ins/K562__rep*.npz (the primary feature's signal)."""
    from cdtest import atac, features
    val = jl("results/atac_local/k562_ENCSR868FGK.validation_vs_encode.json")
    enc_bw = ROOT / val["encode_bw"]
    loc_bw = ROOT / "data/atac_local/k562_ENCSR868FGK/k562_ENCSR868FGK.fc.bigWig"
    g = pd.read_csv(ROOT / "data/annot/genes_gencode_v50.tsv", sep="\t",
                    usecols=["gene_id", "gene_name", "gene_type", "chrom", "tss0", "blacklist_overlap_tss1kb"])
    g = g[~g["blacklist_overlap_tss1kb"].astype(bool)].reset_index(drop=True)
    s, e = g["tss0"] - 1000, g["tss0"] + 1000
    print("S2: ENCODE FC ...", flush=True)
    g["encode_fc"] = atac.bigwig_mean(enc_bw, g["chrom"], s, e)
    print("S2: local FC ...", flush=True)
    g["local_fc"] = atac.bigwig_mean(loc_bw, g["chrom"], s, e)
    print("S2: insertion CPM ...", flush=True)
    sub = g[g.chrom.isin(features.CHROMS)]
    iv = (sub.chrom.values, (sub.tss0 - 1000).clip(lower=0).values, (sub.tss0 + 1000).values)
    m = features.cpm_matrices({"K562": sorted((ROOT / "data/atac_ins").glob("K562__rep*.npz"))}, [iv], procs=3)
    g["k562_insertion_cpm"] = np.nan
    g.loc[sub.index, "k562_insertion_cpm"] = m["cpm"][0]["K562"].values
    g = g.drop(columns="blacklist_overlap_tss1kb")
    # reproduce the stored validation metric before writing
    k = np.isfinite(g.encode_fc) & np.isfinite(g.local_fc)
    rho = stats.spearmanr(np.log2(1 + g.encode_fc[k]), np.log2(1 + g.local_fc[k]))[0]
    assert int(k.sum()) == val["tss1kb_all_genes"]["n"] and abs(rho - val["tss1kb_all_genes"]["spearman"]) < 1e-9, \
        (int(k.sum()), rho)
    VAL_TSV.parent.mkdir(parents=True, exist_ok=True)
    with open(VAL_TSV, "w") as fh:
        fh.write("# Per-gene K562 promoter signal, TSS +/- 1 kb, genes_gencode_v50 minus blacklisted TSS windows.\n"
                 "# encode_fc: ENCODE FC bigWig ENCFF102ARJ (mean0); local_fc: local pipeline FC bigWig "
                 "(data/atac_local/k562_ENCSR868FGK); k562_insertion_cpm: pooled Tn5 insertion CPM from "
                 "data/atac_ins/K562__rep*.npz (chr1-22,X only). Written by scripts/make_figures.py (Fig. S2).\n")
        g.to_csv(fh, sep="\t", index=False, float_format="%.6g")
    print(f"wrote {VAL_TSV.relative_to(ROOT)}", flush=True)


def figS2():
    if not VAL_TSV.exists():
        compute_validation_table()
    d = pd.read_csv(VAL_TSV, sep="\t", comment="#")
    val = jl("results/atac_local/k562_ENCSR868FGK.validation_vs_encode.json")
    fig, axes = plt.subplots(1, 2, figsize=(W2, 2.9))
    fig.subplots_adjust(left=0.07, right=0.88, bottom=0.15, top=0.84, wspace=0.3)
    for ax, (col, ylab, L) in zip(axes, (("local_fc", "local pipeline, log2(1 + FC)", "A"),
                                         ("k562_insertion_cpm", "insertion CPM (data/atac_ins), log2(1 + CPM)", "B"))):
        k = np.isfinite(d.encode_fc) & np.isfinite(d[col])
        xv, yv = np.log2(1 + d.encode_fc[k]), np.log2(1 + d[col][k])
        hb = ax.hexbin(xv, yv, gridsize=70, bins="log", cmap="Greys", mincnt=1, lw=0)
        rs, rp = stats.spearmanr(xv, yv)[0], stats.pearsonr(xv, yv)[0]
        ax.text(0.03, 0.97, f"n = {k.sum():,} genes\nSpearman {rs:.3f}\nPearson {rp:.3f}", transform=ax.transAxes,
                va="top", fontsize=6)
        if col == "local_fc":
            lim = [0, max(xv.max(), yv.max()) * 1.02]
            ax.plot(lim, lim, color=OI["vermillion"], lw=0.6, ls="--")
            # matches the stored validation metric (TSV keeps 6 significant digits, so allow rounding)
            assert abs(rs - val["tss1kb_all_genes"]["spearman"]) < 1e-3, rs
            ax.set_title("Local pipeline vs ENCODE, same reads (ENCSR868FGK)", fontsize=6)
        else:
            ax.set_title("Primary feature signal (insertion CPM) vs ENCODE FC", fontsize=6)
        ax.set_xlabel("ENCODE K562 (ENCFF102ARJ), log2(1 + FC)")
        ax.set_ylabel(ylab)
        letter(ax, L, dx=-28, dy=12)
    cb = fig.colorbar(hb, cax=fig.add_axes([0.9, 0.2, 0.012, 0.6]))
    cb.set_label("genes per hexagon")
    fig.suptitle("K562 promoter signal, TSS ± 1 kb", fontsize=7, x=0.45)
    save(fig, "FigS2")


# ====================================================================================== Fig S3
def figS3():
    from cdtest.k_rule_sim import DESIGNS, FPR_MAX05, REFS, RULES
    pc = pd.read_csv(RES / "aim3/k_rule_selection_by_cell.csv")
    pr = pd.read_csv(RES / "aim3/k_rule_selection.csv").set_index("rule")
    order = list(RULES) + list(REFS)
    chosen = "x4"
    cmap = [OI["vermillion"], OI["orange"], OI["yellow"], OI["sky"], OI["black"], OI["blue"], OI["purple"],
            "#8C564B", OI["green"]]
    cols = dict(zip(order, cmap))
    fig, ax = plt.subplots(1, 3, figsize=(W2, 2.6), gridspec_kw=dict(width_ratios=[1.5, 1, 1]))
    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.17, top=0.88, wspace=0.38)
    designs = list(DESIGNS)
    xs = np.arange(len(designs))
    w = 0.8 / len(order)
    for i, rule in enumerate(order):
        gg = pc[pc.rule == rule].groupby("design").fpr05
        ax[0].bar(xs + i * w, gg.max().reindex(designs), w, color=cols[rule],
                  label=rule + (" (adopted)" if rule == chosen else ""))
        ax[0].scatter(xs + i * w, gg.mean().reindex(designs), color="k", s=3, zorder=3, lw=0)
    ax[0].axhline(0.05, color="k", lw=0.6)
    ax[0].axhline(FPR_MAX05, color=OI["vermillion"], ls="--", lw=0.7)
    ax[0].set_xticks(xs + 0.4 - w / 2, designs)
    ax[0].set_ylabel("false-positive rate, α = 0.05\n(bar: worst cell; dot: mean)")
    ncell = int(pr.n_cells.iloc[0])
    ax[0].set_title(f"Size by module design ({ncell} null cells, {int(pc.n_null.iloc[0]):,} replicates each)", fontsize=6.3)
    ax[0].legend(ncol=3, fontsize=5.4, loc="upper right")
    ax[0].set_ylim(0, pc.fpr05.max() * 1.45)
    letter(ax[0], "A", dx=-30)
    for rule in order:
        r = pr.loc[rule]
        ax[1].errorbar(r.fpr05_worst, r.power_matched_small, xerr=1.96 * r.fpr05_worst_se,
                       yerr=1.96 * r.power_matched_small_se, fmt="o", color=cols[rule], ms=5 if rule == chosen else 3.5,
                       elinewidth=0.6)
        ax[1].annotate(rule, (r.fpr05_worst, r.power_matched_small), fontsize=5.5,
                       xytext=(-3, 2) if rule == "x3" else (3, 2), ha="right" if rule == "x3" else "left",
                       textcoords="offset points")
    ax[1].axvline(FPR_MAX05, color=OI["vermillion"], ls="--", lw=0.7)
    ax[1].set_xlabel("worst-cell false-positive rate (α = 0.05)")
    ax[1].set_ylabel("matched power, partial R² 0.5%")
    ax[1].set_title("Eligibility: left of the red line", fontsize=6.3)
    letter(ax[1], "B", dx=-30)
    for rule in RULES:
        km = [pr.loc[rule, f"k_median_{d}"] for d in designs]
        ax[2].plot(xs, km, marker="o", ms=2.8, color=cols[rule], label=rule, lw=1.4 if rule == chosen else 0.8)
    ax[2].set_xticks(xs, designs)
    ax[2].set_ylabel("median block size k chosen")
    ax[2].set_title("Block size chosen", fontsize=6.3)
    letter(ax[2], "C", dx=-30)
    save(fig, "FigS3")


# ====================================================================================== Fig S4
def figS4():
    fig, axes = plt.subplots(2, 4, figsize=(W2, 3.9))
    fig.subplots_adjust(left=0.08, right=0.99, bottom=0.1, top=0.82, hspace=1.15, wspace=0.3)
    for row, (suf, name) in enumerate((("", "plan v1.2 (no arm terms)"), ("_v13", "plan v1.3 (arm terms in f(C))"))):
        res = jl(f"results/aim2/RPE1_vs_K562_test{suf}.json")
        nulls = np.load(RES / f"aim2/RPE1_vs_K562_null{suf}.npz")
        for k, lab in enumerate(res["features"]):
            ax = axes[row, k]
            d = res["features"][lab]
            m = d["main_block_x4"]
            Ts, T = nulls[lab], d["theta_hat"]
            lim = max(np.abs(Ts).max(), abs(T)) * 1.1
            ax.hist(Ts, bins=70, range=(-lim, lim), color=OI["sky"], lw=0)
            ax.axvline(T, color=OI["vermillion"], lw=1.4)
            ax.axvline(0, color="k", lw=0.5, ls=":")
            ex = d["exchangeability_x4"]
            kind, sig = lab.split("_")
            ax.set_title(f"{kind} ({sig.upper()})\n{TH} = {T:.4f}; two-sided p = {m['p_two_sided']:.4f}\n"
                         f"k = {d['k_x4']}; null mean {m['null_mean']:.4f}\n{ex['frac_genes_movable']:.0%} of genes "
                         "exchangeable", fontsize=5.6)
            ax.set_xlabel("T* (block-permutation null)")
            if k == 0:
                ax.set_ylabel(f"{name}\npermutations (B = {m['B']:,})", fontsize=6)
        letter(axes[row, 0], "AB"[row], dx=-30)
    fig.suptitle(f"Development pair RPE1 {MINUS} K562 (mixed; descriptive): observed {TH} (red) against the "
                 "block-permutation null", fontsize=7, y=0.99)
    save(fig, "FigS4")


# ====================================================================================== Fig S5
def figS5():
    fig, ax = plt.subplots(1, 2, figsize=(W2 * 0.7, 2.5))
    fig.subplots_adjust(left=0.1, right=0.98, bottom=0.35, top=0.82, wspace=0.35)
    tabs = planted_from_csv(RES / "aim3/planted_RPE1_vs_K562.csv")
    planted_axes(ax[0], ax[1], tabs, "")
    fig.suptitle(f"Planted signal in real M, RPE1 {MINUS} K562, plan v1.2 (no arm terms)", fontsize=7)
    h, l = ax[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=4, fontsize=5.6)
    letter(ax[0], "A", dx=-28); letter(ax[1], "B", dx=-28)
    save(fig, "FigS5")


# ====================================================================================== Fig S6
def figS6():
    h2 = jl("results/aim4/h2_location.json")
    fig, axes = plt.subplots(1, 2, figsize=(W2 * 0.85, 2.6), sharey=True)
    fig.subplots_adjust(left=0.17, right=0.98, bottom=0.17, top=0.88, wspace=0.08)
    n = len(PAIRS)
    for ax, key, xl in ((axes[0], "per_sd", f"{TH} × SD(x), 95% CI (main test)"),
                        (axes[1], "pr2", "partial R² / reliability (%), bootstrap 95%")):
        for i, t in enumerate(PAIRS):
            y0 = n - 1 - i
            for j, (f, col, lab) in enumerate((("promoter_cpm", C_PROM, "promoter"), ("distal_cpm", C_DIST, "distal"))):
                o = h2[t][f]
                y = y0 + (0.5 - j) * 0.3
                if key == "per_sd":
                    ax.plot(o["ci95_per_sd"], [y, y], color=col, lw=1.1)
                    ax.plot(o["theta_per_sd"], y, "o", color=col, ms=3.5, label=lab if i == 0 else None,
                            mfc=col if o["p_two_sided"] <= 0.05 else "white")
                else:
                    ax.plot(np.array(o["pr2_rel_boot_ci95"]) * 100, [y, y], color=col, lw=1.1)
                    ax.plot(o["pr2_rel"] * 100, y, "o", color=col, ms=3.5, label=lab if i == 0 else None)
        ax.set_xlabel(xl)
        if key == "per_sd":
            ax.axvline(0, color="k", lw=0.6, ls="--")
        ax.set_yticks(range(n)[::-1], [pl(t) + (" *" if t in MATCHED else "") for t in PAIRS])
        for tl, t in zip(ax.get_yticklabels(), PAIRS):
            tl.set_color(pcol(t))
        ax.tick_params(axis="y", length=0)
    axes[0].legend(handles=[Line2D([], [], color=C_PROM, marker="o", ms=3.5, label="promoter"),
                            Line2D([], [], color=C_DIST, marker="o", ms=3.5, label="distal")], loc="lower left")
    axes[0].text(0, 1.02, "* matched pair; † contains RPE1; filled: two-sided p ≤ 0.05", transform=axes[0].transAxes,
                 fontsize=5.6)
    letter(axes[0], "A", dx=-75); letter(axes[1], "B", dx=-10)
    save(fig, "FigS6")


# ====================================================================================== Fig S7
def figS7():
    h3 = jl("results/aim4/h3_comparators.json")
    fig, axes = plt.subplots(1, 2, figsize=(W2, 2.7))
    fig.subplots_adjust(left=0.08, right=0.82, bottom=0.27, top=0.88, wspace=0.3)
    x = np.arange(len(PAIRS))
    axes[0].bar(x - 0.2, [h3[t]["theta_main"] for t in PAIRS], 0.4, color=[pcol(t) for t in PAIRS],
                label="main test (f(C) + arm terms)")
    axes[0].bar(x + 0.2, [h3[t]["naive_step7"]["coef"] for t in PAIRS], 0.4, color="white",
                edgecolor=[pcol(t) for t in PAIRS], hatch="/////", lw=0.8, label="naive step 9 (1 + expression difference)")
    axes[0].axhline(0, color="k", lw=0.6)
    axes[0].set_ylabel(f"{TH}, combined CPM")
    axes[0].legend(loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2, fontsize=5.6)
    comps = [("naive_ols", "OLS, genes independent"), ("unrestricted", "unrestricted permutation"),
             ("expr_bin", "expression-bin permutation"), ("stratified", "stratified permutation"),
             ("camera", "CAMERA-style VIF")]
    for k, lab in comps:
        st = TEST_STYLE[k]
        axes[1].plot(x, [-np.log10(h3[t]["comparators_p_two_sided"][k]) for t in PAIRS], marker=st[3], color=st[1],
                     ms=3, lw=0.8, label=lab)
    axes[1].plot(x, [-np.log10(h3[t]["p_main_two"]) for t in PAIRS], "s-", color="k", ms=3.5, lw=1.4,
                 label="main test (block, x4)")
    axes[1].plot(x, [-np.log10(h3[t]["naive_step7"]["p_two_sided"]) for t in PAIRS], "^--", color="#777777", ms=3,
                 lw=0.8, label="naive step 9")
    axes[1].axhline(-np.log10(0.05), color="grey", ls=":", lw=0.7)
    axes[1].set_ylabel("−log10 two-sided p")
    axes[1].legend(loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=5.6)
    for ax in axes:
        ax.set_xticks(x, [pl(t) for t in PAIRS], rotation=35, ha="right")
        for tl, t in zip(ax.get_xticklabels(), PAIRS):
            tl.set_color(pcol(t))
    letter(axes[0], "A", dx=-30); letter(axes[1], "B", dx=-30)
    save(fig, "FigS7")


# ====================================================================================== Fig S8
def figS8():
    H4 = pd.read_csv(RES / "aim4/extras/h4_similarity_gradient.csv")
    F = jl("results/aim4/h4_gradient.json")["fits"]
    fig, axes = plt.subplots(2, 3, figsize=(W2, 4.4), sharex=True)
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.1, top=0.86, hspace=0.45, wspace=0.3)
    mk = {"own": ("*", "k", 50, "own x"), "other_pair": ("D", OI["purple"], 12, "other pairs' x"),
          "decoy": ("o", "#888888", 9, "decoy differences")}
    for k, (ax, t) in enumerate(zip(axes.ravel(), PAIRS)):
        d = H4[H4.pair == t]
        ax.axvspan(-0.05, 0.05, color=OI["yellow"], alpha=0.35, lw=0)
        for grp, (m, c, s, lab) in mk.items():
            q = d[d.group == grp]
            ax.scatter(q.similarity, q.theta_hat, marker=m, c=c, s=s, label=lab, lw=0, zorder=3)
        f = F[t]["all"]
        xx = np.linspace(-1, 1, 10)
        ax.plot(xx, f["intercept"] + f["slope"] * xx, color=pcol(t), lw=1.1)
        ax.axhline(0, color="grey", lw=0.5); ax.axvline(0, color="grey", lw=0.5)
        ax.set_title(f"{pl(t)} ({ROLE[t]})\nslope {f['slope']:.4f} (p = {f['slope_p']:.2g}); n = {f['n']}",
                     fontsize=6, color=pcol(t))
        if k >= 3:
            ax.set_xlabel("similarity to own x (Pearson, C-residualised)")
        if k % 3 == 0:
            ax.set_ylabel(TH)
        letter(ax, "ABCDEF"[k], dx=-26)
    axes[0, 0].legend(loc="lower right", fontsize=5.4)
    p = F["pooled_all_pairs"]["all"]
    fig.suptitle(f"H4 similarity gradient (combined CPM; yellow band: decoy negative-control region |similarity| < 0.05)\n"
                 f"pooled n = {p['n']} (descriptive): slope {p['slope']:.4f} (95% CI {p['slope_ci95'][0]:.4f} to "
                 f"{p['slope_ci95'][1]:.4f}), intercept {p['intercept']:.4f} (p = {p['intercept_p']:.3f})", fontsize=6.5)
    save(fig, "FigS8")


# ====================================================================================== Fig S9
def figS9():
    P = jl("results/aim4/p53_control.json")
    fig, ax = plt.subplots(figsize=(W1 * 1.1, 2.6))
    fig.subplots_adjust(left=0.17, right=0.98, bottom=0.3, top=0.88)
    x = np.arange(len(PAIRS))
    z = [P[t]["z_oriented"] for t in PAIRS]
    cols = [C_MIXED if P[t]["p53_differs"] else C_MATCHED for t in PAIRS]
    ax.bar(x, z, color=cols, width=0.65)
    for xi, t in zip(x, PAIRS):
        ax.text(xi, P[t]["z_oriented"] + (0.15 if P[t]["z_oriented"] > 0 else -0.15),
                f"{P[t]['z_oriented']:.2f}".replace("-", MINUS) + f"\n(n = {P[t]['n_p53_targets']})", ha="center",
                va="bottom" if P[t]["z_oriented"] > 0 else "top", fontsize=5.4)
    ax.axhline(0, color="k", lw=0.6)
    ax.axhline(1.645, color="grey", ls=":", lw=0.7)
    ax.axhline(-1.96, color="grey", ls=":", lw=0.7)
    ax.set_xticks(x, [pl(t) for t in PAIRS], rotation=35, ha="right")
    ax.set_ylabel("z, mean M at Fischer core p53 targets\n(+ = larger response in p53-WT line)")
    ax.set_ylim(min(min(z) - 1.2, -2.3), max(z) + 1.5)
    ax.legend(handles=[Patch(color=C_MIXED, label="pair differs in p53 status"),
                       Patch(color=C_MATCHED, label="matched pair (A − B)")], loc="upper left", fontsize=5.6)
    ax.set_title(f"p53 control ({P[PAIRS[0]]['n_sets']:,} expression-matched random sets)", fontsize=6.5)
    save(fig, "FigS9")


# ===================================================================================== Fig S10
JURKAT_QC = RES / "atac_local/jurkat_GSE302221.qc.json"


def figS10():
    """Jurkat CRISPRi-background ATAC sensitivity (GSE302221). Skipped until the QC file exists.
    A: QC gate metrics per replicate, GSE273644 (primary) vs GSE302221 (CRISPRi background).
    B: promoter CPM agreement (results/atac_local/jurkat_GSE302221.concordance.json), if present.
    C: K562 - Jurkat H1 estimate with the alternative Jurkat features, if a result file
       results/atac_local/jurkat_GSE302221.h1_sensitivity.json (keys 'primary'/'alternative', each with
       theta_hat, ci95, p_two_sided) exists."""
    if not JURKAT_QC.exists():
        print("FigS10 skipped: results/atac_local/jurkat_GSE302221.qc.json absent", flush=True)
        return False
    qc = {"GSE273644 (primary)": jl("results/atac_local/jurkat_GSE273644.qc.json"),
          "GSE302221 (CRISPRi bg)": json.loads(JURKAT_QC.read_text())}
    conc_p = RES / "atac_local/jurkat_GSE302221.concordance.json"
    h1_p = RES / "atac_local/jurkat_GSE302221.h1_sensitivity.json"
    fig, axes = plt.subplots(1, 3, figsize=(W2, 2.6))
    fig.subplots_adjust(left=0.07, right=0.98, bottom=0.25, top=0.86, wspace=0.45)
    metrics = ["tss_enrichment", "NRF", "PBC1"]
    sub = axes[0].inset_axes([0, 0, 1, 1]); axes[0].axis("off")
    xs = np.arange(len(metrics))
    for k, (name, q) in enumerate(qc.items()):
        col = OI["blue"] if k == 0 else OI["vermillion"]
        for i, m in enumerate(metrics):
            reps = q["gate"].get(m, {})
            v = [vv[0] for vv in reps.values()] if isinstance(reps, dict) else []
            thr = q.get("gate_thresholds", {}).get(f"{m.lower()}_min" if m != "tss_enrichment" else "tss_enrichment_min")
            vn = np.array(v, float) / (thr if thr else 1)
            sub.scatter(np.full(len(vn), i + (k - 0.5) * 0.3), vn, color=col, s=10, label=name if i == 0 else None)
    sub.axhline(1, color="k", lw=0.6, ls="--")
    sub.set_xticks(xs, ["TSS enrichment", "NRF", "PBC1"])
    sub.set_ylabel("per-replicate value / gate threshold")
    sub.legend(fontsize=5.4, loc="upper right")
    letter(sub, "A", dx=-28)
    if conc_p.exists():
        c = json.loads(conc_p.read_text())
        items = [("primary vs CRISPRi\n(all, no BL)", c["primary_vs_crispri"]["all_genes_no_blacklist"]["spearman"])]
        rs = c["context"]["replicate_spearman_summary"]
        items += [("within primary reps", rs["within_primary_mean"]), ("within CRISPRi reps", rs["within_crispri_mean"]),
                  ("across sets (reps)", rs["across_sets_mean"]),
                  ("Jurkat vs K562", c["context"]["Jurkat_vs_K562"]["spearman"]),
                  ("Jurkat vs GM12878", c["context"]["Jurkat_vs_GM12878"]["spearman"])]
        axes[1].barh(range(len(items))[::-1], [v for _, v in items], color="#777777")
        axes[1].set_yticks(range(len(items))[::-1], [k for k, _ in items], fontsize=5.4)
        axes[1].set_xlabel("Spearman ρ, promoter insertion CPM")
        axes[1].set_xlim(min(v for _, v in items) - 0.05, 1)
    else:
        axes[1].text(0.5, 0.5, "concordance file absent", ha="center", transform=axes[1].transAxes)
        axes[1].axis("off")
    letter(axes[1], "B", dx=-60)
    if h1_p.exists():
        h = json.loads(h1_p.read_text())
        rows = [dict(label=f"Jurkat = {k}", est=h[k]["theta_hat"], ci=h[k]["ci95"], color=C_MATCHED,
                     filled=h[k]["p_two_sided"] <= 0.05) for k in ("primary", "alternative") if k in h]
        forest(axes[2], rows, xlabel=f"{TH}, K562 {MINUS} Jurkat (combined CPM)")
    else:
        axes[2].text(0.5, 0.5, "H1 sensitivity result absent", ha="center", transform=axes[2].transAxes)
        axes[2].axis("off")
    letter(axes[2], "C", dx=-60)
    save(fig, "FigS10")
    return True


# ===================================================================================== Fig S11
def figS11():
    fig = plt.figure(figsize=(W2, 6.8))
    g = fig.add_gridspec(1, 1, left=0.3, right=0.985, top=0.88, bottom=0.07)
    axes = diag_panel(fig, g[0], DIAG_FULL, lambda key, grp: DIAG_GCOL[grp[0]], label_size=5.6)
    D = jl("results/aim4_exploratory/rpe1_diagnostic.json")
    for ax, p in zip(axes, RPE1_PAIRS):
        lj = D["loco_jackknife"][p]
        ax.text(0.02, -0.075, f"LOCO: {lj['n_chrom_positive']}/{lj['n']} drops > 0; jackknife SE {lj['jackknife_se']:.4f}",
                transform=ax.transAxes, fontsize=5.2, color="#444444")
    fig.legend(handles=[Line2D([], [], color=c, lw=1.4, label=k) for k, c in DIAG_GCOL.items()]
               + [Line2D([], [], color="k", marker="o", ls="", ms=3.5, mfc="white", label="two-sided p > 0.05")],
               loc="upper center", ncol=4, fontsize=5.6, bbox_to_anchor=(0.6, 0.962))
    tss = D["tss_enrichment_insertions"]
    others = [v for k, v in tss.items() if k != "RPE1"]
    fig.suptitle("EXPLORATORY (post hoc): RPE1 diagnostic, all variants (filled: two-sided p ≤ 0.05; dotted: primary "
                 f"{TH})\ninsertion TSS enrichment RPE1 {tss['RPE1']:.1f}, other lines {min(others):.1f}–{max(others):.1f}",
                 fontsize=6.8, x=0.6, y=0.995)
    save(fig, "FigS11")


# ===================================================================================== Fig S12
def figS12():
    h5 = jl("results/aim4/h5_tf_cases.json")
    res = pd.DataFrame(h5["cases"] + h5["supplementary_non_lambert"]).sort_values("rho_own", ascending=False)
    rng = np.random.default_rng(h5["seed"])
    fig, ax = plt.subplots(figsize=(W2, 2.9))
    fig.subplots_adjust(left=0.08, right=0.99, bottom=0.2, top=0.86)
    for x, r in enumerate(res.itertuples()):
        ax.scatter(np.full(len(r.rho_swap), x) + rng.uniform(-0.15, 0.15, len(r.rho_swap)), r.rho_swap, s=5,
                   c="#BBBBBB", lw=0, zorder=1)
        c = OI["vermillion"] if r.expected_pattern else (OI["orange"] if r.rho_own > 0 else OI["blue"])
        ax.scatter(x, r.rho_own, s=26, marker="D", c=c, zorder=3, edgecolor="k" if r.chip_source == "both" else "none",
                   lw=0.5)
    ax.axhline(0, c="k", lw=0.6)
    ax.set_xticks(range(len(res)), [f"{t}{'' if s == 'both' else '*'}{'' if lt else '+'}"
                                    for t, s, lt in zip(res.tf, res.chip_source, res.lambert_tf)], rotation=70, fontsize=5.6)
    ax.set_xlim(-0.7, len(res) - 0.3)
    ax.set_ylabel(r"Spearman ρ, $d_g$ vs Δ accessible" + "\nbound sites within 50 kb")
    ax.legend(handles=[Line2D([], [], marker="D", ls="", color=OI["vermillion"], ms=4, label="own ρ > 0, rank 1 of 21"),
                       Line2D([], [], marker="D", ls="", color=OI["orange"], ms=4, label="own ρ > 0"),
                       Line2D([], [], marker="D", ls="", color=OI["blue"], ms=4, label="own ρ ≤ 0"),
                       Line2D([], [], marker="D", ls="", color="w", mec="k", ms=4, label="black edge: ChIP in both lines"),
                       Line2D([], [], marker="o", ls="", color="#BBBBBB", ms=3, label=f"{h5['n_swap']} swap TFs")],
              loc="lower left", ncol=5, fontsize=5.6)
    ax.set_title(f"H5 case studies, K562 {MINUS} HepG2 (mixed pair; illustrative): {h5['n_tfs']} Lambert TFs + "
                 f"{len(h5['supplementary_non_lambert'])} other factors. * ChIP in one line only; + not a Lambert TF",
                 fontsize=6.5)
    save(fig, "FigS12")


FIGS = {"Fig1": fig1, "Fig2": fig2, "Fig3": fig3, "Fig4": fig4, "FigS1": figS1, "FigS2": figS2, "FigS3": figS3,
        "FigS4": figS4, "FigS5": figS5, "FigS6": figS6, "FigS7": figS7, "FigS8": figS8, "FigS9": figS9,
        "FigS10": figS10, "FigS11": figS11, "FigS12": figS12}


def main(argv):
    which = argv or list(FIGS)
    for name in which:
        FIGS[name]()


if __name__ == "__main__":
    main(sys.argv[1:])
