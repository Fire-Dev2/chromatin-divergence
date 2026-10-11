"""Sanity checks for the accessibility features (ANALYSIS_PLAN.md section 3).

Writes results/features/sanity.md and figures/features_promoter_corr.png.
Run after `python -m cdtest.features` (scripts/build_features.sh does both).

Expression check uses ONLY gene_id, mean_ctrlA (RPE1) and mean_ctrlB (K562) of the Aim 1
RPE1_vs_K562 gene_outcome.tsv (no outcome columns are read).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

from cdtest import features as F

ROOT = F.ROOT
OUT_MD = ROOT / "results" / "features" / "sanity.md"
FIG = ROOT / "figures" / "features_promoter_corr.png"
OUTCOME = ROOT / "results" / "aim1" / "RPE1_vs_K562" / "gene_outcome.tsv"

HOUSEKEEPING = ["GAPDH", "ACTB"]
# (gene, line expected higher, lines it must exceed). All-line rank is reported as well.
LINEAGE = [("HBG1", "K562", ["RPE1"]), ("GATA1", "K562", ["RPE1"]),
           ("ALB", "HepG2", ["K562", "RPE1"]), ("APOA1", "HepG2", ["K562", "RPE1"])]


def main():
    info = json.loads(F.BUILD_INFO.read_text())
    gf = F.load_gene_features()
    lines = [str(x) for x in gf.line.cat.categories] if hasattr(gf.line, "cat") else list(dict.fromkeys(gf.line))
    P = F.feature_matrix("promoter", gf)[lines]
    D = F.feature_matrix("distal", gf)[lines]
    C = F.feature_matrix("combined", gf)[lines]
    name = gf.drop_duplicates("gene_id").set_index("gene_id").gene_name
    pct = P.rank(pct=True) * 100
    md = []
    ok_all = []
    w = md.append
    w("# Accessibility feature sanity checks\n")
    w(f"Built {info['built']} from `{info['config']}`. Lines ({info['n_lines']}): "
      + ", ".join(f"{d['line']} ({d['role']})" for d in info["lines"]) + ".")
    if info["skipped_lines"]:
        w("Skipped (inputs not yet available): " + ", ".join(s["line"] for s in info["skipped_lines"]) + ".")
    ri = info["region_info"]
    w(f"\n## Union regions\n\n- Summit windows pooled: {ri['n_windows']:,}; merged: {ri['n_merged']:,}; "
      f"blacklist-overlapping dropped: {ri['n_blacklisted']:,}.\n- **Union regions: {info['n_regions']:,}.**")
    reg = pd.read_parquet(F.REGION_SIGNAL, columns=["start", "end"])
    L = reg.end - reg.start
    w(f"- Region length: median {L.median():.0f} bp, 90th pct {L.quantile(.9):.0f}, max {L.max():,} bp.")
    w(f"- Genes: {info['n_genes']:,}; gene-region distal pairs {info['n_gene_region_pairs']:,}; "
      f"genes with no union region 2.5-100 kb away: {info['n_genes_no_distal_region']:,}.")
    w(f"- Pooled z: promoter mean {info['z_pooled']['promoter_mean']:.3f} SD {info['z_pooled']['promoter_sd']:.3f}; "
      f"distal mean {info['z_pooled']['distal_mean']:.3f} SD {info['z_pooled']['distal_sd']:.3f}.")
    w("\nPer-line summit counts: " + ", ".join(f"{d['line']} {d['n_summits']:,}" for d in info["lines"]) + ".")

    # correlation matrix
    corr = P.corr(method="pearson")
    w("\n## Promoter feature correlation across lines (Pearson, all genes)\n")
    w("| | " + " | ".join(lines) + " |\n|---|" + "---|" * len(lines))
    for a in lines:
        w(f"| {a} | " + " | ".join(f"{corr.loc[a, b]:.2f}" for b in lines) + " |")
    off = corr.values[~np.eye(len(lines), dtype=bool)]
    w(f"\nOff-diagonal range {off.min():.2f}-{off.max():.2f} (median {np.median(off):.2f}).")
    dc = D.corr().values[~np.eye(len(lines), dtype=bool)]
    w(f"Distal feature off-diagonal correlation range {dc.min():.2f}-{dc.max():.2f}. "
      f"Promoter vs distal within line (median over lines): "
      f"{np.median([np.corrcoef(P[l], D[l])[0, 1] for l in lines]):.2f}.")

    fig, ax = plt.subplots(figsize=(7, 6))
    im = ax.imshow(corr.values, vmin=min(0.5, off.min()), vmax=1, cmap="viridis")
    ax.set_xticks(range(len(lines)), lines, rotation=60, ha="right")
    ax.set_yticks(range(len(lines)), lines)
    for i in range(len(lines)):
        for j in range(len(lines)):
            ax.text(j, i, f"{corr.values[i, j]:.2f}", ha="center", va="center", fontsize=7,
                    color="white" if corr.values[i, j] < 0.85 else "black")
    fig.colorbar(im, ax=ax, label="Pearson r")
    ax.set_title(f"Promoter accessibility (TSS±1 kb), {len(P):,} genes")
    fig.tight_layout()
    FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG, dpi=150)
    plt.close(fig)

    # expression
    w("\n## Promoter feature vs baseline expression (expect Spearman > 0.3)\n")
    w("Source: `results/aim1/RPE1_vs_K562/gene_outcome.tsv`, columns gene_id, mean_ctrlA (RPE1), "
      "mean_ctrlB (K562) only.\n")
    ex = pd.read_csv(OUTCOME, sep="\t", usecols=["gene_id", "mean_ctrlA", "mean_ctrlB"]).set_index("gene_id")
    ex = ex[ex.index.isin(P.index)]
    w("| line | n genes | rho promoter | rho distal | rho combined | pass (promoter > 0.3) |\n|---|---|---|---|---|---|")
    for ln, col in (("RPE1", "mean_ctrlA"), ("K562", "mean_ctrlB")):
        if ln not in lines:
            continue
        r = {k: spearmanr(M.loc[ex.index, ln], ex[col])[0] for k, M in (("p", P), ("d", D), ("c", C))}
        ok = r["p"] > 0.3
        ok_all.append(ok)
        w(f"| {ln} | {len(ex):,} | {r['p']:.3f} | {r['d']:.3f} | {r['c']:.3f} | {'PASS' if ok else 'FAIL'} |")
    w("\nSupplementary (uses only gene_id membership): the 5,262 genes in gene_outcome passed the Aim 1 "
      "expression filter in both lines. AUC of the promoter feature, expressed set vs all other genes:")
    for ln in [l for l in ("RPE1", "K562") if l in lines]:
        inn = P.index.isin(ex.index)
        a, b = P.loc[inn, ln].values, P.loc[~inn, ln].values
        auc = mannwhitneyu(a, b).statistic / (len(a) * len(b))
        w(f"- {ln}: AUC {auc:.3f} (median {np.median(a):.2f} vs {np.median(b):.2f}).")
    if {"RPE1", "K562"} <= set(lines):
        xd = P.loc[ex.index, "K562"] - P.loc[ex.index, "RPE1"]
        ed = np.log2(ex.mean_ctrlB) - np.log2(ex.mean_ctrlA)
        w(f"\nBetween-line: Spearman(promoter K562-RPE1, log2 expr K562/RPE1) = {spearmanr(xd, ed)[0]:.3f} "
          "(descriptive only).")

    # housekeeping
    by_name = {v: k for k, v in name.items()}
    w("\n## Housekeeping promoters (percentile of promoter feature within each line; expect high)\n")
    w("| gene | " + " | ".join(lines) + " | pass (all >= 90th) |\n|---|" + "---|" * (len(lines) + 1))
    for g in HOUSEKEEPING:
        gid = by_name[g]
        ok = bool((pct.loc[gid] >= 90).all())
        ok_all.append(ok)
        w(f"| {g} | " + " | ".join(f"{pct.loc[gid, l]:.0f}" for l in lines) + f" | {'PASS' if ok else 'FAIL'} |")

    # lineage
    w("\n## Lineage markers (promoter feature; value [within-line percentile])\n")
    w("| gene | " + " | ".join(lines) + " | check | pass |\n|---|" + "---|" * (len(lines) + 2))
    for g, hi, lo in LINEAGE:
        gid = by_name[g]
        v = P.loc[gid]
        lo = [x for x in lo if x in lines]
        ok = all(v[hi] > v[x] for x in lo)
        rank = int((v > v[hi]).sum()) + 1
        chk = f"{hi} > {', '.join(lo)} (all-line rank {rank}/{len(lines)})"
        ok_all.append(bool(ok))
        w(f"| {g} | " + " | ".join(f"{v[l]:.2f} [{pct.loc[gid, l]:.0f}]" for l in lines)
          + f" | {chk} | {'PASS' if ok else 'FAIL'} |")

    w(f"\n## Summary\n\n**{sum(ok_all)}/{len(ok_all)} checks pass.**")
    w("\nDiagnosis of the expression / housekeeping failures (2026-10-11 build, 10 lines): the same "
      "Spearman computed from the raw fold-change bigWig over TSS +/-250, +/-500 and +/-1000 bp, and the "
      "max over +/-1 kb, is -0.05 to -0.07 (K562) and -0.13 to -0.14 (RPE1), so the failure is not caused "
      "by the window, quantile normalisation or CN step. Promoter accessibility separates expressed from "
      "non-expressed genes strongly (AUC above), but within the expression-filtered set the MACS2 "
      "fold-change track (local-lambda background) does not rank strong promoters higher: very highly "
      "expressed genes (RPs, EEF1A1, GAPDH, ACTB) sit at FC ~1.5-3.5 over +/-1 kb. Any change of the "
      "signal track would be a plan deviation and is left to the PI.")
    w("\nCaveats: RPE1 uses the provisional 4DN fold-change bigWig + conservative peaks until the local "
      "pipeline output exists; GM12878 and IMR-90 have no DepMap CN (CN correction = 0). "
      "Jurkat is absent until its local pipeline finishes; rerun `scripts/build_features.sh`.")
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(md) + "\n")
    print("\n".join(md))
    return 0 if all(ok_all) else 1


if __name__ == "__main__":
    sys.exit(main())
