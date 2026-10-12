"""Sanity checks for the accessibility features (ANALYSIS_PLAN.md section 3, signal per plan v1.2).

Checks BOTH versions: cpm (primary; Tn5 insertion CPM) and fc (sensitivity; MACS2 fold change).
Writes results/features/sanity.md, results/features/replicate_reproducibility.tsv,
results/features/cpm_vs_fc_pairs.tsv, figures/features_promoter_corr.png (cpm) and
figures/features_promoter_corr_fc.png (fc). Run after `python -m cdtest.features --signal both`
(scripts/build_features.sh does both).

Expression check uses ONLY gene_id, mean_ctrlA (RPE1) and mean_ctrlB (K562) of the Aim 1
RPE1_vs_K562 gene_outcome.tsv (no outcome columns are read). The CPM-vs-FC comparison of x uses
features only.
"""
from __future__ import annotations

import itertools
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, pearsonr, spearmanr

from cdtest import features as F

ROOT = F.ROOT
OUT_DIR = ROOT / "results" / "features"
OUT_MD = OUT_DIR / "sanity.md"
REP_TSV = OUT_DIR / "replicate_reproducibility.tsv"
PAIR_TSV = OUT_DIR / "cpm_vs_fc_pairs.tsv"
FIG = {"cpm": ROOT / "figures" / "features_promoter_corr.png",
       "fc": ROOT / "figures" / "features_promoter_corr_fc.png"}
OUTCOME = ROOT / "results" / "aim1" / "RPE1_vs_K562" / "gene_outcome.tsv"

HOUSEKEEPING = ["GAPDH", "ACTB", "EEF1A1", "RPL13A"]
# (gene, lines expected higher, lines expected lower): pass iff min(higher) > max(lower).
LINEAGE = [("HBG1", ["K562"], ["RPE1"]), ("GATA1", ["K562"], ["RPE1"]),
           ("ALB", ["HepG2"], ["K562", "RPE1"]), ("APOA1", ["HepG2"], ["K562", "RPE1"]),
           ("CD3E", ["Jurkat"], ["RPE1", "HepG2"]), ("LCK", ["Jurkat"], ["RPE1", "HepG2"]),
           ("KRT8", ["RPE1", "HepG2"], ["K562", "Jurkat"]), ("KRT18", ["RPE1", "HepG2"], ["K562", "Jurkat"])]
TARGETS = ["RPE1", "K562", "Jurkat", "HepG2"]
NAMED_PAIRS = [("RPE1", "K562"), ("K562", "Jurkat"), ("RPE1", "HepG2")]

EARLIER_FC = """\
## Earlier FC-only result (2026-10-11 build, superseded)

10 lines (Jurkat absent; RPE1 = provisional 4DN FC bigWig + conservative peaks), 338,483 union
regions, MACS2 fold-change track only. Lineage markers passed, but promoter FC vs baseline
expression over the 5,262 Aim 1 genes was Spearman -0.14 (RPE1) and -0.09 (K562), and GAPDH /
ACTB sat at the 53rd-93rd percentile. Same sign from the raw FC bigWig at any window, so it was
not a processing bug; diagnosis = MACS2 local-lambda background compresses strong promoters.
This led to plan amendment v1.2 (primary signal = Tn5 insertion CPM; FC kept as sensitivity)."""


def line_order(gf):
    return [str(x) for x in gf.line.cat.categories] if hasattr(gf.line, "cat") else list(dict.fromkeys(gf.line))


def heatmap(corr, lines, path, title):
    off = corr.values[~np.eye(len(lines), dtype=bool)]
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    im = ax.imshow(corr.values, vmin=min(0.5, off.min()), vmax=1, cmap="viridis")
    ax.set_xticks(range(len(lines)), lines, rotation=60, ha="right")
    ax.set_yticks(range(len(lines)), lines)
    for i in range(len(lines)):
        for j in range(len(lines)):
            ax.text(j, i, f"{corr.values[i, j]:.2f}", ha="center", va="center", fontsize=7,
                    color="white" if corr.values[i, j] < 0.85 else "black")
    fig.colorbar(im, ax=ax, label="Pearson r")
    ax.set_title(title, fontsize=10)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def section(signal: str, ex: pd.DataFrame, w) -> list[tuple[str, bool]]:
    paths = F.output_paths(signal)
    info = json.loads(paths["build_info"].read_text())
    gf = F.load_gene_features(paths["gene_features"])
    lines = line_order(gf)
    P = F.feature_matrix("promoter", gf)[lines]
    D = F.feature_matrix("distal", gf)[lines]
    C = F.feature_matrix("combined", gf)[lines]
    name = gf.drop_duplicates("gene_id").set_index("gene_id").gene_name
    by_name = {v: k for k, v in name.items()}
    pct = P.rank(pct=True) * 100
    checks = []
    label = "CPM (primary)" if signal == "cpm" else "FC (sensitivity)"
    w(f"\n# Version: {label}\n")
    w(f"Signal: {info['signal_desc']}. Built {info['built']}; `{paths['gene_features'].relative_to(ROOT)}`. "
      f"Lines ({info['n_lines']}): " + ", ".join(f"{d['line']} ({d['role']})" for d in info["lines"]) + ".")
    if info["skipped_lines"]:
        w("Skipped: " + ", ".join(s["line"] for s in info["skipped_lines"]) + ".")
    if signal == "cpm":
        w("\nInsertions (chr1-22,X, replicates pooled): " + ", ".join(
            f"{d['line']} {d['n_insertions_total'] / 1e6:.1f}M ({len(d['insertions'])} reps)" for d in info["lines"]) + ".")
    ri = info["region_info"]
    w(f"\n## Union regions ({signal})\n\n- Summit windows pooled: {ri['n_windows']:,}; merged: {ri['n_merged']:,}; "
      f"blacklist-overlapping dropped: {ri['n_blacklisted']:,}.\n- **Union regions: {info['n_regions']:,}.**")
    reg = pd.read_parquet(paths["region_signal"], columns=["start", "end"])
    L = reg.end - reg.start
    w(f"- Region length: median {L.median():.0f} bp, 90th pct {L.quantile(.9):.0f}, max {L.max():,} bp.")
    w(f"- Genes: {info['n_genes']:,}; gene-region distal pairs {info['n_gene_region_pairs']:,}; "
      f"genes with no union region 2.5-100 kb away: {info['n_genes_no_distal_region']:,}.")
    z = info["z_pooled"]
    w(f"- Pooled z: promoter mean {z['promoter_mean']:.3f} SD {z['promoter_sd']:.3f}; "
      f"distal mean {z['distal_mean']:.3f} SD {z['distal_sd']:.3f}.")
    w("- Per-line summit counts: " + ", ".join(f"{d['line']} {d['n_summits']:,}" for d in info["lines"]) + ".")

    corr = P.corr(method="pearson")
    w(f"\n## Promoter feature correlation across lines ({signal}; Pearson, all {len(P):,} genes)\n")
    w("| | " + " | ".join(lines) + " |\n|---|" + "---|" * len(lines))
    for a in lines:
        w(f"| {a} | " + " | ".join(f"{corr.loc[a, b]:.2f}" for b in lines) + " |")
    off = corr.values[~np.eye(len(lines), dtype=bool)]
    dc = D.corr().values[~np.eye(len(lines), dtype=bool)]
    w(f"\nOff-diagonal range {off.min():.2f}-{off.max():.2f} (median {np.median(off):.2f}). "
      f"Distal off-diagonal range {dc.min():.2f}-{dc.max():.2f}. Promoter vs distal within line "
      f"(median over lines): {np.median([np.corrcoef(P[l], D[l])[0, 1] for l in lines]):.2f}. "
      f"Heatmap: `{FIG[signal].relative_to(ROOT)}`.")
    heatmap(corr, lines, FIG[signal], f"Promoter accessibility ({'insertion CPM' if signal == 'cpm' else 'fold change'}, "
            f"TSS±1 kb), {len(P):,} genes")

    w(f"\n## Promoter feature vs baseline expression ({signal}; expect Spearman > 0)\n")
    w("Source: Aim 1 RPE1_vs_K562 `gene_outcome.tsv`, columns gene_id, mean_ctrlA (RPE1), mean_ctrlB (K562) only.\n")
    exs = ex[ex.index.isin(P.index)]
    R = gf.pivot(index="gene_id", columns="line", values="promoter_raw")
    R.columns = [str(c) for c in R.columns]
    w("| line | n genes | rho promoter | rho raw promoter signal (pre log/QN/CN) | rho distal | rho combined "
      "| pass (promoter > 0) | original > 0.3 |\n|---|---|---|---|---|---|---|---|")
    for ln, col in (("RPE1", "mean_ctrlA"), ("K562", "mean_ctrlB")):
        r = {k: spearmanr(M.loc[exs.index, ln], exs[col])[0] for k, M in (("p", P), ("d", D), ("c", C))}
        ok = r["p"] > 0
        checks.append((f"{signal}: expression rho > 0 ({ln})", bool(ok)))
        rr = spearmanr(R.loc[exs.index, ln], exs[col])[0]
        w(f"| {ln} | {len(exs):,} | {r['p']:.3f} | {rr:.3f} | {r['d']:.3f} | {r['c']:.3f} | {'PASS' if ok else 'FAIL'} | "
          f"{'yes' if r['p'] > 0.3 else 'no'} |")
    aucs = []
    for ln in ("RPE1", "K562"):
        inn = P.index.isin(exs.index)
        a, b = P.loc[inn, ln].values, P.loc[~inn, ln].values
        aucs.append(f"{ln} {mannwhitneyu(a, b).statistic / (len(a) * len(b)):.3f}")
    w("\nAUC of the promoter feature, Aim 1 expressed set vs all other genes: " + "; ".join(aucs) + ".")
    xd = P.loc[exs.index, "K562"] - P.loc[exs.index, "RPE1"]
    ed = np.log2(exs.mean_ctrlB) - np.log2(exs.mean_ctrlA)
    w(f"Between-line (descriptive): Spearman(promoter K562-RPE1, log2 expr K562/RPE1) = {spearmanr(xd, ed)[0]:.3f}.")

    w(f"\n## Housekeeping promoters ({signal}; within-line percentile; expect >= 90th)\n")
    w("| gene | " + " | ".join(lines) + " | lines >= 90th | pass (all) |\n|---|" + "---|" * (len(lines) + 2))
    for g in HOUSEKEEPING:
        gid = by_name[g]
        n90 = int((pct.loc[gid] >= 90).sum())
        ok = n90 == len(lines)
        checks.append((f"{signal}: housekeeping {g} >= 90th in all lines", ok))
        w(f"| {g} | " + " | ".join(f"{pct.loc[gid, l]:.0f}" for l in lines) + f" | {n90}/{len(lines)} | {'PASS' if ok else 'FAIL'} |")

    w(f"\n## Lineage markers ({signal}; promoter feature value [within-line percentile])\n")
    w("| gene | " + " | ".join(lines) + " | check | pass |\n|---|" + "---|" * (len(lines) + 2))
    for g, hi, lo in LINEAGE:
        gid = by_name[g]
        v = P.loc[gid]
        ok = bool(min(v[h] for h in hi) > max(v[x] for x in lo))
        checks.append((f"{signal}: {g} {'/'.join(hi)} > {'/'.join(lo)}", ok))
        w(f"| {g} | " + " | ".join(f"{v[l]:.2f} [{pct.loc[gid, l]:.0f}]" for l in lines)
          + f" | {', '.join(hi)} > {', '.join(lo)} | {'PASS' if ok else 'FAIL'} |")
    w("\nNote: KRT8/KRT18 are expressed in K562 too (Aim 1 baseline mean_ctrlB 1.0 / 4.3 vs RPE1 mean_ctrlA "
      "3.9 / 7.9), so K562 promoter accessibility at these genes is plausible; the check is kept as pre-stated.")
    npass = sum(ok for _, ok in checks)
    w(f"\n**{signal}: {npass}/{len(checks)} checks pass.**")
    return checks


def replicate_section(w, ex_ids):
    rep = pd.read_parquet(F.output_paths("cpm")["promoter_reps"]).set_index("gene_id")
    lines = list(dict.fromkeys(c.split("__")[0] for c in rep.columns))
    rows = []
    for ln in lines:
        cols = [c for c in rep.columns if c.split("__")[0] == ln]
        for a, b in itertools.combinations(cols, 2):
            sub = rep.loc[rep.index.isin(ex_ids)]
            rows.append(dict(line=ln, rep_a=a.split("__")[1], rep_b=b.split("__")[1],
                             spearman_all=spearmanr(rep[a], rep[b])[0],
                             spearman_aim1=spearmanr(sub[a], sub[b])[0],
                             pearson_log2cpm_all=pearsonr(np.log2(rep[a] + 1), np.log2(rep[b] + 1))[0]))
    t = pd.DataFrame(rows)
    t.to_csv(REP_TSV, sep="\t", index=False, float_format="%.4f")
    w("\n# Replicate reproducibility of promoter insertion CPM (TSS±1 kb)\n")
    w(f"Per-replicate CPM (each replicate normalised by its own chr1-22,X total). Spearman over all "
      f"{len(rep):,} genes and over the {len(ex_ids):,} Aim 1 genes (membership only); Pearson on log2(CPM+1). "
      f"Table: `{REP_TSV.relative_to(ROOT)}`.\n")
    w("| line | reps | Spearman (all genes) | Spearman (Aim 1 genes) | Pearson log2 (all) |\n|---|---|---|---|---|")
    for _, r in t.iterrows():
        w(f"| {r.line} | {r.rep_a} vs {r.rep_b} | {r.spearman_all:.3f} | {r.spearman_aim1:.3f} | {r.pearson_log2cpm_all:.3f} |")
    w(f"\nRange over all pairs: Spearman (all genes) {t.spearman_all.min():.3f}-{t.spearman_all.max():.3f}; "
      f"Aim 1 genes {t.spearman_aim1.min():.3f}-{t.spearman_aim1.max():.3f}.")
    return t


def pair_section(w, ex_ids):
    rows = []
    gc, gfc = F.load_gene_features(signal="cpm"), F.load_gene_features(signal="fc")
    pairs = NAMED_PAIRS + [p for p in itertools.combinations(TARGETS, 2)
                           if p not in NAMED_PAIRS and p[::-1] not in NAMED_PAIRS]
    for a, b in pairs:
        for kind in F.KINDS:
            xc = F.pair_feature(a, b, kind, gc)
            xf = F.pair_feature(a, b, kind, gfc).reindex(xc.index)
            s = xc.index.isin(ex_ids)
            rows.append(dict(pair=f"{a}-{b}", named=(a, b) in NAMED_PAIRS, kind=kind,
                             pearson_all=pearsonr(xc, xf)[0], spearman_all=spearmanr(xc, xf)[0],
                             pearson_aim1=pearsonr(xc[s], xf[s])[0], spearman_aim1=spearmanr(xc[s], xf[s])[0]))
    t = pd.DataFrame(rows)
    t.to_csv(PAIR_TSV, sep="\t", index=False, float_format="%.4f")
    w("\n# Agreement of x = feature_A - feature_B between CPM and FC versions (features only)\n")
    w(f"Correlation of x_cpm with x_fc over all genes and over the Aim 1 gene set (membership only). "
      f"First three pairs are the named target pairs; the remaining target-line pairs are listed for completeness. "
      f"Table: `{PAIR_TSV.relative_to(ROOT)}`.\n")
    w("| pair | kind | Pearson (all) | Spearman (all) | Pearson (Aim 1) | Spearman (Aim 1) |\n|---|---|---|---|---|---|")
    for _, r in t.iterrows():
        w(f"| {'**' + r.pair + '**' if r.named else r.pair} | {r.kind} | {r.pearson_all:.3f} | {r.spearman_all:.3f} | "
          f"{r.pearson_aim1:.3f} | {r.spearman_aim1:.3f} |")
    return t


def main():
    md = []
    w = md.append
    ex = pd.read_csv(OUTCOME, sep="\t", usecols=["gene_id", "mean_ctrlA", "mean_ctrlB"]).set_index("gene_id")
    w("# Accessibility feature sanity checks (plan v1.2: CPM primary, FC sensitivity)\n")
    w("Generated by `scripts/features_sanity.py`. Expression checks read only gene_id, mean_ctrlA (RPE1) "
      "and mean_ctrlB (K562) of results/aim1/RPE1_vs_K562/gene_outcome.tsv. Pass criteria: expression "
      "Spearman > 0 (v1.2; the original > 0.3 is also shown), housekeeping >= 90th percentile in every line, "
      "lineage marker ordering as listed.")
    summary_at = len(md)
    checks = {}
    for sig in F.SIGNALS:
        if not F.output_paths(sig)["gene_features"].exists():
            w(f"\n# Version {sig}: not built\n")
            continue
        checks[sig] = section(sig, ex, w)
    rep = replicate_section(w, set(ex.index)) if F.output_paths("cpm")["promoter_reps"].exists() else None
    if set(checks) == set(F.SIGNALS):
        pair_section(w, set(ex.index))
    w("\n" + EARLIER_FC)
    summ = ["\n## Summary\n"]
    for sig, ch in checks.items():
        fails = [n for n, ok in ch if not ok]
        summ.append(f"- **{sig}: {sum(ok for _, ok in ch)}/{len(ch)} checks pass.**"
                    + (" Failing: " + "; ".join(n.split(': ', 1)[1] for n in fails) + "." if fails else ""))
    if rep is not None:
        summ.append(f"- Replicate promoter-CPM Spearman (all genes): {rep.spearman_all.min():.3f}-{rep.spearman_all.max():.3f}.")
    md[summary_at:summary_at] = summ
    OUT_MD.parent.mkdir(parents=True, exist_ok=True)
    OUT_MD.write_text("\n".join(md) + "\n")
    print("\n".join(md))
    return 0 if checks and all(ok for ch in checks.values() for _, ok in ch) else 1


if __name__ == "__main__":
    sys.exit(main())
