"""EXPLORATORY (post hoc) report for scripts/aim4x_rpe1_diagnostic.py.

Reads results/aim4_exploratory/rpe1_diagnostic_raw.json and writes
results/aim4_exploratory/rpe1_diagnostic.{md,json} and figures/aim4x_rpe1_diagnostic.png.
"""
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/aim4_exploratory"
RAW = json.loads((OUT / "rpe1_diagnostic_raw.json").read_text())
PAIRS = ("RPE1_vs_K562", "Jurkat_vs_RPE1", "RPE1_vs_HepG2")

# (key, source, display label, group)
ROWS = [
    ("primary", "p", "primary (frozen main test; reproduced)", "0 primary"),
    ("cn_rpe1_depmap", "p", "CN: RPE1 = DepMap RPE1-clone consensus (features + cn_diff)", "1 copy number"),
    ("cn_rpe1_depmap__features_only", "p", "  … features only", "1 copy number"),
    ("cn_rpe1_depmap__covariate_only", "p", "  … cn_diff covariate only", "1 copy number"),
    ("cn_rpe1_none", "p", "CN: no CN correction for RPE1 (features + cn_diff)", "1 copy number"),
    ("rpe1_rep1_only", "p", "RPE1 CPM from rep1 only", "2 accessibility source"),
    ("rpe1_rep2_only", "p", "RPE1 CPM from rep2 only", "2 accessibility source"),
    ("rpe1_4dn_fc_hybrid", "p", "RPE1 = 4DN-processed FC bigWig (others CPM)", "2 accessibility source"),
    ("rpe1_local_fc_hybrid", "p", "RPE1 = local-pipeline FC bigWig (others CPM; control)", "2 accessibility source"),
    ("rpe1_ds_jurkat_depth_116M", "ds", "RPE1 CPM thinned to Jurkat depth (115.8M ins.), mean of 3 seeds", "2 accessibility source"),
    ("rpe1_ds_58M_insertions", "ds", "RPE1 CPM thinned to 58M insertions, mean of 3 seeds", "2 accessibility source"),
    ("fc_all_lines_primary", "p", "all lines FC (frozen FC sensitivity)", "2 accessibility source"),
    ("fc_all_lines_rpe1_4dn", "p", "all lines FC, RPE1 = 4DN FC + 4DN peaks", "2 accessibility source"),
    ("partner_fc_hybrid", "x", "partner = its FC, RPE1 kept CPM (control)", "2 accessibility source"),
    ("both_fc_hybrid", "x", "partner and RPE1 = FC (others CPM)", "2 accessibility source"),
    ("rpe1_cpm_shared_with_fc", "x", "RPE1 CPM part shared with RPE1 FC (RPE1 side only)", "2 accessibility source"),
    ("rpe1_cpm_specific_part", "x", "RPE1 CPM part NOT shared with RPE1 FC (RPE1 side only)", "2 accessibility source"),
    ("promoter", "p", "promoter feature", "3 promoter vs distal"),
    ("distal", "p", "distal feature", "3 promoter vs distal"),
    ("excl_chr10", "p", "exclude chr10", "4 chromosomes"),
    ("excl_chrX", "p", "exclude chrX", "4 chromosomes"),
    ("excl_chr10_and_chrX", "p", "exclude chr10 and chrX", "4 chromosomes"),
    ("excl_fischer_p53", "p", "exclude Fischer p53 core targets", "6 p53 / expression"),
    ("excl_hallmark_p53", "p", "exclude HALLMARK_P53_PATHWAY", "6 p53 / expression"),
    ("excl_p53_union", "p", "exclude both p53 sets", "6 p53 / expression"),
    ("rpe1_expr_bottom_tercile", "p", "bottom tercile of RPE1 expression", "6 p53 / expression"),
    ("rpe1_expr_middle_tercile", "p", "middle tercile of RPE1 expression", "6 p53 / expression"),
    ("rpe1_expr_top_tercile", "p", "top tercile of RPE1 expression", "6 p53 / expression"),
    ("side_rpe1_alone", "p", "RPE1 feature alone (oriented, resid. on C)", "7 which side"),
    ("side_partner_alone", "p", "partner feature alone (oriented, resid. on C)", "7 which side"),
]


def get(pair, key, src):
    if src == "p":
        return RAW["pairs"][pair][key]
    if src == "x":
        return RAW["extra"][pair][key]
    runs = [RAW["pairs"][pair][f"{key}_s{s}"] for s in (1, 2, 3)]
    th = np.array([r["theta_hat"] for r in runs])
    return {"theta_hat": float(th.mean()), "theta_seed_range": [float(th.min()), float(th.max())],
            "ci95": [float(np.mean([r["ci95"][0] for r in runs])), float(np.mean([r["ci95"][1] for r in runs]))],
            "p_two_sided": float(max(r["p_two_sided"] for r in runs)), "pr2_rel": float(np.mean([r["pr2_rel"] for r in runs])),
            "G": runs[0]["G"], "k_x4": runs[0]["k_x4"], "seeds": runs, "note": "CI = mean of seed CIs; p = max over seeds"}


def fmt(v, d=4):
    return f"{v:.{d}f}" if abs(v) >= 1e-3 or v == 0 else f"{v:.2e}"


def main():
    tab = {p: {k: get(p, k, s) for k, s, *_ in ROWS} for p in PAIRS}
    lab = {k: l for k, _, l, _ in ROWS}
    grp = {k: g for k, _, _, g in ROWS}
    ex = RAW["extra"]
    # ---------------------------------------------------------------- json
    js = {"label": "EXPLORATORY (post hoc) -- H1 was run once on pipeline-v1.1; nothing here is confirmatory",
          "date": "2026-10-12", "pipeline": "pipeline-v1.1 (plan v1.3) frozen settings: block permutation, x4 rule, "
          "arm terms in f(C), B = 9,999, aim2 seeds (same seed for every variant), refined test-inversion CI",
          "pairs": {p: {k: {kk: vv for kk, vv in v.items() if kk != "seeds"} for k, v in tab[p].items()} for p in PAIRS},
          "loco_jackknife": {p: RAW["pairs"][p]["loco_jackknife"] for p in PAIRS},
          "p53_sets": {p: RAW["pairs"][p]["p53_sets"] for p in PAIRS},
          "p53_numerator_share": {p: RAW["pairs"][p]["p53_numerator_share"] for p in PAIRS},
          "joint_two_line_ols": ex["joint_two_line_ols_all_pairs"],
          "side_tests_pairs_without_rpe1": {p: ex[p] for p in ("K562_vs_Jurkat", "K562_vs_HepG2", "Jurkat_vs_HepG2")},
          "tss_enrichment_insertions": ex["tss_enrichment_insertions"],
          "rpe1_cpm_vs_fc_feature_r": {p: ex[p]["rpe1_cpm_vs_fc_r"] for p in PAIRS},
          "feature_meta": RAW["feature_meta"]}
    # ---------------------------------------------------------------- figure
    keys = [k for k, *_ in ROWS]
    fig, axes = plt.subplots(1, 3, figsize=(17, 11), sharey=True)
    y = np.arange(len(keys))[::-1]
    cols = {"0": "k", "1": "tab:blue", "2": "tab:red", "3": "tab:green", "4": "tab:purple", "6": "tab:orange",
            "7": "tab:brown"}
    for ax, p in zip(axes, PAIRS):
        prim = tab[p]["primary"]["theta_hat"]
        for yi, k in zip(y, keys):
            r = tab[p][k]
            c = cols[grp[k][0]]
            lo, hi = r["ci95"]
            ax.plot([lo, hi], [yi, yi], color=c, lw=1.6)
            ax.plot(r["theta_hat"], yi, "o", color=c, ms=5, mfc=c if r["p_two_sided"] <= 0.05 else "white")
        ax.axvline(0, color="grey", lw=0.8)
        ax.axvline(prim, color="k", lw=0.6, ls=":")
        ax.set_title(p.replace("_vs_", " − "))
        ax.set_xlabel("θ̂ (95% CI, block permutation)")
    axes[0].set_yticks(y)
    axes[0].set_yticklabels([lab[k].strip() for k in keys], fontsize=7.5)
    fig.suptitle("EXPLORATORY (post hoc): RPE1 diagnostic — θ̂ of the combined accessibility feature across variants\n"
                 "filled = two-sided p ≤ 0.05; dotted = primary θ̂; colour = diagnostic group (1 CN, 2 source, 3 location, "
                 "4 chromosomes, 6 p53/expression, 7 side)", fontsize=10)
    fig.tight_layout(rect=(0.04, 0, 1, 0.95))
    fig.savefig(ROOT / "figures/aim4x_rpe1_diagnostic.png", dpi=150)
    # ---------------------------------------------------------------- markdown
    L = ["# RPE1 diagnostic — EXPLORATORY (post hoc)", "",
         "**EXPLORATORY (post hoc).** H1 was tested once on the frozen pipeline (`pipeline-v1.1`, plan v1.3; "
         "`results/aim4/`). Everything below was decided after seeing that result, is not confirmatory, and "
         "changes no primary file or tag. Question: is the negative accessibility-response association that "
         "appears only in pairs containing RPE1 an RPE1 processing artefact or biology?", "",
         "Method: for each RPE1 pair the frozen main test (block permutation, x4 k rule, chromosome-arm terms "
         "in f(C), B = 9,999, the aim2 seed of the primary feature for every variant, refined test-inversion "
         "CI) is rerun with one input changed. Feature variants are rebuilt from the saved raw signal with the "
         "frozen feature steps (log2 → quantile normalisation across the 11 lines → − log2 CN → distal sum → "
         "combined z); only the RPE1 column (or RPE1 CN track) changes. The rebuild reproduces the primary "
         "feature to ≤ 3e-8 and the primary θ̂, CI and p exactly. pR²/rel = partial R² / the pair's full-data "
         "reliability (no bootstrap). Gene-subset variants refit f(C) (splines, arm terms) and k on the subset. "
         "Code: `scripts/aim4x_rpe1_diagnostic.py` (stages features, fc4dn, tests, extra, joint), "
         "`scripts/aim4x_rpe1_report.py`; logs in `results/aim4_exploratory/logs/`.", ""]
    for p in PAIRS:
        L += [f"## {p.replace('_vs_', ' − ')} (x = A − B)", "",
              "| group | variant | genes | θ̂ | 95% CI | p two-sided | pR²/rel |", "|---|---|---|---|---|---|---|"]
        for k in keys:
            r = tab[p][k]
            extra = f" (seed range {fmt(r['theta_seed_range'][0])} … {fmt(r['theta_seed_range'][1])})" if "theta_seed_range" in r else ""
            L.append(f"| {grp[k][2:]} | {lab[k].strip()}{extra} | {r['G']:,} | {fmt(r['theta_hat'])} | "
                     f"[{fmt(r['ci95'][0])}, {fmt(r['ci95'][1])}] | {r['p_two_sided']:.4g} | {r['pr2_rel']:.3%} |")
        j = RAW["pairs"][p]["loco_jackknife"]
        L += ["", f"Leave-one-chromosome-out (θ̂ only): range {fmt(j['min'][1])} (drop {j['min'][0]}) to "
              f"{fmt(j['max'][1])} (drop {j['max'][0]}); jackknife SE {fmt(j['jackknife_se'])}; "
              f"{j['n_chrom_positive']}/23 drops give θ̂ > 0.", ""]
    hj = RAW["pairs"]["RPE1_vs_HepG2"]["loco_jackknife"]["theta_drop"]
    L += ["## (5) Leave-one-chromosome-out jackknife, RPE1 − HepG2", "",
          "| dropped | " + " | ".join(hj) + " |", "|---|" + "---|" * len(hj),
          "| θ̂ | " + " | ".join(f"{v:.4f}" for v in hj.values()) + " |", ""]
    L += ["## (6) p53 gene sets", "", "| pair | Fischer core in pair | Hallmark in pair | union | share of θ̂ numerator from union | share of genes |",
          "|---|---|---|---|---|---|"]
    for p in PAIRS:
        s, n = RAW["pairs"][p]["p53_sets"], RAW["pairs"][p]["p53_numerator_share"]
        L.append(f"| {p} | {s['fischer_in_pair']} | {s['hallmark_in_pair']} | {s['union_in_pair']} | {n['union']:.1%} | {n['expected_if_average']:.1%} |")
    L += ["", "Fischer: 116 core targets listed (`fischer2017_p53_core_mapped.tsv`); Hallmark: 200 symbols "
          "(`HALLMARK_P53_PATHWAY.v2025.1`), matched by gene name.", ""]
    L += ["## (7) Which side — joint two-line model and pairs without RPE1", "",
          "Joint OLS of M (residualised on f(C)) on both lines' residualised features, each oriented as in x "
          "(A feature, −B feature); descriptive, no permutation. The two residualised features correlate "
          "−0.79 to −0.85 (shared accessibility), so the marginal one-line tests above are the cleaner read.", "",
          "| pair | coef A | coef B | r (resid., oriented) |", "|---|---|---|---|"]
    for p, v in ex["joint_two_line_ols_all_pairs"].items():
        c = [vv for kk, vv in v.items() if kk.startswith("coef_")]
        L.append(f"| {p} | {c[0]:.4f} | {c[1]:.4f} | {v['corr_resid_oriented']:.2f} |")
    L += ["", "One-line-alone main tests in the pairs without RPE1 (context):", "",
          "| pair | variant | θ̂ | p two-sided |", "|---|---|---|---|"]
    for p in ("K562_vs_Jurkat", "K562_vs_HepG2", "Jurkat_vs_HepG2"):
        for k, v in ex[p].items():
            L.append(f"| {p} | {k} | {fmt(v['theta_hat'])} | {v['p_two_sided']:.4g} |")
    te = ex["tss_enrichment_insertions"]
    from cdtest import features as F
    import pandas as pd
    a_, b_ = F.feature_matrix("combined"), F.feature_matrix("combined", signal="fc")
    gg = pd.read_csv(ROOT / "results/aim1/RPE1_vs_HepG2/gene_outcome.tsv", sep="\t").gene_id
    rline = {c: float(np.corrcoef(a_.loc[gg, c], b_.loc[gg, c])[0, 1]) for c in a_.columns}
    js["cpm_vs_fc_combined_r_per_line_RPE1_vs_HepG2_genes"] = rline
    L += ["", "## Data-quality context", "",
          "Insertion TSS enrichment (TSS ± 50 bp vs flanks 1.9–2.0 kb, protein-coding, replicates pooled): " +
          ", ".join(f"{k} {v:.1f}" for k, v in sorted(te.items(), key=lambda t: t[1])) + ". "
          "CPM vs FC combined feature (the two frozen builds), per line over the RPE1 − HepG2 genes: " +
          ", ".join(f"{k} {v:.2f}" for k, v in sorted(rline.items(), key=lambda t: t[1])) + " — RPE1 is not an "
          "outlier; FC is a substantially different summary for every line. Within the hybrid (RPE1 local FC "
          "quantile-normalised against the CPM lines) the RPE1 CPM vs FC correlation over each pair's genes is " +
          ", ".join(f"{p} {ex[p]['rpe1_cpm_vs_fc_r']:.2f}" for p in PAIRS) + ".", ""]
    sl = ex["same_line_link"]
    js["same_line_link"] = sl
    L += ["## Link to the planned same-line control (`results/aim4/same_line_control.md`)", "",
          "M = K562 essential − K562 genome-wide (two screens of the same line: technical divergence only). Frozen "
          "main test, seed keyed on the row label (so p differs slightly from `same_line_control.json`; θ̂ "
          "reproduces). partial r = θ̂·SD(resid. x)/SD(resid. M). For comparison, the cross-line partial r of the "
          "pair x in the RPE1 pairs is −0.074 / −0.071 / −0.087, and of the RPE1 CPM-specific part "
          "−0.073 / −0.054 / −0.083 (oriented RPE1 side).", "",
          "| feature tested against same-line M | θ̂ | 95% CI | p two-sided | partial r |", "|---|---|---|---|---|"]
    for k, r in sl["tests"].items():
        L.append(f"| {k} | {fmt(r['theta_hat'])} | [{fmt(r['ci95'][0])}, {fmt(r['ci95'][1])}] | "
                 f"{r['p_two_sided']:.4g} | {r['theta_std']:.3f} |")
    L += [""]
    L += ["## Reading", ""] + READING + ["", "## Conclusion (plain language)", "", CONCLUSION, ""]
    (OUT / "rpe1_diagnostic.md").write_text("\n".join(L))
    js["reading"] = READING
    js["conclusion"] = CONCLUSION
    (OUT / "rpe1_diagnostic.json").write_text(json.dumps(js, indent=1))
    print("written")


READING = [
    "- **Copy number (1): not the cause.** DepMap RPE1-clone consensus CN or no RPE1 CN correction changes θ̂ by "
    "< 0.0004 in every pair (features, covariate, or both). Arm-level CN is already absorbed by the arm terms in f(C).",
    "- **Replicates and depth (2a, 2c): not the cause.** Each RPE1 replicate alone, and the pooled insertions thinned "
    "to Jurkat's depth (115.8M insertions = ~58M fragments) or to 58M insertions (3 seeds each), give θ̂ within "
    "0.0007 of the primary, all p₂ ≤ 0.0002. (RPE1 at 152M insertions is not deep relative to the other lines anyway.)",
    "- **Signal type (2b): the effect disappears when RPE1 accessibility is measured as fold-change.** Swapping "
    "only the RPE1 column for the 4DN-processed FC bigWig gives θ̂ ≈ +0.002 to +0.004 (p₂ 0.16–0.67) in all three "
    "pairs; the local-pipeline FC bigWig from the same reads does the same — so it is not 4DN-vs-local processing "
    "but CPM-vs-FC for RPE1. All-FC features (RPE1 local or 4DN FC/peaks) are also n.s. Swapping only the partner "
    "to FC keeps the effect (θ̂ −0.010 to −0.017, p₂ ≤ 0.004). With both RPE1 and the partner as FC, a smaller "
    "residual remains in RPE1 − K562 (θ̂ −0.007, p₂ 0.03) but not in the other two pairs (−0.006, −0.003, n.s.). "
    "Splitting the RPE1 CPM feature into the part shared with RPE1 FC and the remainder: the shared part gives "
    "θ̂ ≈ +0.007 to +0.008 (n.s.); the CPM-specific remainder carries all of the negative association "
    "(θ̂ −0.011 to −0.015, p₂ 0.0002–0.012).",
    "- **Location (3):** promoter-driven in all three pairs; distal is n.s. in RPE1 − K562 and Jurkat − RPE1 and "
    "negative (smaller pR²) in RPE1 − HepG2.",
    "- **Chromosomes (4, 5): not driven by chr10 (10q gain) or chrX, nor by any single chromosome.** Exclusions "
    "slightly strengthen θ̂; every leave-one-chromosome-out θ̂ stays negative (RPE1 − HepG2: −0.0154 to −0.0187).",
    "- **p53 (6): not driven by p53 genes.** Excluding Fischer core targets, Hallmark p53 pathway genes or both "
    "changes θ̂ by < 0.001 (all p₂ ≤ 0.0001); in two of three pairs the p53 genes contribute less than their "
    "gene share to the θ̂ numerator.",
    "- **Expression (6):** the association grows with RPE1 expression (bottom tercile −0.004 to −0.009, n.s. or "
    "marginal; top tercile −0.028 to −0.034, p₂ ≤ 0.0011).",
    "- **Side (7): it is RPE1's own accessibility.** RPE1 feature alone: θ̂ −0.008 to −0.013, p₂ ≤ 0.003 in every "
    "pair; partner alone: ≈ 0 (p₂ ≥ 0.5). The same one-line tests in the three pairs without RPE1 are all ≈ 0, "
    "and the joint two-line coefficients are ≈ 0 or positive there but negative for both terms in the RPE1 pairs "
    "(the partner term picks up RPE1's deviation from the shared accessibility profile).",
    "- **Data quality:** RPE1 has the lowest insertion TSS enrichment of the 11 lines (10.6 vs 13.7–28.3).",
    "- **Same-line control (planned, added after this diagnostic was drafted):** accessibility does associate "
    "with purely technical between-screen divergence. Every line's own CPM feature correlates positively with "
    "K562 essential − genome-wide M (partial r +0.046 to +0.079; K562 strongest). The RPE1_vs_K562 x reaches a "
    "partial r of −0.046, ~60% of its cross-line size (−0.074), mostly through its −K562 term. But this technical "
    "channel does not have the structure of the RPE1 effect. Replacing RPE1 CPM with FC removes the cross-line "
    "RPE1 association, yet makes the same-line association stronger (partial r about −0.10, p₂ 1e-4). The RPE1 "
    "CPM-specific part, which carries the whole cross-line effect, also relates to same-line M (+0.080), and "
    "the part shared with FC has the opposite sign (−0.046); the pattern does not match the cross-line one. So the same-line result is a general warning "
    "that non-biological between-screen divergence tracks accessibility at a magnitude comparable to the "
    "cross-line effects. It is not a mechanistic explanation of the RPE1-specific pattern.",
]

CONCLUSION = (
    "EXPLORATORY (post hoc). The evidence moderately favours a non-biological explanation of the negative "
    "association in the RPE1 pairs over RPE1 biology. The planned same-line control adds weight against biology "
    "but makes the exact artefact less certain. Within the RPE1 pairs, the association is not explained by "
    "copy-number handling, either replicate, sequencing depth, the 10q gain, chrX, any single chromosome, or "
    "p53 target / p53-pathway genes. It belongs to RPE1's own accessibility, is strongest in highly expressed "
    "genes, and exists only in the insertion-count (CPM) summary of the RPE1 reads. When the same reads are "
    "summarised as fold-change over local background (by the local pipeline or by 4DN's), it vanishes in all "
    "three pairs. It sits entirely in the part of the RPE1 CPM signal that FC does not share, and RPE1 has the "
    "lowest TSS enrichment of the 11 libraries. A real accessibility-response relation should survive both "
    "summaries of the same reads, so this points to a CPM-specific background/normalisation component of the "
    "RPE1 track. The same-line control (two K562 screens, no biological difference) shows that accessibility "
    "differences involving RPE1 also associate with purely technical between-screen divergence. In standardised "
    "terms that is about 60% of the cross-line size, so cross-line θ̂ values of this magnitude cannot be read as "
    "biology by default. That technical association is not the same signal, though: it is mostly carried by the "
    "partner/K562 accessibility and by the FC-shared part of RPE1, and it gets stronger, not weaker, when RPE1 is "
    "switched to FC. The two non-biological channels therefore look distinct, and these data cannot assign the "
    "RPE1-pair effect wholly to one of them. Biology is not excluded: FC and CPM agree only modestly for every "
    "line (r ≈ 0.4–0.7; RPE1 0.55), and no independent RPE1 ATAC dataset was available. The RPE1-pair negatives "
    "and the matched-pair disagreement in H1 should be reported as most likely technical. They are not evidence "
    "of a biological negative accessibility effect."
)

if __name__ == "__main__":
    main()
