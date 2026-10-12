"""Markdown reports for the development pair (Aim 2 test, Aim 3 planted benchmark).

Run after scripts/run_aim2_dev.py and scripts/run_planted_dev.py:
    PYTHONPATH=src python scripts/report_dev.py
"""
import json
from pathlib import Path

import pandas as pd

from cdtest import planted

ROOT = Path(__file__).resolve().parents[1]
TAG = "RPE1_vs_K562"


def f(v, n=4):
    return "NA" if v is None else (f"{v:.{n}g}" if isinstance(v, float) else str(v))


def ci(v):
    return f"[{v[0]:.4g}, {v[1]:.4g}]"


def aim2_md():
    r = json.loads((ROOT / f"results/aim2/{TAG}_test.json").read_text())
    naive = json.loads((ROOT / f"results/aim1/naive_{TAG}.json").read_text())
    L = [f"# Aim 2 test, development pair RPE1 − K562 (A = RPE1, B = K562)", "",
         "**Status: descriptive.** RPE1 and K562 differ in p53 status, tissue and growth, so this is a",
         "mixed pair; θ̂ here is not a test of H1. It is the pipeline's development run, made before",
         "tag `pipeline-v1.0`. Code: `cdtest.aim2` (`scripts/run_aim2_dev.py`).", "",
         f"- Genes {r['n_genes']:,}; knockdowns {r['n_knockdowns']:,}; reliability of M "
         f"{r['reliability']:.4f} (knockdown-bootstrap 95% {ci(r['reliability_boot_ci95'])}).",
         f"- Covariates C: {', '.join(r['covariates'])}; f(C) has {r['design_columns']} columns "
         "(intercept, ns(df 4) on each line's log2 mean, 5 linear; rank one less because the "
         "expression difference is collinear).",
         f"- Permutations B = {r['B_perm']:,} per test; partial-R² bootstrap {r['n_boot']} resamples of knockdowns.",
         f"- CAMERA ρ̄ (control cells, normalized counts, mean of the two lines) = {r['rho_bar_used']:.5f} "
         f"(RPE1 {r['rho_bar_controls']['RPE1']:.5f}, K562 {r['rho_bar_controls']['K562']:.5f}; on log2(x+1): "
         f"{r['rho_bar_controls_log2p1']['RPE1']:.5f} / {r['rho_bar_controls_log2p1']['K562']:.5f}).", "",
         "## Main test (block permutation, k rule x4 = plan v1.1)", "",
         "| feature | θ̂ | k (x4 / v1.0) | one-sided p | two-sided p | 95% CI (refined) | 95% CI (plan grid) | partial R² / reliability [95% boot] |",
         "|---|---|---|---|---|---|---|---|"]
    for lab, d in r["features"].items():
        m = d["main_block_x4"]
        L.append(f"| {lab} | {d['theta_hat']:.4g} | {d['k_x4']} / {d['k_v1.0']} | {m['p_one_sided']:.4g} | "
                 f"{m['p_two_sided']:.4g} | {ci(m['ci95'])} | {ci(m['ci95_plan_grid'])} | "
                 f"{d['partial_r2_share_reproducible']:.3%} {ci([100 * x for x in d['partial_r2_share_boot_ci95']])}% |")
    L += ["", "Block exchangeability (main x4 rule; v1.0 in brackets): share of genes in blocks that have at",
          "least one other block of the same size, i.e. that the permutation can move.", "",
          "| feature | blocks | movable blocks | genes movable | null mean of T* | null SD of T* |", "|---|---|---|---|---|---|"]
    for lab, d in r["features"].items():
        e, e1, m = d["exchangeability_x4"], d["exchangeability_v1.0"], d["main_block_x4"]
        L.append(f"| {lab} | {e['n_blocks']} [{e1['n_blocks']}] | {e['n_movable_blocks']} [{e1['n_movable_blocks']}] | "
                 f"{e['frac_genes_movable']:.0%} [{e1['frac_genes_movable']:.0%}] | {m['null_mean']:.4g} | {m['null_sd']:.3g} |")
    L += ["", "## Sensitivity analyses and comparators (one-sided p for θ > 0; two-sided p; 95% CI)", ""]
    names = ["block_v1.0", "naive_ols", "unrestricted", "expr_bin", "stratified", "camera",
             "block_coaccess_clusters", "decoy_neg_control"]
    L += ["| feature | " + " | ".join(names) + " |", "|---" * (len(names) + 1) + "|"]
    for lab, d in r["features"].items():
        cells = []
        for n in names:
            c = d["comparators"][n]
            if n == "decoy_neg_control":
                cells.append(f"p={f(c.get('p_one_sided'))} / {f(c.get('p_two_sided'))} (n={c['n_neg_control']})")
            else:
                cells.append(f"{c['p_one_sided']:.3g} / {c['p_two_sided']:.3g} {ci(c['ci95'])}")
        L.append(f"| {lab} | " + " | ".join(cells) + " |")
    L += ["", "## Naive analysis (Aim 1 step 7: M ~ 1 + expression difference + x, genes independent)", "",
          "| feature | coefficient | SE | partial R² | one-sided p | two-sided p |", "|---|---|---|---|---|---|"]
    for lab in r["features"]:
        n = naive[lab]
        L.append(f"| {lab} | {n['coef']:.4g} | {n['se']:.3g} | {n['partial_r2']:.3%} | {n['p_one_sided']:.3g} | "
                 f"{n['p_two_sided']:.3g} |")
    c = r["features"]["combined_cpm"]
    L += ["", "## Reading", "",
          f"- The primary feature (combined, CPM) has a **negative** coefficient, θ̂ = {c['theta_hat']:.4f}, "
          f"so the one-sided test of θ > 0 gives p = {c['main_block_x4']['p_one_sided']:.4f}. Genes more accessible "
          "in RPE1 than in K562 have a *smaller* RPE1-minus-K562 response difference. Two-sided, the block null "
          f"puts it at p = {c['main_block_x4']['p_two_sided']:.4f}; every comparator agrees on the sign. The "
          "promoter feature carries it; distal and the FC version are near zero.",
          f"- Size: partial R² {c['partial_r2_share_reproducible']:.2%} of reproducible variance for combined, "
          f"{r['features']['promoter_cpm']['partial_r2_share_reproducible']:.2%} for promoter.",
          "- This is the development pair, which differs in p53, tissue and growth; it is not evidence for or "
          "against H1.",
          "- **Pipeline issue found here (for the author, before the matched pairs).** At G ≈ 5,300 genes on 41 "
          "arms the x residuals stay autocorrelated past lag 50, so k_v1.0 hits its cap (50) and the x4 rule its "
          "cap (100). Blocks are then 100–199 genes, and because the plan permutes blocks only within equal-size "
          "classes, only ~41% of genes (21 of 49 blocks) can move. The fixed blocks keep their real x, so the "
          "block null is centred near the observed association (null mean of T* ≈ −0.017 for combined), not at 0. "
          "The test stays valid (exchangeable units), but it tests only the movable part and loses power; see the "
          "planted benchmark (results/aim3/planted_RPE1_vs_K562.md). No parameter was changed.",
          "- CI: the plan's 201-point grid spans ±max(4|θ̂|, 0.5/SD(x)); with SD(x) ≈ 0.3 its step (≈0.017) is "
          "about the size of θ̂, so the plan-grid CI collapses to one or two grid points. The table also gives the "
          "same inversion on a 2,001-point grid inside the plan span (docs/methods_aim2_aim3.md §6.2).",
          "- Partial-R² bootstrap: approximate (β̃ and s̃² held at full-data values; see methods §6.3). It "
          "reflects knockdown sampling only, not gene-level dependence.", ""]
    (ROOT / f"results/aim2/{TAG}_test.md").write_text("\n".join(L))


def aim3_md():
    s = pd.read_csv(ROOT / f"results/aim3/planted_{TAG}.csv")
    info = json.loads((ROOT / f"results/aim3/planted_{TAG}_info.json").read_text())
    d = pd.read_parquet(ROOT / f"results/aim3/planted_{TAG}_replicates.parquet")
    L = [f"# Aim 3 planted-signal benchmark on real M, RPE1 − K562", "",
         "Code: `cdtest.planted` (`scripts/run_planted_dev.py`). Plan §6: M' = M + θ·v, v = (random sign) × a "
         "negative-control decoy feature, θ set so the planted partial R² (raw, given f(C)) is the level; "
         f"{info['reps_per_level']} replicates per level, B = {info['B']}, every test applied to v.", "",
         f"- Decoys: 21 decoy-line pairs (combined CPM feature). Rule: {info['decoy_rule']}. "
         f"{info['n_decoys_used']} pass: {', '.join(info['decoys_used'])}.",
         f"- k per decoy (x4 / v1.0): " + ", ".join(f"{k} {v['x4']}/{v['v1.0']}" for k, v in info['k_by_decoy'].items()) + ".",
         f"- Only {2 * info['n_decoys_used']} distinct (decoy, sign) datasets exist; replicates within one differ "
         "only in permutation draws. MC SEs are therefore given two ways: binomial over 500 (`_se`) and "
         "cluster-robust over configurations (`_se_cluster`, the honest one).", "",
         "## Rejection rate at nominal 5% (level 0 = false-positive rate); ± cluster SE", ""]
    piv = s.pivot(index="test", columns="level_partial_r2", values="rej05")
    pse = s.pivot(index="test", columns="level_partial_r2", values="rej05_se_cluster")
    cols = list(piv.columns)
    L += ["| test | " + " | ".join(f"R²={c:.2%}" for c in cols) + " |", "|---" * (len(cols) + 1) + "|"]
    for t in planted.TESTS:
        L.append(f"| {t} | " + " | ".join(f"{piv.loc[t, c]:.3f} ± {pse.loc[t, c]:.3f}" for c in cols) + " |")
    for title, col in (("FPR at 1% / rejection at nominal 1%", "rej01"), ("Power at matched 5% FPR", "power_matched"),
                       ("95% CI coverage of θ", "coverage")):
        pv = s.pivot(index="test", columns="level_partial_r2", values=col)
        L += ["", f"## {title}", "", "| test | " + " | ".join(f"R²={c:.2%}" for c in cols) + " |",
              "|---" * (len(cols) + 1) + "|"]
        for t in planted.TESTS:
            L.append(f"| {t} | " + " | ".join(f"{pv.loc[t, c]:.3f}" for c in cols) + " |")
    b = s[s.test == "block_x4"].set_index("level_partial_r2")
    L += ["", "## Bias of θ̂ (same estimator for every test)", "", "| level | θ planted | bias | SE (cluster) |",
          "|---|---|---|---|"]
    for lv, row in b.iterrows():
        L.append(f"| {lv:.2%} | {row.theta_mean:.4g} | {row.bias:.3g} | {row.bias_se_cluster:.3g} |")
    z = d[d.level == 0]
    g = z.groupby(["decoy", "sign"]).agg(theta0=("theta0", "first"), n=("rep", "size"),
                                         rej_x4=("p_block_x4", lambda v: (v <= .05).mean()),
                                         rej_v10=("p_block_v1.0", lambda v: (v <= .05).mean()),
                                         rej_naive=("p_naive_ols", lambda v: (v <= .05).mean()))
    L += ["", "## Level 0 by configuration (θ₀ = the decoy's own coefficient in the real M)", "",
          "| decoy | sign | n | θ₀ | rej block_x4 | rej block_v1.0 | rej naive |", "|---|---|---|---|---|---|---|"]
    for (dec, sg), row in g.iterrows():
        L.append(f"| {dec} | {sg:+d} | {row.n} | {row.theta0:.4f} | {row.rej_x4:.2f} | {row.rej_v10:.2f} | {row.rej_naive:.2f} |")
    dl4, dl1 = planted.detection_limit(s, "block_x4"), planted.detection_limit(s, "block_v1.0")
    dl4m = planted.detection_limit(s, "block_x4", "power_matched")
    L += ["", "## Detection limit", "",
          f"- Main test (x4 rule): smallest planted partial R² with ≥ 80% power at nominal 5% = "
          f"**{dl4:.2%}**" if dl4 is not None else "- Main test (x4): 80% power not reached",
          f"  (at matched 5% FPR: {dl4m:.2%})." if dl4m is not None else "",
          f"- v1.0 rule (sensitivity): {dl1:.2%}." if dl1 is not None else "- v1.0: not reached.",
          "", "## Reading", "",
          "- Level 0 is not a clean null: the decoy features are not guaranteed to be unrelated to the real M. "
          "Rejections at level 0 come from two configurations with θ₀ > 0 (GM12878−Panc1 +, MCF-7−Panc1 +); for "
          "the other ten the block tests never reject. Naive OLS and plain permutations reject whenever "
          "|θ₀| is a few thousandths; the block tests need a larger θ₀. With only 12 configurations the FPR "
          "estimates are very uncertain (cluster SE ≈ 0.07).",
          "- The x4 main test has coverage 1.00 at every level (conservative CIs) and loses power relative to "
          "v1.0, as the synthetic benchmark predicted. On real M part of that loss comes from the immobile blocks "
          "(see results/aim2/RPE1_vs_K562_test.md: only ~41% of genes are in movable blocks at k = 100).",
          "- CAMERA (VIF ≈ 17 from ρ̄ ≈ 0.003) is far too conservative here.",
          "- Bias of θ̂ is negligible relative to θ at every level (random signs cancel θ₀).", ""]
    (ROOT / f"results/aim3/planted_{TAG}.md").write_text("\n".join(x for x in L if x is not None))


if __name__ == "__main__":
    aim2_md()
    aim3_md()
