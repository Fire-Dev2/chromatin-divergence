"""Markdown reports for the development pair under plan v1.3 (arm fixed effects in C), with the
plan-v1.2 results alongside.

Run after scripts/run_aim2_dev_v13.py and scripts/run_planted_dev_v13.py:
    PYTHONPATH=src python scripts/report_dev_v13.py
Writes results/aim2/RPE1_vs_K562_test_v13.md, results/aim3/planted_RPE1_vs_K562_v13.md and
results/aim3/planted_RPE1_vs_K562_v13_excl.csv (summary without decoys whose own theta_0 is
significant).
"""
import json
from pathlib import Path

import pandas as pd

from cdtest import planted

ROOT = Path(__file__).resolve().parents[1]
TAG = "RPE1_vs_K562"
LAGS = ["1", "5", "10", "25", "50", "100"]
COMPS = ["block_v1.0", "naive_ols", "unrestricted", "expr_bin", "stratified", "camera",
         "block_coaccess_clusters", "decoy_neg_control"]


def f(v, n=4):
    return "NA" if v is None else (f"{v:.{n}g}" if isinstance(v, float) else str(v))


def ci(v):
    return f"[{v[0]:.4g}, {v[1]:.4g}]"


def aim2_md():
    r = json.loads((ROOT / f"results/aim2/{TAG}_test_v13.json").read_text())
    o = json.loads((ROOT / f"results/aim2/{TAG}_test.json").read_text())
    rep = json.loads((ROOT / "results/aim2/v12_reproduction.json").read_text())
    L = ["# Aim 2 test, development pair RPE1 − K562 — plan v1.3 (arm fixed effects in C)", "",
         "**Status: descriptive** (mixed pair: p53, tissue, growth; not an H1 test). Plan v1.3 adds one "
         "indicator per chromosome arm to f(C) (reference arm dropped); the arm terms enter the model "
         "M ~ f(C) + x, the residualisation of x that sets k and is block-permuted, every comparator and "
         "the decoy similarities. Everything else as plan v1.1/v1.2. Code: `cdtest.aim2.run(arm_effects=True)` "
         "(`scripts/run_aim2_dev_v13.py`). The v1.2 results (`RPE1_vs_K562_test.{json,md}`) are kept "
         "unchanged and shown in brackets.", "",
         f"- Genes {r['n_genes']:,} on {r['n_arms']} arms; knockdowns {r['n_knockdowns']:,}; reliability "
         f"{r['reliability']:.4f}.",
         f"- f(C): {r['design_columns']} columns ({o['design_columns']} of v1.2 + {r['n_arms'] - 1} arm "
         f"indicators), rank {r['design_rank']}.",
         f"- B = {r['B_perm']:,} per test; partial-R² bootstrap {r['n_boot']} knockdown resamples; CAMERA ρ̄ "
         f"{r['rho_bar_used']:.5f}.",
         f"- v1.2 reproduction: `arm_effects=False` reruns reproduce `RPE1_vs_K562_test.json` exactly (max "
         f"|Δ| = {rep['aim2_json_max_abs_diff']:g}, null draws identical) and the planted benchmark "
         f"replicates exactly (`results/aim2/v12_reproduction.json`).", "",
         "## Main test (block permutation, x4 rule) — v1.3 [v1.2]", "",
         "| feature | θ̂ | k x4 | k v1.0 | one-sided p | two-sided p | 95% CI (refined) | partial R² / reliability [95% boot] |",
         "|---|---|---|---|---|---|---|---|"]
    for lab, d in r["features"].items():
        m, dd = d["main_block_x4"], o["features"][lab]
        mo = dd["main_block_x4"]
        L.append(f"| {lab} | {d['theta_hat']:.4g} [{dd['theta_hat']:.4g}] | {d['k_x4']} [{dd['k_x4']}] | "
                 f"{d['k_v1.0']} [{dd['k_v1.0']}] | {m['p_one_sided']:.4g} [{mo['p_one_sided']:.4g}] | "
                 f"{m['p_two_sided']:.4g} [{mo['p_two_sided']:.4g}] | {ci(m['ci95'])} [{ci(mo['ci95'])}] | "
                 f"{d['partial_r2_share_reproducible']:.3%} {ci([100 * x for x in d['partial_r2_share_boot_ci95']])}% "
                 f"[{dd['partial_r2_share_reproducible']:.3%}] |")
    L += ["", "Plan-grid CIs (201 points, as specified): " + "; ".join(
        f"{lab} {ci(d['main_block_x4']['ci95_plan_grid'])}" for lab, d in r["features"].items()) + ".", "",
          "## Block diagnostics — v1.3 [v1.2]", "",
          "Exchangeability = share of genes in blocks with at least one other block of the same size (only "
          "those move, plan §5 step 5). Null mean/SD = mean/SD of the B block-permutation T*.", "",
          "| feature | blocks (x4) | genes movable, x4 | genes movable, v1.0 | null mean T* | null SD T* | observed T | (T − null mean)/null SD |",
          "|---|---|---|---|---|---|---|---|"]
    for lab, d in r["features"].items():
        e, e1, m = d["exchangeability_x4"], d["exchangeability_v1.0"], d["main_block_x4"]
        dd = o["features"][lab]
        eo, e1o, mo = dd["exchangeability_x4"], dd["exchangeability_v1.0"], dd["main_block_x4"]
        z, zo = (m["T"] - m["null_mean"]) / m["null_sd"], (mo["T"] - mo["null_mean"]) / mo["null_sd"]
        L.append(f"| {lab} | {e['n_blocks']} [{eo['n_blocks']}] | {e['frac_genes_movable']:.1%} "
                 f"[{eo['frac_genes_movable']:.1%}] | {e1['frac_genes_movable']:.1%} [{e1o['frac_genes_movable']:.1%}] | "
                 f"{m['null_mean']:.4g} [{mo['null_mean']:.4g}] | {m['null_sd']:.3g} [{mo['null_sd']:.3g}] | "
                 f"{m['T']:.4g} [{mo['T']:.4g}] | {z:.1f} [{zo:.1f}] |")
    L += ["", "Within-arm ACF of the x residuals (after f(C); v1.3 includes the arm terms) — v1.3 [v1.2]:", "",
          "| series | " + " | ".join(f"lag {l}" for l in LAGS) + " |", "|---" * (len(LAGS) + 1) + "|"]
    for lab, d in r["features"].items():
        a, b = d["acf_x_resid_within_arm"], o["features"][lab]["acf_x_resid_within_arm"]
        L.append(f"| x: {lab} | " + " | ".join(f"{a[l]:.3f} [{b[l]:.3f}]" for l in LAGS) + " |")
    a, b = r["acf_M_resid_within_arm_diagnostic"], o["acf_M_resid_within_arm_diagnostic"]
    L.append("| M (diagnostic only) | " + " | ".join(f"{a[l]:.3f} [{b[l]:.3f}]" for l in LAGS) + " |")
    L += ["", "## Sensitivity analyses and comparators, v1.3 (one-sided p / two-sided p [95% CI])", "",
          "| feature | " + " | ".join(COMPS) + " |", "|---" * (len(COMPS) + 1) + "|"]
    for lab, d in r["features"].items():
        cells = []
        for n in COMPS:
            c = d["comparators"][n]
            if n == "decoy_neg_control":
                cells.append(f"p={f(c.get('p_one_sided'))} / {f(c.get('p_two_sided'))} (n={c['n_neg_control']})")
            else:
                cells.append(f"{c['p_one_sided']:.3g} / {c['p_two_sided']:.3g} {ci(c['ci95'])}")
        L.append(f"| {lab} | " + " | ".join(cells) + " |")
    L += ["", "Two-sided p under v1.2 for the same comparators:", "",
          "| feature | " + " | ".join(COMPS) + " |", "|---" * (len(COMPS) + 1) + "|"]
    for lab, d in o["features"].items():
        L.append(f"| {lab} | " + " | ".join(f(d["comparators"][n].get("p_two_sided"), 3) for n in COMPS) + " |")
    c, co = r["features"]["combined_cpm"], o["features"]["combined_cpm"]
    p = r["features"]["promoter_cpm"]
    L += ["", "## Reading", "",
          f"- Arm terms fix the block problem. k for combined CPM falls from {co['k_x4']} to {c['k_x4']} (x4) and "
          f"from {co['k_v1.0']} to {c['k_v1.0']} (v1.0); genes in movable blocks rise from "
          f"{co['exchangeability_x4']['frac_genes_movable']:.0%} to {c['exchangeability_x4']['frac_genes_movable']:.0%}; "
          f"the block null is now centred near 0 (mean {c['main_block_x4']['null_mean']:.4f}, v1.2 "
          f"{co['main_block_x4']['null_mean']:.4f}) and wider (SD {c['main_block_x4']['null_sd']:.4f} vs "
          f"{co['main_block_x4']['null_sd']:.4f}).",
          f"- The estimate barely moves: θ̂ = {c['theta_hat']:.4f} (v1.2 {co['theta_hat']:.4f}); the association is "
          "within arms, not an arm-scale artefact. One-sided p (θ > 0) = "
          f"{c['main_block_x4']['p_one_sided']:.4g}; two-sided {c['main_block_x4']['p_two_sided']:.4g} "
          f"(the minimum at B = 9,999). Partial R² / reliability {c['partial_r2_share_reproducible']:.2%} "
          f"(v1.2 {co['partial_r2_share_reproducible']:.2%}); promoter {p['partial_r2_share_reproducible']:.2%}.",
          "- Sign and pattern as in v1.2: genes more accessible in RPE1 than K562 have a smaller RPE1-minus-K562 "
          "response difference, carried by the promoter feature; distal and FC near zero. All comparators "
          "agree on the sign.",
          "- Not evidence for or against H1 (mixed pair).", ""]
    (ROOT / f"results/aim2/{TAG}_test_v13.md").write_text("\n".join(L))


def _piv(s, col, fmt="{:.3f}", se=None):
    pv = s.pivot(index="test", columns="level_partial_r2", values=col)
    ps = s.pivot(index="test", columns="level_partial_r2", values=se) if se else None
    cols = list(pv.columns)
    L = ["| test | " + " | ".join(f"R²={c:.2%}" for c in cols) + " |", "|---" * (len(cols) + 1) + "|"]
    for t in planted.TESTS:
        L.append(f"| {t} | " + " | ".join(fmt.format(pv.loc[t, c]) + (f" ± {ps.loc[t, c]:.3f}" if se else "")
                                          for c in cols) + " |")
    return L


def aim3_md():
    s = pd.read_csv(ROOT / f"results/aim3/planted_{TAG}_v13.csv")
    so = pd.read_csv(ROOT / f"results/aim3/planted_{TAG}.csv")
    info = json.loads((ROOT / f"results/aim3/planted_{TAG}_v13_info.json").read_text())
    d = pd.read_parquet(ROOT / f"results/aim3/planted_{TAG}_v13_replicates.parquet")
    ex = planted.theta0_significant(info, "p_block_x4_two_sided")
    ex_n = planted.theta0_significant(info, "p_naive_two_sided")
    se_ = planted.summarize(d[~d.decoy.isin(ex)])
    se_n = planted.summarize(d[~d.decoy.isin(ex_n)])
    out = pd.concat([se_.assign(excluded="theta0_block_x4_p2<=0.05:" + ";".join(ex)),
                     se_n.assign(excluded="theta0_naive_p2<=0.05:" + ";".join(ex_n))])
    out.to_csv(ROOT / f"results/aim3/planted_{TAG}_v13_excl.csv", index=False)
    nconf = 2 * info["n_decoys_used"]
    L = ["# Aim 3 planted-signal benchmark on real M, RPE1 − K562 — plan v1.3 (arm fixed effects in C)", "",
         "Code: `cdtest.planted.run(arm_effects=True)` (`scripts/run_planted_dev_v13.py`). Protocol as v1.2 "
         "(plan §6): M' = M + θ·v, v = random sign × a negative-control decoy, θ set so the raw planted partial "
         f"R² given f(C) is the level; {info['reps_per_level']} replicates/level, B = {info['B']}; every test applied "
         "to v. Under v1.3, f(C) includes the arm indicators everywhere: decoy similarities (negative-control "
         "set), k per decoy, θ calibration and all tests. v1.2 results: `planted_RPE1_vs_K562.{csv,md}` "
         "(unchanged).", "",
         f"- Decoys (|similarity| < 0.05 recomputed under v1.3 C): {info['n_decoys_used']} of 21 — "
         f"{', '.join(info['decoys_used'])} (v1.2: 6). {nconf} (decoy, sign) configurations.",
         "- k per decoy (x4 / v1.0): " + ", ".join(f"{k} {v['x4']}/{v['v1.0']}" for k, v in info["k_by_decoy"].items())
         + ".",
         "- MC SEs: binomial over 500 and cluster-robust over (decoy, sign) configurations (± below, the honest one).",
         "", "## Each decoy's own θ₀ in the real M (v1.3 f(C); two-sided)", "",
         "| decoy | θ₀ | k x4 | p block x4 (B = 9,999) | p naive OLS |", "|---|---|---|---|---|"]
    for k, v in info["theta0_tests"].items():
        L.append(f"| {k} | {v['theta0']:.4f} | {v['k_x4']} | {v['p_block_x4_two_sided']:.4g} | {v['p_naive_two_sided']:.3g} |")
    L += ["", f"Significantly nonzero θ₀ by the main test (two-sided p ≤ 0.05): {', '.join(ex) or 'none'}; "
          f"by naive OLS: {', '.join(ex_n) or 'none'}. Exclusion drops both signs of such a decoy.", "",
          "## Rejection at nominal 5% (level 0 = FPR) ± cluster SE — all configurations", ""]
    L += _piv(s, "rej05", se="rej05_se_cluster")
    L += ["", f"## Rejection at nominal 5% ± cluster SE — excluding decoys with significant θ₀ (main test: {', '.join(ex)})", ""]
    L += _piv(se_, "rej05", se="rej05_se_cluster")
    L += ["", f"## Same, excluding by the naive-OLS θ₀ test ({', '.join(ex_n)})", ""]
    L += _piv(se_n, "rej05", se="rej05_se_cluster")
    L += ["", "## Main test vs v1.2 (rejection at nominal 5%, all configurations)", "",
          "| test | plan | " + " | ".join(f"R²={c:.2%}" for c in planted.LEVELS) + " |", "|---" * (len(planted.LEVELS) + 2) + "|"]
    for t in ("block_x4", "block_v1.0"):
        for lab, ss in (("v1.3", s), ("v1.2", so)):
            g = ss[ss.test == t].set_index("level_partial_r2").rej05
            L.append(f"| {t} | {lab} | " + " | ".join(f"{g.loc[c]:.3f}" for c in planted.LEVELS) + " |")
    for title, col in (("Rejection at nominal 1% (all configurations)", "rej01"),
                       ("Power at matched 5% FPR (all configurations)", "power_matched"),
                       ("95% CI coverage of the planted θ (all configurations)", "coverage")):
        L += ["", f"## {title}", ""] + _piv(s, col)
    L += ["", "## 95% CI coverage of the planted θ, excluding decoys with significant θ₀ (main test)", ""] + _piv(se_, "coverage")
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
    L += ["", "## Level 0 by configuration (θ₀ = the signed decoy's own coefficient in the real M)", "",
          "| decoy | sign | n | θ₀ | rej block_x4 | rej block_v1.0 | rej naive |", "|---|---|---|---|---|---|---|"]
    for (dec, sg), row in g.iterrows():
        L.append(f"| {dec} | {sg:+d} | {row.n} | {row.theta0:.4f} | {row.rej_x4:.2f} | {row.rej_v10:.2f} | {row.rej_naive:.2f} |")
    dl = {k: planted.detection_limit(v, "block_x4") for k, v in (("all", s), ("excl", se_), ("excl_naive", se_n))}
    dl1 = {k: planted.detection_limit(v, "block_v1.0") for k, v in (("all", s), ("excl", se_), ("excl_naive", se_n))}
    dlm = planted.detection_limit(s, "block_x4", "power_matched")
    pc = lambda v: "not reached" if v is None else f"{v:.2%}"
    f0 = lambda ss, t="block_x4": ss[(ss.test == t) & (ss.level_partial_r2 == 0)].iloc[0]
    a0, e0, n0 = f0(s), f0(se_), f0(se_n)
    L += ["", "## Headline (main test = block permutation, x4 rule)", "",
          f"- **Detection limit (80% power at nominal 5%): {pc(dl['all'])}** with all configurations; "
          f"{pc(dl['excl'])} excluding decoys with significant θ₀ (main-test criterion); {pc(dl['excl_naive'])} "
          f"with the naive-OLS criterion. At matched 5% FPR: {pc(dlm)}. v1.2: 1.00%.",
          f"  Power at 0.5%: {b.loc[0.005, 'rej05']:.3f} (all) / "
          f"{se_[(se_.test == 'block_x4') & (se_.level_partial_r2 == 0.005)].rej05.iloc[0]:.3f} (excluded); "
          f"v1.2 {so[(so.test == 'block_x4') & (so.level_partial_r2 == 0.005)].rej05.iloc[0]:.3f}.",
          f"- **Level-0 FPR at nominal 5%: {a0.rej05:.3f} ± {a0.rej05_se_cluster:.3f}** (all {nconf} configurations); "
          f"**{e0.rej05:.3f} ± {e0.rej05_se_cluster:.3f}** excluding {', '.join(ex)} "
          f"({e0.n_configs} configurations, n = {e0.n}); {n0.rej05:.3f} with the naive criterion. v1.2: 0.096.",
          f"- v1.0 rule (sensitivity): detection limit {pc(dl1['all'])} (all) / {pc(dl1['excl'])} (excluded); "
          f"level-0 FPR {f0(s, 'block_v1.0').rej05:.3f} / {f0(se_, 'block_v1.0').rej05:.3f}.",
          "", "## Reading", "",
          "- All level-0 rejections of the main test come from configurations whose decoy is itself associated "
          "with the real M with θ₀ > 0 (A549−Panc1 +, IMR-90−PC-3 −). Both decoys have two-sided main-test "
          "p ≈ 0.02 for θ₀ by themselves; without them the main test never rejects at level 0. Level 0 is "
          "therefore not a clean null; the decoys are not guaranteed signal-free.",
          "- Under v1.3 the negative-control set changed (7 decoys: 3 new, 2 of the v1.2 six dropped) and three of them carry |θ₀| ≈ 0.016–0.018, "
          "comparable to θ at the 0.25–0.5% levels. Because every test estimates θ + θ₀, CI coverage of the planted θ "
          "falls (main test 0.71–0.77 over all configurations); excluding the significant-θ₀ decoys the main test "
          "covers 1.00 again. Coverage here measures decoy contamination as much as calibration.",
          "- The main test gains power over v1.2 at low levels (0.25%: 0.65 vs 0.57; 0.5%: 0.78 vs 0.60) but just "
          "misses 80% at 0.5%, so the detection limit stays 1%. Decoy residuals keep long-range dependence even "
          "with arm terms (k x4 = 44–100), unlike the real features (k 20–24).",
          "- Only 14 distinct datasets: cluster SEs (≈ 0.1 at level 0) are the honest uncertainty.",
          "- CAMERA (VIF from ρ̄ ≈ 0.003) remains far too conservative.", ""]
    (ROOT / f"results/aim3/planted_{TAG}_v13.md").write_text("\n".join(L))


if __name__ == "__main__":
    aim2_md()
    aim3_md()
