"""Aim 4 reports and figures (after scripts/run_aim4_aim1.sh and scripts/run_aim4.sh).

Reads the stage outputs under results/aim4/{tests,sens_v12,sens_groupmean,planted,extras} (and the
development-pair files results/aim2/RPE1_vs_K562_test_v13.json, reused), computes the H3 naive-vs-
main gap (knockdown bootstrap) and writes one JSON + MD per analysis, results/aim4/SUMMARY.md and
figures/aim4_*.png.
    PYTHONPATH=src python scripts/report_aim4.py
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from cdtest import aim2, aim4, planted  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
R4 = ROOT / "results/aim4"
FIG = ROOT / "figures"
PAIRS = aim4.PAIRS
TAGS = [aim4.tag(a, b) for a, b, _ in PAIRS]
ROLE = {aim4.tag(a, b): r for a, b, r in PAIRS}
MATCHED = [t for t in TAGS if ROLE[t] == "matched"]
MIXED = [t for t in TAGS if ROLE[t] == "mixed"]
FEATS = ["combined_cpm", "promoter_cpm", "distal_cpm", "combined_fc"]


def jl(p):
    return json.loads(Path(p).read_text())


def jw(p, d):
    Path(p).write_text(json.dumps(d, indent=1, default=aim2._json_default))


def test_json(t):
    p = R4 / f"tests/{t}_test_v13.json"
    return jl(p if p.exists() else ROOT / f"results/aim2/{t}_test_v13.json")


def g(v, n=3):
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return "NA"
    return f"{v:.{n}g}"


def ci(v, n=3):
    return f"[{g(v[0], n)}, {g(v[1], n)}]"


def pc(v, n=3):
    return "NA" if v is None else f"{100 * v:.{n}g}%"


def holm(ps):
    ps = np.asarray(ps, float)
    o = np.argsort(ps)
    adj = np.empty_like(ps)
    run = 0.0
    for i, j in enumerate(o):
        run = max(run, min(1.0, (len(ps) - i) * ps[j]))
        adj[j] = run
    return adj


T = {t: test_json(t) for t in TAGS}
EX = R4 / "extras"
POOL = jl(EX / "pooled_contrast.json")
H4F = jl(EX / "h4_gradient_fits.json")
H4 = pd.read_csv(EX / "h4_similarity_gradient.csv")
P53 = jl(EX / "p53_control.json")
TP = jl(EX / "target_promoter_check.json")
SUMM = {t: jl(ROOT / f"results/aim1/{t}/summary.json") for t in TAGS}
SAME = jl(ROOT / "results/aim1/same_line_baseline.json")


# =========================================================================== 1. Aim 1 per pair
def aim1_report():
    rows = []
    for t in TAGS:
        s = SUMM[t]
        kd = pd.read_csv(ROOT / f"results/aim1/{t}/knockdowns.tsv", sep="\t")
        rows.append({"pair": t, "role": ROLE[t], "A": s["A"], "B": s["B"], "n_genes": s["n_genes"],
                     "n_kd_kept": s["n_kd_kept"], "n_kd_shared": s["n_kd_total"],
                     "n_few_cells": int((~kd.cells_ok).sum()), "n_target_not_expr": int((~kd.target_expr_ok).sum()),
                     "n_weak_repression": int((~kd.repressed).sum()),
                     "split_r_median": s["split_r_median"], "reliability": s["reliability"],
                     "generic_share_A": s["generic_share_median_A"], "generic_share_B": s["generic_share_median_B"],
                     "M_mean": s["M_mean"], "M_sd": s["M_sd"],
                     "reproducible_var_vs_same_line": s["M_sd"] ** 2 * s["reliability"]})
    df = pd.DataFrame(rows)
    sl = SAME
    jw(R4 / "aim1_pairs.json", {"pairs": rows, "same_line_baseline": sl})
    L = ["# Aim 4 step 1 — Aim 1 outcome M for all six pairs (frozen pipeline pipeline-v1.1)", "",
         "Streaming pseudobulk pass (`screens.run`) on Jurkat and HepG2, then `aim1.run_pair(A, B)` with plan §0 "
         "orientation (M = A − B). RPE1 − K562 is the development-pair run (reused). Logs: `results/aim1/screens_*.log`, "
         "`results/aim4/logs/aim1_*.log`.", "",
         "- Jurkat: 262,956 cells, 257,157 kept, 11,892 controls, 2,393 knockdowns, 8,882 genes.",
         "- HepG2: 145,473 cells, 142,813 kept, 4,926 controls, 2,393 knockdowns, 9,624 genes.", "",
         "| pair | role | genes | knockdowns kept / shared | few cells | target not expr. | weak repression | split-half r | reliability | generic share A / B | mean M | SD M |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        L.append(f"| {r['pair']} | {r['role']} | {r['n_genes']:,} | {r['n_kd_kept']:,} / {r['n_kd_shared']:,} | "
                 f"{r['n_few_cells']} | {r['n_target_not_expr']} | {r['n_weak_repression']} | {r['split_r_median']:.3f} | "
                 f"{r['reliability']:.3f} | {r['generic_share_A']:.1%} / {r['generic_share_B']:.1%} | "
                 f"{r['M_mean']:.4f} | {r['M_sd']:.4f} |")
    sd_sl = sl.get("same_line", sl).get("M_sd", None) if isinstance(sl, dict) else None
    L += ["", "Same-line baseline (K562 essential − K562 genome-wide; `results/aim1/same_line_baseline.json`): "
          "SD M 0.0095, reliability 0.72, generic share 2.5% / 2.9%.", "",
          "Reproducible variance of M (SD² × reliability) relative to the same-line baseline (0.0095² × 0.72 = 6.5e-5):", ""]
    for r in rows:
        L.append(f"- {r['pair']}: {r['reproducible_var_vs_same_line']:.2e} ({r['reproducible_var_vs_same_line'] / (0.0095 ** 2 * 0.72):.1f}× baseline)")
    L += ["", "Notes: K562 − Jurkat has the smallest and least reliable M (SD 0.017, reliability 0.86), about 4× "
          "the same-line baseline in reproducible variance; the other pairs are 40–150×. The generic share is "
          "below 13% in every line (largest in RPE1, 12–13%)."]
    (R4 / "aim1_pairs.md").write_text("\n".join(L) + "\n")
    return df


# =========================================================================== 2. H1 primary
def planted_summary(t):
    s = pd.read_csv(R4 / f"planted/planted_{t}_v13.csv")
    e = pd.read_csv(R4 / f"planted/planted_{t}_v13_excl.csv")
    info = jl(R4 / f"planted/planted_{t}_v13_info.json")
    ex_b = e[e.excluded.str.startswith("theta0_p_block")]
    ex_n = e[e.excluded.str.startswith("theta0_p_naive")]
    out = {"decoy_rule": info["decoy_rule"], "decoys_used": info["decoys_used"],
           "theta0_tests": info["theta0_tests"],
           "excluded_main": planted.theta0_significant(info, "p_block_x4_two_sided"),
           "excluded_naive": planted.theta0_significant(info, "p_naive_two_sided"),
           "detection_limit_all": planted.detection_limit(s, "block_x4"),
           "detection_limit_excl_main": planted.detection_limit(ex_b, "block_x4"),
           "detection_limit_excl_naive": planted.detection_limit(ex_n, "block_x4"),
           "detection_limit_matched_fpr": planted.detection_limit(s, "block_x4", col="power_matched"),
           "detection_limit_v1.0_all": planted.detection_limit(s, "block_v1.0"),
           "detection_limit_v1.0_excl_main": planted.detection_limit(ex_b, "block_v1.0")}
    for name, d in (("all", s), ("excl_main", ex_b), ("excl_naive", ex_n)):
        tab = {}
        for test in planted.TESTS:
            q = d[d.test == test].sort_values("level_partial_r2")
            tab[test] = {f"{lv:.4f}": {"rej05": float(r.rej05), "se_cluster": float(r.rej05_se_cluster),
                                         "rej01": float(r.rej01), "coverage": float(r.coverage),
                                         "power_matched": float(r.power_matched), "n": int(r.n),
                                         "n_configs": int(r.n_configs), "bias": float(r.bias)}
                         for lv, r in zip(q.level_partial_r2, q.itertuples())}
        out[f"table_{name}"] = tab
    return out


def planted_md(t, ps):
    L = [f"# Planted-signal benchmark (plan §6) on matched pair {t} — v1.3", "",
         "Code `planted.run(arm_effects=True)`, 500 replicates/level, B = 999, θ₀ tests B = 9,999 "
         "(`scripts/run_aim4.py planted`). Decoys: " + ps["decoy_rule"] + f" under this pair's C: "
         f"{len(ps['decoys_used'])} of 21 — " + ", ".join(ps["decoys_used"]) + ".", "",
         "## Each decoy's own θ₀ in the real M (two-sided)", "", "| decoy | θ₀ | k x4 | p block x4 | p naive |",
         "|---|---|---|---|---|"]
    for d, r in ps["theta0_tests"].items():
        L.append(f"| {d} | {r['theta0']:.4f} | {r['k_x4']} | {r['p_block_x4_two_sided']:.4g} | {r['p_naive_two_sided']:.3g} |")
    for name, lab in (("all", "all configurations"),
                      ("excl_main", "excluding decoys with significant θ₀ (main test: " + ", ".join(ps["excluded_main"]) + ")"),
                      ("excl_naive", "excluding by naive θ₀ test (" + ", ".join(ps["excluded_naive"]) + ")")):
        tab = ps[f"table_{name}"]
        lv = list(next(iter(tab.values())).keys())
        L += ["", f"## Rejection at nominal 5% ± cluster SE — {lab}", "",
              "| test | " + " | ".join(f"R²={float(x):.2%}" for x in lv) + " |", "|---" * (len(lv) + 1) + "|"]
        for test, row in tab.items():
            L.append(f"| {test} | " + " | ".join(f"{row[x]['rej05']:.3f} ± {row[x]['se_cluster']:.3f}" for x in lv) + " |")
        if name == "all":
            L += ["", "Main-test CI coverage of planted θ: " + ", ".join(
                f"{float(x):.2%}: {tab['block_x4'][x]['coverage']:.3f}" for x in lv) +
                  "; power at matched 5% FPR: " + ", ".join(f"{float(x):.2%}: {tab['block_x4'][x]['power_matched']:.3f}" for x in lv)]
    L += ["", "## Headline", "",
          f"- **Detection limit (main test, 80% power at nominal 5%): {pc(ps['detection_limit_all'])}** with all "
          f"configurations; {pc(ps['detection_limit_excl_main'])} excluding decoys with significant θ₀ (main test); "
          f"{pc(ps['detection_limit_excl_naive'])} excluding by naive θ₀; {pc(ps['detection_limit_matched_fpr'])} at matched FPR.",
          f"- Level-0 FPR (main test): {ps['table_all']['block_x4']['0.0000']['rej05']:.3f} ± "
          f"{ps['table_all']['block_x4']['0.0000']['se_cluster']:.3f} (all); "
          f"{ps['table_excl_main']['block_x4']['0.0000']['rej05']:.3f} excluding significant-θ₀ decoys "
          f"(n configs {ps['table_excl_main']['block_x4']['0.0000']['n_configs']}).",
          f"- v1.0 rule: detection limit {pc(ps['detection_limit_v1.0_all'])} (all) / {pc(ps['detection_limit_v1.0_excl_main'])} (excluded)."]
    (R4 / f"planted_{t}.md").write_text("\n".join(L) + "\n")


def _p_less(t, label="combined_cpm"):
    """p for theta < 0 from the stored main-test null draws (same B = 9,999 draws)."""
    from cdtest import crt
    d = ROOT / "results/aim4/tests" / f"{t}_null_v13.npz"
    d = d if d.exists() else ROOT / f"results/aim2/{t}_null_v13.npz"
    Ts = np.load(d)[label]
    return crt.pvalue(T[t]["features"][label]["theta_hat"], Ts, "less")


def h1_report():
    ps = {t: planted_summary(t) for t in MATCHED}
    for t in MATCHED:
        planted_md(t, ps[t])
        jw(R4 / f"planted_{t}.json", ps[t])
    p1 = [T[t]["features"]["combined_cpm"]["main_block_x4"]["p_one_sided"] for t in MATCHED]
    adj = holm(p1)
    out = {"family": "H1, combined CPM feature, two matched pairs, one-sided theta > 0, Holm at FWER 0.05",
           "B": 9999, "pairs": {}}
    for t, pa in zip(MATCHED, adj):
        f = T[t]["features"]["combined_cpm"]
        m = f["main_block_x4"]
        out["pairs"][t] = {"theta_hat": f["theta_hat"], "k_x4": f["k_x4"], "p_one_sided": m["p_one_sided"],
                           "p_two_sided": m["p_two_sided"], "p_holm": float(pa),
                           "p_less": float(_p_less(t)), "p_two_sided_equal_tailed": float(min(1, 2 * min(m["p_one_sided"], _p_less(t)))), "ci95_refined": m["ci95"],
                           "ci95_plan_grid": m["ci95_plan_grid"], "null_mean": m["null_mean"], "null_sd": m["null_sd"],
                           "partial_r2_raw": f["partial_r2_raw"], "reliability": f["reliability"],
                           "partial_r2_share_reproducible": f["partial_r2_share_reproducible"],
                           "partial_r2_share_boot_ci95": f["partial_r2_share_boot_ci95"],
                           "theta_per_sd_x": f["theta_hat"] * f["sd_x"], "sd_x": f["sd_x"],
                           "exchangeability_x4": f["exchangeability_x4"],
                           "detection_limit_all": ps[t]["detection_limit_all"],
                           "detection_limit_excl_sig_theta0": ps[t]["detection_limit_excl_main"],
                           "detection_limit_excl_naive_theta0": ps[t]["detection_limit_excl_naive"],
                           "rejects_at_holm_0.05": bool(pa <= 0.05)}
    sup = [out["pairs"][t]["rejects_at_holm_0.05"] and out["pairs"][t]["theta_hat"] > 0 for t in MATCHED]
    out["verdict"] = ("supported in both matched pairs" if all(sup) else
                      "supported in one matched pair only" if any(sup) else
                      "not supported: no matched pair shows a positive effect at 5% (Holm) under the main test")
    jw(R4 / "h1_primary.json", out)
    L = ["# H1 (primary) — matched pairs, combined CPM feature, main test", "",
         "Run once on the frozen pipeline (tag `pipeline-v1.1`, plan v1.3): block permutation, x4 rule, chromosome-arm "
         "terms in f(C), B = 9,999, one-sided θ > 0, Holm across the two matched pairs. 95% CI by test inversion "
         "(refined grid inside the plan span; plan-grid interval also given). Partial R² / reliability with 500-resample "
         "knockdown bootstrap. Detection limit from the planted-signal benchmark on each pair (`planted_<pair>.md`).", "",
         "| pair (A − B) | θ̂ | k | one-sided p (θ > 0) | two-sided p [equal-tailed] | Holm p | 95% CI (refined) | plan-grid CI | partial R²/rel [boot 95%] | detection limit (all / excl. sig. θ₀ decoys) |",
         "|---|---|---|---|---|---|---|---|---|---|"]
    for t in MATCHED:
        r = out["pairs"][t]
        L.append(f"| {t} | {r['theta_hat']:.4g} | {r['k_x4']} | {r['p_one_sided']:.4g} | {r['p_two_sided']:.4g} [{r['p_two_sided_equal_tailed']:.4g}] | "
                 f"{r['p_holm']:.4g} | {ci(r['ci95_refined'], 4)} | {ci(r['ci95_plan_grid'], 3)} | "
                 f"{r['partial_r2_share_reproducible']:.3%} [{r['partial_r2_share_boot_ci95'][0]:.3%}, "
                 f"{r['partial_r2_share_boot_ci95'][1]:.3%}] | {pc(r['detection_limit_all'])} / "
                 f"{pc(r['detection_limit_excl_sig_theta0'])} |")
    L += ["", f"**Verdict: H1 {out['verdict']}.**", "",
          "Proposal 'Supported if': the effect of the signed accessibility difference on the signed response difference "
          "is positive beyond expression, at 5% under the main test. Neither matched pair meets it.", "",
          "- **K562 − Jurkat:** θ̂ ≈ 0 (2e-4; 95% CI [−0.0021, 0.0014]; partial R²/rel 0.001%, bootstrap upper bound "
          "0.023%). The planted benchmark on this pair detects 0.5% partial R² with ≥ 80% power (level-0 FPR 0.08 all "
          "decoys, 0.00 without the two decoys whose own θ₀ is significant), so an effect of the size the study can "
          "see is absent. 'If not' column: accessibility adds nothing detectable here once p53, tissue and growth mode "
          "are held constant.",
          "- **RPE1 − HepG2:** a clear association of the *opposite* sign (θ̂ = −0.017, CI [−0.026, −0.014], two-sided "
          "p 1e-4 at the B = 9,999 floor, partial R²/rel 0.77%): genes more accessible in one line respond *less* in "
          "that line. The one-sided H1 test cannot reject (p = 1). The p53 control shows RPE1 − HepG2 still diverges at "
          "p53 targets (z = 3.3), so p53 activity is not held constant in this pair. The planted benchmark here is "
          "contaminated by decoys with their own θ₀ (4 of 6 significant): detection limit 2% with all decoys, 0.5% "
          "with the 2 clean decoys (4 configurations only).",
          "- Two-sided p in the frozen code is (1 + #|T*| ≥ |T|)/(1 + B); the equal-tailed version 2·min(p>, p<) is "
          "given in brackets. They differ when the block null is not centred at 0."]
    (R4 / "h1_primary.md").write_text("\n".join(L) + "\n")
    return out, ps


# =========================================================================== 3. H2
def h2_report():
    out = {}
    L = ["# H2 (secondary) — promoter vs distal CPM, matched pairs", "",
         "Main test per feature (v1.3, B = 9,999). θ is in M units per unit of the feature difference; the two "
         "features are on different scales, so the scale-free comparison is θ per SD of x and partial R² / "
         "reliability (500-resample knockdown bootstrap). Supported only if distal is larger with non-overlapping "
         "intervals.", "",
         "| pair | feature | θ̂ | 95% CI | θ̂ × SD(x) [CI × SD(x)] | one-/two-sided p | partial R²/rel [boot 95%] |",
         "|---|---|---|---|---|---|---|"]
    for t in TAGS:
        out[t] = {}
        for lab in ("promoter_cpm", "distal_cpm"):
            f = T[t]["features"][lab]
            m = f["main_block_x4"]
            sd = f["sd_x"]
            out[t][lab] = {"theta_hat": f["theta_hat"], "ci95": m["ci95"], "sd_x": sd,
                           "theta_per_sd": f["theta_hat"] * sd, "ci95_per_sd": [m["ci95"][0] * sd, m["ci95"][1] * sd],
                           "p_one_sided": m["p_one_sided"], "p_two_sided": m["p_two_sided"],
                           "pr2_rel": f["partial_r2_share_reproducible"], "pr2_rel_boot_ci95": f["partial_r2_share_boot_ci95"]}
            if t in MATCHED:
                o = out[t][lab]
                L.append(f"| {t} | {lab} | {o['theta_hat']:.4g} | {ci(o['ci95'], 3)} | {o['theta_per_sd']:.4g} "
                         f"{ci(o['ci95_per_sd'], 3)} | {o['p_one_sided']:.3g} / {o['p_two_sided']:.3g} | "
                         f"{o['pr2_rel']:.3%} [{o['pr2_rel_boot_ci95'][0]:.3%}, {o['pr2_rel_boot_ci95'][1]:.3%}] |")
        p, d = out[t]["promoter_cpm"], out[t]["distal_cpm"]
        out[t]["distal_larger_nonoverlap_pr2"] = bool(d["pr2_rel_boot_ci95"][0] > p["pr2_rel_boot_ci95"][1])
        out[t]["distal_larger_nonoverlap_theta_per_sd"] = bool(
            d["theta_per_sd"] > 0 and d["ci95_per_sd"][0] > max(p["ci95_per_sd"][1], 0))
        out[t]["promoter_larger_nonoverlap_pr2"] = bool(p["pr2_rel_boot_ci95"][0] > d["pr2_rel_boot_ci95"][1])
    out["verdict"] = {t: ("supported" if out[t]["distal_larger_nonoverlap_pr2"] and out[t]["distal_larger_nonoverlap_theta_per_sd"]
                          else "not supported") for t in MATCHED}
    L += ["", "Mixed pairs (descriptive):", "", "| pair | promoter θ̂ [CI] p₂ | distal θ̂ [CI] p₂ | pR²/rel promoter / distal |", "|---|---|---|---|"]
    for t in MIXED:
        p, d = out[t]["promoter_cpm"], out[t]["distal_cpm"]
        L.append(f"| {t} | {p['theta_hat']:.4g} {ci(p['ci95'])} {p['p_two_sided']:.3g} | {d['theta_hat']:.4g} "
                 f"{ci(d['ci95'])} {d['p_two_sided']:.3g} | {p['pr2_rel']:.3%} / {d['pr2_rel']:.3%} |")
    L += ["", "**Verdict (matched pairs):** " + "; ".join(f"{t}: {v}" for t, v in out["verdict"].items()) + ".", "",
          "- K562 − Jurkat: distal θ̂ > 0 (one-sided p 0.021) and promoter slightly negative, but the partial-R² "
          "bootstrap intervals overlap (distal 0.092% [0.036, 0.169] vs promoter 0.041% [0.008, 0.108]) and both "
          "shares are below 0.1%.",
          "- RPE1 − HepG2: both features negative (opposite to H1); promoter carries the larger share "
          "(0.59% vs 0.23%, non-overlapping intervals) — the reverse of H2.",
          "- Across all six pairs: the promoter feature carries a negative association in every RPE1-containing "
          "pair; in the three pairs without RPE1 the distal feature is weakly positive (K562 − Jurkat p₂ 0.028, "
          "K562 − HepG2 0.029, Jurkat − HepG2 0.12) and the promoter ≈ 0. Descriptive (secondary)."]
    jw(R4 / "h2_location.json", out)
    (R4 / "h2_location.md").write_text("\n".join(L) + "\n")
    return out


# =========================================================================== 4. H3
def h3_report():
    out = {}
    for a, b, _ in PAIRS:
        t = aim4.tag(a, b)
        P = aim2.load_pair(a, b, arm_effects=True)
        naive7 = aim2.naive_analysis(P)["combined_cpm"]
        Mb, relb, _ = aim2.bootstrap_reliability_M(t, 500)
        x = P["X"]["combined_cpm"]
        Q = P["Q"]
        rx = x - Q @ (Q.T @ x)
        th_main_b = ((Mb - (Mb @ Q) @ Q.T) @ rx) / (rx @ rx)
        Z7 = np.column_stack([np.ones(len(x)), P["C"].expr_diff.to_numpy()])
        import cdtest.crt as crt
        Q7 = crt.orthonormal_basis(Z7)
        rx7 = x - Q7 @ (Q7.T @ x)
        th_n7_b = ((Mb - (Mb @ Q7) @ Q7.T) @ rx7) / (rx7 @ rx7)
        f = T[t]["features"]["combined_cpm"]
        gap = naive7["coef"] - f["theta_hat"]
        gb = th_n7_b - th_main_b
        pg = float(min(1.0, 2 * min((gb <= 0).mean(), (gb >= 0).mean())))
        comps = f["comparators"]
        m = f["main_block_x4"]
        nols = comps["naive_ols"]
        out[t] = {"role": ROLE[t], "theta_main": f["theta_hat"], "p_main_one": m["p_one_sided"],
                  "p_main_two": m["p_two_sided"], "ci_main": m["ci95"], "null_sd_main": m["null_sd"],
                  "naive_step7": naive7,
                  "gap_naive7_minus_main": gap, "gap_boot_ci95": [float(np.quantile(gb, .025)), float(np.quantile(gb, .975))],
                  "gap_boot_p_two_sided": pg,
                  "naive_ols_fC": {"theta": nols["T"], "se": nols["se"], "p_two": nols["p_two_sided"], "ci": nols["ci95"]},
                  "se_ratio_naive_vs_block_null_sd": nols["se"] / m["null_sd"],
                  "comparators_p_two_sided": {k: comps[k]["p_two_sided"] for k in
                                              ("naive_ols", "unrestricted", "expr_bin", "stratified", "camera", "block_v1.0",
                                               "block_coaccess_clusters", "decoy_neg_control")},
                  "comparators_p_one_sided": {k: comps[k]["p_one_sided"] for k in
                                              ("naive_ols", "unrestricted", "expr_bin", "stratified", "camera", "block_v1.0",
                                               "block_coaccess_clusters", "decoy_neg_control")},
                  "comparators_ci95": {k: comps[k].get("ci95") for k in
                                       ("naive_ols", "unrestricted", "expr_bin", "stratified", "camera", "block_v1.0",
                                        "block_coaccess_clusters")}}
    jw(R4 / "h3_comparators.json", out)
    L = ["# H3 (primary in the proposal) — naive tests vs the main test on real data, all six pairs", "",
         "Combined CPM feature. Comparators as in plan §5 (all with the v1.3 f(C)); 'naive step 7' is the "
         "Aim 1 naive analysis M ~ 1 + expr_diff + x (OLS, genes independent, no f(C)). The naive-vs-main gap = "
         "θ̂(naive step 7) − θ̂(main) with a 500-resample knockdown bootstrap (β̃, s̃² held fixed; same "
         "approximation as the partial-R² bootstrap). The OLS-with-f(C) comparator has the same θ̂ as the main "
         "test; its inflation is in the p-value (SE ratio = naive SE / SD of the block-permutation null).", "",
         "| pair | role | θ̂ main | p₂ main | θ̂ naive step 7 | p₂ naive step 7 | gap [boot 95%] p | naive SE / block-null SD |",
         "|---|---|---|---|---|---|---|---|"]
    for t, r in out.items():
        n7 = r["naive_step7"]
        L.append(f"| {t} | {r['role']} | {r['theta_main']:.4g} | {r['p_main_two']:.3g} | {n7['coef']:.4g} | "
                 f"{n7['p_two_sided']:.3g} | {r['gap_naive7_minus_main']:.4g} {ci(r['gap_boot_ci95'])} p={('<0.002' if r['gap_boot_p_two_sided'] == 0 else format(r['gap_boot_p_two_sided'], '.3g'))} | "
                 f"{r['se_ratio_naive_vs_block_null_sd']:.2f} |")
    L += ["", "Two-sided p by comparator (combined CPM):", "",
          "| pair | main (block x4) | naive OLS | unrestricted | expr-bin | stratified | CAMERA | block v1.0 | co-access. clusters | decoy (n) |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for t, r in out.items():
        c = r["comparators_p_two_sided"]
        nd = T[t]["features"]["combined_cpm"]["comparators"]["decoy_neg_control"].get("n_neg_control")
        L.append(f"| {t} | {r['p_main_two']:.3g} | {c['naive_ols']:.3g} | {c['unrestricted']:.3g} | {c['expr_bin']:.3g} | "
                 f"{c['stratified']:.3g} | {c['camera']:.3g} | {c['block_v1.0']:.3g} | {c['block_coaccess_clusters']:.3g} | "
                 f"{g(c['decoy_neg_control'])} ({nd}) |")
    L += ["", "Reading:",
          "- Simulation half of H3 (Aim 3, unchanged): naive tests reject a true null well above 5% under shared "
          "structure; on the matched-pair planted benchmarks the naive OLS level-0 FPR is " +
          "; ".join(f"{t} {jl(R4 / f'planted_{t}.json')['table_all']['naive_ols']['0.0000']['rej05']:.2f}" for t in MATCHED)
          + " vs main test " + "; ".join(f"{jl(R4 / f'planted_{t}.json')['table_all']['block_x4']['0.0000']['rej05']:.2f}" for t in MATCHED) + ".",
          "- Real-data half: the naive step-7 estimate is *not* larger than the main-test estimate. In absolute value "
          "it is smaller in 4 of 6 pairs (RPE1 − HepG2 −0.0097 vs −0.0173; Jurkat − RPE1 −0.0129 vs −0.0219; "
          "RPE1 − K562 −0.0194 vs −0.0224; Jurkat − HepG2 0.0013 vs 0.0024), has the opposite sign in K562 − HepG2 "
          "(−0.0026 vs +0.0071) and is equal within noise in K562 − Jurkat. With the full f(C) the naive OLS SE is "
          "within 3–18% of the SD of the block-permutation null, so naive OLS p-values are smaller mainly because "
          "a t-test has no 1/(B+1) floor; the only pair where naive-type tests reject at 5% and the main test does "
          "not is K562 − HepG2 (mixed; naive OLS p₂ 0.022, unrestricted 0.008, expr-bin 0.010, stratified 0.016 vs "
          "main 0.055). CAMERA is non-significant everywhere (ρ̄ 0.0025–0.0035, VIF ≈ 15–20, SE × ~4).",
          "- H3 verdict: the simulation half holds (Aim 3; and on these pairs' planted benchmarks naive OLS level-0 "
          "FPR is 0.27–0.31 even after excluding decoys with significant θ₀, vs 0.00 for the main test); the "
          "real-data half ('a larger real-data effect than the main test') is not supported."]
    (R4 / "h3_comparators.md").write_text("\n".join(L) + "\n")
    return out


# =========================================================================== 5. mixed pooled
def mixed_report():
    pm, pmat, gap, six, pu = POOL["mixed"], POOL["matched"], POOL["gap"], POOL["all6"], POOL["mixed_union_rank_coupling"]
    L = ["# Mixed pairs — pooled contrast (plan §7) and matched-vs-mixed gap", "",
         "Statistic: mean of the four mixed-pair θ̂ (combined CPM; RPE1 − K562, K562 − HepG2, Jurkat − RPE1, "
         "Jurkat − HepG2; each pair's own signed x = A − B and M = A − B, so under the effect model every θ is "
         "positive regardless of orientation).", "",
         "**Shared permutation (documented choice).** One block permutation is applied identically to all four "
         "pairs. Exact sharing requires the same genes, so the statistic is computed on the genes kept in all four "
         f"pairs ({pm['n_common_genes']:,} genes = genes kept in all four lines; the same set for the two matched "
         "pairs and for all six). Genes in genome order, blocks cut within arms with k = the largest of the four "
         f"pairs' x4 k on these genes (k = {pm['k_shared']}; {pm['exchangeability']['frac_genes_movable']:.1%} of genes "
         "in movable blocks), one permutation of equal-size blocks per draw applied to every pair's x residuals "
         "(each pair keeps its own f(C) with arm terms), B = 9,999. CI by inverting the pooled statistic with one "
         "common θ. The union-order alternative requested (each pair on its own genes, blocks on the union order) "
         "cannot move the same gene to the same gene in every pair; an implementation that couples the pairs "
         "through a shared ordering of union blocks gave null correlations ≈ 0 even for a 90% subset of one pair "
         "against itself, i.e. it behaves like independent permutations. It is reported as secondary.", "",
         "| contrast | genes | θ̂ pooled | one-sided p | two-sided p | 95% CI | null mean / SD |", "|---|---|---|---|---|---|---|"]
    for lab, d in (("mixed (4 pairs), primary", pm), ("matched (2 pairs)", pmat), ("all six", six)):
        L.append(f"| {lab} | {d['n_common_genes']:,} | {d['T_pooled']:.4g} | {d['p_one_sided']:.4g} | {d['p_two_sided']:.4g} | "
                 f"{ci(d['ci95'], 3)} | {d['null_mean']:.3g} / {d['null_sd']:.3g} |")
    L.append(f"| mixed, union-rank coupling (secondary; full gene sets) | {pu['n_union_genes']:,} union | {pu['T_pooled']:.4g} | "
             f"{pu['p_one_sided']:.4g} | {pu['p_two_sided']:.4g} | {ci(pu['ci95'], 3)} | {pu['null_mean']:.3g} / {pu['null_sd']:.3g} |")
    L += ["", "Per-pair θ̂ on the common genes (shared draws) vs on the pair's full gene set (own main test):", "",
          "| pair | θ̂ common genes | p₂ shared draws | θ̂ full set | p₂ own main test |", "|---|---|---|---|---|"]
    for t in MATCHED + MIXED:
        r = six["per_pair"][t]
        m = T[t]["features"]["combined_cpm"]
        L.append(f"| {t} | {r['T_common_genes']:.4g} | {r['p_two_sided']:.3g} | {m['theta_hat']:.4g} | {m['main_block_x4']['p_two_sided']:.3g} |")
    L += ["", "Null correlation of the four mixed-pair T* under the shared permutation: " + str(pm["null_corr_between_pairs"]) + ".", "",
          f"**Gap (mixed − matched)** on the common genes, one permutation shared by all six pairs: "
          f"{gap['gap_mixed_minus_matched']:.4g} (mixed {gap['T_mixed_mean_common6']:.4g}, matched {gap['T_matched_mean_common6']:.4g}); "
          f"null SD {gap['null_sd_gap']:.3g}; approximate 95% interval {ci(gap['approx_ci95'], 3)}; two-sided p (joint "
          f"null) {gap['p_two_sided_joint_null']:.3g}.", "",
          "Reading: the mixed pooled θ̂ is negative (opposite to H1) and clearly nonzero; the matched pooled θ̂ is "
          "also negative, driven by RPE1 − HepG2. The proposal expected a larger positive effect in the mixed pairs; "
          "instead both are negative and the gap is modest and not clearly different from 0. Within the mixed pairs the "
          "sign splits by line: the two RPE1-containing pairs are strongly negative, K562 − HepG2 is positive "
          "(p₂ ≈ 0.05), Jurkat − HepG2 ≈ 0."]
    jw(R4 / "mixed_pooled.json", {"mixed": pm, "matched": pmat, "all6": six, "gap": gap, "mixed_union_rank_coupling": pu})
    (R4 / "mixed_pooled.md").write_text("\n".join(L) + "\n")


# =========================================================================== 6. H4
def h4_report():
    fits, neg = H4F["fits"], H4F["decoy_negative_control"]
    L = ["# H4 (secondary) — similarity gradient (plan §7)", "",
         "For each pair: θ̂ (M ~ f(C) + v, v3 f(C) with arm terms) for every available combined-CPM feature "
         "difference v — the pair's own x, the other five pairs' x (each in its plan orientation) and the 21 "
         "decoy-line differences — against similarity = Pearson correlation of v and the pair's own x, both "
         "residualised on C. OLS of θ̂ on similarity per pair (n = 27) and pooled over pairs (n = 162). H4 "
         "'supported if' effect increases with similarity across other pairs and decoys and is near zero at low "
         "similarity: slope > 0 and intercept ≈ 0. 'θ_std' = θ̂ × SD(r_v)/SD(r_x), which equals θ_own × similarity "
         "exactly under the effect model (diagnostic).", "",
         "| pair | slope [95% CI] p | intercept [95% CI] p | R² | slope excl. own [p] | slope θ_std [p] |", "|---|---|---|---|---|---|"]
    for t, f in fits.items():
        a, e, s = f["all"], f["excluding_own"], f["standardised_theta"]
        L.append(f"| {t} | {a['slope']:.4g} {ci(a['slope_ci95'])} {a['slope_p']:.3g} | {a['intercept']:.4g} "
                 f"{ci(a['intercept_ci95'])} {a['intercept_p']:.3g} | {a['r2']:.2f} | {e['slope']:.4g} [{e['slope_p']:.3g}] | "
                 f"{s['slope']:.4g} [{s['slope_p']:.3g}] |")
    L += ["", "Pooled OLS treats the 162 points as independent (they are not: features and Ms are shared across "
          "pairs); descriptive.", "", "## Decoy negative control (|similarity| < 0.05; each decoy's θ̂ vs its own block-permutation null, x4, B = 9,999, two-sided)", "",
          "| pair | n decoys | n with p ≤ 0.05 | decoys significant (θ̂) |", "|---|---|---|---|"]
    negsum = {}
    for t, d in neg.items():
        sig = {k: v for k, v in d.items() if v["p_block_x4_two_sided"] <= 0.05}
        negsum[t] = {"n": len(d), "n_sig": len(sig)}
        L.append(f"| {t} | {len(d)} | {len(sig)} | " + ", ".join(f"{k} ({v['theta0']:.4f})" for k, v in sig.items()) + " |")
    tot = sum(v["n"] for v in negsum.values()); ts = sum(v["n_sig"] for v in negsum.values())
    L += ["", f"Overall {ts} of {tot} low-similarity decoy tests reject at 5% (expected ≈ {0.05 * tot:.1f} if decoys "
          "carried no signal and the null were exact). Low similarity to the pair's x does not make a feature "
          "signal-free in M.", "",
          "Reading: the pooled slope is negative (θ̂ falls as similarity rises) and the intercept is ≈ 0. A negative "
          "slope is what the effect model predicts when the pair's own θ is negative: the gradient follows the own "
          "effect, whose sign is opposite to H1. H4 as stated (effect increases with similarity) is not supported; "
          "the near-zero intercept says the features do not carry a generic, similarity-independent association."]
    jw(R4 / "h4_gradient.json", {"fits": fits, "decoy_negative_control": neg, "decoy_negative_summary": negsum,
                                  "points": H4.to_dict(orient="records")})
    (R4 / "h4_gradient.md").write_text("\n".join(L) + "\n")
    return negsum


# =========================================================================== 7. p53
def p53_report():
    L = ["# p53 control — all six pairs", "",
         "Aim 1 p53 gate statistic: mean M over the Fischer 2017 core targets kept for the pair vs 1,000 random gene "
         "sets matched on expression deciles in both lines (`aim1.p53_gate` algorithm and seed; RPE1 − K562 "
         "reproduces the stored gate exactly). Oriented so that positive = more response in the p53-WT line "
         "(K562 and Jurkat p53-mutant/null; RPE1 and HepG2 WT). Matched pairs keep A − B and are read two-sided. "
         "Expectation (proposal): divergence in the four mixed pairs, none in the two matched pairs.", "",
         "| pair | p53 differs | n targets | stat (oriented) | null mean | z | one-sided p (oriented) | two-sided p | as expected? |",
         "|---|---|---|---|---|---|---|---|---|"]
    for t, r in P53.items():
        exp = (r["p_one_sided_oriented"] < 0.05) if r["p53_differs"] else (r["p_two_sided"] >= 0.05)
        r["as_expected"] = bool(exp)
        L.append(f"| {t} | {r['p53_differs']} | {r['n_p53_targets']} | {r['stat_oriented']:.4f} | {r['null_mean_oriented']:.4f} | "
                 f"{r['z_oriented']:.2f} | {r['p_one_sided_oriented']:.3g} | {r['p_two_sided']:.3g} | {'yes' if exp else 'NO'} |")
    L += ["", "Reading: all four mixed pairs diverge at p53 targets in the expected direction (z 4–7). K562 − Jurkat "
          "(both p53-deficient) shows none, as expected. **RPE1 − HepG2 (both p53 WT) diverges strongly (z = 3.3, "
          "RPE1 > HepG2)** — the matched pair is not matched in p53 pathway activity, so the RPE1 − HepG2 "
          "comparison does not fully hold the p53 axis constant."]
    jw(R4 / "p53_control.json", P53)
    (R4 / "p53_control.md").write_text("\n".join(L) + "\n")


# =========================================================================== 8. sensitivity
def sensitivity_report():
    out = {}
    for t in MATCHED:
        f = T[t]["features"]["combined_cpm"]
        c = f["comparators"]
        v12 = jl(R4 / f"sens_v12/{t}_test.json")["features"]["combined_cpm"]
        gm = jl(R4 / f"sens_groupmean/{t}_groupmean_test.json")
        fc = T[t]["features"]["combined_fc"]
        rows = {
            "main (v1.3, x4, CPM)": (f["theta_hat"], f["main_block_x4"]["p_one_sided"], f["main_block_x4"]["p_two_sided"], f["main_block_x4"]["ci95"]),
            "v1.2 (no arm terms)": (v12["theta_hat"], v12["main_block_x4"]["p_one_sided"], v12["main_block_x4"]["p_two_sided"], v12["main_block_x4"]["ci95"]),
            "FC features (combined_fc)": (fc["theta_hat"], fc["main_block_x4"]["p_one_sided"], fc["main_block_x4"]["p_two_sided"], fc["main_block_x4"]["ci95"]),
            "v1.0 k rule": (f["theta_hat"], c["block_v1.0"]["p_one_sided"], c["block_v1.0"]["p_two_sided"], c["block_v1.0"]["ci95"]),
            "stratified null": (f["theta_hat"], c["stratified"]["p_one_sided"], c["stratified"]["p_two_sided"], c["stratified"]["ci95"]),
            f"decoy null (n = {c['decoy_neg_control'].get('n_neg_control')})": (f["theta_hat"], c["decoy_neg_control"].get("p_one_sided"), c["decoy_neg_control"].get("p_two_sided"), None),
            "co-accessibility-cluster blocks": (f["theta_hat"], c["block_coaccess_clusters"]["p_one_sided"], c["block_coaccess_clusters"]["p_two_sided"], c["block_coaccess_clusters"]["ci95"]),
            "group-mean M (v1.3)": (gm["v1.3"]["theta_hat"], gm["v1.3"]["main_block_x4"]["p_one_sided"], gm["v1.3"]["main_block_x4"]["p_two_sided"], gm["v1.3"]["main_block_x4"]["ci95"]),
            "group-mean M (v1.2)": (gm["v1.2"]["theta_hat"], gm["v1.2"]["main_block_x4"]["p_one_sided"], gm["v1.2"]["main_block_x4"]["p_two_sided"], gm["v1.2"]["main_block_x4"]["ci95"]),
        }
        out[t] = {k: {"theta_hat": v[0], "p_one_sided": v[1], "p_two_sided": v[2], "ci95": v[3]} for k, v in rows.items()}
        out[t]["_groupmean_reliability"] = gm["reliability"]
        out[t]["target_promoter_check"] = TP[t]
    jw(R4 / "sensitivity_h1.json", out)
    L = ["# H1 sensitivity analyses (labelled sensitivity; not part of the confirmatory test)", "",
         "Combined feature unless stated; one-sided p is for θ > 0.", ""]
    for t in MATCHED:
        L += [f"## {t}", "", "| analysis | θ̂ | one-sided p | two-sided p | 95% CI |", "|---|---|---|---|---|"]
        for k, v in out[t].items():
            if k.startswith("_") or k == "target_promoter_check":
                continue
            L.append(f"| {k} | {v['theta_hat']:.4g} | {g(v['p_one_sided'], 4)} | {g(v['p_two_sided'], 4)} | "
                     f"{ci(v['ci95'], 3) if v['ci95'] else '—'} |")
        tp = TP[t]
        L += ["", f"Target-promoter accessibility, knockdown level (exploratory, not pre-specified in detail): Spearman "
              f"between each knockdown's response divergence (mean over genes of (β̃²−s̃²)_A − (β̃²−s̃²)_B) and its "
              f"target's promoter CPM difference (A − B), both residualised on the target's expression difference "
              f"and on-target log2FC difference: ρ = {tp['spearman_resid']:.3f} over {tp['n_knockdowns']:,} knockdowns, "
              f"one-sided p {tp['p_one_sided_perm']:.3g}, two-sided p {tp['p_two_sided_perm']:.3g} (9,999 knockdown "
              f"permutations; knockdowns treated as exchangeable).", ""]
    L += ["## Target-promoter check, all pairs (exploratory)", "", "| pair | n kd | ρ (resid.) | p₁ | p₂ |", "|---|---|---|---|---|"]
    for t, tp in TP.items():
        L.append(f"| {t} | {tp['n_knockdowns']} | {tp['spearman_resid']:.3f} | {tp['p_one_sided_perm']:.3g} | {tp['p_two_sided_perm']:.3g} |")
    L += ["", "Reading: no sensitivity analysis turns either matched pair positive. K562 − Jurkat stays ≈ 0 under every "
          "variant; RPE1 − HepG2 stays negative (two-sided p ≤ 0.01 under every variant that permutes x; the decoy "
          "null, with few decoys, has a coarse p)."]
    (R4 / "sensitivity_h1.md").write_text("\n".join(L) + "\n")
    return out


# =========================================================================== figures
def fig_forest(h1):
    rows = []
    for t in MATCHED + MIXED:
        m = T[t]["features"]["combined_cpm"]
        rows.append((t + (" (matched)" if t in MATCHED else ""), m["theta_hat"], m["main_block_x4"]["ci95"], "C0" if t in MATCHED else "C1"))
    rows.append(("pooled matched (common genes)", POOL["matched"]["T_pooled"], POOL["matched"]["ci95"], "C0"))
    rows.append(("pooled mixed (common genes)", POOL["mixed"]["T_pooled"], POOL["mixed"]["ci95"], "C1"))
    fig, ax = plt.subplots(figsize=(7, 4.2))
    for i, (lab, th, c, col) in enumerate(rows[::-1]):
        ax.plot(c, [i, i], color=col, lw=2)
        ax.plot(th, i, "s" if "pooled" in lab else "o", color=col, ms=7)
    ax.axvline(0, color="k", lw=0.8, ls="--")
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows[::-1]], fontsize=8)
    ax.set_xlabel("θ̂ (M per unit combined-CPM difference), 95% CI by test inversion")
    ax.set_title("H1: main test (block permutation, x4, arm terms), combined CPM", fontsize=9)
    fig.tight_layout(); fig.savefig(FIG / "aim4_h1_forest.png", dpi=150); plt.close(fig)


def fig_h2(h2):
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6))
    for ax, key, lab in ((axes[0], "per_sd", "θ̂ × SD(x) [95% CI]"), (axes[1], "pr2", "partial R² / reliability [boot 95%]")):
        for i, t in enumerate(TAGS):
            for j, (f, col) in enumerate((("promoter_cpm", "C2"), ("distal_cpm", "C4"))):
                o = h2[t][f]
                y = i + (j - 0.5) * 0.3
                if key == "per_sd":
                    ax.plot(o["ci95_per_sd"], [y, y], color=col); ax.plot(o["theta_per_sd"], y, "o", color=col,
                                                                        label=f if i == 0 else None)
                else:
                    ax.plot(np.array(o["pr2_rel_boot_ci95"]) * 100, [y, y], color=col)
                    ax.plot(o["pr2_rel"] * 100, y, "o", color=col, label=f if i == 0 else None)
        ax.set_yticks(range(len(TAGS))); ax.set_yticklabels([t + (" *" if t in MATCHED else "") for t in TAGS], fontsize=7)
        ax.set_xlabel(lab + (" (%)" if key == "pr2" else ""), fontsize=8)
        if key == "per_sd":
            ax.axvline(0, color="k", lw=0.8, ls="--")
        ax.legend(fontsize=7)
    fig.suptitle("H2: promoter vs distal CPM (* matched pairs)", fontsize=9)
    fig.tight_layout(); fig.savefig(FIG / "aim4_h2_promoter_distal.png", dpi=150); plt.close(fig)


def fig_h3(h3):
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6))
    x = np.arange(len(TAGS))
    th_m = [h3[t]["theta_main"] for t in TAGS]; th_n = [h3[t]["naive_step7"]["coef"] for t in TAGS]
    axes[0].bar(x - 0.2, th_m, 0.4, label="main (f(C) + arms)")
    axes[0].bar(x + 0.2, th_n, 0.4, label="naive step 7 (1 + expr_diff)")
    axes[0].axhline(0, color="k", lw=0.8); axes[0].set_ylabel("θ̂ combined CPM"); axes[0].legend(fontsize=7)
    comps = ["naive_ols", "unrestricted", "expr_bin", "stratified", "camera"]
    for k in comps:
        axes[1].plot(x, [-np.log10(h3[t]["comparators_p_two_sided"][k]) for t in TAGS], "o-", ms=4, label=k, alpha=0.8)
    axes[1].plot(x, [-np.log10(h3[t]["p_main_two"]) for t in TAGS], "ks-", ms=6, label="main (block x4)")
    axes[1].plot(x, [-np.log10(h3[t]["naive_step7"]["p_two_sided"]) for t in TAGS], "r^--", ms=5, label="naive step 7")
    axes[1].axhline(-np.log10(0.05), color="grey", ls=":")
    axes[1].set_ylabel("−log10 two-sided p"); axes[1].legend(fontsize=6)
    for ax in axes:
        ax.set_xticks(x); ax.set_xticklabels(TAGS, rotation=30, ha="right", fontsize=7)
    fig.suptitle("H3: naive analyses vs main test on real data", fontsize=9)
    fig.tight_layout(); fig.savefig(FIG / "aim4_h3_naive_vs_main.png", dpi=150); plt.close(fig)


def fig_h4():
    fig, axes = plt.subplots(2, 3, figsize=(11, 6.5), sharex=True)
    mk = {"own": ("*", "k", 120), "other_pair": ("D", "C3", 35), "decoy": ("o", "C0", 20)}
    for ax, t in zip(axes.ravel(), TAGS):
        d = H4[H4.pair == t]
        for grp, (m, c, s) in mk.items():
            q = d[d.group == grp]
            ax.scatter(q.similarity, q.theta_hat, marker=m, c=c, s=s, label=grp)
        f = H4F["fits"][t]["all"]
        xx = np.linspace(-1, 1, 10)
        ax.plot(xx, f["intercept"] + f["slope"] * xx, "k-", lw=1)
        ax.axhline(0, color="grey", lw=0.6); ax.axvline(0, color="grey", lw=0.6)
        ax.axvspan(-0.05, 0.05, color="yellow", alpha=0.25)
        ax.set_title(f"{t} ({ROLE[t]})\nslope {f['slope']:.3g} (p {f['slope_p']:.2g}), int. {f['intercept']:.2g}", fontsize=8)
    axes[0, 0].legend(fontsize=7)
    for ax in axes[1]:
        ax.set_xlabel("similarity to own x (Pearson, C-residualised)")
    for ax in axes[:, 0]:
        ax.set_ylabel("θ̂")
    fig.suptitle("H4 similarity gradient (combined CPM; yellow = decoy negative-control band |sim| < 0.05)", fontsize=9)
    fig.tight_layout(); fig.savefig(FIG / "aim4_similarity_gradient.png", dpi=150); plt.close(fig)


def fig_p53():
    fig, ax = plt.subplots(figsize=(6.5, 3.4))
    for i, t in enumerate(TAGS):
        r = P53[t]
        col = "C1" if r["p53_differs"] else "C0"
        ax.bar(i, r["z_oriented"], color=col)
    ax.axhline(0, color="k", lw=0.8); ax.axhline(1.645, color="grey", ls=":"); ax.axhline(-1.96, color="grey", ls=":")
    ax.set_xticks(range(len(TAGS))); ax.set_xticklabels(TAGS, rotation=30, ha="right", fontsize=7)
    ax.set_ylabel("z of mean M at p53 targets\n(oriented: + = p53-WT line responds more)", fontsize=8)
    ax.set_title("p53 control (orange: pair differs in p53; blue: matched)", fontsize=9)
    fig.tight_layout(); fig.savefig(FIG / "aim4_p53_control.png", dpi=150); plt.close(fig)


def summary(a1, h1, h2, h3, negsum):
    pm, pmat, gap = POOL["mixed"], POOL["matched"], POOL["gap"]
    fits = H4F["fits"]
    oc = [
        ("Positive effect in both matched pairs, rising with similarity", "no",
         "K562 − Jurkat θ̂ ≈ 0; RPE1 − HepG2 θ̂ < 0; H4 slopes ≤ 0 where the own effect is nonzero."),
        ("Effect in the mixed pairs only", "no (as stated)",
         f"Mixed pooled θ̂ {pm['T_pooled']:.4f} (CI {ci(pm['ci95'])}) is nonzero but negative, and a matched pair "
         f"(RPE1 − HepG2) shows the same negative sign; gap mixed − matched {gap['gap_mixed_minus_matched']:.4f} "
         f"(approx. CI {ci(gap['approx_ci95'])}, p {gap['p_two_sided_joint_null']:.2g})."),
        ("The generic response accounts for most of the raw difference", "no",
         "Median generic share 2–13% per line in every pair."),
        ("Naive tests find an effect; the main test does not", "no (one borderline mixed pair)",
         "Naive and main agree in sign and significance in 5 of 6 pairs; K562 − HepG2 (mixed): permutation-free "
         "naive tests p₂ 0.008–0.022 vs main 0.055. The naive step-7 estimate is not larger than the main one."),
        ("No effect under any test, with adequate power", "yes, for K562 − Jurkat",
         "θ̂ ≈ 0 under every test and sensitivity; planted detection limit 0.5% partial R²; pR²/rel upper bound 0.023%."),
        ("Cross-line numbers sit near the same-line baseline", "partly, for K562 − Jurkat",
         f"K562 − Jurkat reproducible variance of M is 3.6× the same-line baseline (other pairs 43–154×); "
         f"reliability 0.86."),
        ("The two matched pairs disagree", "YES — lead",
         "K562 − Jurkat null; RPE1 − HepG2 significantly negative (opposite to H1). What differs: RPE1 − HepG2 still "
         "diverges at p53 targets (z 3.3) and every RPE1-containing pair shows the same promoter-driven negative "
         "association, while the three pairs without RPE1 do not."),
        ("Signal in distal accessibility only (secondary)", "no (weak hint only)",
         "Distal θ̂ weakly positive in the pairs without RPE1 (K562 − Jurkat one-sided p 0.021) but H2's "
         "non-overlap criterion fails in both matched pairs."),
    ]
    jw(R4 / "outcome_classification.json", {"rows": [{"outcome": a, "supported": b, "evidence": c} for a, b, c in oc],
                                            "lead_with": "pair-by-pair comparison of the two matched pairs (they "
                                            "disagree) with K562 − Jurkat's detection limit (null with adequate power) "
                                            "and RPE1 − HepG2's opposite-sign association traced to RPE1-containing pairs "
                                            "and residual p53 divergence"})
    L = ["# Aim 4 — summary (frozen pipeline `pipeline-v1.1`, plan v1.3; H1 test run once)", "",
         "Date: 2026-10-12. All analyses: `scripts/run_aim4_aim1.sh`, `scripts/run_aim4.sh` (`scripts/run_aim4.py`), "
         "`scripts/report_aim4.py`; new code `src/cdtest/aim4.py` (pooled contrast, H4, p53 control, target-promoter "
         "check); every test, rule and parameter as frozen (B = 9,999; x4 block rule; arm terms in f(C); CPM "
         "features; 500 bootstrap resamples; planted 500 reps/level, B = 999).", "",
         "## Pairs (Aim 1; M = A − B)", "",
         "| pair | role | genes | knockdowns kept | reliability | generic share A / B | SD M |", "|---|---|---|---|---|---|---|"]
    for r in a1.to_dict(orient="records"):
        L.append(f"| {r['pair']} | {r['role']} | {r['n_genes']:,} | {r['n_kd_kept']:,} / {r['n_kd_shared']:,} | "
                 f"{r['reliability']:.3f} | {r['generic_share_A']:.1%} / {r['generic_share_B']:.1%} | {r['M_sd']:.4f} |")
    L += ["", "## Headline table (combined CPM feature, main test)", "",
          "| pair | role | θ̂ | 95% CI | p (θ > 0) | p two-sided | Holm p | pR²/rel [boot] | detection limit | p53 control z (oriented) |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for t in MATCHED + MIXED:
        f = T[t]["features"]["combined_cpm"]; m = f["main_block_x4"]
        hp = h1["pairs"][t]["p_holm"] if t in MATCHED else None
        dl = (f"{pc(h1['pairs'][t]['detection_limit_all'])} / {pc(h1['pairs'][t]['detection_limit_excl_sig_theta0'])}"
              if t in MATCHED else "—")
        L.append(f"| {t} | {ROLE[t]} | {f['theta_hat']:.4f} | {ci(m['ci95'])} | {m['p_one_sided']:.4g} | "
                 f"{m['p_two_sided']:.4g} | {g(hp, 3) if hp is not None else '—'} | {f['partial_r2_share_reproducible']:.3%} "
                 f"[{f['partial_r2_share_boot_ci95'][0]:.3%}, {f['partial_r2_share_boot_ci95'][1]:.3%}] | {dl} | "
                 f"{P53[t]['z_oriented']:.1f} |")
    L += [f"| pooled mixed (4 pairs, {pm['n_common_genes']:,} common genes, one shared block permutation) | mixed | "
          f"{pm['T_pooled']:.4f} | {ci(pm['ci95'])} | {pm['p_one_sided']:.4g} | {pm['p_two_sided']:.4g} | — | — | — | — |",
          f"| pooled matched (same genes) | matched | {pmat['T_pooled']:.4f} | {ci(pmat['ci95'])} | {pmat['p_one_sided']:.4g} | "
          f"{pmat['p_two_sided']:.4g} | — | — | — | — |", "",
          "Detection limit = smallest planted partial R² recovered with 80% power at nominal 5% (all decoys / "
          "excluding decoys whose own θ₀ is significant). p53 z oriented to the p53-WT line (matched: A − B).", "",
          "## Hypotheses", "",
          f"- **H1 (primary): not supported** — no matched pair shows a positive effect at 5% (Holm) under the main test. K562 − Jurkat: θ̂ ≈ 0 with adequate power "
          f"(detection limit 0.5%). RPE1 − HepG2: θ̂ = {h1['pairs']['RPE1_vs_HepG2']['theta_hat']:.4f}, opposite sign, "
          "two-sided p 1e-4. Holm p = 1 for both.",
          "- **H2 (secondary): not supported.** Distal is not larger with non-overlapping intervals in either matched "
          "pair; in RPE1 − HepG2 the promoter carries more (0.59% vs 0.23%).",
          "- **H3: simulation half supported (Aim 3, planted naive FPR 0.27–0.43), real-data half not supported** — "
          "the naive step-7 estimate is smaller in magnitude than the main estimate in 4/6 pairs and flips sign in "
          "K562 − HepG2; naive and main tests agree except borderline K562 − HepG2.",
          f"- **H4 (secondary): not supported as stated.** Pooled slope {fits['pooled_all_pairs']['all']['slope']:.4f} "
          f"(p {fits['pooled_all_pairs']['all']['slope_p']:.2g}), intercept {fits['pooled_all_pairs']['all']['intercept']:.4f} "
          f"(p {fits['pooled_all_pairs']['all']['intercept_p']:.2g}); per-pair slopes are negative where the own θ is "
          "negative and ≈ 0 where it is ≈ 0 — the gradient tracks the own effect (consistent with the effect model, "
          "but with negative θ). Decoy negative control: "
          f"{sum(v['n_sig'] for v in negsum.values())} of {sum(v['n'] for v in negsum.values())} low-similarity decoys "
          "reject at 5% — low similarity does not imply no association.",
          f"- **Mixed vs matched:** mixed pooled {pm['T_pooled']:.4f} vs matched pooled {pmat['T_pooled']:.4f}; gap "
          f"{gap['gap_mixed_minus_matched']:.4f} (approx. 95% {ci(gap['approx_ci95'])}, p {gap['p_two_sided_joint_null']:.2g}).",
          "- **p53 control:** all four mixed pairs diverge (z 4.1–7.0, oriented as expected); K562 − Jurkat none "
          "(z −0.7); **RPE1 − HepG2 diverges (z 3.3)** although both are p53 WT.",
          "- **Sensitivity (H1):** v1.2, FC features, v1.0 k rule, stratified, decoy null, co-accessibility blocks and "
          "group-mean M: no variant gives a positive effect in either matched pair (`sensitivity_h1.md`). The "
          "RPE1 − HepG2 negative association holds under every CPM variant but not with FC features (θ̂ −0.004, "
          "p₂ 0.27), and is weaker with co-accessibility blocks (p₂ 0.04) and the 6-decoy null (p₂ 0.14). Target-promoter knockdown-level check: no positive association (ρ −0.05 to 0.00).", "",
          "## Outcome classification (proposal 'Interpreting outcomes')", "",
          "| outcome | supported? | evidence |", "|---|---|---|"]
    for a, b, c in oc:
        L.append(f"| {a} | {b} | {c} |")
    L += ["", "**The report should lead with:** the pair-by-pair comparison — *the two matched pairs disagree*. "
          "K562 − Jurkat gives a null with adequate power (state its detection limit, 0.5% partial R²); RPE1 − HepG2 "
          "gives a significant association of the opposite sign to H1. Then 'what differs': the negative, "
          "promoter-carried association appears in every pair containing RPE1 (RPE1 − K562, Jurkat − RPE1, "
          "RPE1 − HepG2) and in none without it, and RPE1 − HepG2 is not p53-matched in practice. Candidate "
          "explanations to state, not test: RPE1-specific feature/CN processing (local reprocessing of 4DN ATAC, "
          "karyotype-only CN track), RPE1 biology (strong p53 response), cross-study batch. H1 as pre-specified is "
          "not supported.", "",
          "## Files", "",
          "- `aim1_pairs.{json,md}`, `h1_primary.{json,md}`, `planted_<pair>.{json,md}`, `h2_location.{json,md}`, "
          "`h3_comparators.{json,md}`, `mixed_pooled.{json,md}`, `h4_gradient.{json,md}`, `p53_control.{json,md}`, "
          "`sensitivity_h1.{json,md}`, `outcome_classification.json`; raw stage outputs in `tests/`, `sens_v12/`, "
          "`sens_groupmean/`, `planted/`, `extras/`, logs in `logs/`.",
          "- Figures: `figures/aim4_h1_forest.png`, `aim4_h2_promoter_distal.png`, `aim4_h3_naive_vs_main.png`, "
          "`aim4_similarity_gradient.png`, `aim4_p53_control.png`."]
    (R4 / "SUMMARY.md").write_text("\n".join(L) + "\n")


if __name__ == "__main__":
    a1 = aim1_report()
    h1, ps = h1_report()
    h2 = h2_report()
    h3 = h3_report()
    mixed_report()
    negsum = h4_report()
    p53_report()
    sens = sensitivity_report()
    fig_forest(h1); fig_h2(h2); fig_h3(h3); fig_h4(); fig_p53()
    summary(a1, h1, h2, h3, negsum)
    print("reports written")
