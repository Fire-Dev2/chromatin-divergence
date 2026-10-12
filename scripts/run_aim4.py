"""Aim 4 driver (frozen pipeline pipeline-v1.1, plan v1.3). One stage per call:

  python scripts/run_aim4.py tests A B        main test + comparators + sensitivities, v1.3 (aim2.run)
  python scripts/run_aim4.py tests_v12 A B    same under plan v1.2 (arm_effects=False), sensitivity
  python scripts/run_aim4.py planted A B      planted-signal benchmark (plan §6) on a matched pair
  python scripts/run_aim4.py groupmean A B    group-mean M variant, combined CPM, main test (sensitivity)
  python scripts/run_aim4.py extras           pooled mixed contrast, H4 gradient, decoy negative
                                              control, p53 control, target-promoter check

No analysis parameter is set here that differs from the frozen defaults (B = 9,999, 500 bootstrap
resamples, x4 rule, arm effects on; planted 500 reps/level, B = 999, theta0 B = 9,999).
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from cdtest import aim2, aim4, crt, planted

ROOT = Path(__file__).resolve().parents[1]
R4 = ROOT / "results/aim4"
STATUS = {"matched": "Aim 4 matched pair (primary H1 test, run once on pipeline-v1.1)",
          "mixed": "Aim 4 mixed pair (p53/tissue/growth confounded; enters the pooled contrast; descriptive)"}
ROLE = {aim4.tag(a, b): r for a, b, r in aim4.PAIRS}


def _set_status(path: Path, t: str):
    d = json.loads(path.read_text())
    d["status"] = STATUS[ROLE[t]]
    path.write_text(json.dumps(d, indent=1, default=aim2._json_default))


def stage_tests(a, b, v12=False):
    out = R4 / ("sens_v12" if v12 else "tests")
    aim2.run(a, b, B=9999, n_boot=500, out_dir=out, arm_effects=not v12)
    t = aim4.tag(a, b)
    _set_status(out / f"{t}_test{'' if v12 else '_v13'}.json", t)


def stage_planted(a, b):
    out = R4 / "planted"
    df = planted.run(a, b, reps=500, B=999, n_jobs=10, arm_effects=True, out_dir=out)
    t = aim4.tag(a, b)
    s = planted.summarize(df)
    s.to_csv(out / f"planted_{t}_v13.csv", index=False)
    info = json.loads((out / f"planted_{t}_v13_info.json").read_text())
    parts = []
    for col in ("p_block_x4_two_sided", "p_naive_two_sided"):
        ex = planted.theta0_significant(info, col)
        keep = df[~df.decoy.isin(ex)]
        se_ = planted.summarize(keep) if len(keep) else s.iloc[0:0]
        parts.append(se_.assign(excluded=f"theta0_{col}<=0.05:" + ";".join(ex)))
    pd.concat(parts).to_csv(out / f"planted_{t}_v13_excl.csv", index=False)
    print("detection limit x4:", planted.detection_limit(s, "block_x4"), flush=True)


def stage_groupmean(a, b):
    t = f"{a}_vs_{b}_groupmean"
    out = R4 / "sens_groupmean"
    out.mkdir(parents=True, exist_ok=True)
    res = {}
    for arm in (True, False):
        P = aim2.load_pair(a, b, tag=t, arm_effects=arm)
        r = aim2.test_feature(P, aim4.PRIMARY, B=9999)
        r["comparators"] = {k: v for k, v in r["comparators"].items() if k in ("block_v1.0", "naive_ols")}
        res["v1.3" if arm else "v1.2"] = r
    res.update(pair=t, A=a, B=b, reliability=P["summary"]["reliability"], n_genes=len(P["genes"]),
               n_knockdowns=P["summary"]["n_kd_kept"])
    (out / f"{t}_test.json").write_text(json.dumps(res, indent=1, default=aim2._json_default))


def stage_extras(only=None):
    out = R4 / "extras"
    out.mkdir(parents=True, exist_ok=True)
    Ps = {aim4.tag(a, b): aim2.load_pair(a, b, arm_effects=True) for a, b, _ in aim4.PAIRS}
    # ---- pooled mixed contrast (+ matched pooled, for the gap)
    mixed = [Ps[aim4.tag(a, b)] for a, b, r in aim4.PAIRS if r == "mixed"]
    matched = [Ps[aim4.tag(a, b)] for a, b, r in aim4.PAIRS if r == "matched"]
    # primary: one block permutation applied identically to all pairs, on their common genes
    pm = aim4.pooled_shared_common(mixed, seed_key="pooled_mixed_common")
    pmat = aim4.pooled_shared_common(matched, seed_key="pooled_matched_common")
    six = aim4.pooled_shared_common(matched + mixed, seed_key="pooled_all6_common")
    Ts6 = six.pop("_Tstar_pairs"); six.pop("_cstar_pairs")
    Te = np.array([six["per_pair"][P["tag"]]["T_common_genes"] for P in matched + mixed])
    gap = float(Te[2:].mean() - Te[:2].mean())
    gnull = Ts6[2:].mean(0) - Ts6[:2].mean(0)
    gsd = float(gnull.std())
    gap_d = {"pairs_common": six["pairs"], "n_common_genes_all6": six["n_common_genes"], "k_shared": six["k_shared"],
             "T_mixed_mean_common6": float(Te[2:].mean()), "T_matched_mean_common6": float(Te[:2].mean()),
             "gap_mixed_minus_matched": gap, "null_sd_gap": gsd, "null_mean_gap": float(gnull.mean()),
             "approx_ci95": [gap - 1.96 * gsd, gap + 1.96 * gsd],
             "p_two_sided_joint_null": float((1 + (np.abs(gnull - gnull.mean()) >= abs(gap - gnull.mean())).sum()) / (1 + gnull.size)),
             "gap_full_gene_sets": pm["T_pooled"] - pmat["T_pooled"],
             "note": "Gap on the genes common to all six pairs; null = the same single block permutation applied "
                     "to all six pairs (B = 9,999); interval = gap +/- 1.96 x null SD of the gap (approximate, "
                     "descriptive). gap_full_gene_sets uses the 4-pair and 2-pair common gene sets."}
    # secondary (documented): union-order coupling, each pair on its full gene set
    pu = aim4.pooled_mixed_contrast(mixed)
    pu.pop("_Tstar_pairs"); pu.pop("_cstar_pairs")
    np.savez_compressed(out / "pooled_nulls.npz", mixed=pm.pop("_Tstar_pairs"), mixed_c=pm.pop("_cstar_pairs"),
                        matched=pmat.pop("_Tstar_pairs"), matched_c=pmat.pop("_cstar_pairs"), all6=Ts6)
    (out / "pooled_contrast.json").write_text(json.dumps({"mixed": pm, "matched": pmat, "gap": gap_d,
                                                          "all6": six, "mixed_union_rank_coupling": pu},
                                                         indent=1, default=aim2._json_default))
    print("pooled mixed", pm["T_pooled"], pm["p_one_sided"], "matched", pmat["T_pooled"], flush=True)
    if only == "pooled":
        return
    # ---- H4 gradient
    grads = pd.concat([aim4.similarity_gradient(P) for P in Ps.values()], ignore_index=True)
    grads.to_csv(out / "h4_similarity_gradient.csv", index=False)
    fits = {}
    for t, g in grads.groupby("pair", sort=False):
        fits[t] = {"all": aim4.ols_line(g.similarity, g.theta_hat),
                   "excluding_own": aim4.ols_line(g[g.group != "own"].similarity, g[g.group != "own"].theta_hat),
                   "standardised_theta": aim4.ols_line(g.similarity, g.theta_std)}
    fits["pooled_all_pairs"] = {"all": aim4.ols_line(grads.similarity, grads.theta_hat),
                                "excluding_own": aim4.ols_line(grads[grads.group != "own"].similarity,
                                                               grads[grads.group != "own"].theta_hat),
                                "standardised_theta": aim4.ols_line(grads.similarity, grads.theta_std)}
    # ---- decoy negative control: |sim| < 0.05 decoys vs their own block-permutation null
    neg = {}
    for t, P in Ps.items():
        D, labels = aim2.decoy_matrix(P["genes"])
        sim = crt.decoy_similarity(P["X"][aim4.PRIMARY], D, P["Q"])
        idx = np.flatnonzero(np.abs(sim) < aim2.NEG_CONTROL_SIM)
        r = planted.decoy_theta0_tests(P, D, labels, idx, B=9999)
        for l in r:
            r[l]["similarity"] = float(sim[labels.index(l)])
        neg[t] = r
    (out / "h4_gradient_fits.json").write_text(json.dumps({"fits": fits, "decoy_negative_control": neg},
                                                          indent=1, default=aim2._json_default))
    # ---- p53 control
    p53 = {}
    for a, b, _ in aim4.PAIRS:
        r = aim4.p53_control(a, b)
        r.pop("_null")
        p53[aim4.tag(a, b)] = r
    (out / "p53_control.json").write_text(json.dumps(p53, indent=1, default=aim2._json_default))
    # ---- target-promoter knockdown-level check
    tp = {aim4.tag(a, b): aim4.target_promoter_check(a, b) for a, b, _ in aim4.PAIRS}
    (out / "target_promoter_check.json").write_text(json.dumps(tp, indent=1, default=aim2._json_default))


if __name__ == "__main__":
    t0 = time.time()
    st = sys.argv[1]
    if st == "tests":
        stage_tests(*sys.argv[2:4])
    elif st == "tests_v12":
        stage_tests(*sys.argv[2:4], v12=True)
    elif st == "planted":
        stage_planted(*sys.argv[2:4])
    elif st == "groupmean":
        stage_groupmean(*sys.argv[2:4])
    elif st == "extras":
        stage_extras(*sys.argv[2:3])
    else:
        raise SystemExit(__doc__)
    print(f"{st} done {time.time() - t0:.0f}s", flush=True)
