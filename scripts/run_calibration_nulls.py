"""Calibration of the main test on the real outcomes without decoy features (added in v0.2.0,
in response to review; not part of the frozen confirmatory analysis).

1. circular   Decoy-free real-data null. For every pair (six pairs and the same-line pair), the
              pair's own accessibility difference x is circularly shifted within each chromosome
              arm by an independent random offset (genes in genome order). The shifted feature keeps
              x's values, its within-arm autocorrelation and its arm composition (which the arm terms
              in f(C) absorb) but is no longer aligned with M, so its true effect is 0 by
              construction. The frozen main test (x4 rule, arm terms, B = 999) is run on each shifted
              feature; naive OLS on the same features for contrast. Rejection rate at two-sided 5%
              should be ~5% if the test is calibrated on that M. For the same-line M each replicate
              shifts one of the six pair features (cycled).
2. sameline   Planted signal on the same-line M (K562 essential - K562 genome-wide). Null features:
              circular shifts of the six pair features (decoy-free) and, separately, the 21 decoy-line
              differences with random sign. Signal planted at partial R2 in {0, 0.25, 0.5, 1, 2}%
              exactly as in cdtest.planted; main test and naive OLS. Level 0 is the false-positive
              rate on this M; the decoys' own association with the same-line M (B = 9,999) is also
              reported, to compare with the 3 of 6 pair features that reject (Figure 4A).
3. report     results/calibration/{circular_shift_null,sameline_planted}.{csv,md} and summary.json

    python scripts/run_calibration_nulls.py circular [reps=2000] [B=999] [pair tags ...]
    python scripts/run_calibration_nulls.py sameline [reps=500] [B=999]
    python scripts/run_calibration_nulls.py report
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from cdtest import aim2, aim4, crt, features, planted  # noqa: E402
from cdtest import covariates as cov  # noqa: E402

OUT = ROOT / "results/calibration"
SAME_TAG = "K562_vs_K562gw"
N_JOBS = int(os.environ.get("CDTEST_JOBS", max(1, (os.cpu_count() or 2) - 1)))


def same_line_pair() -> dict:
    gs = pd.read_csv(ROOT / f"results/aim1/{SAME_TAG}/gene_outcome.tsv", sep="\t")
    genes = gs.gene_id.tolist()
    C = cov.build_covariates("K562", "K562gw", genes, SAME_TAG)
    Z = cov.design(C)
    P = dict(a="K562", b="K562gw", tag=SAME_TAG, genes=genes, M=gs.M.to_numpy(float), C=C, Z=Z,
             gen=cov.genome_info(genes), Q=crt.orthonormal_basis(Z),
             summary=json.loads((ROOT / f"results/aim1/{SAME_TAG}/summary.json").read_text()), X={})
    for a, b, _ in aim4.PAIRS:
        x = features.pair_feature(a, b, "combined", signal="cpm").reindex(genes).to_numpy(float)
        P["X"][aim4.tag(a, b)] = np.where(np.isfinite(x), x, np.nanmean(x))
    return P


def arm_segments(P: dict):
    g = P["gen"]
    order = crt.genome_order(g.chrom.to_numpy(), g.pos.to_numpy())
    arm_o = g.arm.to_numpy()[order]
    cuts = np.flatnonzero(arm_o[1:] != arm_o[:-1]) + 1
    return order, np.split(np.arange(len(order)), cuts)


def circular_shift(x: np.ndarray, order: np.ndarray, segs, rng) -> np.ndarray:
    xo = x[order].copy()
    for s in segs:
        if len(s) > 1:
            xo[s] = np.roll(xo[s], int(rng.integers(1, len(s))))
    out = np.empty_like(x)
    out[order] = xo
    return out


def _test(M, v, Q, chrom, pos, arm, order, B, seed, theta_plant=0.0):
    """Main test (x4 k from v's residuals) and naive OLS of v on M + theta_plant * v."""
    rng = np.random.default_rng(seed)
    Mp = M + theta_plant * v
    prep = crt.prepare(Mp, v, Q=Q)
    k = crt.choose_k(prep.r_x[order], arm[order], rule="x4")
    r = crt.block_permutation_test(Mp, v, chrom=chrom, pos=pos, arm=arm, k=k, B=B, rng=rng, prep=prep, order=order)
    n = crt.naive_ols_test(Mp, v, prep=prep)
    return {"theta_hat": prep.theta_hat, "k_x4": int(k), "p_block_two_sided": float(r.p_two_sided),
            "p_block_one_sided": float(r.p), "p_naive_two_sided": float(n["p_two_sided"])}


def _theta_for(M, v, Q, R):
    """theta giving planted partial R2 = R (as cdtest.planted._one)."""
    if R <= 0:
        return 0.0
    rv = v - Q @ (Q.T @ v)
    ssv = float(rv @ rv)
    rM = M - Q @ (Q.T @ M)
    th0 = float(rv @ rM) / ssv
    rss0 = float(rM @ rM) - th0 ** 2 * ssv
    return float(np.sqrt(R * rss0 / ((1 - R) * ssv)))


def stage_circular(reps: int = 2000, B: int = 999, *only: str):
    """``only``: optional pair tags to (re)run; rows of the other pairs already in the CSV are kept."""
    OUT.mkdir(parents=True, exist_ok=True)
    pairs = [(aim4.tag(a, b), lambda a=a, b=b: aim2.load_pair(a, b)) for a, b, _ in aim4.PAIRS]
    pairs.append((SAME_TAG, same_line_pair))
    if only:
        pairs = [p for p in pairs if p[0] in only]
    f = OUT / "circular_shift_null.csv"
    keep = pd.read_csv(f) if f.exists() else pd.DataFrame()
    if len(keep):
        keep = keep[~keep.pair.isin([t for t, _ in pairs])]
    rows = [keep] if len(keep) else []
    for tag, loader in pairs:
        t0 = time.time()
        P = loader()
        g = P["gen"]
        chrom, pos, arm = g.chrom.to_numpy(), g.pos.to_numpy(), g.arm.to_numpy()
        order, segs = arm_segments(P)
        labs = list(P["X"]) if tag == SAME_TAG else ["combined_cpm"]
        rng = aim2.job_seed("calibration_circular", tag)
        tasks, meta = [], []
        for rep in range(reps):
            lab = labs[rep % len(labs)]
            v = circular_shift(P["X"][lab], order, segs, rng)
            tasks.append((P["M"], v, P["Q"], chrom, pos, arm, order, B, int(rng.integers(2 ** 63 - 1))))
            meta.append({"pair": tag, "rep": rep, "feature": lab})
        res = Parallel(n_jobs=N_JOBS, batch_size=16)(delayed(_test)(*t) for t in tasks)
        df = pd.concat([pd.DataFrame(meta), pd.DataFrame(res)], axis=1)
        rows.append(df)
        print(f"{tag}: block {np.mean(df.p_block_two_sided <= 0.05):.3f}, naive {np.mean(df.p_naive_two_sided <= 0.05):.3f} "
              f"({time.time() - t0:.0f}s)", flush=True)
        out = pd.concat(rows)
        out["B"] = out["B"].fillna(B) if "B" in out else B
        order = [aim4.tag(a, b) for a, b, _ in aim4.PAIRS] + [SAME_TAG]
        out = out.sort_values(["pair", "rep"], key=lambda c: c.map(order.index) if c.name == "pair" else c)
        out.to_csv(f, index=False)


def stage_sameline(reps: int = 500, B: int = 999):
    OUT.mkdir(parents=True, exist_ok=True)
    P = same_line_pair()
    g = P["gen"]
    chrom, pos, arm = g.chrom.to_numpy(), g.pos.to_numpy(), g.arm.to_numpy()
    order, segs = arm_segments(P)
    D, labels = aim2.decoy_matrix(P["genes"], "combined", "cpm")
    labs = list(P["X"])
    rng = aim2.job_seed("calibration_sameline", SAME_TAG)
    tasks, meta = [], []
    for null in ("circular", "decoy"):
        for R in planted.LEVELS:
            for rep in range(reps):
                if null == "circular":
                    lab = labs[rep % len(labs)]
                    v = circular_shift(P["X"][lab], order, segs, rng)
                else:
                    j = int(rng.integers(len(labels)))
                    lab = labels[j]
                    v = int(rng.choice([-1, 1])) * D[:, j]
                th = _theta_for(P["M"], v, P["Q"], R)
                tasks.append((P["M"], v, P["Q"], chrom, pos, arm, order, B, int(rng.integers(2 ** 63 - 1)), th))
                meta.append({"null": null, "level": R, "rep": rep, "feature": lab, "theta": th})
    t0 = time.time()
    res = Parallel(n_jobs=N_JOBS, batch_size=16)(delayed(_test)(*t) for t in tasks)
    df = pd.concat([pd.DataFrame(meta), pd.DataFrame(res)], axis=1).assign(B=B)
    df.to_csv(OUT / "sameline_planted.csv", index=False)
    print(f"planted on same-line M: {len(df)} tests in {time.time() - t0:.0f}s", flush=True)
    # each decoy's own association with the same-line M, frozen settings (B = 9,999)
    own = {}
    for j, lab in enumerate(labels):
        r = _test(P["M"], D[:, j], P["Q"], chrom, pos, arm, order, 9999,
                  int(aim2.job_seed("theta0_sameline", lab).integers(2 ** 63 - 1)))
        own[lab] = r
    (OUT / "sameline_decoy_own.json").write_text(json.dumps(own, indent=1))


def _rate(p, a=0.05):
    p = np.asarray(p, float)
    n, k = len(p), int((p <= a).sum())
    r = k / n
    se = np.sqrt(r * (1 - r) / n) if n else np.nan
    # Wilson 95% interval
    z = 1.959964
    c = (r + z * z / (2 * n)) / (1 + z * z / n)
    h = z * np.sqrt(r * (1 - r) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return {"n": n, "k": k, "rate": r, "se": float(se), "wilson95": [float(c - h), float(c + h)]}


def stage_report():
    summ = {}
    L = ["# Calibration of the main test on the real outcomes (v0.2.0; not part of the frozen analysis)", ""]
    f = OUT / "circular_shift_null.csv"
    if f.exists():
        df = pd.read_csv(f)
        L += ["## 1. Decoy-free null: circular shifts of the pair's own feature within chromosome arms", "",
              f"Main test: block permutation, x4 rule, arm terms, B = {int(df.B.iloc[0]):,}; "
              f"{df.groupby('pair').size().iloc[0]:,} shifted features per pair. Rejection at two-sided 5% "
              "(Wilson 95% interval).", "",
              "The last column uses the shifted features as a decoy-free reference distribution for the observed θ̂ "
              "(two-sided: share of shifted |θ̂| at least the observed |θ̂|).", "",
              "| pair | main test | naive OLS | observed θ̂ | circular-shift p |", "|---|---|---|---|---|"]
        summ["circular"] = {}
        for tag, d in df.groupby("pair", sort=False):
            rb, rn = _rate(d.p_block_two_sided), _rate(d.p_naive_two_sided)
            summ["circular"][tag] = {"block_x4": rb, "naive_ols": rn}
            obs, pc = "", ""
            if tag != SAME_TAG:
                f = ROOT / f"results/aim4/tests/{tag}_test_v13.json"
                f = f if f.exists() else ROOT / f"results/aim2/{tag}_test_v13.json"
                th = json.loads(f.read_text())["features"]["combined_cpm"]["theta_hat"]
                p_c = float((1 + (np.abs(d.theta_hat) >= abs(th)).sum()) / (1 + len(d)))
                summ["circular"][tag].update(theta_obs=th, p_circular_two_sided=p_c,
                                             null_sd=float(d.theta_hat.std()), null_mean=float(d.theta_hat.mean()))
                obs, pc = f"{th:.4f}", f"{p_c:.4g}"
            L.append(f"| {tag} | {100 * rb['rate']:.1f}% [{100 * rb['wilson95'][0]:.1f}, {100 * rb['wilson95'][1]:.1f}] | "
                     f"{100 * rn['rate']:.1f}% [{100 * rn['wilson95'][0]:.1f}, {100 * rn['wilson95'][1]:.1f}] | {obs} | {pc} |")
        six = df[df.pair != SAME_TAG]
        rb, rn = _rate(six.p_block_two_sided), _rate(six.p_naive_two_sided)
        summ["circular"]["six_pairs_pooled"] = {"block_x4": rb, "naive_ols": rn}
        L += [f"| six pairs pooled | {100 * rb['rate']:.1f}% [{100 * rb['wilson95'][0]:.1f}, {100 * rb['wilson95'][1]:.1f}] | "
              f"{100 * rn['rate']:.1f}% [{100 * rn['wilson95'][0]:.1f}, {100 * rn['wilson95'][1]:.1f}] | | |", ""]
    f = OUT / "sameline_planted.csv"
    if f.exists():
        df = pd.read_csv(f)
        L += ["## 2. Planted signal on the same-line M (K562 essential − K562 genome-wide)", "",
              "Rejection rate of the main test (and naive OLS) at two-sided 5% by planted partial R²; level 0 is the "
              "false-positive rate on this M.", "",
              "| null features | planted partial R² | main test | naive OLS |", "|---|---|---|---|"]
        summ["sameline_planted"] = {}
        for (null, R), d in df.groupby(["null", "level"], sort=False):
            rb, rn = _rate(d.p_block_two_sided), _rate(d.p_naive_two_sided)
            summ["sameline_planted"][f"{null}_{R}"] = {"block_x4": rb, "naive_ols": rn}
            L.append(f"| {null} | {100 * R:g}% | {100 * rb['rate']:.1f}% [{100 * rb['wilson95'][0]:.1f}, "
                     f"{100 * rb['wilson95'][1]:.1f}] | {100 * rn['rate']:.1f}% |")
        own = OUT / "sameline_decoy_own.json"
        if own.exists():
            o = json.loads(own.read_text())
            kb = sum(v["p_block_two_sided"] <= 0.05 for v in o.values())
            kn = sum(v["p_naive_two_sided"] <= 0.05 for v in o.values())
            summ["sameline_decoy_own"] = {"n": len(o), "reject_block_x4": kb, "reject_naive": kn}
            L += ["", f"Decoy-line differences tested against the same-line M (B = 9,999): {kb} of {len(o)} reject "
                  f"with the main test ({kn} of {len(o)} with naive OLS); the six pair features: 3 of 6 (Figure 4A)."]
    (OUT / "calibration.md").write_text("\n".join(L) + "\n")
    (OUT / "summary.json").write_text(json.dumps(summ, indent=1))
    print("\n".join(L))


if __name__ == "__main__":
    st = sys.argv[1] if len(sys.argv) > 1 else "report"
    args = [int(x) if x.isdigit() else x for x in sys.argv[2:]]
    {"circular": stage_circular, "sameline": stage_sameline, "report": stage_report}[st](*args)
