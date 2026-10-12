"""Aim 3 planted-signal benchmark on a real response difference M (ANALYSIS_PLAN §6).

For level R in {0, 0.25%, 0.5%, 1%, 2%} and each replicate:
  1. draw a decoy feature x_d at random from the negative-control set (decoy-line pair
     differences with |similarity| < 0.05 to the pair's own x; similarity = Pearson
     correlation of features residualized on C) and a random sign s; the tested feature is
     v = s * x_d;
  2. plant M' = M + theta * v with theta >= 0 chosen so that the planted partial R^2,
     theta^2 |r_v|^2 / (theta^2 |r_v|^2 + RSS(M ~ f(C) + v)), equals R (raw scale, not divided
     by reliability);
  3. run every test of v on M' (B = 999) and record p, theta-hat and CI coverage of theta.

Every test assumes the decoy itself carries no signal in the real M, so theta-hat estimates
theta + theta_0(v), with theta_0 the decoy's own coefficient in the real M (recorded). The
random sign makes theta_0 average to zero over replicates.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from . import aim2, crt, screens
from . import covariates as cov

ROOT = screens.ROOT
LEVELS = (0.0, 0.0025, 0.005, 0.01, 0.02)
TESTS = ("block_x4", "block_v1.0", "naive_ols", "camera", "unrestricted", "expr_bin", "stratified")
MIN_DECOYS = 5


def select_decoys(P: dict, label: str = "combined_cpm", kind: str = "combined", signal: str = "cpm"):
    D, labels = aim2.decoy_matrix(P["genes"], kind, signal)
    sim = crt.decoy_similarity(P["X"][label], D, P["Q"])
    ok = np.abs(sim) < aim2.NEG_CONTROL_SIM
    rule = "|similarity| < 0.05"
    if ok.sum() < MIN_DECOYS:
        ok = np.zeros(len(sim), bool)
        ok[np.argsort(np.abs(sim))[:MIN_DECOYS]] = True
        rule = f"fewer than {MIN_DECOYS} with |similarity| < 0.05: the {MIN_DECOYS} lowest |similarity| used"
    return D, labels, sim, ok, rule


def _one(args):
    (M, Z_Q, chrom, pos, arm, C_expr, C_cn, rho_bar, v, k_x4, k_v10, R, B, seed, order) = args
    rng = np.random.default_rng(seed)
    Q = Z_Q
    rv = v - Q @ (Q.T @ v)
    ssv = float(rv @ rv)
    rM = M - Q @ (Q.T @ M)
    theta0 = float(rv @ rM) / ssv
    rss0 = float(rM @ rM) - theta0 ** 2 * ssv
    theta = float(np.sqrt(R * rss0 / ((1 - R) * ssv))) if R > 0 else 0.0
    Mp = M + theta * v
    prep = crt.prepare(Mp, v, Q=Q)
    row = {"theta": theta, "theta0": theta0, "theta_hat": prep.theta_hat,
           "partial_r2": crt.partial_r2(Mp, v, prep=prep), "k_x4": k_x4, "k_v1.0": k_v10}
    res = {
        "block_x4": crt.block_permutation_test(Mp, v, chrom=chrom, pos=pos, arm=arm, k=k_x4, B=B, rng=rng,
                                               prep=prep, order=order),
        "block_v1.0": crt.block_permutation_test(Mp, v, chrom=chrom, pos=pos, arm=arm, k=k_v10, B=B,
                                                 rng=rng, prep=prep, order=order),
        "unrestricted": crt.unrestricted_permutation_test(Mp, v, B=B, rng=rng, prep=prep),
        "expr_bin": crt.expression_bin_permutation_test(Mp, v, expr_diff=C_expr, B=B, rng=rng, prep=prep),
        "stratified": crt.stratified_permutation_test(Mp, v, strata=crt.make_strata(C_expr, C_cn), B=B,
                                                      rng=rng, prep=prep),
    }
    for name, r in res.items():
        row[f"p_{name}"] = r.p
        row[f"cover_{name}"] = r.covers(theta)
    for name, r in (("naive_ols", crt.naive_ols_test(Mp, v, prep=prep)),
                    ("camera", crt.camera_vif_test(Mp, v, prep=prep, rho_bar=rho_bar))):
        row[f"p_{name}"] = r["p"]
        row[f"cover_{name}"] = bool(r["ci"][0] <= theta <= r["ci"][1])
    return row


def decoy_theta0_tests(P: dict, D: np.ndarray, labels, idx, B: int = 9999) -> dict:
    """Each used decoy's own association with the real M (theta_0): main block test (x4 rule,
    k from the decoy's residuals, B draws) and naive OLS, both two-sided. Used to report the
    level-0 FPR also without configurations whose decoy is itself associated with M."""
    gen = P["gen"]
    chrom, pos, arm = gen.chrom.to_numpy(), gen.pos.to_numpy(), gen.arm.to_numpy()
    order = crt.genome_order(chrom, pos)
    out = {}
    for j in idx:
        v = D[:, j]
        prep = crt.prepare(P["M"], v, Q=P["Q"])
        k = crt.choose_k(prep.r_x[order], arm[order], rule="x4")
        r = crt.block_permutation_test(P["M"], v, chrom=chrom, pos=pos, arm=arm, k=k, B=B,
                                       rng=aim2.job_seed("theta0", P["tag"], labels[j]), prep=prep, order=order)
        n = crt.naive_ols_test(P["M"], v, prep=prep)
        out[labels[j]] = {"theta0": prep.theta_hat, "k_x4": int(k), "p_block_x4_two_sided": float(r.p_two_sided),
                          "p_naive_two_sided": float(n["p_two_sided"]), "B": B}
    return out


def run(a: str = "RPE1", b: str = "K562", reps: int = 500, B: int = 999, n_jobs: int = 10,
        rho_bar: float | None = None, out_dir: Path = ROOT / "results/aim3",
        arm_effects: bool = cov.ARM_EFFECTS_DEFAULT, suffix: str | None = None,
        theta0_B: int = 9999) -> pd.DataFrame:
    """Writes planted_<tag><suffix>_replicates.parquet and _info.json; suffix defaults to
    aim2.plan_suffix(arm_effects) ('_v13' under plan v1.3, '' = the v1.2 files)."""
    from joblib import Parallel, delayed
    P = aim2.load_pair(a, b, arm_effects=arm_effects)
    tag = P["tag"]
    suffix = aim2.plan_suffix(arm_effects) if suffix is None else suffix
    if rho_bar is None:
        rho_bar = float(np.mean([aim2.rho_bar_controls(ln, P["genes"]) for ln in (a, b)]))
    D, labels, sim, ok, rule = select_decoys(P)
    idx = np.flatnonzero(ok)
    gen = P["gen"]
    chrom, pos, arm = gen.chrom.to_numpy(), gen.pos.to_numpy(), gen.arm.to_numpy()
    order = crt.genome_order(chrom, pos)
    Q = P["Q"]
    ks = {}
    for j in idx:
        r = D[:, j] - Q @ (Q.T @ D[:, j])
        ks[j] = (crt.choose_k(r[order], arm[order], rule="x4"), crt.choose_k(r[order], arm[order], rule="v1.0"))
    rng = aim2.job_seed("planted", tag)
    tasks, meta = [], []
    for li, R in enumerate(LEVELS):
        for rep in range(reps):
            j = int(rng.choice(idx))
            s = int(rng.choice([-1, 1]))
            seed = int(rng.integers(2 ** 63 - 1))
            tasks.append((P["M"], Q, chrom, pos, arm, P["C"].expr_diff.to_numpy(), P["C"].cn_diff.to_numpy(),
                          rho_bar, s * D[:, j], ks[j][0], ks[j][1], R, B, seed, order))
            meta.append({"level": R, "rep": rep, "decoy": labels[j], "sign": s, "similarity": float(sim[j])})
    rows = Parallel(n_jobs=n_jobs, verbose=2, batch_size=8)(delayed(_one)(t) for t in tasks)
    df = pd.concat([pd.DataFrame(meta), pd.DataFrame(rows)], axis=1)
    df["B"] = B
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out_dir / f"planted_{tag}{suffix}_replicates.parquet", index=False)
    info = {"pair": tag, "decoy_rule": rule, "n_decoys_total": len(labels), "n_decoys_used": int(ok.sum()),
            "decoys_used": [labels[j] for j in idx], "similarity": dict(zip(labels, np.round(sim, 4).tolist())),
            "k_by_decoy": {labels[j]: {"x4": ks[j][0], "v1.0": ks[j][1]} for j in idx},
            "rho_bar": rho_bar, "reps_per_level": reps, "B": B, "levels": list(LEVELS)}
    if arm_effects:
        info["plan"] = "v1.3"
        info["arm_effects"] = True
        info["theta0_tests"] = decoy_theta0_tests(P, D, labels, idx, B=theta0_B)
    (out_dir / f"planted_{tag}{suffix}_info.json").write_text(json.dumps(info, indent=1))
    return df


def _rate(v):
    v = np.asarray(v, float)
    p = v.mean()
    return p, np.sqrt(p * (1 - p) / v.size)


def _cluster_se(v, cl):
    """SE of a mean over replicates that share a (decoy, sign) configuration: the replicates
    within a configuration differ only in permutation draws, so the configuration is the unit."""
    s = pd.Series(np.asarray(v, float)).groupby(np.asarray(cl))
    m, n = s.mean(), s.size()
    if len(m) < 2:
        return np.nan
    w = n / n.sum()
    mu = float((w * m).sum())
    # ratio-estimator (cluster-robust) variance of a weighted mean
    return float(np.sqrt((w ** 2 * (m - mu) ** 2).sum() * len(m) / (len(m) - 1)))


def summarize(df: pd.DataFrame, tests=TESTS) -> pd.DataFrame:
    df = df.assign(config=df.decoy + ":" + df.sign.astype(str))
    null = df[df.level == 0]
    rows = []
    for R, g in df.groupby("level"):
        bias = g.theta_hat - g.theta
        for t in tests:
            p = g[f"p_{t}"].to_numpy()
            row = {"level_partial_r2": R, "test": t, "n": len(g), "n_configs": g.config.nunique(),
                   "theta_mean": g.theta.mean(), "bias": bias.mean(),
                   "bias_se": bias.std(ddof=1) / np.sqrt(len(g)), "bias_se_cluster": _cluster_se(bias, g.config),
                   "mean_partial_r2_obs": g.partial_r2.mean()}
            row["rej05"], row["rej05_se"] = _rate(p <= 0.05)
            row["rej05_se_cluster"] = _cluster_se(p <= 0.05, g.config)
            row["rej01"], row["rej01_se"] = _rate(p <= 0.01)
            row["rej01_se_cluster"] = _cluster_se(p <= 0.01, g.config)
            row["coverage"], row["coverage_se"] = _rate(g[f"cover_{t}"].to_numpy())
            row["coverage_se_cluster"] = _cluster_se(g[f"cover_{t}"].to_numpy(), g.config)
            thr = np.quantile(null[f"p_{t}"].to_numpy(), 0.05, method="lower")
            row["matched_thr05"] = thr
            row["power_matched"], row["power_matched_se"] = _rate(p <= thr)
            rows.append(row)
    return pd.DataFrame(rows)


def theta0_significant(info: dict, col: str = "p_block_x4_two_sided", alpha: float = 0.05) -> list:
    """Decoys whose own theta_0 in the real M is significantly nonzero (two-sided, both signs)."""
    return sorted(d for d, r in info.get("theta0_tests", {}).items() if r[col] <= alpha)


def detection_limit(summ: pd.DataFrame, test: str = "block_x4", col: str = "rej05", power: float = 0.8):
    s = summ[(summ.test == test) & (summ.level_partial_r2 > 0)].sort_values("level_partial_r2")
    hit = s[s[col] >= power]
    return float(hit.level_partial_r2.iloc[0]) if len(hit) else None
