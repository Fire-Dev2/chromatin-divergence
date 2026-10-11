"""Aim 3 synthetic benchmark (ANALYSIS_PLAN.md section 6).

Generative model (one replicate)
--------------------------------
Genome: G genes on 22 autosomes x 2 arms (arm sizes random), in genome order.
Modules: contiguous runs of genes within an arm (sizes 1 + Geometric, mean
``module_mean``), standing in for co-regulated / co-accessible neighbourhoods.

Covariates C (per gene): log2 mean expression in line A and B (e_A, e_B,
correlated 0.85), log2 variance in each line, expression difference d = e_A - e_B,
copy-number difference cn (arm-level gains/losses).

Feature  x = kappa * z(d) + sqrt(1 - kappa^2) * s,
         s = sqrt(rho_x) * u_module + sqrt(1 - rho_x) * a_gene         (Var s = 1)
Response M = f(C) + theta * x + sigma_g * (sqrt(rho_M) * v_module + sqrt(1 - rho_M) * e_gene)
         f(C) = 0.4 tanh(e_A - 2) - 0.4 tanh(e_B - 2) + 0.5 z(d) + 0.3 cn
u, v, a, e iid N(0, 1); u and v independent (x shares M's modules, not its values).
kappa = feature-expression confounding (corr of x with d).
sigma_g = 1 (homoscedastic) or sigma_g^2 proportional to 0.3 + 100 / (n_eff_g * 2^{mean expr}),
          n_eff_g ~ LogNormal(log 100, 0.5) per-gene effective cell count; scaled to mean 1.
theta gives partial R^2 = R2 for x: theta = sqrt(R2 / (1 - R2) / (1 - kappa^2)).

Decoys: independent draws of x (same kappa, same modules, new u and a) - i.e.
features exchangeable with x under the null.

Tests per replicate: naive OLS, CAMERA VIF (rho_bar = true mean inter-gene
correlation of the M noise), unrestricted permutation, expression-bin permutation,
stratified permutation (expr bins x CN state), decoy null, block permutation with
k = choose_k (main), wrong blocks k/4, 4k, random gene order, and oracle blocks
(the true modules).

CLI
---
python -m cdtest.simulate --reps 1000 --B 999 --n-jobs 14 \
       --out results/aim3/synthetic_grid.parquet
"""
from __future__ import annotations

import argparse
import itertools
import os
import time
from dataclasses import asdict, dataclass

import numpy as np
import pandas as pd

from . import crt

MASTER_SEED = 20261011
R2_LEVELS = {"zero": 0.0, "small": 0.005, "moderate": 0.02, "large": 0.05}
RHO_M = (0.0, 0.1, 0.3, 0.5)
RHO_X = (0.0, 0.3, 0.6)
CONFOUND = (0.0, 0.3)
NOISE = ("homo", "hetero")

PERM_TESTS = ("unrestricted", "expr_bin", "stratified", "decoy", "block_main",
              "block_k_div4", "block_k_x4", "block_random_order", "block_oracle")
PARAM_TESTS = ("naive_ols", "camera")
ALL_TESTS = PARAM_TESTS + PERM_TESTS


@dataclass(frozen=True)
class Setting:
    rho_M: float
    rho_x: float
    confound: float
    noise: str
    effect: str
    G: int = 2000
    module_mean: float = 10.0

    @property
    def r2(self) -> float:
        return R2_LEVELS[self.effect]

    @property
    def theta(self) -> float:
        r2 = self.r2
        return float(np.sqrt(r2 / (1 - r2) / (1 - self.confound ** 2)))


def grid_settings(G: int = 2000, module_mean: float = 10.0) -> list[Setting]:
    return [Setting(a, b, c, n, e, G, module_mean)
            for a, b, c, n, e in itertools.product(RHO_M, RHO_X, CONFOUND, NOISE, R2_LEVELS)]


# --------------------------------------------------------------------------
# Generator
# --------------------------------------------------------------------------

def make_genome(G: int, rng, module_mean: float = 10.0) -> pd.DataFrame:
    n_arms = 44
    w = rng.gamma(4.0, size=n_arms)
    sizes = np.maximum(1, np.floor(w / w.sum() * G).astype(int))
    sizes[-1] += G - sizes.sum()
    if sizes[-1] < 1:  # pragma: no cover - only for tiny G
        sizes = np.full(n_arms, G // n_arms)
        sizes[-1] += G - sizes.sum()
    arm_id = np.repeat(np.arange(n_arms), sizes)
    chrom = arm_id // 2 + 1
    arm = np.array([f"chr{c}{'pq'[a % 2]}" for c, a in zip(chrom, arm_id)])
    pos = np.cumsum(rng.exponential(1.0, size=G))  # increasing within each chromosome
    module = np.empty(G, dtype=np.int64)
    p = 1.0 / module_mean
    nxt, start = 0, 0
    for s in sizes:
        i = 0
        while i < s:
            L = min(int(rng.geometric(p)), s - i)
            module[start + i:start + i + L] = nxt
            nxt += 1
            i += L
        start += s
    return pd.DataFrame({"chrom": [f"chr{c}" for c in chrom], "arm": arm, "pos": pos,
                         "module": module})


def _z(v):
    return (v - v.mean()) / v.std()


def simulate_replicate(s: Setting, rng, n_decoys: int = 99) -> dict:
    G = s.G
    gen = make_genome(G, rng, s.module_mean)
    mod = gen["module"].to_numpy()
    n_mod = mod.max() + 1

    # covariates
    eA = rng.normal(2.0, 1.2, G)
    eB = 2.0 + 0.85 * (eA - 2.0) + np.sqrt(1 - 0.85 ** 2) * rng.normal(0, 1.2, G)
    vA = 1.1 * eA + rng.normal(0, 0.3, G)
    vB = 1.1 * eB + rng.normal(0, 0.3, G)
    d = eA - eB
    arm_codes = pd.factorize(gen["arm"])[0]
    arm_state = rng.choice([-1.0, 0.0, 0.585], size=arm_codes.max() + 1, p=[0.1, 0.8, 0.1])
    cn = arm_state[arm_codes]
    C = pd.DataFrame({"eA": eA, "eB": eB, "vA": vA, "vB": vB, "d": d, "cn": cn})
    zd = _z(d)

    def draw_feature():
        u = rng.normal(size=n_mod)[mod]
        a = rng.normal(size=G)
        sx = np.sqrt(s.rho_x) * u + np.sqrt(1 - s.rho_x) * a
        return s.confound * zd + np.sqrt(1 - s.confound ** 2) * sx

    x = draw_feature()
    decoys = np.column_stack([draw_feature() for _ in range(n_decoys)]) if n_decoys else None

    # noise
    if s.noise == "homo":
        sig2 = np.ones(G)
    else:
        n_eff = rng.lognormal(np.log(100), 0.5, G)
        ebar = 0.5 * (eA + eB)
        sig2 = 0.3 + 100.0 / (n_eff * 2.0 ** ebar)
        sig2 = sig2 / sig2.mean()
    v = rng.normal(size=n_mod)[mod]
    e = rng.normal(size=G)
    eta = np.sqrt(sig2) * (np.sqrt(s.rho_M) * v + np.sqrt(1 - s.rho_M) * e)
    fC = 0.4 * np.tanh(eA - 2) - 0.4 * np.tanh(eB - 2) + 0.5 * zd + 0.3 * cn
    M = fC + s.theta * x + eta

    sizes = np.bincount(mod)
    rho_bar = s.rho_M * float((sizes * (sizes - 1)).sum()) / (G * (G - 1))
    return {"M": M, "x": x, "C": C, "genome": gen, "decoys": decoys, "theta": s.theta,
            "rho_bar": rho_bar}


# --------------------------------------------------------------------------
# Analysis of one replicate
# --------------------------------------------------------------------------

SPLINE_COLS = ("eA", "eB")


def analyze_replicate(sim: dict, B: int, rng) -> dict:
    M, x, C, gen = sim["M"], sim["x"], sim["C"], sim["genome"]
    theta = sim["theta"]
    Z = crt.build_design(C, SPLINE_COLS, df=4)
    prep = crt.prepare(M, x, Z)
    G = prep.G
    out = {"theta_hat": prep.theta_hat, "partial_r2": crt.partial_r2(M, x, prep=prep)}

    # parametric comparators
    for name, res in (("naive_ols", crt.naive_ols_test(M, x, prep=prep)),
                      ("camera", crt.camera_vif_test(M, x, prep=prep, rho_bar=sim["rho_bar"]))):
        out[f"p_{name}"] = res["p"]
        out[f"p2_{name}"] = res["p_two_sided"]
        lo, hi = res["ci"]
        out[f"cover_{name}"] = bool(lo <= theta <= hi)

    chrom, arm, pos = gen["chrom"].to_numpy(), gen["arm"].to_numpy(), gen["pos"].to_numpy()
    order = crt.genome_order(chrom, pos)
    k = crt.choose_k(prep.r_x[order], arm[order])
    out["k"] = k
    strata = crt.make_strata(C["d"].to_numpy(), C["cn"].to_numpy())
    rand_order = rng.permutation(G)

    results = {
        "unrestricted": crt.unrestricted_permutation_test(M, x, B=B, rng=rng, prep=prep),
        "expr_bin": crt.expression_bin_permutation_test(M, x, expr_diff=C["d"].to_numpy(), B=B,
                                                        rng=rng, prep=prep),
        "stratified": crt.stratified_permutation_test(M, x, strata=strata, B=B, rng=rng, prep=prep),
        "decoy": crt.decoy_test(M, x, decoys=sim["decoys"], prep=prep),
        "block_main": crt.block_permutation_test(M, x, chrom=chrom, pos=pos, arm=arm, k=k, B=B,
                                                 rng=rng, prep=prep, order=order),
        "block_k_div4": crt.block_permutation_test(M, x, chrom=chrom, pos=pos, arm=arm,
                                                   k=max(1, int(round(k / 4))), B=B, rng=rng,
                                                   prep=prep, order=order),
        "block_k_x4": crt.block_permutation_test(M, x, chrom=chrom, pos=pos, arm=arm, k=4 * k,
                                                 B=B, rng=rng, prep=prep, order=order),
        "block_random_order": crt.block_permutation_test(M, x, k=k, B=B, rng=rng, prep=prep,
                                                         order=rand_order),
        "block_oracle": crt.block_permutation_test(M, x, B=B, rng=rng, prep=prep, order=order,
                                                   blocks=gen["module"].to_numpy()),
    }
    for name, r in results.items():
        out[f"p_{name}"] = r.p
        out[f"p2_{name}"] = r.p_two_sided
        out[f"cover_{name}"] = r.covers(theta)
    # plan's grid CI for the main test (width and grid-based coverage)
    lo, hi, info = results["block_main"].ci(sd_x=float(np.std(x)))
    out["ci_lo_block_main"], out["ci_hi_block_main"] = lo, hi
    out["ci_edge_block_main"] = info.get("at_edge", False)
    return out


def _run_batch(setting: Setting, rep_ids, B: int, n_decoys: int, seed: int,
               part_dir: str | None = None) -> list[dict]:
    part = None
    if part_dir:
        part = os.path.join(part_dir, f"s{_setting_key(setting):03d}_G{setting.G}_B{B}_"
                                      f"r{rep_ids[0]:05d}-{rep_ids[-1]:05d}.parquet")
        if os.path.exists(part):
            return pd.read_parquet(part).to_dict("records")
    rows = []
    for r in rep_ids:
        ss = np.random.SeedSequence([seed, _setting_key(setting), r])
        rng = np.random.default_rng(ss)
        sim = simulate_replicate(setting, rng, n_decoys)
        row = analyze_replicate(sim, B, rng)
        row.update(asdict(setting))
        row.update({"rep": r, "theta": setting.theta, "B": B})
        rows.append(row)
    if part:
        tmp = part + ".tmp"
        pd.DataFrame(rows).to_parquet(tmp, index=False)
        os.replace(tmp, part)
    return rows


def _setting_key(s: Setting) -> int:
    key = (RHO_M.index(s.rho_M), RHO_X.index(s.rho_x), CONFOUND.index(s.confound),
           NOISE.index(s.noise), list(R2_LEVELS).index(s.effect))
    return int(np.ravel_multi_index(key, (4, 3, 2, 2, 4)))


def run_grid(reps: int = 200, B: int = 999, G: int = 2000, n_jobs: int = 14, n_decoys: int = 99,
             seed: int = MASTER_SEED, batch: int = 25, settings=None, verbose: int = 5,
             part_dir: str | None = None) -> pd.DataFrame:
    """Run the grid in parallel. With ``part_dir`` every batch is checkpointed and a
    rerun resumes (replicate r of a setting always uses the same seed, so a
    200-replicate run is the first 200 replicates of the 1000-replicate run)."""
    from joblib import Parallel, delayed

    if part_dir:
        os.makedirs(part_dir, exist_ok=True)
    settings = settings or grid_settings(G)
    tasks = [(s, list(range(i, min(i + batch, reps))))
             for s in settings for i in range(0, reps, batch)]
    out = Parallel(n_jobs=n_jobs, verbose=verbose)(
        delayed(_run_batch)(s, ids, B, n_decoys, seed, part_dir) for s, ids in tasks)
    return pd.DataFrame([row for rows in out for row in rows])


# --------------------------------------------------------------------------
# Scoring
# --------------------------------------------------------------------------

def _rate(v):
    v = np.asarray(v, float)
    n = v.size
    p = v.mean() if n else np.nan
    return p, np.sqrt(p * (1 - p) / n) if n else np.nan


def summarize(df: pd.DataFrame, tests=ALL_TESTS) -> pd.DataFrame:
    """Per setting x test: FPR@5%/1% (null), power at 5% and at matched FPR, bias, coverage; MC SEs."""
    keys = ["rho_M", "rho_x", "confound", "noise"]
    df = df.assign(effect=df.effect.replace({"null": "zero"}))  # old label ("null" reads as NA in CSV)
    rows = []
    for kv, g in df.groupby(keys):
        null = g[g.effect == "zero"]
        for eff, ge in g.groupby("effect"):
            bias = ge.theta_hat - ge.theta
            for t in tests:
                p = ge[f"p_{t}"].to_numpy()
                row = dict(zip(keys, kv))
                row.update({"effect": eff, "test": t, "n": len(ge), "theta": ge.theta.iloc[0],
                            "mean_k": ge.k.mean(),
                            "bias": bias.mean(), "bias_se": bias.std(ddof=1) / np.sqrt(len(ge)),
                            "mean_partial_r2": ge.partial_r2.mean()})
                row["rej05"], row["rej05_se"] = _rate(p <= 0.05)
                row["rej01"], row["rej01_se"] = _rate(p <= 0.01)
                row["coverage"], row["coverage_se"] = _rate(ge[f"cover_{t}"].to_numpy())
                # power at matched false-positive rate: threshold = 5% quantile of null p's
                if len(null):
                    thr = np.quantile(null[f"p_{t}"].to_numpy(), 0.05, method="lower")
                    row["matched_thr"] = thr
                    row["power_matched"], row["power_matched_se"] = _rate(p <= thr)
                rows.append(row)
    return pd.DataFrame(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reps", type=int, default=200, help="replicates per setting (plan: 1000)")
    ap.add_argument("--B", type=int, default=999, help="permutations per test")
    ap.add_argument("--G", type=int, default=2000, help="genes per replicate")
    ap.add_argument("--n-decoys", type=int, default=99)
    ap.add_argument("--n-jobs", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--seed", type=int, default=MASTER_SEED)
    ap.add_argument("--batch", type=int, default=25)
    ap.add_argument("--out", default="results/aim3/synthetic_grid.parquet")
    ap.add_argument("--summary", default="results/aim3/synthetic_grid_summary.csv")
    ap.add_argument("--part-dir", default="results/aim3/_parts",
                    help="checkpoint directory (resumable); '' to disable")
    a = ap.parse_args(argv)
    t0 = time.time()
    df = run_grid(a.reps, a.B, a.G, a.n_jobs, a.n_decoys, a.seed, a.batch,
                  part_dir=a.part_dir or None)
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    if a.out.endswith(".parquet"):
        df.to_parquet(a.out, index=False)
    else:
        df.to_csv(a.out, index=False)
    summarize(df).to_csv(a.summary, index=False)
    print(f"{len(df)} rows in {time.time() - t0:.0f} s -> {a.out}, {a.summary}")


if __name__ == "__main__":
    main()
