"""Aim 2 on one pair: the block-permutation test of M ~ f(C) + x and every comparator.

Inputs (all produced earlier in the pipeline):
  results/aim1/<tag>/gene_outcome.tsv   M per gene (A minus B)
  results/aim1/<tag>/terms.npz          full-data per-(knockdown, gene) terms (bootstrap)
  data/aim1_cache/<tag>_halfD.npy       half-split per-(knockdown, gene) terms (bootstrap)
  data/features/gene_features[_fc].parquet
  covariates.build_covariates(A, B)     C and f(C)

Only call this on pairs the plan allows at the current stage (before tag pipeline-v1.0:
the development pair RPE1-K562 only).
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from . import covariates as cov
from . import crt, features, screens

ROOT = screens.ROOT
AIM1 = ROOT / "results/aim1"
CACHE = ROOT / "data/aim1_cache"
MASTER_SEED = screens.MASTER_SEED
DECOY_LINES = ("GM12878", "HCT116", "MCF-7", "A549", "IMR-90", "PC-3", "Panc1")
#: (label, kind, signal); the first is the primary (H1) feature.
FEATURES = (("combined_cpm", "combined", "cpm"), ("promoter_cpm", "promoter", "cpm"),
            ("distal_cpm", "distal", "cpm"), ("combined_fc", "combined", "fc"))
NEG_CONTROL_SIM = 0.05
ACF_LAGS = (1, 2, 5, 10, 25, 50, 100)


def job_seed(*keys) -> np.random.Generator:
    """Generator derived from the master seed and string/int keys (plan §8)."""
    ints = [k if isinstance(k, int) else sum((i + 1) * ord(c) for i, c in enumerate(str(k))) for k in keys]
    return np.random.default_rng(np.random.SeedSequence([MASTER_SEED, *ints]))


# --------------------------------------------------------------------------- data
def plan_suffix(arm_effects: bool) -> str:
    """File suffix of a result: '_v13' for plan v1.3 (arm effects), '' for v1.2 files."""
    return "_v13" if arm_effects else ""


def load_pair(a: str, b: str, tag: str | None = None, arm_effects: bool = cov.ARM_EFFECTS_DEFAULT) -> dict:
    """Everything one pair's tests need. ``arm_effects`` True = plan v1.3 f(C) (arm indicators
    in the model and in the residualisation of x, hence in every test, comparator, decoy
    similarity and k choice that uses Q); False = plan v1.2."""
    tag = tag or f"{a}_vs_{b}"
    gs = pd.read_csv(AIM1 / tag / "gene_outcome.tsv", sep="\t")
    genes = gs.gene_id.tolist()
    C = cov.build_covariates(a, b, genes, tag)
    gen = cov.genome_info(genes)
    Z = cov.design(C, arm_effects=arm_effects, arm=gen.arm.to_numpy())
    X = {}
    for label, kind, signal in FEATURES:
        X[label] = features.pair_feature(a, b, kind, signal=signal).reindex(genes).to_numpy(float)
        if not np.isfinite(X[label]).all():
            raise ValueError(f"{label}: missing feature values")
    return dict(a=a, b=b, tag=tag, genes=genes, M=gs.M.to_numpy(float), gs=gs, C=C, Z=Z, gen=gen, X=X,
                arm_effects=bool(arm_effects),
                Q=crt.orthonormal_basis(Z), summary=json.loads((AIM1 / tag / "summary.json").read_text()))


def decoy_matrix(genes, kind: str = "combined", signal: str = "cpm", lines=DECOY_LINES):
    """Feature differences between all pairs of decoy lines (G, 21) and their labels."""
    m = features.feature_matrix(kind, signal=signal).reindex(genes)
    labels, cols = [], []
    for u, v in itertools.combinations(lines, 2):
        labels.append(f"{u}-{v}")
        cols.append((m[u] - m[v]).to_numpy(float))
    return np.column_stack(cols), labels


def rho_bar_controls(line: str, genes, log: bool = False, max_cells: int | None = None) -> float:
    """Mean off-diagonal inter-gene Pearson correlation over the line's control cells.

    sum_{g,h} R_gh = (1/(n-1)) sum_cells (sum_g z_cg)^2, so the mean off-diagonal entry is
    (that - G) / (G (G - 1)); no G x G matrix is formed.
    """
    d = screens.OUT / line
    gt = pd.read_csv(d / "genes.tsv", sep="\t", index_col=0)
    pos = pd.Series(np.arange(len(gt)), index=gt.gene_id.values).loc[list(genes)].to_numpy()
    X = np.asarray(np.load(d / "ctrl_norm.npy", mmap_mode="r")[:, pos], dtype=np.float64)
    if max_cells:
        X = X[:max_cells]
    if log:
        X = np.log2(X + 1)
    sd = X.std(0, ddof=1)
    Zs = (X - X.mean(0)) / np.where(sd > 0, sd, 1)
    n, G = Zs.shape
    tot = float((Zs.sum(1) ** 2).sum() / (n - 1))
    return (tot - G) / (G * (G - 1))


def coaccess_clusters(genes, signal: str = "cpm", cut: float = 0.5) -> np.ndarray:
    """Secondary block definition (plan §5): average-linkage clusters of genes on their
    promoter + distal profile across the 11 lines, distance 1 - Pearson, cut at 0.5.
    Profile = [z(promoter) over 11 lines, z(distal) over 11 lines] with the pooled z of the
    combined feature [choice]."""
    from scipy.cluster.hierarchy import fcluster, linkage
    prom = features.feature_matrix("promoter", signal=signal).reindex(genes)
    dist = features.feature_matrix("distal", signal=signal).reindex(genes)
    zp = (prom - np.nanmean(prom.values)) / np.nanstd(prom.values)
    zd = (dist - np.nanmean(dist.values)) / np.nanstd(dist.values)
    P = np.hstack([zp.values, zd.values])
    flat = P.std(1) == 0
    P[flat] += np.random.default_rng(0).normal(0, 1e-9, P[flat].shape)
    L = linkage(P, method="average", metric="correlation")
    return fcluster(L, t=cut, criterion="distance")


# --------------------------------------------------------------------------- naive (Aim 1 step 7)
def naive_expr_only(M, x, expr_diff) -> dict:
    """OLS M ~ 1 + expression difference + x, genes independent."""
    Z = np.column_stack([np.ones_like(M), expr_diff])
    prep = crt.prepare(M, x, Z)
    r = crt.naive_ols_test(M, x, prep=prep)
    return {"coef": r["theta"], "se": r["se"], "t": r["t"], "p_one_sided": r["p"],
            "p_two_sided": r["p_two_sided"], "partial_r2": crt.partial_r2(M, x, prep=prep),
            "ci95": list(r["ci"]), "n_genes": len(M)}


def naive_analysis(P: dict) -> dict:
    out = {"pair": P["tag"], "A": P["a"], "B": P["b"], "model": "M ~ 1 + expr_diff + x (OLS, genes independent)",
           "note": "Aim 1 step 7; descriptive only. Development pair = mixed pair (p53, tissue, growth)."}
    for label, x in P["X"].items():
        out[label] = naive_expr_only(P["M"], x, P["C"].expr_diff.to_numpy())
    return out


# --------------------------------------------------------------------------- bootstrap
def bootstrap_reliability_M(tag: str, n_boot: int = 500, rng=None, half_path: Path | None = None):
    """Knockdown bootstrap (approximation, docs/methods_aim2_aim3.md §6.3).

    Resamples kept knockdowns with replacement and recomputes M (full data) and the
    split-half reliability (Spearman-Brown of the median split r over the stored splits)
    from the stored per-(knockdown, gene) terms. beta~ and s~^2 are held at their
    full-data values (the generic mean and the fake-screen null are not recomputed).
    Returns (M_boot (n_boot, G), rel_boot (n_boot,), W).
    """
    rng = rng if rng is not None else job_seed("boot", tag)
    t = np.load(AIM1 / tag / "terms.npz", allow_pickle=True)
    mask = t["mask"]
    D = (t["btA"].astype(np.float64) ** 2 - t["s2A"]) - (t["btB"].astype(np.float64) ** 2 - t["s2B"])
    D = np.where(mask, 0.0, D)
    keep = (~mask).astype(np.float64)
    P = D.shape[0]
    W = rng.multinomial(P, np.full(P, 1.0 / P), size=n_boot).astype(np.float64)
    Mb = (W @ D) / (W @ keep)
    half_path = half_path or CACHE / f"{tag}_halfD.npy"
    H = np.load(half_path, mmap_mode="r")
    S = H.shape[0]
    R = np.empty((n_boot, S))
    den = W @ keep
    for s in range(S):
        m = []
        for h in range(2):
            Dh = np.nan_to_num(np.asarray(H[s, h], dtype=np.float64), nan=0.0)
            m.append((W @ Dh) / den)
        a, b = m
        a = a - a.mean(1, keepdims=True); b = b - b.mean(1, keepdims=True)
        R[:, s] = (a * b).sum(1) / np.sqrt((a * a).sum(1) * (b * b).sum(1))
    rmed = np.median(R, axis=1)
    rel = 2 * rmed / (1 + rmed)
    # full-data check: unweighted halves must reproduce the stored split r's
    return Mb, rel, W


def partial_r2_many(Mb: np.ndarray, x: np.ndarray, Q: np.ndarray) -> np.ndarray:
    rx = x - Q @ (Q.T @ x)
    RM = Mb - (Mb @ Q) @ Q.T
    num = (RM @ rx) ** 2 / (rx @ rx)
    return num / (RM * RM).sum(1)


# --------------------------------------------------------------------------- test
def ci_refined(res: crt.TestResult, grid, n_fine: int = 2001):
    """Test-inversion CI on a refined grid inside the plan grid's span [choice, docs §6.2].

    The plan grid (201 points over +/- max(4|theta|, 0.5/SD(x))) is too coarse when theta is
    small relative to 0.5/SD(x): the step can exceed the interval width. The same acceptance
    rule (each one-sided p > 0.025, same B draws) is evaluated on n_fine points spanning one
    plan step beyond the plan-grid interval (or +/- 2 steps around theta-hat if empty).
    """
    grid = np.asarray(grid, float)
    step = grid[1] - grid[0]
    lo, hi, _ = res.ci(grid=grid)
    if not np.isfinite(lo):
        lo = hi = res.T
        pad = 2 * step
    else:
        pad = step
    fine = np.linspace(max(grid[0], lo - pad), min(grid[-1], hi + pad), n_fine)
    flo, fhi, info = res.ci(grid=fine)
    info["step"] = float(fine[1] - fine[0])
    return flo, fhi, info


def block_exchangeability(gen: pd.DataFrame, k: int) -> dict:
    """How much of the genome the block permutation can actually move (size classes >= 2 blocks)."""
    order = crt.genome_order(gen.chrom.to_numpy(), gen.pos.to_numpy())
    bid = crt.make_blocks(gen.arm.to_numpy()[order], int(k))
    classes = crt._block_classes(bid)
    movable = [(nb, s) for nb, s, _ in classes if nb >= 2]
    return {"n_blocks": int(bid.max() + 1), "n_size_classes": len(classes),
            "n_movable_blocks": int(sum(nb for nb, _ in movable)),
            "n_movable_classes": len(movable),
            "frac_genes_movable": float(sum(nb * s for nb, s in movable) / len(bid))}


def _summ(res: crt.TestResult, grid) -> dict:
    lo, hi, info = res.ci(grid=grid)
    flo, fhi, finfo = ci_refined(res, grid)
    return {"T": res.T, "p_one_sided": res.p, "p_two_sided": res.p_two_sided, "B": res.B,
            "ci95_plan_grid": [lo, hi], "ci_info": info, "ci95": [flo, fhi], "ci_fine_info": finfo,
            "null_sd": float(np.nanstd(res.Tstar)), "null_mean": float(np.nanmean(res.Tstar)),
            **{k: v for k, v in res.info.items() if k != "sd_x"}}


def _param(r: dict) -> dict:
    return {"T": r["theta"], "se": r["se"], "t": r["t"], "p_one_sided": r["p"],
            "p_two_sided": r["p_two_sided"], "ci95": list(r["ci"]), **({"vif": r["vif"]} if "vif" in r else {})}


def test_feature(P: dict, label: str, B: int = 9999, rho_bar: float = 0.0, decoys=None,
                 clusters=None, boot=None, keep_null: bool = False) -> dict:
    M, x, Z, Q, gen, C = P["M"], P["X"][label], P["Z"], P["Q"], P["gen"], P["C"]
    prep = crt.prepare(M, x, Q=Q)
    chrom, pos, arm = gen.chrom.to_numpy(), gen.pos.to_numpy(), gen.arm.to_numpy()
    order = crt.genome_order(chrom, pos)
    k_x4 = crt.choose_k(prep.r_x[order], arm[order], rule="x4")
    k_v10 = crt.choose_k(prep.r_x[order], arm[order], rule="v1.0")
    grid = crt.theta_grid(prep.theta_hat, float(np.std(x)))
    rng = job_seed("aim2", P["tag"], label)
    main = crt.block_permutation_test(M, x, chrom=chrom, pos=pos, arm=arm, k=k_x4, B=B, rng=rng,
                                      prep=prep, order=order)
    rel = P["summary"]["reliability"]
    pr2 = crt.partial_r2(M, x, prep=prep)
    out = {"feature": label, "G": prep.G, "design_rank": prep.rank, "sd_x": float(np.std(x)),
           "theta_hat": prep.theta_hat, "k_x4": k_x4, "k_v1.0": k_v10,
           "main_block_x4": _summ(main, grid),
           "acf_x_resid_within_arm": dict(zip(map(str, ACF_LAGS), np.round(
               crt.acf(prep.r_x[order], max(ACF_LAGS), arm[order])[np.array(ACF_LAGS) - 1], 4).tolist())),
           "exchangeability_x4": block_exchangeability(gen, k_x4),
           "exchangeability_v1.0": block_exchangeability(gen, k_v10),
           "theta_grid": [float(grid[0]), float(grid[-1]), len(grid)],
           "partial_r2_raw": pr2, "reliability": rel, "partial_r2_share_reproducible": pr2 / rel}
    if boot is not None:
        Mb, relb = boot
        b = partial_r2_many(Mb, x, Q) / relb
        out["partial_r2_share_boot_ci95"] = [float(np.quantile(b, 0.025)), float(np.quantile(b, 0.975))]
        out["partial_r2_share_boot_median"] = float(np.median(b))
        rx = prep.r_x
        th_b = ((Mb - (Mb @ Q) @ Q.T) @ rx) / (rx @ rx)
        out["theta_boot_se"] = float(th_b.std(ddof=1))
    comps = {}
    comps["block_v1.0"] = _summ(crt.block_permutation_test(M, x, chrom=chrom, pos=pos, arm=arm, k=k_v10, B=B,
                                                           rng=rng, prep=prep, order=order), grid)
    comps["naive_ols"] = _param(crt.naive_ols_test(M, x, prep=prep))
    comps["unrestricted"] = _summ(crt.unrestricted_permutation_test(M, x, B=B, rng=rng, prep=prep), grid)
    comps["expr_bin"] = _summ(crt.expression_bin_permutation_test(
        M, x, expr_diff=C.expr_diff.to_numpy(), B=B, rng=rng, prep=prep), grid)
    comps["stratified"] = _summ(crt.stratified_permutation_test(
        M, x, strata=crt.make_strata(C.expr_diff.to_numpy(), C.cn_diff.to_numpy()), B=B, rng=rng,
        prep=prep), grid)
    comps["camera"] = _param(crt.camera_vif_test(M, x, prep=prep, rho_bar=rho_bar))
    comps["camera"]["rho_bar"] = rho_bar
    if decoys is not None:
        D, labels = decoys
        sim = crt.decoy_similarity(x, D, Q)
        ok = np.abs(sim) < NEG_CONTROL_SIM
        entry = {"similarity": dict(zip(labels, np.round(sim, 4).tolist())), "n_neg_control": int(ok.sum())}
        if ok.sum():
            r = crt.decoy_test(M, x, decoys=D[:, ok], prep=prep)
            entry.update({"T": r.T, "p_one_sided": r.p, "p_two_sided": r.p_two_sided,
                          "T_decoys": dict(zip(np.array(labels)[ok].tolist(), np.round(r.Tstar, 5).tolist()))})
        comps["decoy_neg_control"] = entry
    if clusters is not None:
        r = crt.block_permutation_test(M, x, blocks=clusters, B=B, rng=rng, prep=prep, order=order)
        comps["block_coaccess_clusters"] = _summ(r, grid)
        comps["block_coaccess_clusters"]["n_clusters"] = int(np.unique(clusters).size)
    out["comparators"] = comps
    if keep_null:
        out["_null"] = main.Tstar
    return out


def run(a: str = "RPE1", b: str = "K562", B: int = 9999, n_boot: int = 500,
        out_dir: Path = ROOT / "results/aim2", arm_effects: bool = cov.ARM_EFFECTS_DEFAULT,
        suffix: str | None = None) -> dict:
    """Writes <tag>_test<suffix>.json and <tag>_null<suffix>.npz; suffix defaults to
    plan_suffix(arm_effects) so v1.3 results never overwrite the v1.2 files."""
    P = load_pair(a, b, arm_effects=arm_effects)
    tag = P["tag"]
    suffix = plan_suffix(arm_effects) if suffix is None else suffix
    rho = {ln: rho_bar_controls(ln, P["genes"]) for ln in (a, b)}
    rho_log = {ln: rho_bar_controls(ln, P["genes"], log=True) for ln in (a, b)}
    rho_bar = float(np.mean(list(rho.values())))
    Mb, relb, _ = bootstrap_reliability_M(tag, n_boot)
    clusters = coaccess_clusters(P["genes"])
    res = {"pair": tag, "A": a, "B": b, "status": "development pair (mixed: p53, tissue, growth) -- descriptive, not H1",
           "B_perm": B, "n_boot": n_boot, "n_genes": len(P["genes"]),
           "n_knockdowns": P["summary"]["n_kd_kept"], "reliability": P["summary"]["reliability"],
           "reliability_boot_ci95": [float(np.quantile(relb, .025)), float(np.quantile(relb, .975))],
           "rho_bar_controls": rho, "rho_bar_controls_log2p1": rho_log, "rho_bar_used": rho_bar,
           "acf_M_resid_within_arm_diagnostic": _acf_M(P),
           "covariates": cov.COLUMNS, "design_columns": int(P["Z"].shape[1]), "features": {}}
    if arm_effects:
        res["plan"] = "v1.3"
        res["arm_effects"] = True
        res["n_arms"] = int(P["gen"].arm.nunique())
        res["design_rank"] = int(P["Q"].shape[1])
    nulls = {}
    for label, kind, signal in FEATURES:
        dec = decoy_matrix(P["genes"], kind, signal)
        r = test_feature(P, label, B, rho_bar, dec, clusters, (Mb, relb), keep_null=True)
        nulls[label] = r.pop("_null")
        res["features"][label] = r
        print(f"{label}: T={r['theta_hat']:.4g} p={r['main_block_x4']['p_one_sided']:.4g} "
              f"k={r['k_x4']} pR2/rel={r['partial_r2_share_reproducible']:.4%}", flush=True)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / f"{tag}_test{suffix}.json").write_text(json.dumps(res, indent=1, default=_json_default))
    np.savez_compressed(out_dir / f"{tag}_null{suffix}.npz", **nulls)
    return res


def _acf_M(P):
    """Diagnostic only (never used to choose k): within-arm ACF of M residualized on f(C)."""
    g = P["gen"]
    o = crt.genome_order(g.chrom.to_numpy(), g.pos.to_numpy())
    r = P["M"] - P["Q"] @ (P["Q"].T @ P["M"])
    a = crt.acf(r[o], max(ACF_LAGS), g.arm.to_numpy()[o])
    return dict(zip(map(str, ACF_LAGS), np.round(a[np.array(ACF_LAGS) - 1], 4).tolist()))


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(type(o))
