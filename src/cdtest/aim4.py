"""Aim 4 extras on the frozen pipeline (pipeline-v1.1, plan v1.3): new code only for analyses the
plan names in §7 / the proposal's Aim 4 that have no earlier implementation. Every statistic is
built from the frozen pieces (aim2.load_pair, crt engine, block rules, planted.decoy_theta0_tests);
no analysis parameter of the pipeline is changed.

* pooled_mixed_contrast  -- mean of the four mixed-pair T statistics under ONE shared block
  permutation (plan §7).
* similarity_gradient    -- H4: theta-hat of every available feature difference vs its similarity
  to the pair's own x (plan §5 sensitivity 2, §7).
* p53_control            -- the Aim 1 p53 gate statistic on every pair, oriented to the p53-WT line.
* target_promoter_check  -- exploratory knockdown-level check (sensitivity, not pre-specified).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from . import aim1, aim2, crt, features, screens

ROOT = screens.ROOT
PAIRS = (("K562", "Jurkat", "matched"), ("RPE1", "HepG2", "matched"),
         ("RPE1", "K562", "mixed"), ("K562", "HepG2", "mixed"),
         ("Jurkat", "RPE1", "mixed"), ("Jurkat", "HepG2", "mixed"))
#: p53 status (proposal: K562 and Jurkat p53-mutant/null; RPE1 and HepG2 p53 wild-type)
P53_WT = {"K562": False, "Jurkat": False, "RPE1": True, "HepG2": True}
PRIMARY = "combined_cpm"


def tag(a, b):
    return f"{a}_vs_{b}"


# --------------------------------------------------------------------------- pooled mixed contrast
def _pair_blocks(P: dict, label: str = PRIMARY):
    gen = P["gen"]
    chrom, pos, arm = gen.chrom.to_numpy(), gen.pos.to_numpy(), gen.arm.to_numpy()
    order = crt.genome_order(chrom, pos)
    prep = crt.prepare(P["M"], P["X"][label], Q=P["Q"])
    k = crt.choose_k(prep.r_x[order], arm[order], rule="x4")
    bid = crt.make_blocks(arm[order], int(k))
    return prep, order, int(k), bid


def pooled_mixed_contrast(Ps: list, label: str = PRIMARY, B: int = 9999, seed_key: str = "pooled_mixed",
                          chunk: int = 250) -> dict:
    """SECONDARY pooled statistic (each pair on its full gene set), union-order coupling.

    Limitation found on use: because each pair has its own genes and blocks, the coupling moves
    regions only approximately (a gene does not receive exactly the same source gene in every pair),
    and the pairs' null T* come out ~uncorrelated even for a 90% subset of the same pair (checked).
    It therefore behaves like independent permutations; the primary pooled contrast is
    :func:`pooled_shared_common` (exact shared permutation on the common genes).

    Construction:
    1. Each pair keeps its own main-test blocks (its own kept genes in genome order, x4 k from its
       own x residuals, blocks within arms, leftover joins the previous block, size classes).
    2. The union of the pairs' genes is put in genome order and cut into union blocks within arms
       with k_u = max over pairs of the pair k.  Each pair block is anchored to the union block that
       contains its first gene.
    3. One draw = one random ordering of the union blocks (shared by all pairs).  Within each
       equal-size class of each pair, slot i (in genome order) receives the pair block whose anchor
       has the i-th smallest shared key (ties between pair blocks with the same anchor broken at
       random).  So genomic regions are moved jointly: a region drawn early goes to the first slot
       of its class in every pair.
    Each pair's marginal null is therefore exactly a uniform within-class block permutation (the
    main test's null); the pairs' draws are coupled through the shared union ordering.  T* per draw
    is the mean of the pairs' T*; the CI inverts the pooled statistic assuming one common theta.
    """
    rng = aim2.job_seed(seed_key, *[P["tag"] for P in Ps])
    info = []
    gids = []
    for P in Ps:
        prep, order, k, bid = _pair_blocks(P, label)
        info.append(dict(P=P, prep=prep, order=order, k=k, bid=bid))
        gids.append(pd.Index(P["genes"]))
    union = gids[0]
    for g in gids[1:]:
        union = union.union(g)
    from . import covariates as cov
    ug = cov.genome_info(list(union))
    uorder = crt.genome_order(ug.chrom.to_numpy(), ug.pos.to_numpy())
    k_u = max(d["k"] for d in info)
    ubid_sorted = crt.make_blocks(ug.arm.to_numpy()[uorder], k_u)
    ublock_of_gene = pd.Series(np.empty(len(union), int), index=ug.index)
    ublock_of_gene.iloc[uorder] = ubid_sorted
    NU = int(ubid_sorted.max() + 1)
    for d in info:
        genes_sorted = np.array(d["P"]["genes"])[d["order"]]
        classes = crt._block_classes(d["bid"])
        cl = []
        for nb, s, pos in classes:
            if nb < 2:
                continue
            first = pos[:, 0]  # sorted-coordinate position of each block's first gene
            anchors = ublock_of_gene.loc[genes_sorted[first]].to_numpy()
            cl.append((nb, s, pos, anchors))
        d["classes"] = cl
        T0, c0 = crt._engine(d["prep"], d["prep"].r_x[None, :])
        d["T"], d["c_obs"] = float(T0[0]), float(c0[0])
    G_list = [d["prep"].G for d in info]
    Ts = np.zeros((len(info), B)); Cs = np.zeros((len(info), B))
    done = 0
    while done < B:
        n = min(chunk, B - done)
        keys = np.argsort(rng.random((n, NU)), axis=1).astype(float)  # shared ordering of union blocks
        for j, d in enumerate(info):
            G = G_list[j]
            idx = np.tile(np.arange(G), (n, 1))
            for nb, s, pos, anchors in d["classes"]:
                K = keys[:, anchors] + rng.random((n, nb)) * 0.5
                Pm = np.argsort(K, axis=1)
                idx[:, pos.ravel()] = pos[Pm].reshape(n, nb * s)
            out = np.empty_like(idx)
            out[:, d["order"]] = d["order"][idx]
            t, c = crt._engine(d["prep"], d["prep"].r_x[out])
            Ts[j, done:done + n] = t; Cs[j, done:done + n] = c
        done += n
    T = float(np.mean([d["T"] for d in info]))
    res = crt.TestResult("pooled_block", T, Ts.mean(0), Cs.mean(0), float(np.mean([d["c_obs"] for d in info])),
                         {"sd_x": float(np.mean([np.std(d["P"]["X"][label]) for d in info]))})
    grid = crt.theta_grid(T, res.info["sd_x"])
    lo, hi, ci_info = res.ci(grid=grid)
    flo, fhi, finfo = aim2.ci_refined(res, grid)
    per_pair = {}
    for j, d in enumerate(info):
        r = crt.TestResult("block_shared", d["T"], Ts[j], Cs[j], d["c_obs"])
        per_pair[d["P"]["tag"]] = {"T": d["T"], "k_x4": d["k"], "n_genes": d["prep"].G,
                                   "p_one_sided_shared_draws": r.p, "p_two_sided_shared_draws": r.p_two_sided,
                                   "null_sd_shared_draws": float(np.std(Ts[j]))}
    corr = np.corrcoef(Ts)
    return {"label": label, "pairs": [d["P"]["tag"] for d in info], "B": B, "T_pooled": T,
            "p_one_sided": res.p, "p_two_sided": res.p_two_sided, "ci95_plan_grid": [lo, hi],
            "ci95": [flo, fhi], "ci_info": ci_info, "ci_fine_info": finfo,
            "null_mean": float(res.Tstar.mean()), "null_sd": float(res.Tstar.std()),
            "k_union": k_u, "n_union_genes": len(union), "n_union_blocks": NU,
            "null_corr_between_pairs": np.round(corr, 3).tolist(), "per_pair": per_pair,
            "_Tstar_pairs": Ts, "_cstar_pairs": Cs}


def restrict_pair(P: dict, genes) -> dict:
    """The pair's data on a subset of its genes (in the given order); f(C) rows are taken from
    the pair's full design (same spline knots, same arm indicators) and re-orthonormalised."""
    pos = pd.Series(np.arange(len(P["genes"])), index=P["genes"]).loc[list(genes)].to_numpy()
    Z = P["Z"][pos]
    return dict(P, genes=list(genes), M=P["M"][pos], X={k: v[pos] for k, v in P["X"].items()},
                gen=P["gen"].iloc[pos], C=P["C"].iloc[pos], Z=Z, Q=crt.orthonormal_basis(Z))


def pooled_shared_common(Ps: list, label: str = PRIMARY, B: int = 9999, seed_key: str = "pooled_common",
                         chunk: int = 250) -> dict:
    """Mean of the pairs' T under ONE block permutation applied identically to every pair.

    Exact sharing needs the same genes in every pair, so the pooled statistic is computed on the
    genes kept in all of the given pairs (the intersection), in genome order. Each pair keeps its
    own M, x and f(C) (rows of its full design, arm effects on). k = the largest of the pairs' x4 k
    (each from its own x residuals on these genes); blocks are cut once within arms; one draw is one
    permutation of blocks within equal-size classes, applied to every pair's x residuals.
    """
    common = pd.Index(Ps[0]["genes"])
    for P in Ps[1:]:
        common = common.intersection(pd.Index(P["genes"]))
    from . import covariates as cov
    g = cov.genome_info(list(common))
    o = crt.genome_order(g.chrom.to_numpy(), g.pos.to_numpy())
    common = list(np.array(common)[o])
    subs = [restrict_pair(P, common) for P in Ps]
    arm = subs[0]["gen"].arm.to_numpy()
    preps, ks = [], []
    for S in subs:
        prep = crt.prepare(S["M"], S["X"][label], Q=S["Q"])
        preps.append(prep)
        ks.append(int(crt.choose_k(prep.r_x, arm, rule="x4")))
    k = max(ks)
    bid = crt.make_blocks(arm, k)
    rng = aim2.job_seed(seed_key, *[P["tag"] for P in Ps])
    n = len(Ps)
    Ts, Cs = np.zeros((n, B)), np.zeros((n, B))
    done = 0
    for idx in crt._block_index_gen(bid, B, rng, chunk):
        m = len(idx)
        for j, prep in enumerate(preps):
            t, c = crt._engine(prep, prep.r_x[idx])
            Ts[j, done:done + m] = t; Cs[j, done:done + m] = c
        done += m
    T_each = [p.theta_hat for p in preps]
    T = float(np.mean(T_each))
    res = crt.TestResult("pooled_block_shared", T, Ts.mean(0), Cs.mean(0), 1.0,
                         {"sd_x": float(np.mean([np.std(S["X"][label]) for S in subs]))})
    grid = crt.theta_grid(T, res.info["sd_x"])
    lo, hi, ci_info = res.ci(grid=grid)
    flo, fhi, finfo = aim2.ci_refined(res, grid)
    per = {}
    for j, P in enumerate(Ps):
        r = crt.TestResult("block_shared", T_each[j], Ts[j], Cs[j], 1.0)
        per[P["tag"]] = {"T_common_genes": T_each[j], "k_x4_common": ks[j], "n_genes_full": len(P["genes"]),
                         "p_one_sided": r.p, "p_two_sided": r.p_two_sided, "null_sd": float(Ts[j].std())}
    exch = aim2.block_exchangeability(subs[0]["gen"], k)
    return {"label": label, "pairs": [P["tag"] for P in Ps], "B": B, "n_common_genes": len(common),
            "k_shared": k, "exchangeability": exch, "T_pooled": T, "p_one_sided": res.p,
            "p_two_sided": res.p_two_sided, "ci95_plan_grid": [lo, hi], "ci95": [flo, fhi],
            "ci_info": ci_info, "ci_fine_info": finfo, "null_mean": float(res.Tstar.mean()),
            "null_sd": float(res.Tstar.std()), "null_corr_between_pairs": np.round(np.corrcoef(Ts), 3).tolist(),
            "per_pair": per, "_Tstar_pairs": Ts, "_cstar_pairs": Cs}


# --------------------------------------------------------------------------- H4
def candidate_features(P: dict, kind: str = "combined", signal: str = "cpm"):
    """Own x, the other five pairs' x (each in its own plan orientation) and the 21 decoys."""
    cols, labels, groups = [P["X"][PRIMARY]], [f"own:{P['tag']}"], ["own"]
    for a, b, _ in PAIRS:
        if tag(a, b) == P["tag"]:
            continue
        v = features.pair_feature(a, b, kind, signal=signal).reindex(P["genes"]).to_numpy(float)
        cols.append(v); labels.append(f"pair:{tag(a, b)}"); groups.append("other_pair")
    D, dl = aim2.decoy_matrix(P["genes"], kind, signal)
    for j, l in enumerate(dl):
        cols.append(D[:, j]); labels.append(f"decoy:{l}"); groups.append("decoy")
    return np.column_stack(cols), labels, groups


def similarity_gradient(P: dict) -> pd.DataFrame:
    X, labels, groups = candidate_features(P)
    sim = crt.decoy_similarity(P["X"][PRIMARY], X, P["Q"])
    rows = []
    rx_own = P["X"][PRIMARY] - P["Q"] @ (P["Q"].T @ P["X"][PRIMARY])
    for j, l in enumerate(labels):
        v = X[:, j]
        prep = crt.prepare(P["M"], v, Q=P["Q"])
        n = crt.naive_ols_test(P["M"], v, prep=prep)
        rows.append({"pair": P["tag"], "feature": l, "group": groups[j], "similarity": float(sim[j]),
                     "theta_hat": prep.theta_hat, "se_naive": float(n["se"]),
                     # theta-hat rescaled to the own feature's residual SD: equals theta_own * sim
                     # exactly if M = theta_own * x_own + noise orthogonal to v (diagnostic only)
                     "theta_std": prep.theta_hat * np.sqrt(prep.ssx / float(rx_own @ rx_own))})
    return pd.DataFrame(rows)


def ols_line(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    X = np.column_stack([np.ones_like(x), x])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ beta
    df = len(y) - 2
    s2 = res @ res / df
    cov_ = s2 * np.linalg.inv(X.T @ X)
    se = np.sqrt(np.diag(cov_))
    t = beta / se
    p = 2 * stats.t.sf(np.abs(t), df)
    r2 = 1 - (res @ res) / ((y - y.mean()) @ (y - y.mean()))
    q = stats.t.ppf(0.975, df)
    return {"n": len(y), "intercept": float(beta[0]), "intercept_se": float(se[0]), "intercept_p": float(p[0]),
            "intercept_ci95": [float(beta[0] - q * se[0]), float(beta[0] + q * se[0])],
            "slope": float(beta[1]), "slope_se": float(se[1]), "slope_p": float(p[1]),
            "slope_ci95": [float(beta[1] - q * se[1]), float(beta[1] + q * se[1])], "r2": float(r2)}


# --------------------------------------------------------------------------- p53 control
def p53_orientation(a: str, b: str) -> int:
    """+1 if A is the p53-WT line of a pair differing in p53, -1 if B is; 0 if they do not differ."""
    if P53_WT[a] == P53_WT[b]:
        return 0
    return 1 if P53_WT[a] else -1


def p53_control(a: str, b: str, n_sets: int = 1000, seed: int = screens.MASTER_SEED) -> dict:
    """The Aim 1 p53 gate statistic (aim1.p53_gate, identical algorithm and seed), with the null
    kept so it can be oriented: s = +1/-1 so that positive = more response in the p53-WT line;
    matched pairs (no p53 difference) keep the A - B orientation and are read two-sided."""
    t = tag(a, b)
    gs = pd.read_csv(aim1.OUT / t / "gene_outcome.tsv", sep="\t")
    annot = pd.read_csv(ROOT / "data/annot/genes_gencode_v50.tsv", sep="\t", low_memory=False)
    p53 = set(annot.gene_id[annot.p53_fischer2017_core.astype(str) == "True"])
    gs["p53"] = gs.gene_id.isin(p53)
    gs["binA"] = pd.qcut(np.log2(gs.mean_ctrlA), 10, labels=False, duplicates="drop")
    gs["binB"] = pd.qcut(np.log2(gs.mean_ctrlB), 10, labels=False, duplicates="drop")
    obs = gs.M[gs.p53].mean()
    rng = np.random.default_rng(seed)
    tgt = gs[gs.p53]
    pools = {k: v.index.values for k, v in gs[~gs.p53].groupby(["binA", "binB"])}
    allpool = gs.index[~gs.p53].values
    null = np.empty(n_sets)
    for i in range(n_sets):
        pick = [rng.choice(pools.get((row.binA, row.binB), allpool)) for _, row in tgt.iterrows()]
        null[i] = gs.M.values[pick].mean()
    s = p53_orientation(a, b)
    so = s if s != 0 else 1
    o_obs, o_null = so * obs, so * null
    return {"pair": t, "A": a, "B": b, "p53_differs": s != 0,
            "orientation": "positive = more response in p53-WT line" if s else "A - B (no p53 difference)",
            "sign_applied": so, "n_p53_targets": int(gs.p53.sum()), "mean_M_p53_raw": float(obs),
            "stat_oriented": float(o_obs), "null_mean_oriented": float(o_null.mean()),
            "null_sd": float(null.std()), "z_oriented": float((o_obs - o_null.mean()) / o_null.std()),
            "p_one_sided_oriented": float((1 + (o_null >= o_obs).sum()) / (1 + n_sets)),
            "p_two_sided": float((1 + (np.abs(null - null.mean()) >= abs(obs - null.mean())).sum()) / (1 + n_sets)),
            "n_sets": n_sets, "_null": null}


# --------------------------------------------------------------------------- knockdown-level check
def target_promoter_check(a: str, b: str, n_perm: int = 9999) -> dict:
    """Exploratory (not pre-specified): per kept knockdown, its response divergence
    d_p = mean_g[(b~_A^2 - s~_A^2) - (b~_B^2 - s~_B^2)] (on-target masked) against the promoter
    CPM accessibility difference of its own target gene (A - B), after regressing both on the
    target's expression difference (log2 control mean A - B). Spearman + permutation over knockdowns."""
    t = tag(a, b)
    z = np.load(aim1.OUT / t / "terms.npz", allow_pickle=True)
    mask = z["mask"]
    D = (z["btA"].astype(float) ** 2 - z["s2A"]) - (z["btB"].astype(float) ** 2 - z["s2B"])
    d = np.nanmean(np.where(mask, np.nan, D), axis=1)
    kd = pd.read_csv(aim1.OUT / t / "knockdowns.tsv", sep="\t")
    kept = kd[kd.kept].reset_index(drop=True)
    assert (kept.target.values == z["targets"]).all()
    gs = pd.read_csv(aim1.OUT / t / "gene_outcome.tsv", sep="\t").set_index("gene_id")
    xp = features.pair_feature(a, b, "promoter", signal="cpm").reindex(kept.target_id).to_numpy(float)
    ed = np.log2(gs.mean_ctrlA.reindex(kept.target_id).to_numpy()) - np.log2(gs.mean_ctrlB.reindex(kept.target_id).to_numpy())
    lfc_d = (kept.lfcA - kept.lfcB).to_numpy(float)
    ok = np.isfinite(xp) & np.isfinite(ed) & np.isfinite(d) & np.isfinite(lfc_d)
    Zc = np.column_stack([np.ones(ok.sum()), ed[ok], lfc_d[ok]])
    Qc = crt.orthonormal_basis(Zc)
    rd = d[ok] - Qc @ (Qc.T @ d[ok]); rx = xp[ok] - Qc @ (Qc.T @ xp[ok])
    rho = stats.spearmanr(rd, rx).statistic
    rng = aim2.job_seed("target_promoter", t)
    rr = stats.rankdata(rd); rxr = stats.rankdata(rx)
    rr = (rr - rr.mean()) / rr.std(); rxr = (rxr - rxr.mean()) / rxr.std()
    null = np.array([np.mean(rr * rxr[rng.permutation(rxr.size)]) for _ in range(n_perm)])
    return {"pair": t, "n_knockdowns": int(ok.sum()), "spearman_resid": float(rho),
            "p_one_sided_perm": float((1 + (null >= rho - 1e-12).sum()) / (1 + n_perm)),
            "p_two_sided_perm": float((1 + (np.abs(null) >= abs(rho) - 1e-12).sum()) / (1 + n_perm)),
            "spearman_raw": float(stats.spearmanr(d[ok], xp[ok]).statistic),
            "covariates_partialled": ["target expr_diff", "target on-target lfc difference (A - B)"],
            "mean_d": float(np.mean(d[ok])), "n_perm": n_perm}
