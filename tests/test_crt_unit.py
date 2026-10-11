"""Unit tests for cdtest.crt building blocks."""
import numpy as np
import pytest
import statsmodels.api as sm

from cdtest import crt


def test_spline_basis_spans_linear_and_has_df_columns():
    rng = np.random.default_rng(0)
    x = rng.normal(size=500)
    B = crt.natural_spline_basis(x, df=4)
    assert B.shape == (500, 4)
    Q = crt.orthonormal_basis(np.column_stack([np.ones(500), B]))
    assert Q.shape[1] == 5
    lin = 3 * x + 1
    assert np.allclose(Q @ (Q.T @ lin), lin)
    # a cubic is not in the span (natural spline is linear in the tails)
    cub = x ** 3
    assert not np.allclose(Q @ (Q.T @ cub), cub)


def test_statistic_matches_statsmodels_and_partial_r2():
    rng = np.random.default_rng(1)
    G = 800
    C = rng.normal(size=(G, 3))
    Z = crt.build_design(C, [0])
    x = rng.normal(size=G) + C[:, 1]
    M = 0.3 * x + np.sin(C[:, 0]) + rng.normal(size=G)
    fit = sm.OLS(M, np.column_stack([Z, x])).fit()
    res = crt.naive_ols_test(M, x, Z)
    assert np.isclose(res["theta"], fit.params[-1])
    assert np.isclose(res["se"], fit.bse[-1])
    assert np.isclose(res["p_two_sided"], fit.pvalues[-1])
    red = sm.OLS(M, Z).fit()
    pr2 = (red.ssr - fit.ssr) / red.ssr
    assert np.isclose(crt.partial_r2(M, x, Z), pr2)
    assert np.isclose(crt.partial_r2(M, x, Z, reliability=0.5), 2 * pr2)


def test_engine_equals_refit_for_replacement_features():
    rng = np.random.default_rng(2)
    G = 300
    C = rng.normal(size=(G, 2))
    Z = crt.build_design(C, [0])
    x = rng.normal(size=G)
    M = rng.normal(size=G)
    r = crt.block_permutation_test(M, x, Z, k=7, B=5, rng=3)
    # recompute by explicit refit with x* = fitted + permuted residuals
    order = np.arange(G)
    ids = crt.make_blocks(np.zeros(G), 7)
    Xs = crt.block_permuted_features(x, Z, ids, order, 5, rng=3)
    for b in range(5):
        fit = sm.OLS(M, np.column_stack([Z, Xs[b]])).fit()
        assert np.isclose(fit.params[-1], r.Tstar[b])
    assert np.isclose(r.T, sm.OLS(M, np.column_stack([Z, x])).fit().params[-1])


def test_make_blocks_leftover_merges_into_previous_block():
    arm = np.array(["1p"] * 23 + ["1q"] * 4 + ["2p"] * 10)
    ids = crt.make_blocks(arm, 10)
    sizes = np.bincount(ids)
    assert list(sizes) == [10, 13, 4, 10]
    # blocks never span arms
    for b in np.unique(ids):
        assert len(set(arm[ids == b])) == 1


def test_block_permutation_keeps_blocks_intact_and_sizes_matched():
    rng = np.random.default_rng(4)
    arm = np.array(["a"] * 25 + ["b"] * 20)
    ids = crt.make_blocks(arm, 10)  # sizes 10, 15, 10, 10
    gen = crt._block_index_gen(ids, 50, rng, 50)
    idx = next(gen)
    for row in idx:
        for b in np.unique(ids):
            pos = np.flatnonzero(ids == b)
            src = row[pos]
            # consecutive within the source block, i.e. order kept
            assert np.all(np.diff(src) == 1)
            assert np.bincount(ids)[ids[src[0]]] == pos.size
            assert len(set(ids[src])) == 1
    # the size-15 block is alone in its class and never moves
    pos15 = np.flatnonzero(ids == 1)
    assert np.all(idx[:, pos15] == pos15)


def test_genome_order_natural_chromosome_order():
    chrom = np.array(["chr10", "chr2", "chrX", "chr2", "chr1"])
    pos = np.array([5, 9, 1, 3, 7])
    o = crt.genome_order(chrom, pos)
    assert list(chrom[o]) == ["chr1", "chr2", "chr2", "chr10", "chrX"]
    assert list(pos[o]) == [7, 3, 9, 5, 1]


def test_pvalue_formula():
    Ts = np.array([0.1, 0.5, 0.9, -1.0])
    assert crt.pvalue(0.5, Ts) == (1 + 2) / 5
    assert crt.pvalue(0.5, Ts, "two-sided") == (1 + 3) / 5
    assert crt.pvalue(10, Ts) == 1 / 5


def test_decoy_test_and_similarity():
    rng = np.random.default_rng(5)
    G = 500
    x = rng.normal(size=G)
    D = rng.normal(size=(G, 40))
    D[:, 0] = x + 0.01 * rng.normal(size=G)
    sim = crt.decoy_similarity(x, D)
    assert sim[0] > 0.99 and np.all(np.abs(sim[1:]) < 0.3)
    M = x + rng.normal(size=G)
    r = crt.decoy_test(M, x, None, D, max_abs_similarity=0.05)
    assert r.info["n_kept"] < 40
    assert r.p == 1 / (1 + r.info["n_kept"])


def test_camera_vif_reduces_to_ols_at_zero_rho():
    rng = np.random.default_rng(6)
    x = rng.normal(size=200)
    M = rng.normal(size=200)
    a = crt.naive_ols_test(M, x)
    b = crt.camera_vif_test(M, x, rho_bar=0.0)
    assert np.isclose(a["p"], b["p"])
    c = crt.camera_vif_test(M, x, rho_bar=0.01)
    assert np.isclose(c["vif"], 1 + 199 * 0.01)
    assert np.isclose(c["se"], a["se"] * np.sqrt(c["vif"]))


def test_ci_linear_shortcut_matches_direct_test():
    rng = np.random.default_rng(7)
    G = 400
    x = rng.normal(size=G)
    M = 0.2 * x + rng.normal(size=G)
    r = crt.block_permutation_test(M, x, k=5, B=199, rng=11)
    for th in (-0.1, 0.0, 0.2, 0.35):
        direct = crt.block_permutation_test(M - th * x, x, k=5, B=199, rng=11)
        assert np.isclose(direct.p, r.p_at(th))
    lo, hi, info = r.ci()
    assert lo < 0.2 < hi and info["contiguous"]


def test_vectorized_speed_10k_genes_9999_perms():
    import time
    rng = np.random.default_rng(8)
    G = 10_000
    chrom = np.repeat(np.arange(1, 23), G // 22 + 1)[:G]
    C = rng.normal(size=(G, 3))
    Z = crt.build_design(C, [0, 1])
    x = rng.normal(size=G)
    M = rng.normal(size=G)
    t = time.time()
    r = crt.block_permutation_test(M, x, Z, chrom=chrom, pos=np.arange(G), B=9999, rng=1)
    assert r.B == 9999
    assert time.time() - t < 15
