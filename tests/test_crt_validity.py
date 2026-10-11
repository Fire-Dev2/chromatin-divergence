"""Monte Carlo validity checks (exactness, coverage, naive inflation, choose_k)."""
import numpy as np
import pytest

from cdtest import crt

N_REPS = 2000
MC_TOL = 3.0  # accept within 3 Monte Carlo SE


def _block_exchangeable_x(rng, n_blocks, k, phi=0.8):
    """x built from iid stationary AR(1) segments of length k: exchangeable at block level."""
    e = rng.normal(size=(n_blocks, k))
    out = np.empty_like(e)
    out[:, 0] = e[:, 0]
    for j in range(1, k):
        out[:, j] = phi * out[:, j - 1] + np.sqrt(1 - phi ** 2) * e[:, j]
    return out.ravel()


def _shared_module_noise(rng, n_blocks, k, rho):
    v = np.repeat(rng.normal(size=n_blocks), k)
    return np.sqrt(rho) * v + np.sqrt(1 - rho) * rng.normal(size=n_blocks * k)


def _null_rep(seed, n_blocks=30, k=10, B=199, Z=None):
    rng = np.random.default_rng(seed)
    G = n_blocks * k
    x = _block_exchangeable_x(rng, n_blocks, k)
    # M strongly correlated along the same blocks, heteroscedastic, independent of x
    M = _shared_module_noise(rng, n_blocks, k, 0.6) * rng.uniform(0.5, 2, G)
    blk = crt.block_permutation_test(M, x, Z, k=k, B=B, rng=rng)
    ols = crt.naive_ols_test(M, x, Z)
    return blk.p, blk.p_two_sided, ols["p"]


def test_block_test_exact_under_block_exchangeable_null():
    ps = np.array([_null_rep(s) for s in range(N_REPS)])
    se = np.sqrt(0.05 * 0.95 / N_REPS)
    fpr_one = np.mean(ps[:, 0] <= 0.05)
    fpr_two = np.mean(ps[:, 1] <= 0.05)
    assert abs(fpr_one - 0.05) < MC_TOL * se, fpr_one
    assert abs(fpr_two - 0.05) < MC_TOL * se, fpr_two
    se1 = np.sqrt(0.01 * 0.99 / N_REPS)
    assert abs(np.mean(ps[:, 0] <= 0.01) - 0.01) < MC_TOL * se1
    # naive OLS inflates under the same shared block structure
    assert np.mean(ps[:, 2] <= 0.05) > 0.05 + 5 * se


def test_block_test_near_exact_with_covariates():
    rng = np.random.default_rng(123)
    G = 300
    C = rng.normal(size=(G, 2))
    Z = crt.build_design(C, [0])
    ps = np.array([_null_rep(10_000 + s, Z=Z)[0] for s in range(N_REPS)])
    se = np.sqrt(0.05 * 0.95 / N_REPS)
    assert abs(np.mean(ps <= 0.05) - 0.05) < MC_TOL * se


def test_ci_coverage_under_true_effect():
    theta, n = 0.25, 1000
    cover, grid_cover = [], []
    for s in range(n):
        rng = np.random.default_rng(50_000 + s)
        x = _block_exchangeable_x(rng, 30, 10)
        M = theta * x + _shared_module_noise(rng, 30, 10, 0.5)
        r = crt.block_permutation_test(M, x, k=10, B=199, rng=rng)
        cover.append(r.covers(theta))
        lo, hi, _ = r.ci()
        grid_cover.append(lo <= theta <= hi)
    se = np.sqrt(0.95 * 0.05 / n)
    assert abs(np.mean(cover) - 0.95) < MC_TOL * se, np.mean(cover)
    # grid endpoints round inward by at most one grid step
    assert np.mean(grid_cover) > 0.95 - MC_TOL * se - 0.02


def test_naive_ols_inflates_only_with_shared_modules():
    def fpr(rho_M, rho_x, n=600):
        rej = 0
        for s in range(n):
            rng = np.random.default_rng(90_000 + s)
            nb, k = 40, 10
            u = np.repeat(rng.normal(size=nb), k)
            x = np.sqrt(rho_x) * u + np.sqrt(1 - rho_x) * rng.normal(size=nb * k)
            M = _shared_module_noise(rng, nb, k, rho_M)
            rej += crt.naive_ols_test(M, x)["p"] <= 0.05
        return rej / n
    se = np.sqrt(0.05 * 0.95 / 600)
    assert abs(fpr(0.5, 0.0) - 0.05) < MC_TOL * se
    assert abs(fpr(0.0, 0.6) - 0.05) < MC_TOL * se
    assert fpr(0.5, 0.6) > 0.15


def test_choose_k_iid_hits_lower_bound():
    r = np.random.default_rng(0).normal(size=5000)
    assert crt.choose_k(r, rule="v1.0") == 5


def test_choose_k_ar1_matches_theory():
    rng = np.random.default_rng(1)
    phi, G = 0.85, 50_000
    e = rng.normal(size=G)
    r = np.empty(G)
    r[0] = e[0]
    for i in range(1, G):
        r[i] = phi * r[i - 1] + np.sqrt(1 - phi ** 2) * e[i]
    expected = int(np.ceil(np.log(0.05) / np.log(phi)))  # 19
    assert abs(crt.choose_k(r, rule="v1.0") - expected) <= 2


def test_choose_k_upper_bound_and_groups():
    G = 4000
    r = np.repeat(np.random.default_rng(2).normal(size=G // 100), 100)
    assert crt.choose_k(r, rule="v1.0") == 50
    # alternating constant groups: within-group pairs only
    smooth = np.sin(np.arange(G) / 3.0)
    k_all = crt.choose_k(smooth, rule="v1.0")
    assert 5 <= k_all <= 50
    assert crt.choose_k(np.random.default_rng(3).normal(size=G), groups=np.arange(G) // 10, rule="v1.0") == 5


# --- plan v1.1: k-rule options -------------------------------------------------

def _ar1(phi, G, seed):
    rng = np.random.default_rng(seed)
    e = rng.normal(size=G)
    r = np.empty(G)
    r[0] = e[0]
    for i in range(1, G):
        r[i] = phi * r[i - 1] + np.sqrt(1 - phi ** 2) * e[i]
    return r


def test_k_rules_bounds_and_relations():
    r = _ar1(0.85, 20_000, 4)
    k = {rule: crt.choose_k(r, rule=rule) for rule in crt.K_RULES}
    assert all(isinstance(v, int) and 5 <= v <= 100 for v in k.values())
    assert k["v1.0"] <= 50
    for m in (2, 3, 4):
        assert k[f"x{m}"] == min(m * k["v1.0"], 100)
    assert k["acf01"] >= k["v1.0"]
    assert k["tau2"] >= k["tau1"]
    # explicit v1.0 equals the frozen rule regardless of the default
    assert crt.choose_k(r, rule="v1.0") == crt._v10_k(r, None)
    with pytest.raises(ValueError):
        crt.choose_k(r, rule="nope")


def test_tau_matches_ar1_theory_and_pw_iid():
    phi = 0.8
    tau = crt.integrated_acf_time(_ar1(phi, 50_000, 5))
    L = int(np.ceil(np.log(0.05) / np.log(phi)))
    expected = 1 + 2 * sum(phi ** h for h in range(1, L + 1))
    assert abs(tau - expected) < 0.6
    assert crt.choose_k(np.random.default_rng(6).normal(size=5000), rule="pw") == 5
    # PW grows with dependence
    assert crt.politis_white_block(_ar1(0.9, 5000, 7)) > crt.politis_white_block(_ar1(0.5, 5000, 7))


def test_default_rule_is_plan_v11():
    assert crt.K_RULE_DEFAULT == "x4"
    r = _ar1(0.7, 20_000, 8)
    assert crt.choose_k(r) == min(4 * crt.choose_k(r, rule="v1.0"), 100)
