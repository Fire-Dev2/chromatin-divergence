"""Aim 2: randomization test of M ~ f(C) + x that randomizes the feature x.

Everything follows ANALYSIS_PLAN.md section 5.

Notation
--------
G genes. ``M`` (G,) response difference, held fixed. ``x`` (G,) accessibility
feature difference, the thing that is randomized. ``Z`` (G, p) covariate design
for f(C) (intercept + natural cubic splines + linear terms), see
:func:`build_design`.

Statistic
---------
T = OLS coefficient on x in M ~ Z + x. By Frisch-Waugh-Lovell,
T = <r_x, M> / <r_x, r_x>, where r_x = (I - H_Z) x. For any replacement feature v
the statistic is T*(v) = <v, r_M> / (|v|^2 - |Q'v|^2), with Q an orthonormal
basis of col(Z) and r_M = (I - H_Z) M. Both r_M and Q are computed once, so each
permutation costs a few dot products, and a whole chunk of permutations is one
matrix product.

Confidence intervals
--------------------
Testing H0: theta = theta0 means testing M - theta0 * x. The statistic is linear
in M, so T(theta0) = T - theta0 and T*_b(theta0) = T*_b - theta0 * c*_b with
c*_b = <v_b, r_x> / (|v_b|^2 - |Q'v_b|^2). Keeping (T*_b, c*_b) for the same B
replacement features gives the p-value at every theta in closed form.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterator, Sequence

import numpy as np
from scipy import stats

__all__ = [
    "natural_spline_basis", "build_design", "orthonormal_basis", "Prepared", "prepare",
    "fit_effect_model", "partial_r2", "TestResult", "pvalue", "choose_k",
    "genome_order", "make_blocks", "block_permutation_test", "block_permuted_features",
    "stratified_permutation_test", "make_strata", "expression_bin_permutation_test",
    "unrestricted_permutation_test", "decoy_test", "decoy_similarity",
    "naive_ols_test", "camera_vif_test", "theta_grid", "acf",
]

# --------------------------------------------------------------------------
# Design matrix
# --------------------------------------------------------------------------

def natural_spline_basis(x, df: int = 4, knots=None, boundary=None) -> np.ndarray:
    """Natural cubic spline basis without intercept (same span as R ``ns(x, df)``).

    df basis columns; interior knots at the (1..df-1)/df quantiles, boundary
    knots at min/max. The basis is linear beyond the boundary knots. Uses the
    truncated-power construction of Hastie, Tibshirani & Friedman (ESL eq. 5.4-5.5)
    on x rescaled to [0, 1] for conditioning.
    """
    x = np.asarray(x, dtype=float)
    if boundary is None:
        boundary = (np.min(x), np.max(x))
    lo, hi = boundary
    scale = hi - lo if hi > lo else 1.0
    u = (x - lo) / scale
    if knots is None:
        q = np.linspace(0, 1, df + 1)[1:-1]
        inner = (np.quantile(x, q) - lo) / scale
    else:
        inner = (np.asarray(knots, float) - lo) / scale
    xi = np.concatenate([[0.0], np.sort(inner), [1.0]])  # K = df + 1 knots
    K = len(xi)

    def d(k):
        return (np.clip(u - xi[k], 0, None) ** 3 - np.clip(u - xi[K - 1], 0, None) ** 3) / (
            xi[K - 1] - xi[k])

    cols = [u]
    dK2 = d(K - 2)
    for k in range(K - 2):
        cols.append(d(k) - dK2)
    B = np.column_stack(cols)
    return B[:, :df]


def build_design(C, spline_cols: Sequence = (), df: int = 4, intercept: bool = True) -> np.ndarray:
    """Design for f(C): natural cubic spline (``df``) on ``spline_cols``, linear for the rest.

    ``C`` is a pandas DataFrame (columns referenced by name) or a 2-D array
    (columns referenced by integer index). Collinear columns are allowed; the
    test uses an orthonormal basis of the column space.
    """
    if hasattr(C, "columns"):
        names = list(C.columns)
        arr = C.to_numpy(dtype=float)
        idx = [names.index(c) for c in spline_cols]
    else:
        arr = np.asarray(C, dtype=float)
        if arr.ndim == 1:
            arr = arr[:, None]
        idx = list(spline_cols)
    parts = [np.ones((arr.shape[0], 1))] if intercept else []
    for j in range(arr.shape[1]):
        if j in idx:
            parts.append(natural_spline_basis(arr[:, j], df=df))
        else:
            parts.append(arr[:, [j]])
    return np.column_stack(parts)


def orthonormal_basis(Z: np.ndarray, rtol: float = 1e-10) -> np.ndarray:
    """Orthonormal basis Q of col(Z) via SVD (handles rank deficiency)."""
    Z = np.asarray(Z, float)
    if Z.ndim == 1:
        Z = Z[:, None]
    # column scaling for conditioning; does not change the column space
    norms = np.linalg.norm(Z, axis=0)
    norms[norms == 0] = 1.0
    U, s, _ = np.linalg.svd(Z / norms, full_matrices=False)
    keep = s > rtol * s[0]
    return U[:, keep]


# --------------------------------------------------------------------------
# Prepared quantities and the vectorized engine
# --------------------------------------------------------------------------

@dataclass
class Prepared:
    M: np.ndarray
    x: np.ndarray
    Q: np.ndarray
    r_M: np.ndarray
    r_x: np.ndarray
    theta_hat: float
    ssx: float
    rank: int

    @property
    def G(self) -> int:
        return self.M.shape[0]

    @property
    def df_resid(self) -> int:
        return self.G - self.rank - 1

    @property
    def fitted_x(self) -> np.ndarray:
        return self.x - self.r_x


def prepare(M, x, Z=None, Q=None) -> Prepared:
    M = np.asarray(M, float)
    x = np.asarray(x, float)
    if Q is None:
        Q = orthonormal_basis(np.ones((len(M), 1)) if Z is None else Z)
    r_M = M - Q @ (Q.T @ M)
    r_x = x - Q @ (Q.T @ x)
    ssx = float(r_x @ r_x)
    if ssx <= 1e-12 * max(1.0, float(x @ x)):
        raise ValueError("x is (numerically) in the span of the covariates")
    return Prepared(M, x, Q, r_M, r_x, float(r_x @ r_M) / ssx, ssx, Q.shape[1])


def _engine(prep: Prepared, V: np.ndarray):
    """T* and c* for replacement features in the rows of V (n, G)."""
    V = np.atleast_2d(V)
    W = np.column_stack([prep.r_M, prep.r_x])
    sc = V @ W                       # (n, 2)
    proj = V @ prep.Q                # (n, p)
    den = np.einsum("ij,ij->i", V, V) - np.einsum("ij,ij->i", proj, proj)
    den = np.where(den > 1e-12, den, np.nan)
    return sc[:, 0] / den, sc[:, 1] / den


def _chunk_size(G: int, target: int = 4_000_000) -> int:
    return max(1, target // max(G, 1))


# --------------------------------------------------------------------------
# Test results, p-values and test inversion
# --------------------------------------------------------------------------

def pvalue(T: float, Tstar: np.ndarray, alternative: str = "greater") -> float:
    """Randomization p-value (1 + #{T* >= T}) / (1 + B) (ties counted, small tolerance).

    alternative: 'greater' (plan default, H1: theta > 0), 'less', 'two-sided'
    (|T*| >= |T|).
    """
    Tstar = np.asarray(Tstar, float)
    Tstar = Tstar[np.isfinite(Tstar)]
    B = Tstar.size
    tol = 1e-10 * max(1.0, abs(T))
    if alternative == "greater":
        n = np.sum(Tstar >= T - tol)
    elif alternative == "less":
        n = np.sum(Tstar <= T + tol)
    elif alternative == "two-sided":
        n = np.sum(np.abs(Tstar) >= abs(T) - tol)
    else:
        raise ValueError(alternative)
    return (1.0 + n) / (1.0 + B)


@dataclass
class TestResult:
    method: str
    T: float              # observed statistic (theta-hat)
    Tstar: np.ndarray     # (B,) statistic for each replacement feature
    cstar: np.ndarray     # (B,) slope of T* in theta (see module docstring)
    c_obs: float = 1.0
    info: dict = field(default_factory=dict)

    @property
    def B(self) -> int:
        return int(np.isfinite(self.Tstar).sum())

    @property
    def p(self) -> float:
        """One-sided p-value for H1: theta > 0 (the plan's primary p)."""
        return pvalue(self.T, self.Tstar, "greater")

    @property
    def p_two_sided(self) -> float:
        return pvalue(self.T, self.Tstar, "two-sided")

    def p_at(self, theta0, alternative: str = "greater") -> np.ndarray:
        """p-value(s) for H0: theta = theta0 (vectorized over theta0)."""
        th = np.atleast_1d(np.asarray(theta0, float))
        ok = np.isfinite(self.Tstar)
        Ts, cs = self.Tstar[ok], self.cstar[ok]
        Tobs = self.T - th * self.c_obs                    # (n,)
        Tst = Ts[None, :] - th[:, None] * cs[None, :]      # (n, B)
        tol = 1e-10 * np.maximum(1.0, np.abs(Tobs))[:, None]
        if alternative == "greater":
            n = (Tst >= Tobs[:, None] - tol).sum(1)
        elif alternative == "less":
            n = (Tst <= Tobs[:, None] + tol).sum(1)
        elif alternative == "two-sided":
            n = (np.abs(Tst) >= np.abs(Tobs)[:, None] - tol).sum(1)
        else:
            raise ValueError(alternative)
        p = (1.0 + n) / (1.0 + ok.sum())
        return p if np.ndim(theta0) else p[0]

    def covers(self, theta0, alpha: float = 0.05) -> bool:
        """True if theta0 is not rejected by the equal-tailed two-sided level-alpha test."""
        return bool(self.p_at(theta0, "greater") > alpha / 2 and self.p_at(theta0, "less") > alpha / 2)

    def ci(self, grid=None, alpha: float = 0.05, sd_x: float | None = None):
        """Test-inversion CI over a theta grid (default: plan's 201-point grid).

        A grid point is kept if the equal-tailed two-sided test (each tail at
        alpha/2) does not reject. Returns (lo, hi, info) where info flags an
        interval touching the grid edge or a non-contiguous accepted set.
        """
        if grid is None:
            grid = theta_grid(self.T, sd_x if sd_x is not None else self.info.get("sd_x", 1.0))
        grid = np.asarray(grid, float)
        keep = (self.p_at(grid, "greater") > alpha / 2) & (self.p_at(grid, "less") > alpha / 2)
        if not keep.any():
            return np.nan, np.nan, {"empty": True}
        ii = np.flatnonzero(keep)
        info = {"at_edge": bool(ii[0] == 0 or ii[-1] == len(grid) - 1),
                "contiguous": bool(ii[-1] - ii[0] + 1 == ii.size), "empty": False}
        return grid[ii[0]], grid[ii[-1]], info


def theta_grid(theta_naive: float, sd_x: float, n: int = 201, min_halfwidth: float = 0.5) -> np.ndarray:
    """Plan grid: n points over +/- max(4 |naive theta|, 0.5 M-units per SD of x)."""
    half = max(4 * abs(theta_naive), min_halfwidth / max(sd_x, 1e-12))
    return np.linspace(-half, half, n)


def _run(prep: Prepared, method: str, base: np.ndarray, index_gen: Iterator[np.ndarray],
         info: dict | None = None) -> TestResult:
    """Apply index arrays (chunk, G) to ``base`` and score through the engine."""
    T_obs, c_obs = _engine(prep, base[None, :])
    Ts, cs = [], []
    for idx in index_gen:
        t, c = _engine(prep, base[idx])
        Ts.append(t)
        cs.append(c)
    info = dict(info or {})
    info.setdefault("sd_x", float(np.std(prep.x)))
    return TestResult(method, float(T_obs[0]), np.concatenate(Ts), np.concatenate(cs),
                      float(c_obs[0]), info)


def _chunks(B: int, chunk: int):
    done = 0
    while done < B:
        n = min(chunk, B - done)
        yield n
        done += n


# --------------------------------------------------------------------------
# Genome order, block size, blocks
# --------------------------------------------------------------------------

def _chrom_key(c) -> tuple:
    s = str(c)
    s = s[3:] if s.lower().startswith("chr") else s
    if s.isdigit():
        return (0, int(s), "")
    return (1, {"X": 23, "Y": 24, "M": 25, "MT": 25}.get(s.upper(), 99), s)


def genome_order(chrom, pos) -> np.ndarray:
    """Indices that sort genes by (chromosome in natural order, position)."""
    chrom = np.asarray(chrom)
    pos = np.asarray(pos, float)
    uniq = sorted(set(chrom.tolist()), key=_chrom_key)
    rank = {c: i for i, c in enumerate(uniq)}
    cr = np.array([rank[c] for c in chrom.tolist()])
    return np.lexsort((pos, cr))


def acf(r: np.ndarray, max_lag: int = 50, groups=None) -> np.ndarray:
    """Autocorrelation of r at lags 1..max_lag along the given order.

    If ``groups`` (e.g. chromosome arm, same order as r) is given, only pairs
    within the same group contribute.
    """
    r = np.asarray(r, float)
    r = r - r.mean()
    v = r @ r / r.size
    out = np.empty(max_lag)
    g = None if groups is None else np.asarray(groups)
    for lag in range(1, max_lag + 1):
        if lag >= r.size:
            out[lag - 1:] = np.nan
            break
        a, b = r[:-lag], r[lag:]
        if g is not None:
            same = g[:-lag] == g[lag:]
            a, b = a[same], b[same]
        out[lag - 1] = (a @ b / a.size) / v if a.size else np.nan
    return out


def choose_k(x_resid_genome_order, groups=None, threshold: float = 0.05,
             kmin: int = 5, kmax: int = 50) -> int:
    """Smallest gene lag with autocorrelation of the x residuals < threshold, bounded to [kmin, kmax].

    Uses the feature only (plan section 5, step 4). If the autocorrelation never
    drops below the threshold up to kmax, returns kmax.
    """
    a = acf(x_resid_genome_order, kmax, groups)
    below = np.flatnonzero(a < threshold)
    lag = int(below[0] + 1) if below.size else kmax
    return int(min(max(lag, kmin), kmax))


def make_blocks(arm_sorted, k: int) -> np.ndarray:
    """Block ids for genes already in genome order.

    Each run of equal ``arm_sorted`` labels (a chromosome arm) is cut into
    consecutive blocks of k genes; a leftover shorter than k joins the previous
    block. An arm with fewer than k genes is a single block.
    """
    arm_sorted = np.asarray(arm_sorted)
    G = arm_sorted.size
    ids = np.empty(G, dtype=np.int64)
    starts = np.flatnonzero(np.r_[True, arm_sorted[1:] != arm_sorted[:-1]])
    ends = np.r_[starts[1:], G]
    nxt = 0
    for s, e in zip(starts, ends):
        n = e - s
        nb = max(1, n // k)
        local = np.minimum(np.arange(n) // k, nb - 1)
        ids[s:e] = nxt + local
        nxt += nb
    return ids


def _block_classes(block_ids: np.ndarray):
    """Group blocks by size: list of (n_blocks, size, pos_mat (n_blocks, size))."""
    order = np.argsort(block_ids, kind="stable")  # positions grouped by block, in order
    bid_sorted = block_ids[order]
    starts = np.flatnonzero(np.r_[True, bid_sorted[1:] != bid_sorted[:-1]])
    ends = np.r_[starts[1:], bid_sorted.size]
    sizes = ends - starts
    classes = []
    for s in np.unique(sizes):
        sel = np.flatnonzero(sizes == s)
        pos = np.stack([order[starts[i]:ends[i]] for i in sel])
        classes.append((len(sel), int(s), pos))
    return classes


def _block_index_gen(block_ids: np.ndarray, B: int, rng, chunk: int):
    """Index arrays that move whole blocks onto other blocks of the same size."""
    G = block_ids.size
    classes = _block_classes(block_ids)
    for n in _chunks(B, chunk):
        idx = np.tile(np.arange(G), (n, 1))
        for nb, s, pos in classes:
            if nb < 2:
                continue
            P = rng.permuted(np.tile(np.arange(nb), (n, 1)), axis=1)  # slot i <- block P[i]
            idx[:, pos.ravel()] = pos[P].reshape(n, nb * s)
        yield idx


def _strata_index_gen(strata: np.ndarray, B: int, rng, chunk: int):
    G = strata.size
    groups = [np.flatnonzero(strata == s) for s in np.unique(strata)]
    for n in _chunks(B, chunk):
        idx = np.tile(np.arange(G), (n, 1))
        for g in groups:
            if g.size < 2:
                continue
            P = rng.permuted(np.tile(np.arange(g.size), (n, 1)), axis=1)
            idx[:, g] = g[P]
        yield idx


def _rng(rng):
    return rng if isinstance(rng, np.random.Generator) else np.random.default_rng(rng)


# --------------------------------------------------------------------------
# Main test: block permutation
# --------------------------------------------------------------------------

def block_permutation_test(M, x, Z=None, chrom=None, pos=None, arm=None, k: int | None = None,
                           B: int = 9999, rng=None, blocks=None, chunk: int | None = None,
                           prep: Prepared | None = None, order: np.ndarray | None = None,
                           k_bounds=(5, 50)) -> TestResult:
    """Main null (plan section 5): block permutation of the x residuals.

    1. Residualize x on C (Z); order genes by (chrom, pos).
    2. Cut each arm (``arm``; defaults to chromosome) into blocks of k consecutive
       genes, leftover joins the previous block.
    3. Permute residual blocks among blocks of equal size, keeping within-block
       order, and add back the fitted part (adding back is a no-op for T because
       the fitted part lies in col(Z)).
    4. k defaults to :func:`choose_k` on the residuals in genome order.

    ``blocks`` (gene-level ids in the original gene order) overrides steps 1-2
    (e.g. co-accessibility clusters, or oracle blocks in simulation); genes in a
    block keep their genome order. ``order`` overrides the genome order (e.g. a
    random order for the wrong-block benchmark).
    """
    rng = _rng(rng)
    prep = prep if prep is not None else prepare(M, x, Z)
    G = prep.G
    chunk = chunk or _chunk_size(G)
    if order is None:
        if chrom is None:
            order = np.arange(G)
        else:
            order = genome_order(chrom, np.zeros(G) if pos is None else pos)
    order = np.asarray(order)
    if blocks is None:
        grp = None
        if arm is not None:
            grp = np.asarray(arm)[order]
        elif chrom is not None:
            grp = np.asarray(chrom)[order]
        r_sorted = prep.r_x[order]
        if k is None:
            k = choose_k(r_sorted, grp, kmin=k_bounds[0], kmax=k_bounds[1])
        bid_sorted = make_blocks(grp if grp is not None else np.zeros(G), int(k))
    else:
        b = np.asarray(blocks)[order]
        _, bid_sorted = np.unique(b, return_inverse=True)
    # work in sorted coordinates: permute r_x[order] by blocks, map back
    inv_gen = _block_index_gen(bid_sorted, B, rng, chunk)

    def gen():
        for idx_sorted in inv_gen:
            # replacement feature in sorted coordinates is r_sorted[idx_sorted];
            # in original coordinates gene order[j] receives r_x[order[idx_sorted[:, j]]]
            out = np.empty_like(idx_sorted)
            out[:, order] = order[idx_sorted]
            yield out

    n_blocks = int(bid_sorted.max() + 1)
    return _run(prep, "block", prep.r_x, gen(),
                {"k": None if k is None else int(k), "n_blocks": n_blocks})


def block_permuted_features(x, Z, block_ids_sorted, order, B, rng=None) -> np.ndarray:
    """Explicit replacement features x* = fitted + block-permuted residuals (B, G).

    For inspection and for passing to other statistics; the test itself uses the
    equivalent shortcut in :func:`block_permutation_test`.
    """
    rng = _rng(rng)
    Q = orthonormal_basis(Z)
    x = np.asarray(x, float)
    fitted = Q @ (Q.T @ x)
    r = x - fitted
    out = np.empty((B, len(x)))
    i = 0
    for idx_sorted in _block_index_gen(np.asarray(block_ids_sorted), B, rng, _chunk_size(len(x))):
        o = np.empty_like(idx_sorted)
        o[:, order] = order[idx_sorted]
        out[i:i + len(o)] = fitted + r[o]
        i += len(o)
    return out


# --------------------------------------------------------------------------
# Sensitivity nulls and comparators
# --------------------------------------------------------------------------

def make_strata(expr_diff, cn_diff=None, n_bins: int = 10, cn_thresh: float = 0.2) -> np.ndarray:
    """Strata: n_bins quantile bins of expression difference x 3 copy-number states."""
    e = np.asarray(expr_diff, float)
    edges = np.quantile(e, np.linspace(0, 1, n_bins + 1)[1:-1])
    eb = np.searchsorted(edges, e, side="right")
    if cn_diff is None:
        return eb
    cn = np.asarray(cn_diff, float)
    cs = np.where(cn < -cn_thresh, 0, np.where(cn > cn_thresh, 2, 1))
    return eb * 3 + cs


def stratified_permutation_test(M, x, Z=None, strata=None, B: int = 9999, rng=None,
                                chunk=None, prep=None, method="stratified") -> TestResult:
    """Sensitivity 1: shuffle x freely within strata (see :func:`make_strata`)."""
    rng = _rng(rng)
    prep = prep if prep is not None else prepare(M, x, Z)
    strata = np.zeros(prep.G, int) if strata is None else np.asarray(strata)
    return _run(prep, method, prep.x,
                _strata_index_gen(strata, B, rng, chunk or _chunk_size(prep.G)),
                {"n_strata": int(np.unique(strata).size)})


def expression_bin_permutation_test(M, x, Z=None, expr_diff=None, n_bins=10, B=9999, rng=None,
                                    chunk=None, prep=None) -> TestResult:
    """Comparator 3: permutation of x within expression-difference bins only."""
    return stratified_permutation_test(M, x, Z, make_strata(expr_diff, None, n_bins), B, rng,
                                       chunk, prep, method="expr_bin")


def unrestricted_permutation_test(M, x, Z=None, B=9999, rng=None, chunk=None, prep=None) -> TestResult:
    """Comparator 2: unrestricted permutation of x across all genes."""
    return stratified_permutation_test(M, x, Z, None, B, rng, chunk, prep, method="unrestricted")


def decoy_similarity(x, decoys, Z=None) -> np.ndarray:
    """Pearson correlation of each decoy (columns of ``decoys``) with x, both residualized on C."""
    Q = orthonormal_basis(np.ones((len(x), 1)) if Z is None else Z)
    D = np.asarray(decoys, float)
    rx = x - Q @ (Q.T @ x)
    rD = D - Q @ (Q.T @ D)
    rx = rx - rx.mean()
    rD = rD - rD.mean(0)
    return (rD.T @ rx) / (np.linalg.norm(rD, axis=0) * np.linalg.norm(rx))


def decoy_test(M, x, Z=None, decoys=None, prep=None, max_abs_similarity: float | None = None) -> TestResult:
    """Sensitivity 2: decoy features (G, n_decoys) as replacement features.

    If ``max_abs_similarity`` is given (plan negative control: 0.05), decoys with
    larger |similarity| to x are dropped.
    """
    prep = prep if prep is not None else prepare(M, x, Z)
    D = np.asarray(decoys, float)
    info = {}
    if max_abs_similarity is not None:
        sim = decoy_similarity(prep.x, D, prep.Q)
        D = D[:, np.abs(sim) < max_abs_similarity]
        info["n_kept"] = int(D.shape[1])
    if D.shape[1] == 0:
        raise ValueError("no decoys left after the similarity filter")
    T_obs, c_obs = _engine(prep, prep.x[None, :])
    t, c = _engine(prep, D.T)
    info = dict(info)
    info.setdefault("sd_x", float(np.std(prep.x)))
    return TestResult("decoy", float(T_obs[0]), t, c, float(c_obs[0]), info)


def naive_ols_test(M, x, Z=None, prep=None, alpha: float = 0.05) -> dict:
    """Comparator 1: OLS t-test treating genes as independent."""
    prep = prep if prep is not None else prepare(M, x, Z)
    rss = float(prep.r_M @ prep.r_M) - prep.theta_hat ** 2 * prep.ssx
    df = prep.df_resid
    se = np.sqrt(rss / df / prep.ssx)
    t = prep.theta_hat / se
    q = stats.t.ppf(1 - alpha / 2, df)
    return {"theta": prep.theta_hat, "se": se, "t": t, "df": df,
            "p": stats.t.sf(t, df), "p_two_sided": 2 * stats.t.sf(abs(t), df),
            "ci": (prep.theta_hat - q * se, prep.theta_hat + q * se)}


def camera_vif_test(M, x, Z=None, rho_bar: float = 0.0, n: int | None = None, prep=None,
                    alpha: float = 0.05) -> dict:
    """Comparator 4: CAMERA-style variance inflation of the OLS SE.

    VIF = 1 + (n - 1) * rho_bar, with rho_bar the mean inter-gene correlation
    (from control cells in real data) and n the number of genes (default G).
    """
    prep = prep if prep is not None else prepare(M, x, Z)
    res = naive_ols_test(M, x, Z, prep, alpha)
    n = prep.G if n is None else n
    vif = max(1.0 + (n - 1) * rho_bar, 1e-12)
    se = res["se"] * np.sqrt(vif)
    t = prep.theta_hat / se
    df = res["df"]
    q = stats.t.ppf(1 - alpha / 2, df)
    return {"theta": prep.theta_hat, "se": se, "t": t, "df": df, "vif": vif,
            "p": stats.t.sf(t, df), "p_two_sided": 2 * stats.t.sf(abs(t), df),
            "ci": (prep.theta_hat - q * se, prep.theta_hat + q * se)}


def partial_r2(M, x, Z=None, reliability: float = 1.0, prep=None) -> float:
    """Partial R^2 of x given C: (RSS_C - RSS_{C,x}) / RSS_C, divided by reliability."""
    prep = prep if prep is not None else prepare(M, x, Z)
    rss_c = float(prep.r_M @ prep.r_M)
    return (prep.theta_hat ** 2 * prep.ssx / rss_c) / reliability


def fit_effect_model(M, x, Z) -> dict:
    """Fit M ~ f(C) + x by OLS. Returns theta, naive SE/p, partial R^2 and residual pieces."""
    prep = prepare(M, x, Z)
    out = naive_ols_test(M, x, Z, prep)
    out["partial_r2"] = partial_r2(M, x, Z, prep=prep)
    out["prep"] = prep
    return out
