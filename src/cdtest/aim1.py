"""Aim 1: the outcome M_g for a pair of lines, with its noise correction and baselines.

Implements ANALYSIS_PLAN §1–2 exactly:

  beta[p,g,k]   = log2(mean_p + eps) - log2(mean_ctrl + eps)                (eps = 0.1)
  beta~         = beta - mean over kept knockdowns (generic response removed)
  s~^2[p,g,k]   = mean over K fake screens of beta~_fake^2, where each fake screen
                  draws n_p control cells as a pseudo-knockdown for every real p and
                  passes through the same generic subtraction
  M_g           = mean_p [(beta~_A^2 - s~_A^2) - (beta~_B^2 - s~_B^2)]
  unsigned_g    = mean_p [(beta~_A - beta~_B)^2 - s~_A^2 - s~_B^2]

The entry for knockdown p and its own target gene is excluded from every average.
Reliability: M recomputed on random cell halves (knockdowns and controls split),
Spearman–Brown 2r/(1+r), median over splits.
"""
from __future__ import annotations

import json
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from . import screens

ROOT = screens.ROOT
PB = screens.OUT
OUT = ROOT / "results/aim1"
EPS = 0.1
K_FAKE = 50
MIN_CELLS = 20
REPRESSION = np.log2(0.7)
EXPR_MEAN = 0.2
EXPR_DETECT = 0.05


# --------------------------------------------------------------------------- loading
@dataclass
class Line:
    name: str
    groups: pd.DataFrame          # group, n_cells (row 0 = controls)
    genes: pd.DataFrame           # gene_id, gene_name
    sums: np.ndarray              # (n_groups, G) float64
    ctrl: np.ndarray              # (Nc, G) float32 memmap
    sums_a: np.ndarray            # (S, n_groups, G) float32
    n_a: np.ndarray               # (S, n_groups)
    halves_ctrl: np.ndarray       # (S, Nc) int8


def load_line(name: str) -> Line:
    d = PB / name
    g = pd.read_csv(d / "groups.tsv", sep="\t")
    genes = pd.read_csv(d / "genes.tsv", sep="\t", index_col=0)
    return Line(name, g, genes, np.load(d / "sums.npy"), np.load(d / "ctrl_norm.npy", mmap_mode="r"),
                np.load(d / "sums_halfA.npy", mmap_mode="r"), np.load(d / "n_halfA.npy"),
                np.load(d / "halves_ctrl.npy"))


# --------------------------------------------------------------------------- core maths
def log_ratio(kd_mean: np.ndarray, ctrl_mean: np.ndarray) -> np.ndarray:
    return np.log2(kd_mean + EPS) - np.log2(ctrl_mean + EPS)


def generic_subtract(beta: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Subtract each gene's mean over knockdowns, ignoring masked (on-target) entries."""
    b = np.where(mask, np.nan, beta)
    return b - np.nanmean(b, axis=0, keepdims=True)


def fake_screen_s2(ctrl: np.ndarray, n_kd: np.ndarray, mask: np.ndarray, rng: np.random.Generator,
                   k: int = K_FAKE, cols: np.ndarray | None = None) -> np.ndarray:
    """Noise variance s~^2 for every (p, g): K fake screens drawn from control cells."""
    nc = ctrl.shape[0]
    C = np.asarray(ctrl if cols is None else ctrl[:, cols], dtype=np.float32)
    total = C.sum(axis=0, dtype=np.float64)
    acc = np.zeros((len(n_kd), C.shape[1]), dtype=np.float64)
    for _ in range(k):
        rows, cols_i = [], []
        for p, n in enumerate(n_kd):
            pick = rng.choice(nc, size=int(n), replace=False)
            rows.append(np.full(len(pick), p)); cols_i.append(pick)
        rows = np.concatenate(rows); cols_i = np.concatenate(cols_i)
        ind = sp.csr_matrix((np.ones(len(rows), np.float32), (rows, cols_i)), shape=(len(n_kd), nc))
        sub = np.asarray(ind @ C, dtype=np.float64)
        kd_mean = sub / n_kd[:, None]
        pc_mean = (total[None, :] - sub) / (nc - n_kd[:, None])
        bt = generic_subtract(log_ratio(kd_mean, pc_mean), mask)
        acc += np.nan_to_num(bt) ** 2
    return acc / k


def line_terms(sums_kd, n_kd, ctrl_sum, n_ctrl, ctrl_cells, mask, rng, k=K_FAKE):
    """Return beta (raw), beta~ and s~^2 for one line on the kept knockdowns/genes."""
    beta = log_ratio(sums_kd / n_kd[:, None], (ctrl_sum / n_ctrl)[None, :])
    bt = generic_subtract(beta, mask)
    s2 = fake_screen_s2(ctrl_cells, n_kd, mask, rng, k=k)
    return beta, bt, s2


def m_from_terms(btA, s2A, btB, s2B, mask):
    dA = np.where(mask, np.nan, btA ** 2 - s2A)
    dB = np.where(mask, np.nan, btB ** 2 - s2B)
    M = np.nanmean(dA - dB, axis=0)
    U = np.nanmean(np.where(mask, np.nan, (btA - btB) ** 2 - s2A - s2B), axis=0)
    return M, U, np.nanmean(dA, axis=0), np.nanmean(dB, axis=0)


def generic_share(beta, mask):
    b = np.where(mask, np.nan, beta)
    mean_sq = np.nanmean(b ** 2, axis=0)
    return np.nanmean(b, axis=0) ** 2 / mean_sq


# --------------------------------------------------------------------------- pair pipeline
def select(A: Line, B: Line):
    """Genes and knockdowns kept for a pair (ANALYSIS_PLAN §1)."""
    def expr_ok(L):
        mean = L.sums[0] / L.groups.n_cells[0]
        det = np.asarray((L.ctrl > 0).mean(axis=0)).ravel()
        return pd.Series((mean >= EXPR_MEAN) & (det >= EXPR_DETECT), index=L.genes.gene_id.values)
    eA, eB = expr_ok(A), expr_ok(B)
    shared = [g for g in A.genes.gene_id if g in set(B.genes.gene_id)]
    genes = [g for g in shared if eA.get(g, False) and eB.get(g, False)]
    annot = pd.read_csv(ROOT / "data/annot/genes_gencode_v50.tsv", sep="\t", low_memory=False)
    annot = annot.set_index("gene_id")
    ok_type = set(annot.index[annot.gene_type.isin(["protein_coding", "lncRNA"])])
    genes = [g for g in genes if g in ok_type]

    ia = pd.Series(np.arange(len(A.genes)), index=A.genes.gene_id.values)
    ib = pd.Series(np.arange(len(B.genes)), index=B.genes.gene_id.values)
    name_to_id = {}
    for L in (A, B):
        name_to_id.update(dict(zip(L.genes.gene_name, L.genes.gene_id)))
    gA = A.groups.set_index("group"); gB = B.groups.set_index("group")
    kds = sorted(set(gA.index[1:]) & set(gB.index[1:]))
    rows = []
    for t in kds:
        nA, nB = int(gA.n_cells[t]), int(gB.n_cells[t])
        tid = name_to_id.get(t)
        r = dict(target=t, target_id=tid, nA=nA, nB=nB, cells_ok=(nA >= MIN_CELLS and nB >= MIN_CELLS),
                 target_expr_ok=tid in set(genes), lfcA=np.nan, lfcB=np.nan)
        if tid in ia.index and tid in ib.index:
            ra, rb = list(A.groups.group).index(t), list(B.groups.group).index(t)
            r["lfcA"] = log_ratio(A.sums[ra, ia[tid]] / nA, A.sums[0, ia[tid]] / A.groups.n_cells[0])
            r["lfcB"] = log_ratio(B.sums[rb, ib[tid]] / nB, B.sums[0, ib[tid]] / B.groups.n_cells[0])
        r["repressed"] = bool(r["lfcA"] <= REPRESSION and r["lfcB"] <= REPRESSION)
        r["kept"] = bool(r["cells_ok"] and r["target_expr_ok"] and r["repressed"])
        rows.append(r)
    kd = pd.DataFrame(rows)
    return genes, kd, ia, ib


def run_pair(a: str, b: str, tag: str | None = None, k: int = K_FAKE, n_splits: int | None = None,
             k_half: int = 10) -> dict:
    t0 = time.time()
    tag = tag or f"{a}_vs_{b}"
    out = OUT / tag
    out.mkdir(parents=True, exist_ok=True)
    A, B = load_line(a), load_line(b)
    genes, kd, ia, ib = select(A, B)
    kd.to_csv(out / "knockdowns.tsv", sep="\t", index=False)
    kept = kd[kd.kept].reset_index(drop=True)
    colsA, colsB = ia[genes].values, ib[genes].values
    rA = np.array([list(A.groups.group).index(t) for t in kept.target])
    rB = np.array([list(B.groups.group).index(t) for t in kept.target])
    gpos = pd.Series(np.arange(len(genes)), index=genes)
    mask = np.zeros((len(kept), len(genes)), bool)
    for p, tid in enumerate(kept.target_id):
        if tid in gpos.index:
            mask[p, gpos[tid]] = True
    print(f"[{tag}] {len(genes):,} genes, {len(kept):,}/{len(kd):,} knockdowns kept "
          f"({(~kd.cells_ok).sum()} few cells, {(~kd.target_expr_ok).sum()} target not expressed, "
          f"{(~kd.repressed).sum()} weak repression)", flush=True)

    seed = screens.MASTER_SEED + 7919 * (sum(map(ord, a)) + 3 * sum(map(ord, b)))
    rng = np.random.default_rng(seed)
    ctrlA = np.asarray(A.ctrl[:, colsA]); ctrlB = np.asarray(B.ctrl[:, colsB])
    betaA, btA, s2A = line_terms(A.sums[rA][:, colsA], kept.nA.values.astype(float), A.sums[0, colsA],
                                 A.groups.n_cells[0], ctrlA, mask, rng, k)
    betaB, btB, s2B = line_terms(B.sums[rB][:, colsB], kept.nB.values.astype(float), B.sums[0, colsB],
                                 B.groups.n_cells[0], ctrlB, mask, rng, k)
    M, U, rawA, rawB = m_from_terms(btA, s2A, btB, s2B, mask)
    print(f"[{tag}] full-data M done {time.time()-t0:.0f}s", flush=True)

    # reliability over random halves
    S = A.sums_a.shape[0] if n_splits is None else n_splits
    halves = []
    for s in range(S):
        ms = []
        for side in (1, 0):
            terms = []
            for L, rows, cols, ctrl in ((A, rA, colsA, ctrlA), (B, rB, colsB, ctrlB)):
                sa = np.asarray(L.sums_a[s][:, cols], dtype=np.float64)
                na = L.n_a[s].astype(float)
                if side == 1:
                    s_kd, n_kd, s_c, n_c = sa[rows], na[rows], sa[0], na[0]
                    cmask = L.halves_ctrl[s] == 1
                else:
                    tot = L.sums[:, cols]; nt = L.groups.n_cells.values.astype(float)
                    s_kd, n_kd = tot[rows] - sa[rows], nt[rows] - na[rows]
                    s_c, n_c = tot[0] - sa[0], nt[0] - na[0]
                    cmask = L.halves_ctrl[s] == 0
                n_kd = np.maximum(n_kd, 1)
                terms.append(line_terms(s_kd, n_kd, s_c, n_c, ctrl[cmask], mask, rng, k_half))
            ms.append(m_from_terms(terms[0][1], terms[0][2], terms[1][1], terms[1][2], mask)[0])
        r = np.corrcoef(ms[0], ms[1])[0, 1]
        halves.append(r)
        print(f"[{tag}] split {s+1}/{S}: r = {r:.3f}  ({time.time()-t0:.0f}s)", flush=True)
    r_med = float(np.median(halves))
    rel = 2 * r_med / (1 + r_med)

    gs = pd.DataFrame({"gene_id": genes, "M": M, "unsigned": U, "specificA": rawA, "specificB": rawB,
                       "generic_shareA": generic_share(betaA, mask), "generic_shareB": generic_share(betaB, mask),
                       "mean_ctrlA": A.sums[0, colsA] / A.groups.n_cells[0],
                       "mean_ctrlB": B.sums[0, colsB] / B.groups.n_cells[0]})
    gs.to_csv(out / "gene_outcome.tsv", sep="\t", index=False)
    np.savez_compressed(out / "terms.npz", btA=btA.astype(np.float32), btB=btB.astype(np.float32),
                        s2A=s2A.astype(np.float32), s2B=s2B.astype(np.float32), mask=mask,
                        targets=kept.target.values, genes=np.array(genes))
    summ = dict(pair=tag, A=a, B=b, n_genes=len(genes), n_kd_total=len(kd), n_kd_kept=len(kept),
                split_r=[round(x, 4) for x in halves], split_r_median=r_med, reliability=rel,
                generic_share_median_A=float(np.nanmedian(gs.generic_shareA)),
                generic_share_median_B=float(np.nanmedian(gs.generic_shareB)),
                M_mean=float(np.nanmean(M)), M_sd=float(np.nanstd(M)),
                K_fake=k, K_fake_halves=k_half, seed=seed, seconds=round(time.time() - t0, 1))
    (out / "summary.json").write_text(json.dumps(summ, indent=1))
    print(json.dumps(summ, indent=1), flush=True)
    return summ


def p53_gate(tag: str = "RPE1_vs_K562", n_sets: int = 1000, seed: int = screens.MASTER_SEED) -> dict:
    """ANALYSIS_PLAN §2 gate: p53 core targets should have M > 0 for A = RPE1, B = K562."""
    gs = pd.read_csv(OUT / tag / "gene_outcome.tsv", sep="\t")
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
        pick = []
        for _, row in tgt.iterrows():
            pool = pools.get((row.binA, row.binB), allpool)
            pick.append(rng.choice(pool))
        null[i] = gs.M.values[pick].mean()
    p = (1 + (null >= obs).sum()) / (1 + n_sets)
    res = dict(pair=tag, n_p53_targets=int(gs.p53.sum()), mean_M_p53=float(obs),
               null_mean=float(null.mean()), null_sd=float(null.std()), z=float((obs - null.mean()) / null.std()),
               p_one_sided=float(p), passed=bool(p < 0.05),
               top_p53_targets=gs[gs.p53].nlargest(10, "M").merge(
                   annot[["gene_id", "gene_name"]], on="gene_id")[["gene_name", "M"]].round(4).values.tolist())
    (OUT / tag / "p53_gate.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))
    return res


if __name__ == "__main__":
    import sys
    a, b = (sys.argv[1:3] if len(sys.argv) >= 3 else ("RPE1", "K562"))
    run_pair(a, b)
    if (a, b) == ("RPE1", "K562"):
        p53_gate()
