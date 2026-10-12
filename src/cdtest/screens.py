"""Aim 1, step 0: one streaming pass over a Perturb-seq screen.

Reads the dense raw-count matrix of a Replogle/Nadig h5ad in row chunks (it never
fits comfortably in memory twice) and writes everything later steps need:

* per-knockdown sums of normalized counts (P x G) and cell counts,
* the normalized control-cell matrix (Nc x G, float32) for the noise null,
* per-knockdown sums for half A of each of S random half-splits (S x P x G);
  half B is total minus half A,
* cell-level metadata (target, gem group, UMI count, mito %, kept flag, half bits).

Normalization (ANALYSIS_PLAN §1): each cell is scaled to the median total UMI count
of that line's kept control cells: norm = raw * median_ctrl_UMI / UMI_count.
No log transform before averaging.

Cell filter: mitochondrial fraction below the line's 99th percentile among control
cells. Every cell carries one target (dual-guide vectors target one gene).
"""
from __future__ import annotations

import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[2]
SCREEN_FILES = {
    "K562": ROOT / "data/perturbseq/K562_essential_raw_singlecell_01.h5ad",
    "RPE1": ROOT / "data/perturbseq/rpe1_raw_singlecell_01.h5ad",
    "Jurkat": ROOT / "data/perturbseq/GSE264667_jurkat_raw_singlecell_01.h5ad",
    "HepG2": ROOT / "data/perturbseq/GSE264667_hepg2_raw_singlecell_01.h5ad",
    # same-line baseline (plan §0): Replogle 2022 K562 genome-wide screen, same obs fields
    "K562gw": ROOT / "data/perturbseq/K562_gwps_raw_singlecell_01.h5ad",
}
#: screens whose knockdown set is restricted to another screen's targets (shared knockdowns
#: only), so the per-split sums stay small: K562gw has ~9,900 targets.
RESTRICT_TO = {"K562gw": "K562"}
OUT = ROOT / "results/aim1/pseudobulk"
CONTROL = "non-targeting"
MASTER_SEED = 20261011
N_SPLITS = 20
CHUNK = 4000


def _read_col(grp: h5py.Group, name: str) -> np.ndarray:
    """Read an obs/var column from either old ('__categories') or new anndata layout."""
    obj = grp[name]
    if isinstance(obj, h5py.Group):  # new-style categorical
        cats = obj["categories"][:].astype(str)
        return cats[obj["codes"][:]]
    vals = obj[:]
    if "__categories" in grp and name in grp["__categories"]:
        cats = grp["__categories"][name][:].astype(str)
        return cats[vals]
    if vals.dtype.kind in "SO":
        return vals.astype(str)
    return vals


def read_obs_var(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    with h5py.File(path, "r") as h:
        obs = pd.DataFrame({c: _read_col(h["obs"], c) for c in
                            ("gene", "gene_id", "gem_group", "UMI_count", "mitopercent", "sgID_AB")
                            if c in h["obs"]})
        var = pd.DataFrame({c: _read_col(h["var"], c) for c in ("gene_name", "gene_id", "chr", "start", "end")
                            if c in h["var"]})
        if "_index" in h["var"]:
            var.index = h["var"]["_index"][:].astype(str)
    return obs, var


def line_seed(line: str) -> int:
    return MASTER_SEED + sum(ord(c) for c in line)


def run(line: str, n_splits: int = N_SPLITS, chunk: int = CHUNK, targets=None) -> Path:
    """Streaming pass. ``targets``: optional set of target names; knockdown cells of other
    targets are not aggregated (controls always are). Defaults to RESTRICT_TO[line]'s groups."""
    path = SCREEN_FILES[line]
    if targets is None and line in RESTRICT_TO:
        ref = pd.read_csv(OUT / RESTRICT_TO[line] / "groups.tsv", sep="\t")
        targets = set(ref.group[ref.group != CONTROL])
    out = OUT / line
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    obs, var = read_obs_var(path)
    n_cells, n_genes = len(obs), len(var)

    is_ctrl = obs["gene"].values == CONTROL
    mito_cut = np.percentile(obs.loc[is_ctrl, "mitopercent"].astype(float), 99)
    keep = obs["mitopercent"].astype(float).values < mito_cut
    keep &= obs["UMI_count"].astype(float).values > 0
    med_umi = float(np.median(obs.loc[keep & is_ctrl, "UMI_count"].astype(float)))
    scale = (med_umi / obs["UMI_count"].astype(float).values).astype(np.float32)

    found = set(obs.loc[keep & ~is_ctrl, "gene"])
    n_targets_all = len(found)
    if targets is not None:
        found &= set(targets)
        keep &= is_ctrl | obs["gene"].isin(found).values
    targets = np.array(sorted(found))
    groups = np.concatenate([[CONTROL], targets])          # row 0 = controls
    gidx = pd.Series(np.arange(len(groups)), index=groups)
    cell_group = np.full(n_cells, -1)
    cell_group[keep] = gidx.loc[obs.loc[keep, "gene"].values].values
    n_per_group = np.bincount(cell_group[keep], minlength=len(groups))

    rng = np.random.default_rng(line_seed(line))
    halves = rng.integers(0, 2, size=(n_splits, n_cells), dtype=np.int8)  # 1 = half A

    ctrl_rows = np.flatnonzero(keep & is_ctrl)
    ctrl_pos = -np.ones(n_cells, dtype=np.int64)
    ctrl_pos[ctrl_rows] = np.arange(len(ctrl_rows))

    sums = np.zeros((len(groups), n_genes), dtype=np.float64)
    sums_a = np.zeros((n_splits, len(groups), n_genes), dtype=np.float32)
    ctrl_mat = np.lib.format.open_memmap(out / "ctrl_norm.npy", mode="w+", dtype=np.float32,
                                         shape=(len(ctrl_rows), n_genes))

    with h5py.File(path, "r") as h:
        X = h["X"]
        for start in range(0, n_cells, chunk):
            stop = min(start + chunk, n_cells)
            block = X[start:stop].astype(np.float32) * scale[start:stop, None]
            g = cell_group[start:stop]
            ok = g >= 0
            rows = np.flatnonzero(ok)
            ind = sp.csr_matrix((np.ones(len(rows), np.float32), (g[ok], rows)),
                                shape=(len(groups), stop - start))
            sums += ind @ block
            for s in range(n_splits):
                a = rows[halves[s, start:stop][rows] == 1]
                ind_a = sp.csr_matrix((np.ones(len(a), np.float32), (g[a], a)),
                                      shape=(len(groups), stop - start))
                sums_a[s] += ind_a @ block
            cp = ctrl_pos[start:stop]
            m = cp >= 0
            if m.any():
                ctrl_mat[cp[m]] = block[m]
            if (start // chunk) % 10 == 0:
                print(f"  {line}: {stop:,}/{n_cells:,} cells  {time.time()-t0:5.0f}s", flush=True)
    ctrl_mat.flush()

    np.save(out / "sums.npy", sums)
    np.save(out / "sums_halfA.npy", sums_a)
    np.save(out / "halves_ctrl.npy", halves[:, ctrl_rows])
    np.save(out / "n_halfA.npy", np.stack([np.bincount(cell_group[keep & (halves[s] == 1)],
                                                       minlength=len(groups)) for s in range(n_splits)]))
    pd.DataFrame({"group": groups, "n_cells": n_per_group}).to_csv(out / "groups.tsv", sep="\t", index=False)
    var.to_csv(out / "genes.tsv", sep="\t")
    meta = obs.assign(kept=keep, group_index=cell_group)
    meta.to_parquet(out / "cells.parquet")
    pd.Series({"line": line, "n_cells": n_cells, "n_kept": int(keep.sum()), "n_ctrl_kept": len(ctrl_rows),
               "n_targets": len(targets), "n_targets_in_screen": n_targets_all,
               "restricted_to": RESTRICT_TO.get(line), "n_genes": n_genes, "median_ctrl_umi": med_umi,
               "mito_cut_pct": float(mito_cut), "n_splits": n_splits, "seed": line_seed(line),
               "seconds": round(time.time() - t0, 1)}).to_json(out / "summary.json", indent=1)
    print(f"{line} done in {time.time()-t0:.0f}s: {keep.sum():,} kept cells, {len(ctrl_rows):,} controls, "
          f"{len(targets):,} knockdowns, {n_genes:,} genes", flush=True)
    return out


if __name__ == "__main__":
    import sys
    for ln in sys.argv[1:] or ["K562", "RPE1"]:
        run(ln)
