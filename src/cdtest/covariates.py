"""Covariate table C and design f(C) for a pair of lines (ANALYSIS_PLAN §4).

C has, per gene (the pair's kept genes, in the order of gene_outcome.tsv):

  log2_meanA, log2_meanB   log2 mean normalized control expression in A and in B
  log2_varA,  log2_varB    log2 variance (ddof = 1) of normalized expression over the
                           line's kept control cells (results/aim1/pseudobulk/<line>/ctrl_norm.npy)
  expr_diff                log2_meanA - log2_meanB
  cn_diff                  relative log2 copy number at the TSS, A - B
                           (data/annot/gene_cn_log2.tsv; RPE1 uses the karyotype column "RPE1")
  cpg_oe                   promoter CpG observed/expected over TSS +/- 500 bp (GENCODE v50 table)

f(C) = intercept + natural cubic spline (df 4) on log2_meanA and on log2_meanB + every
other column linear. expr_diff is an exact linear combination of the two spline bases'
linear parts, so the design is rank-deficient by construction; crt works with an
orthonormal basis of the column space, so this is harmless (docs/methods_aim2_aim3.md §2.1).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from . import crt, screens

ROOT = screens.ROOT
PB = screens.OUT
AIM1 = ROOT / "results/aim1"
ANNOT = ROOT / "data/annot/genes_gencode_v50.tsv"
CN_TSV = ROOT / "data/annot/gene_cn_log2.tsv"

#: column of gene_cn_log2.tsv used for each screen line. RPE1 -> karyotype-based track (plan §3).
CN_COLUMN = {"RPE1": "RPE1", "K562": "K562", "K562gw": "K562", "Jurkat": "Jurkat", "HepG2": "HepG2"}
COLUMNS = ["log2_meanA", "log2_meanB", "log2_varA", "log2_varB", "expr_diff", "cn_diff", "cpg_oe"]
SPLINE_COLS = ("log2_meanA", "log2_meanB")
SPLINE_DF = 4


def control_moments(line: str, gene_ids, pb_dir: Path = PB) -> pd.DataFrame:
    """Mean and variance (ddof 1) of normalized expression over the line's kept control cells."""
    d = Path(pb_dir) / line
    genes = pd.read_csv(d / "genes.tsv", sep="\t", index_col=0)
    pos = pd.Series(np.arange(len(genes)), index=genes.gene_id.values)
    cols = pos.loc[list(gene_ids)].to_numpy()
    X = np.asarray(np.load(d / "ctrl_norm.npy", mmap_mode="r")[:, cols], dtype=np.float64)
    return pd.DataFrame({"mean": X.mean(0), "var": X.var(0, ddof=1)}, index=list(gene_ids))


def build_covariates(a: str, b: str, gene_ids=None, tag: str | None = None, pb_dir: Path = PB,
                     annot: pd.DataFrame | None = None, cn: pd.DataFrame | None = None) -> pd.DataFrame:
    """Covariate table C (index gene_id) for pair A = a, B = b."""
    tag = tag or f"{a}_vs_{b}"
    if gene_ids is None:
        gene_ids = pd.read_csv(AIM1 / tag / "gene_outcome.tsv", sep="\t", usecols=["gene_id"]).gene_id
    gene_ids = list(gene_ids)
    mA, mB = control_moments(a, gene_ids, pb_dir), control_moments(b, gene_ids, pb_dir)
    if annot is None:
        annot = pd.read_csv(ANNOT, sep="\t", low_memory=False, usecols=["gene_id", "cpg_oe_tss500"])
    annot = annot.set_index("gene_id") if "gene_id" in annot.columns else annot
    if cn is None:
        cn = pd.read_csv(CN_TSV, sep="\t")
    cn = cn.set_index("gene_id") if "gene_id" in cn.columns else cn
    C = pd.DataFrame(index=pd.Index(gene_ids, name="gene_id"))
    C["log2_meanA"] = np.log2(mA["mean"].to_numpy())
    C["log2_meanB"] = np.log2(mB["mean"].to_numpy())
    C["log2_varA"] = np.log2(mA["var"].to_numpy())
    C["log2_varB"] = np.log2(mB["var"].to_numpy())
    C["expr_diff"] = C.log2_meanA - C.log2_meanB
    C["cn_diff"] = (cn.loc[gene_ids, CN_COLUMN[a]].to_numpy(float)
                    - cn.loc[gene_ids, CN_COLUMN[b]].to_numpy(float))
    C["cpg_oe"] = annot.loc[gene_ids, "cpg_oe_tss500"].to_numpy(float)
    bad = ~np.isfinite(C.to_numpy()).all(1)
    if bad.any():
        raise ValueError(f"{bad.sum()} genes have non-finite covariates: {list(C.index[bad][:5])}")
    return C[COLUMNS]


def design(C: pd.DataFrame) -> np.ndarray:
    """f(C): intercept + ns(df=4) on each line's log2 mean expression + other columns linear."""
    return crt.build_design(C[COLUMNS], SPLINE_COLS, df=SPLINE_DF)


def genome_info(gene_ids, annot: pd.DataFrame | None = None) -> pd.DataFrame:
    """chrom, TSS (1-based) and arm label (chrom + p/q) for ordering and blocks."""
    if annot is None:
        annot = pd.read_csv(ANNOT, sep="\t", low_memory=False, usecols=["gene_id", "chrom", "tss1", "arm"])
    a = (annot.set_index("gene_id") if "gene_id" in annot.columns else annot).loc[list(gene_ids)]
    return pd.DataFrame({"chrom": a.chrom.astype(str).values, "pos": a.tss1.astype(float).values,
                         "arm": (a.chrom.astype(str) + a.arm.astype(str)).values},
                        index=pd.Index(list(gene_ids), name="gene_id"))
