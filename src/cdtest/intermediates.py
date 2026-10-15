"""Deposited intermediates: the small processed inputs from which every confirmatory number in the
paper can be recomputed without the raw screens (~95 GB) or the local ATAC-seq pipeline.

The repository ships them under ``intermediates/`` with the same relative paths the analysis code
reads (``data/annot``, ``data/features``, ``results/aim1/control_moments``); the per-pair outcome
vectors M already live in ``results/aim1/<pair>/``. :func:`install` copies any file that is missing
from its working location and never overwrites an existing one, so a full local build always wins.

    python -c "from cdtest import intermediates; intermediates.install()"

Contents (see intermediates/MANIFEST.tsv for sizes and md5 sums):

* data/annot/genes_gencode_v50.tsv     gene annotation (chrom, TSS, arm, CpG o/e), genes used
* data/annot/gene_cn_log2.tsv          gene-level relative copy number per line, genes used
* data/features/gene_features*.parquet gene-level accessibility features (CPM primary, FC
                                       sensitivity) for the 4 screened and 7 decoy lines, genes used
* data/features/build_info*.json       provenance of the feature build
* results/aim1/control_moments/        per-line control-cell expression mean and variance (covariates)
"""
from __future__ import annotations

import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "intermediates"


def files(src: Path = SRC) -> list[Path]:
    """Relative paths of the deposited files (MANIFEST and README excluded)."""
    if not src.exists():
        return []
    return sorted(p.relative_to(src) for p in src.rglob("*")
                  if p.is_file() and p.name not in ("MANIFEST.tsv", "README.md"))


def available(src: Path = SRC) -> bool:
    return bool(files(src))


def is_deposited(rel: str | Path, root: Path = ROOT, src: Path = SRC) -> bool:
    """True if the working copy of ``rel`` is the deposited (gene-subset) file, not a full build."""
    w, d = Path(root) / rel, src / rel
    return w.exists() and d.exists() and w.stat().st_size == d.stat().st_size and w.read_bytes() == d.read_bytes()


def install(root: Path = ROOT, src: Path = SRC, verbose: bool = True) -> list[Path]:
    """Copy deposited files that are missing under ``root``; return the paths copied."""
    copied = []
    for rel in files(src):
        dst = Path(root) / rel
        if dst.exists():
            continue
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src / rel, dst)
        copied.append(rel)
    if verbose:
        print(f"intermediates: {len(copied)} file(s) installed, {len(files(src)) - len(copied)} already present",
              flush=True)
    return copied
