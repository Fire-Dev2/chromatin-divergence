"""Write the deposited intermediates (intermediates/ in the repository, and a tarball for GigaDB /
Zenodo) from a full local build.

Keeps only the genes that appear in at least one pair's outcome table (results/aim1/*/gene_outcome.tsv),
which is every gene any confirmatory analysis touches. Run after scripts/fetch_control_moments.py
(or after the streaming pass in cdtest.screens has produced results/aim1/pseudobulk/).

    python scripts/export_intermediates.py            # -> intermediates/ and dist/cdtest_intermediates_v<ver>.tar.gz
"""
from __future__ import annotations

import hashlib
import shutil
import sys
import tarfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import cdtest  # noqa: E402

OUT = ROOT / "intermediates"
LINES_MOMENTS = ("K562", "RPE1", "Jurkat", "HepG2", "K562gw")


def md5(p: Path) -> str:
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main():
    genes = set()
    for f in sorted((ROOT / "results/aim1").glob("*/gene_outcome.tsv")):
        genes |= set(pd.read_csv(f, sep="\t", usecols=["gene_id"]).gene_id)
    print(f"{len(genes):,} genes used by at least one pair", flush=True)
    if OUT.exists():
        shutil.rmtree(OUT)
    (OUT / "data/annot").mkdir(parents=True)
    (OUT / "data/features").mkdir(parents=True)
    (OUT / "results/aim1/control_moments").mkdir(parents=True)

    a = pd.read_csv(ROOT / "data/annot/genes_gencode_v50.tsv", sep="\t", low_memory=False)
    a[a.gene_id.isin(genes)].to_csv(OUT / "data/annot/genes_gencode_v50.tsv", sep="\t", index=False)
    cn = pd.read_csv(ROOT / "data/annot/gene_cn_log2.tsv", sep="\t")
    cn[cn.gene_id.isin(genes)].to_csv(OUT / "data/annot/gene_cn_log2.tsv", sep="\t", index=False)
    for name in ("gene_features.parquet", "gene_features_fc.parquet"):
        g = pd.read_parquet(ROOT / "data/features" / name)
        g[g.gene_id.isin(genes)].reset_index(drop=True).to_parquet(OUT / "data/features" / name, index=False)
    for name in ("build_info.json", "build_info_fc.json"):
        shutil.copy2(ROOT / "data/features" / name, OUT / "data/features" / name)
    for ln in LINES_MOMENTS:
        for ext in (".tsv.gz", ".json"):
            src = ROOT / "results/aim1/control_moments" / f"{ln}{ext}"
            if not src.exists():
                raise FileNotFoundError(f"{src} missing: run scripts/fetch_control_moments.py first")
            shutil.copy2(src, OUT / "results/aim1/control_moments" / src.name)

    readme = (OUT / "README.md")
    readme.write_text(
        "# cdtest deposited intermediates\n\n"
        f"Processed inputs for every confirmatory analysis in the cdtest paper (cdtest {cdtest.__version__}), "
        "restricted to the genes used by at least one pair. With these files and the per-pair outcome tables in "
        "`results/aim1/<pair>/`, `python scripts/reproduce_from_intermediates.py` recomputes Table 2 and Figures 3–4 "
        "without the raw screens or the ATAC-seq pipeline.\n\n"
        "| path | content |\n|---|---|\n"
        "| data/annot/genes_gencode_v50.tsv | GENCODE v50 gene annotation (chromosome, TSS, arm, promoter CpG o/e) |\n"
        "| data/annot/gene_cn_log2.tsv | relative log2 copy number at the TSS per line (DepMap 24Q4; RPE1 karyotype) |\n"
        "| data/features/gene_features.parquet | gene-level accessibility features, insertion CPM (primary), "
        "4 screened + 7 decoy lines |\n"
        "| data/features/gene_features_fc.parquet | the same from fold-change tracks (sensitivity) |\n"
        "| data/features/build_info*.json | provenance of the feature build (inputs, md5 sums, parameters) |\n"
        "| results/aim1/control_moments/<line>.tsv.gz | mean and variance (ddof 1) of normalized expression over "
        "kept control cells (covariates) |\n\n"
        "The outcome vectors M, knockdown lists and reliabilities are in `results/aim1/<pair>/` "
        "(gene_outcome.tsv, knockdowns.tsv, summary.json) in the repository and in the tarball. "
        "Install into a clone with `python -c \"from cdtest import intermediates; intermediates.install()\"` "
        "(copies only missing files).\n")
    rows = []
    for p in sorted(OUT.rglob("*")):
        if p.is_file() and p.name != "MANIFEST.tsv":
            rows.append((str(p.relative_to(OUT)), p.stat().st_size, md5(p)))
    pd.DataFrame(rows, columns=["path", "bytes", "md5"]).to_csv(OUT / "MANIFEST.tsv", sep="\t", index=False)

    dist = ROOT / "dist"
    dist.mkdir(exist_ok=True)
    tar = dist / f"cdtest_intermediates_v{cdtest.__version__}.tar.gz"
    with tarfile.open(tar, "w:gz") as t:
        t.add(OUT, arcname="intermediates")
        for d in sorted((ROOT / "results/aim1").glob("*_vs_*")):
            if d.is_dir():
                for name in ("gene_outcome.tsv", "knockdowns.tsv", "summary.json"):
                    if (d / name).exists():
                        t.add(d / name, arcname=f"results/aim1/{d.name}/{name}")
    total = sum(r[1] for r in rows)
    print(f"intermediates/: {len(rows)} files, {total / 1e6:.1f} MB; tarball {tar.relative_to(ROOT)} "
          f"{tar.stat().st_size / 1e6:.1f} MB (md5 {md5(tar)})", flush=True)
    (dist / f"{tar.name}.md5").write_text(f"{md5(tar)}  {tar.name}\n")


if __name__ == "__main__":
    main()
