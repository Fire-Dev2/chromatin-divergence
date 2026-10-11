"""Build gene-level annotation tables (data/annot/*.tsv). See docs/annotations.md.

Usage:  python scripts/build_annotations.py [--skip-fasta]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from cdtest import annot as A  # noqa: E402


def main(skip_fasta: bool = False, fasta: str = str(A.FASTA_PATH)) -> None:
    # 1. genes ---------------------------------------------------------------
    g = A.parse_gencode_canonical()
    print(f"genes: {len(g)}; genes of these types lacking a canonical tx: {g.attrs['n_genes_without_canonical']}")
    cb = A.load_cytoband()
    g = pd.concat([g, A.annotate_cytoband(g.chrom, g.tss0, cb)], axis=1)
    cgi = A.load_cpg_islands()
    g["cgi_overlap_tss500"] = A.overlaps_intervals(g.chrom, g.tss0 - A.PROMOTER_HALF, g.tss0 + A.PROMOTER_HALF, cgi)
    bl = A.load_blacklist()
    g["blacklist_overlap_tss1kb"] = A.overlaps_intervals(g.chrom, g.tss0 - 1000, g.tss0 + 1000, bl)

    # 2. promoter CpG ----------------------------------------------------------
    if not skip_fasta:
        cpg = A.promoter_cpg(g.chrom.values, g.tss0.values, fasta=Path(fasta))
        g = pd.concat([g, cpg], axis=1)

    # 5. p53 targets -------------------------------------------------------------
    hgnc = A.load_hgnc()
    fis = A.load_fischer2017()
    mp = A.map_symbols(fis.gene_symbol, g, hgnc)
    mp = mp.merge(fis.rename(columns={"gene_symbol": "symbol"}), on="symbol")
    mp.to_csv(A.ANNOT / "fischer2017_p53_core_mapped.tsv", sep="\t", index=False)
    ok = mp.dropna(subset=["gene_id"])
    g["p53_fischer2017_core"] = g.gene_id.isin(ok.gene_id)
    g["p53_fischer2017_ndatasets"] = g.gene_id.map(dict(zip(ok.gene_id, ok.n_genomewide_datasets))).astype("Int64")
    print(f"Fischer 2017: {len(fis)} symbols, mapped {ok.gene_id.nunique()}; methods "
          f"{mp.method.value_counts().to_dict()}; unmapped {mp[mp.gene_id.isna()].symbol.tolist()}")

    # 6. knockdown groups ----------------------------------------------------------
    kg = A.build_knockdown_groups(g)
    kg.to_csv(A.KNOCKDOWN_GROUPS_TSV, sep="\t", index=False)
    g["knockdown_group"] = kg.group.values
    print("knockdown groups:", kg.group.value_counts().to_dict())
    print("ties:", kg[kg.tiebreak != ""][["gene_name", "group", "tiebreak"]].to_string(index=False))

    cols = ["gene_id", "gene_id_version", "gene_name", "gene_type", "hgnc_id", "chrom", "strand", "tss0", "tss1",
            "canonical_transcript_id", "canonical_transcript_type", "gene_start1", "gene_end1", "arm", "cytoband",
            "cpg_oe_tss500", "gc_frac_tss500", "n_cpg_tss500", "n_N_tss500", "cgi_overlap_tss500",
            "blacklist_overlap_tss1kb", "p53_fischer2017_core", "p53_fischer2017_ndatasets", "knockdown_group"]
    g = g[[c for c in cols if c in g.columns]]
    g.to_csv(A.GENES_TSV, sep="\t", index=False, float_format="%.5g")
    print("gene types:", g.gene_type.value_counts().to_dict())

    # 4. CN ------------------------------------------------------------------------
    seg, pm = A.build_cn_segments()
    seg.to_csv(A.CN_SEGMENTS_TSV, sep="\t", index=False, float_format="%.6g")
    pm.to_csv(A.CN_PROFILE_MAP_TSV, sep="\t", index=False)
    print(pm.to_string(index=False))
    cn = A.CNLookup(seg)
    out = pd.DataFrame({"gene_id": g.gene_id, "gene_name": g.gene_name})
    gap = pd.DataFrame({"gene_id": g.gene_id})
    for line in A.CN_LINES:
        out[line] = cn.point(line, g.chrom, g.tss0)
        unc = ~cn.covered(line, g.chrom, g.tss0)
        gap[line] = unc.astype(int)
        print(f"CN {line}: TSS outside any segment (nearest-segment filled): {unc.sum()}")
    gap.to_csv(A.ANNOT / "gene_cn_log2_gapfilled.tsv", sep="\t", index=False)
    absent = pm[pm.status.str.startswith("absent")].line.tolist()
    out.to_csv(A.GENE_CN_TSV, sep="\t", index=False, float_format="%.5f")
    print("CN lines set to 0 (absent from DepMap):", absent)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-fasta", action="store_true")
    ap.add_argument("--fasta", default=str(A.FASTA_PATH), help="uncompressed/bgzipped FASTA (chr-prefixed)")
    main(**vars(ap.parse_args()))
