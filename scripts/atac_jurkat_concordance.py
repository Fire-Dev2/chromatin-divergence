#!/usr/bin/env python3
"""Concordance of the Jurkat CRISPRi-background sensitivity ATAC set (GSE302221) with the primary
Jurkat set (GSE273644).

1. Promoter insertion CPM over TSS +/-1 kb ([tss0-1000, tss0+1000), same windows as
   cdtest.features) for genes in data/annot/genes_gencode_v50.tsv on chr1-22,X; replicates pooled
   (counts and totals summed). Spearman/Pearson(log2(1+CPM)) pooled primary vs pooled CRISPRi, for
   all genes, all genes minus blacklisted TSS windows, and protein-coding. Context: the same
   correlation for each replicate pair within and across sets, and Jurkat vs K562 / GM12878
   (between-line reference).
2. IDR peak sharing: fraction of each set's IDR peaks overlapping >= 1 bp of the other set's IDR
   peaks (bedtools intersect -u), and bp Jaccard of the merged sets.
Writes results/atac_local/jurkat_GSE302221.concordance.json.
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from cdtest import features  # noqa: E402

INS = ROOT / "data" / "atac_ins"
OUT = ROOT / "results" / "atac_local" / "jurkat_GSE302221.concordance.json"
PRIM_PEAKS = ROOT / "data/atac_local/jurkat_GSE273644/jurkat_GSE273644.idr.narrowPeak.gz"
CRI_PEAKS = ROOT / "data/atac_local/jurkat_GSE302221/jurkat_GSE302221.idr.narrowPeak.gz"


def sh(cmd):
    return subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True,
                          executable="/bin/bash").stdout


def corr(a, b):
    la, lb = np.log2(1 + a), np.log2(1 + b)
    return dict(n=int(len(a)), spearman=float(spearmanr(a, b)[0]), pearson_log2=float(pearsonr(la, lb)[0]))


def main():
    g = pd.read_csv(ROOT / "data/annot/genes_gencode_v50.tsv", sep="\t")
    g = g[g.chrom.isin(features.CHROMS)].reset_index(drop=True)
    iv = (g.chrom.values, (g.tss0 - 1000).clip(lower=0).values, (g.tss0 + 1000).values)
    ins = {"Jurkat": sorted(INS.glob("Jurkat__rep*.npz")),
           "JurkatCRISPRi": sorted(INS.glob("JurkatCRISPRi__rep*.npz")),
           "K562": sorted(INS.glob("K562__rep*.npz")),
           "GM12878": sorted(INS.glob("GM12878__rep*.npz"))}
    assert len(ins["Jurkat"]) == 3 and len(ins["JurkatCRISPRi"]) == 3, ins
    m = features.cpm_matrices(ins, [iv], procs=6)
    cpm, rep = m["cpm"][0], m["rep_cpm"][0]
    bl = g.blacklist_overlap_tss1kb.astype(str).str.lower().eq("true").values
    pc = g.gene_type.eq("protein_coding").values
    subsets = {"all_genes": np.ones(len(g), bool), "all_genes_no_blacklist": ~bl,
               "protein_coding_no_blacklist": pc & ~bl}
    res = {"windows": "TSS +/-1 kb [tss0-1000, tss0+1000), genes_gencode_v50.tsv, chr1-22,X",
           "insertion_totals": {k: int(v) for k, v in m["totals"].items()},
           "rep_insertion_totals": {k: int(v) for k, v in m["rep_totals"].items()},
           "primary_vs_crispri": {}, "context": {}}
    for s, msk in subsets.items():
        res["primary_vs_crispri"][s] = corr(cpm.Jurkat.values[msk], cpm.JurkatCRISPRi.values[msk])
    msk = subsets["all_genes_no_blacklist"]
    for other in ("K562", "GM12878"):
        res["context"][f"Jurkat_vs_{other}"] = corr(cpm.Jurkat.values[msk], cpm[other].values[msk])
        res["context"][f"JurkatCRISPRi_vs_{other}"] = corr(cpm.JurkatCRISPRi.values[msk], cpm[other].values[msk])
    cols = [c for c in rep.columns if c.startswith("Jurkat")]
    pairs = {}
    for i, a in enumerate(cols):
        for b in cols[i + 1:]:
            pairs[f"{a}|{b}"] = float(spearmanr(rep[a].values[msk], rep[b].values[msk])[0])
    within_p = [v for k, v in pairs.items() if k.count("CRISPRi") == 0]
    within_c = [v for k, v in pairs.items() if k.count("CRISPRi") == 2]
    across = [v for k, v in pairs.items() if k.count("CRISPRi") == 1]
    res["context"]["replicate_spearman_all_genes_no_blacklist"] = pairs
    res["context"]["replicate_spearman_summary"] = dict(
        within_primary_mean=float(np.mean(within_p)), within_crispri_mean=float(np.mean(within_c)),
        across_sets_mean=float(np.mean(across)), across_sets_min=float(np.min(across)))

    # IDR peaks
    tmp = ROOT / "data/atac_local/jurkat_GSE302221/tmp"; tmp.mkdir(parents=True, exist_ok=True)
    for name, f in (("prim", PRIM_PEAKS), ("cri", CRI_PEAKS)):
        sh(f"zcat < {f} | cut -f1-3 | sort -k1,1 -k2,2n > {tmp}/{name}.bed")
        sh(f"bedtools merge -i {tmp}/{name}.bed > {tmp}/{name}.merged.bed")
    n_p = int(sh(f"wc -l < {tmp}/prim.bed")); n_c = int(sh(f"wc -l < {tmp}/cri.bed"))
    c_in_p = int(sh(f"bedtools intersect -u -a {tmp}/cri.bed -b {tmp}/prim.bed | wc -l"))
    p_in_c = int(sh(f"bedtools intersect -u -a {tmp}/prim.bed -b {tmp}/cri.bed | wc -l"))
    jac = sh(f"bedtools jaccard -a {tmp}/prim.merged.bed -b {tmp}/cri.merged.bed").splitlines()
    jd = dict(zip(jac[0].split(), jac[1].split()))
    res["idr_peaks"] = dict(primary_n=n_p, crispri_n=n_c,
                            frac_crispri_overlapping_primary=c_in_p / n_c,
                            frac_primary_overlapping_crispri=p_in_c / n_p,
                            bp_jaccard_merged=float(jd["jaccard"]),
                            rule=">= 1 bp overlap, bedtools intersect -u; IDR bfilt peaks of each pooled set")
    OUT.write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
