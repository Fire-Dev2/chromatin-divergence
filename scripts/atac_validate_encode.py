#!/usr/bin/env python
"""Validate the local ATAC pipeline against ENCODE's own processing of the same experiment.

Compares data/atac_local/<sample>/<sample>.fc.bigWig and .idr.narrowPeak.gz with ENCODE's
pooled FC bigWig and pooled-pseudoreplicate IDR peaks (K562 ENCSR868FGK: ENCFF102ARJ,
ENCFF948AFM):
  * Pearson and Spearman of log2(1 + mean FC) over (a) the union of both IDR peak sets
    (summit +/- 250 bp, merged) and (b) TSS +/- 1 kb of genes_gencode_v50.tsv (all genes
    and protein-coding only); windows overlapping blacklist v2 are dropped.
  * Peak overlap: Jaccard (bp, bedtools jaccard on merged peaks) and fraction of ENCODE
    IDR peaks overlapped by >=1 bp by a local IDR peak (and vice versa).
Pass criterion (pre-specified): Spearman >= 0.9 over TSS windows (all genes).

  python scripts/atac_validate_encode.py --sample k562_ENCSR868FGK \
      --encode-bw data/atac/K562__ENCFF102ARJ.bigWig --encode-peaks data/atac/K562__ENCFF948AFM.bed.gz
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from cdtest import atac  # noqa: E402

BL = ROOT / "data/atac_local/hg38-blacklist.v2.bed"
GENES = ROOT / "data/annot/genes_gencode_v50.tsv"


def run(cmd):
    return subprocess.run(["bash", "-o", "pipefail", "-c", cmd], check=True, capture_output=True, text=True,
                          env={"LC_ALL": "C", "PATH": "/opt/homebrew/bin:/usr/bin:/bin"}).stdout


def summit_windows(peaks: pd.DataFrame, half=250) -> pd.DataFrame:
    s = peaks["start"] + peaks["summit"]
    return pd.DataFrame({"chrom": peaks["chrom"], "start": (s - half).clip(lower=0), "end": s + half})


def corr(a, b):
    k = np.isfinite(a) & np.isfinite(b)
    x, y = np.log2(1 + a[k]), np.log2(1 + b[k])
    return {"n": int(k.sum()), "pearson": float(stats.pearsonr(x, y)[0]),
            "spearman": float(stats.spearmanr(x, y)[0])}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--encode-bw", required=True)
    ap.add_argument("--encode-peaks", required=True)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    sd = ROOT / "data/atac_local" / a.sample
    tmp = sd / "tmp_validate"; tmp.mkdir(exist_ok=True)
    loc_bw, loc_pk = sd / f"{a.sample}.fc.bigWig", sd / f"{a.sample}.idr.narrowPeak.gz"
    enc_pk = atac.read_narrowpeak(a.encode_peaks)
    our_pk = atac.read_narrowpeak(loc_pk)
    res = {"sample": a.sample, "encode_bw": a.encode_bw, "encode_peaks": a.encode_peaks,
           "n_peaks_encode_idr": len(enc_pk), "n_peaks_local_idr": len(our_pk)}

    # union of summit +/-250 windows, merged, blacklist removed
    w = pd.concat([summit_windows(enc_pk), summit_windows(our_pk)])
    w.to_csv(tmp / "w.bed", sep="\t", header=False, index=False)
    run(f"sort -k1,1 -k2,2n {tmp}/w.bed | bedtools merge -i - | bedtools intersect -v -a - -b {BL} > {tmp}/union.bed")
    u = pd.read_csv(tmp / "union.bed", sep="\t", header=None, names=["chrom", "start", "end"])
    print(f"union regions {len(u)}", flush=True)
    e = atac.bigwig_mean(a.encode_bw, u["chrom"], u["start"], u["end"])
    o = atac.bigwig_mean(loc_bw, u["chrom"], u["start"], u["end"])
    res["union_peaks_summit250"] = corr(e, o)
    print(res["union_peaks_summit250"], flush=True)

    g = pd.read_csv(GENES, sep="\t", usecols=["gene_id", "gene_type", "chrom", "tss0", "blacklist_overlap_tss1kb"])
    g = g[~g["blacklist_overlap_tss1kb"].astype(bool)]
    e = atac.bigwig_mean(a.encode_bw, g["chrom"], g["tss0"] - 1000, g["tss0"] + 1000)
    o = atac.bigwig_mean(loc_bw, g["chrom"], g["tss0"] - 1000, g["tss0"] + 1000)
    res["tss1kb_all_genes"] = corr(e, o)
    pc = (g["gene_type"] == "protein_coding").to_numpy()
    res["tss1kb_protein_coding"] = corr(e[pc], o[pc])
    print(res["tss1kb_all_genes"], res["tss1kb_protein_coding"], flush=True)

    # peak overlap
    enc_pk.iloc[:, :3].to_csv(tmp / "enc.bed", sep="\t", header=False, index=False)
    our_pk.iloc[:, :3].to_csv(tmp / "our.bed", sep="\t", header=False, index=False)
    run(f"sort -k1,1 -k2,2n {tmp}/enc.bed > {tmp}/enc.s.bed; sort -k1,1 -k2,2n {tmp}/our.bed > {tmp}/our.s.bed")
    j = run(f"bedtools jaccard -a <(bedtools merge -i {tmp}/enc.s.bed) -b <(bedtools merge -i {tmp}/our.s.bed)").splitlines()
    hdr, val = j[0].split("\t"), j[1].split("\t")
    res["jaccard_bp"] = float(dict(zip(hdr, val))["jaccard"])
    n_enc_rec = int(run(f"bedtools intersect -u -a {tmp}/enc.s.bed -b {tmp}/our.s.bed | wc -l").split()[0])
    n_our_in = int(run(f"bedtools intersect -u -a {tmp}/our.s.bed -b {tmp}/enc.s.bed | wc -l").split()[0])
    res["frac_encode_idr_recovered"] = n_enc_rec / len(enc_pk)
    res["frac_local_idr_in_encode"] = n_our_in / len(our_pk)
    # recovery of ENCODE peaks by local pooled (top300k) peaks too
    pooled = sd / f"{a.sample}.pooled.pval0.01.300K.bfilt.narrowPeak.gz"
    if pooled.exists():
        n2 = int(run(f"bedtools intersect -u -a {tmp}/enc.s.bed -b <(gzip -dc {pooled}) | wc -l").split()[0])
        res["frac_encode_idr_in_local_pooled300k"] = n2 / len(enc_pk)
    res["pass_spearman_tss_ge_0.9"] = bool(res["tss1kb_all_genes"]["spearman"] >= 0.9)
    for f in tmp.iterdir():
        f.unlink()
    tmp.rmdir()
    out = Path(a.out) if a.out else sd / "validation_vs_encode.json"
    out.write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
