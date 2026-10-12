#!/usr/bin/env python
"""Local ENCODE-equivalent ATAC-seq pipeline (provisional stand-in for ENCODE ATAC v1.9.1).

Usage (from the repo root, after `unset PYTHONPATH; source env/venv/bin/activate`):

  # one replicate (R1/R2 may be comma-separated lists of runs, processed per run and merged)
  python scripts/atac_pipeline.py rep  --sample jurkat_GSE273644 --rep rep1 \
      --r1 data/atac_fastq/.../X_1.fastq.gz --r2 data/atac_fastq/.../X_2.fastq.gz
  # pooled peaks / FC bigWig / IDR / overlap / FRiP / qc.json
  python scripts/atac_pipeline.py pool --sample jurkat_GSE273644 --reps rep1 rep2 rep3
  # both, from a sample sheet (TSV: sample rep r1 r2)
  python scripts/atac_pipeline.py sheet --sheet scripts/atac_samples.tsv --sample jurkat_GSE273644

Steps (see docs/atac_local_pipeline.md for the ENCODE correspondence):
  fastp (Nextera adapter, both reads) | bowtie2 -X 2000 --mm  -> raw BAM
  samtools view -F 1804 -f 2 -q 30 | fixmate -r -m | view -F 1804 -f 2   (= ENCODE filter)
  PBC/NRF from the filtered, pre-dedup, name-sorted BAM (ENCODE bedpe definition, chrM excluded)
  sort | markdup -r  -> nodup; drop chrM -> <rep>.nodup.no_chrM.bam
  bamtobed -bedpe -> drop fragments touching blacklist v2 -> tagAlign (2 reads/fragment)
  Tn5 shift (+4 on '+', -5 on '-' end; ENCODE awk) -> <rep>.tn5.tagAlign.gz
  pooled: MACS3 callpeak (-p 0.01 --shift -75 --extsize 150 --nomodel -B --SPMR --keep-dup all
  --call-summits), top 300k by p, bdgcmp -m FE -> bigWig, pooled pseudoreplicates (each rep
  split 50/50 by fragment, seed-derived) -> IDR 0.05 (--use-best-multisummit-IDR) and naive
  overlap; blacklist + chr[0-9XY]+ filter on final peaks; FRiP; TSS enrichment.
Resumable: each step writes a .done marker in the replicate/pool directory.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from cdtest import atac  # noqa: E402

IDX = ROOT / "data/genome/bowtie2_idx/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
BLACKLIST_GZ = ROOT / "data/annot/blacklist/hg38-blacklist.v2.bed.gz"
GENES = ROOT / "data/annot/genes_gencode_v50.tsv"
OUT_ROOT = ROOT / "data/atac_local"
IDR_BIN = ROOT / "env/idr-venv/bin/idr"
MACS3 = Path(sys.executable).parent / "macs3"
ADAPTER = "CTGTCTCTTATACACATCT"   # Nextera / Tn5
MASTER_SEED = 20261011            # ANALYSIS_PLAN.md §8
CAP_NUM_PEAK = 300_000
IDR_THRESH = 0.05
CHR_REGEX = r"^chr[0-9XY]+$"      # ENCODE regex_bfilt_peak_chr_name 'chr[\dXY]+'
MITO = "chrM"

GATE = {  # log/RESEARCH_LOG.md 'Jurkat ATAC selection'
    "tss_enrichment_min": 5.0, "nrf_min": 0.8, "pbc1_min": 0.7,
    "nondup_nonmito_fragments_merged_min": 25_000_000, "idr_peaks_min": 70_000,
}


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def sh(cmd: str, **kw):
    log("$ " + cmd if len(cmd) < 600 else "$ " + cmd[:600] + " ...")
    env = dict(os.environ, LC_ALL="C", LC_COLLATE="C")
    r = subprocess.run(["bash", "-o", "pipefail", "-c", cmd], env=env, **kw)
    if r.returncode != 0:
        raise RuntimeError(f"command failed ({r.returncode}): {cmd}")
    return r


def sh_out(cmd: str) -> str:
    env = dict(os.environ, LC_ALL="C", LC_COLLATE="C")
    return subprocess.run(["bash", "-o", "pipefail", "-c", cmd], env=env, check=True,
                          capture_output=True, text=True).stdout


def step(d: Path, name: str, fn):
    m = d / f".{name}.done"
    if m.exists():
        log(f"skip {name} (done)")
        return
    t = time.time()
    fn()
    m.write_text(f"{time.time() - t:.0f}s\n")
    log(f"done {name} in {time.time() - t:.0f}s")


def derive_seed(*keys) -> int:
    h = hashlib.sha256(("|".join(map(str, (MASTER_SEED,) + keys))).encode()).hexdigest()
    return int(h[:8], 16)


def common_files(tmp: Path):
    """Uncompressed blacklist + chrom sizes (from the bowtie2 index) in data/atac_local."""
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    bl = OUT_ROOT / "hg38-blacklist.v2.bed"
    if not bl.exists():
        sh(f"gzip -dc {BLACKLIST_GZ} | sort -k1,1 -k2,2n > {bl}")
    cs = OUT_ROOT / "GRCh38_no_alt.chrom.sizes"
    if not cs.exists():
        sh(f"bowtie2-inspect -s {IDX} 2>/dev/null | awk -F'\\t' '$1 ~ /^Sequence-/ {{split($2,a,\" \"); print a[1]\"\\t\"$3}}' > {cs}")
        if os.path.getsize(cs) == 0:
            raise RuntimeError("could not derive chrom sizes")
    return bl, cs


def tss_table(protein_coding_only=True) -> pd.DataFrame:
    g = pd.read_csv(GENES, sep="\t", usecols=["gene_id", "gene_type", "chrom", "strand", "tss0"])
    if protein_coding_only:
        g = g[g["gene_type"] == "protein_coding"]
    return g.drop_duplicates(["chrom", "tss0", "strand"])


# --------------------------------------------------------------------------- replicate
def run_rep(sample: str, rep: str, r1s: list[str], r2s: list[str], threads: int = 10):
    d = OUT_ROOT / sample / rep
    tmp = OUT_ROOT / sample / "tmp" / rep
    d.mkdir(parents=True, exist_ok=True); tmp.mkdir(parents=True, exist_ok=True)
    bl, cs = common_files(tmp)
    p = f"{d}/{rep}"
    st = max(2, threads // 2)
    assert len(r1s) == len(r2s) and r1s

    def align():
        bams = []
        for i, (r1, r2) in enumerate(zip(r1s, r2s), 1):
            for f in (r1, r2):
                if not Path(f).exists():
                    raise FileNotFoundError(f)
            raw = f"{tmp}/{rep}.run{i}.raw.bam"
            sh(f"fastp -i {r1} -I {r2} --stdout -w 4 "
               f"--adapter_sequence {ADAPTER} --adapter_sequence_r2 {ADAPTER} "
               f"--disable_quality_filtering --disable_trim_poly_g --length_required 5 "
               f"-j {p}.run{i}.fastp.json -h {tmp}/{rep}.run{i}.fastp.html 2> {p}.run{i}.fastp.log "
               f"| bowtie2 -X 2000 --mm -p {threads} -x {IDX} --interleaved - 2> {p}.run{i}.bowtie2.log "
               f"| samtools view -@ 2 -b -o {raw} -")
            bams.append(raw)
        if len(bams) == 1:
            os.replace(bams[0], f"{tmp}/{rep}.raw.bam")
        else:
            sh(f"samtools cat -o {tmp}/{rep}.raw.bam {' '.join(bams)} && rm -f {' '.join(bams)}")
    step(d, "align", align)

    def raw_qc():
        sh(f"samtools flagstat -@ {st} {tmp}/{rep}.raw.bam > {p}.raw.flagstat.txt")
        out = sh_out(f"samtools view -@ {st} -F 0x904 {tmp}/{rep}.raw.bam | cut -f3 | "
                     f"awk '{{n++; if($1==\"{MITO}\" || $1==\"MT\") m++}} END {{print n+0\"\\t\"m+0}}'")
        n, m = map(int, out.split())
        (d / f"{rep}.raw_mito.json").write_text(json.dumps({"mapped_primary_reads": n, "mito_reads": m}))
    step(d, "raw_qc", raw_qc)

    def filt():
        # ENCODE rm_unmapped_lowq_reads_pe (bowtie2 output is mate-collated, so no name sort needed)
        sh(f"samtools view -@ {st} -F 1804 -f 2 -q 30 -u {tmp}/{rep}.raw.bam "
           f"| samtools fixmate -@ {st} -r -m - - "
           f"| samtools view -@ {st} -F 1804 -f 2 -b -o {tmp}/{rep}.nmsrt.filt.bam -")
        r1, r2 = sh_out(f"samtools view -c -f 64 {tmp}/{rep}.nmsrt.filt.bam"), sh_out(f"samtools view -c -f 128 {tmp}/{rep}.nmsrt.filt.bam")
        if int(r1) != int(r2):
            raise RuntimeError(f"orphans after filtering: R1={r1} R2={r2}")
    step(d, "filter", filt)

    def pbc():
        # ENCODE pbc_qc_pe: bedpe of name-sorted filtered (pre-dedup) BAM, chrM excluded
        sh(f"bedtools bamtobed -bedpe -i {tmp}/{rep}.nmsrt.filt.bam 2> {tmp}/bamtobed.err "
           f"| awk 'BEGIN{{OFS=\"\\t\"}}{{print $1,$2,$4,$6,$9,$10}}' | grep -v '^{MITO}\\s' "
           f"| sort -S 6G --parallel={st} -T {tmp} | uniq -c "
           f"| awk 'BEGIN{{mt=0;m0=0;m1=0;m2=0}} ($1==1){{m1=m1+1}} ($1==2){{m2=m2+1}} {{m0=m0+1}} {{mt=mt+$1}} "
           f"END{{m1_m2=-1.0; if(m2>0) m1_m2=m1/m2; m0_mt=0; if (mt>0) m0_mt=m0/mt; m1_m0=0; if (m0>0) m1_m0=m1/m0; "
           f"printf \"%d\\t%d\\t%d\\t%d\\t%f\\t%f\\t%f\\n\",mt,m0,m1,m2,m0_mt,m1_m0,m1_m2}}' > {p}.lib_complexity.qc")
    step(d, "pbc", pbc)

    def dedup():
        sh(f"samtools sort -@ {st} -m 1500M -T {tmp}/{rep}.srt -o {tmp}/{rep}.filt.bam {tmp}/{rep}.nmsrt.filt.bam")
        sh(f"samtools flagstat -@ {st} {tmp}/{rep}.filt.bam > {p}.filt.flagstat.txt")
        sh(f"samtools markdup -@ {st} -r -s -f {p}.markdup.txt -T {tmp}/{rep}.md {tmp}/{rep}.filt.bam {tmp}/{rep}.nodup.bam")
        sh(f"samtools view -@ {st} -b -e 'rname!=\"{MITO}\" && rname!=\"MT\"' -o {p}.nodup.no_chrM.bam {tmp}/{rep}.nodup.bam")
        sh(f"samtools index -@ {st} {p}.nodup.no_chrM.bam")
        sh(f"samtools flagstat -@ {st} {tmp}/{rep}.nodup.bam > {p}.nodup.flagstat.txt")
        sh(f"samtools flagstat -@ {st} {p}.nodup.no_chrM.bam > {p}.nodup.no_chrM.flagstat.txt")
        for f in ("nmsrt.filt.bam", "filt.bam", "nodup.bam"):
            Path(f"{tmp}/{rep}.{f}").unlink(missing_ok=True)
    step(d, "dedup", dedup)

    def tagalign():
        sh(f"samtools sort -n -@ {st} -m 1500M -T {tmp}/{rep}.nsrt -o {tmp}/{rep}.nmsrt.final.bam {p}.nodup.no_chrM.bam")
        sh(f"bedtools bamtobed -bedpe -mate1 -i {tmp}/{rep}.nmsrt.final.bam 2> {tmp}/bamtobed2.err "
           f"| tee >(wc -l > {tmp}/n_frag_pre_bl.txt) "
           f"| bedtools pairtobed -type neither -a stdin -b {bl} "
           f"| awk 'BEGIN{{OFS=\"\\t\"}}{{printf \"%s\\t%s\\t%s\\tN\\t1000\\t%s\\n%s\\t%s\\t%s\\tN\\t1000\\t%s\\n\",$1,$2,$3,$9,$4,$5,$6,$10}}' "
           f"| awk 'BEGIN {{OFS = \"\\t\"}}{{ if ($6 == \"+\") {{$2 = $2 + 4}} else if ($6 == \"-\") {{$3 = $3 - 5}} print $0}}' "
           f"| bgzip -@ {st} -l 6 > {p}.tn5.tagAlign.gz")
        time.sleep(2)
        n_pre = int(Path(f"{tmp}/n_frag_pre_bl.txt").read_text().split()[0])
        n_ta = atac.count_lines(f"{p}.tn5.tagAlign.gz")
        (d / f"{rep}.tagalign_counts.json").write_text(json.dumps(
            {"fragments_nodup_nonmito": n_pre, "fragments_after_blacklist": n_ta // 2, "tagalign_reads": n_ta}))
        Path(f"{tmp}/{rep}.nmsrt.final.bam").unlink(missing_ok=True)
    step(d, "tagalign", tagalign)

    def tss():
        prof = atac.tss_profile([f"{p}.tn5.tagAlign.gz"], tss_table(), atac.read_chrom_sizes(cs))
        te = atac.tss_enrichment(prof)
        if not np.isfinite(te["tss_enrichment"]):
            raise RuntimeError(f"TSS enrichment not finite (profile sum {prof.sum()}); check chrom names/sizes")
        np.save(f"{p}.tss_profile.npy", prof)
    step(d, "tss", tss)

    # cleanup (always): raw BAM and the replicate tmp dir
    Path(f"{tmp}/{rep}.raw.bam").unlink(missing_ok=True)
    shutil.rmtree(tmp, ignore_errors=True)
    rq = rep_qc(sample, rep)
    (d / f"{rep}.qc.json").write_text(json.dumps(rq, indent=1))
    log(json.dumps(rq))
    return rq


def _flagstat(path) -> dict:
    out = {}
    for line in Path(path).read_text().splitlines():
        n = int(line.split()[0])
        if "in total" in line:
            out["total"] = n
        elif "primary mapped" in line:
            out["primary_mapped"] = n
        elif "properly paired" in line:
            out["properly_paired"] = n
    return out


def rep_qc(sample, rep) -> dict:
    d = OUT_ROOT / sample / rep
    p = f"{d}/{rep}"
    q: dict = {"sample": sample, "rep": rep}
    fj = sorted(d.glob(f"{rep}.run*.fastp.json"))
    tot = sum(json.loads(f.read_text())["summary"]["before_filtering"]["total_reads"] for f in fj)
    aft = sum(json.loads(f.read_text())["summary"]["after_filtering"]["total_reads"] for f in fj)
    adapt = sum(json.loads(f.read_text()).get("adapter_cutting", {}).get("adapter_trimmed_reads", 0) for f in fj)
    q["read_pairs_input"] = tot // 2
    q["read_pairs_after_trim"] = aft // 2
    q["frac_reads_adapter_trimmed"] = adapt / tot if tot else None
    bt = [f.read_text() for f in sorted(d.glob(f"{rep}.run*.bowtie2.log"))]
    q["bowtie2_overall_alignment_rate"] = [l.split()[0] for t in bt for l in t.splitlines() if "overall alignment rate" in l]
    rm = json.loads((d / f"{rep}.raw_mito.json").read_text())
    q["frac_mito_reads_raw"] = rm["mito_reads"] / rm["mapped_primary_reads"]
    mt, m0, m1, m2, nrf, pbc1, pbc2 = Path(f"{p}.lib_complexity.qc").read_text().split()
    q.update({"pbc_total_fragments": int(mt), "pbc_distinct_fragments": int(m0),
              "NRF": float(nrf), "PBC1": float(pbc1), "PBC2": float(pbc2)})
    fs = _flagstat(f"{p}.filt.flagstat.txt")
    q["filtered_fragments_with_mito"] = fs["total"] // 2
    nd = _flagstat(f"{p}.nodup.flagstat.txt")
    q["frac_duplicate_fragments"] = 1 - nd["total"] / fs["total"]
    nm = _flagstat(f"{p}.nodup.no_chrM.flagstat.txt")
    q["nondup_nonmito_fragments"] = nm["total"] // 2
    q["frac_mito_filtered_nodup"] = 1 - nm["total"] / nd["total"]
    q.update(json.loads((d / f"{rep}.tagalign_counts.json").read_text()))
    q.update(atac.tss_enrichment(np.load(f"{p}.tss_profile.npy")))
    return q


# --------------------------------------------------------------------------- pooled
def postprocess_macs(np_in: str, out_gz: str, cs: Path):
    """ENCODE macs2(): sort by p desc, rename, fix summit, cap 300k, clip to chrom sizes."""
    sh(f"sort -k 8gr,8gr {np_in} "
       f"| awk 'BEGIN{{OFS=\"\\t\"}}{{$4=\"Peak_\"NR; if ($2<0) $2=0; if ($3<0) $3=0; if ($10==-1) $10=$2+int(($3-$2+1)/2.0); print $0}}' "
       f"| awk 'NR<={CAP_NUM_PEAK}' "
       f"| awk 'BEGIN{{OFS=\"\\t\"}} NR==FNR {{s[$1]=$2; next}} ($1 in s) && $3<=s[$1] && $2<$3' {cs} - "
       f"| gzip -nc > {out_gz}")


def bfilt(peak_gz: str, out_gz: str, bl: Path):
    sh(f"bedtools intersect -nonamecheck -v -a {peak_gz} -b {bl} "
       f"| awk 'BEGIN{{OFS=\"\\t\"}} $1 ~ /{CHR_REGEX}/' | gzip -nc > {out_gz}")


def frip_counts(ta: str, peak_gz: str) -> tuple[int, int]:
    """(reads in peaks, total reads) for a tagAlign (ENCODE frip: intersect -wa -u)."""
    n = atac.count_lines(ta)
    k = int(sh_out(f"bedtools intersect -nonamecheck -a <(gzip -dc {ta}) -b <(gzip -dc {peak_gz}) -wa -u | wc -l").split()[0])
    return k, n


def split_pseudoreps(ta_gz: str, out1, out2, seed: int):
    """ENCODE spr_pe: shuffle fragments (line pairs) and split in two exact halves."""
    n = atac.count_lines(ta_gz)
    assert n % 2 == 0
    npair = n // 2
    rng = np.random.default_rng(seed)
    first = np.zeros(npair, dtype=bool)
    first[rng.permutation(npair)[: (npair + 1) // 2]] = True
    proc = subprocess.Popen(["gzip", "-dc", ta_gz], stdout=subprocess.PIPE)
    i = 0
    buf1, buf2 = [], []
    for line in proc.stdout:
        (buf1 if first[i // 2] else buf2).append(line)
        i += 1
        if len(buf1) > 2_000_000:
            out1.write(b"".join(buf1)); buf1 = []
        if len(buf2) > 2_000_000:
            out2.write(b"".join(buf2)); buf2 = []
    out1.write(b"".join(buf1)); out2.write(b"".join(buf2))
    proc.wait()
    assert i == n


def naive_overlap(p1: str, p2: str, pooled: str, out_gz: str, tmp: Path):
    """ENCODE naive_overlap(): pooled peaks overlapping a peak in BOTH p1 and p2 by >=50% of
    either peak."""
    t1, t2, tp = tmp / "ov1.np", tmp / "ov2.np", tmp / "ovp.np"
    for a, b in ((p1, t1), (p2, t2), (pooled, tp)):
        sh(f"gzip -dc {a} > {b}")
    aw = "{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}"
    sh(f"intersectBed -nonamecheck -wo -a {tp} -b {t1} | awk 'BEGIN{{FS=\"\\t\";OFS=\"\\t\"}} {aw}' "
       f"| cut -f 1-10 | sort | uniq "
       f"| intersectBed -nonamecheck -wo -a stdin -b {t2} | awk 'BEGIN{{FS=\"\\t\";OFS=\"\\t\"}} {aw}' "
       f"| cut -f 1-10 | sort | uniq | gzip -nc > {out_gz}")
    for t in (t1, t2, tp):
        t.unlink()


def run_pool(sample: str, reps: list[str], threads: int = 10):
    sd = OUT_ROOT / sample
    pd_ = sd / "pooled"; tmp = sd / "tmp" / "pooled"
    pd_.mkdir(parents=True, exist_ok=True); tmp.mkdir(parents=True, exist_ok=True)
    bl, cs = common_files(tmp)
    tas = [f"{sd}/{r}/{r}.tn5.tagAlign.gz" for r in reps]
    # pooled tagAlign = the replicate tagAligns given together to MACS3 (-t a b c); not
    # written to disk to avoid duplicating them.
    pooled_ta = " ".join(tas)
    pr1, pr2 = f"{tmp}/pooled.pr1.tagAlign", f"{tmp}/pooled.pr2.tagAlign"

    def spr():
        with open(pr1, "wb") as o1, open(pr2, "wb") as o2:
            for r, ta in zip(reps, tas):
                split_pseudoreps(ta, o1, o2, derive_seed(sample, r, "pseudorep"))
    step(pd_, "spr", spr)

    def macs():
        cmds = [
            f"{MACS3} callpeak -t {pooled_ta} -f BED -n pooled --outdir {tmp} -g hs -p 0.01 "
            f"--shift -75 --extsize 150 --nomodel -B --SPMR --keep-dup all --call-summits --tempdir {tmp} "
            f"> {pd_}/macs3.pooled.log 2>&1",
            f"{MACS3} callpeak -t {pr1} -f BED -n pr1 --outdir {tmp} -g hs -p 0.01 "
            f"--shift -75 --extsize 150 --nomodel --keep-dup all --call-summits --tempdir {tmp} "
            f"> {pd_}/macs3.pr1.log 2>&1",
            f"{MACS3} callpeak -t {pr2} -f BED -n pr2 --outdir {tmp} -g hs -p 0.01 "
            f"--shift -75 --extsize 150 --nomodel --keep-dup all --call-summits --tempdir {tmp} "
            f"> {pd_}/macs3.pr2.log 2>&1",
        ]
        # MACS3 needs ~35 bytes RAM per tag; run the three calls concurrently only when the
        # pooled data are small enough (<= 250M tags), else one after another.
        n_tags = sum(json.loads((sd / r / f"{r}.tagalign_counts.json").read_text())["tagalign_reads"] for r in reps)
        if n_tags <= 250_000_000:
            procs = [subprocess.Popen(["bash", "-o", "pipefail", "-c", c]) for c in cmds]
            for c, pr in zip(cmds, procs):
                if pr.wait() != 0:
                    raise RuntimeError("macs3 failed: " + c)
        else:
            for c in cmds:
                sh(c)
        for n in ("pooled", "pr1", "pr2"):
            postprocess_macs(f"{tmp}/{n}_peaks.narrowPeak", f"{pd_}/{n}.pval0.01.300K.narrowPeak.gz", cs)
        for f in tmp.glob("*_summits.bed"):
            f.unlink()
    step(pd_, "macs3", macs)

    def fc():
        sh(f"{MACS3} bdgcmp -t {tmp}/pooled_treat_pileup.bdg -c {tmp}/pooled_control_lambda.bdg "
           f"--outdir {tmp} --o-prefix pooled -m FE > {pd_}/macs3.bdgcmp.log 2>&1")
        n = atac.bedgraph_to_bigwig(f"{tmp}/pooled_FE.bdg", atac.read_chrom_sizes(cs), f"{sd}/{sample}.fc.bigWig")
        log(f"bigWig intervals: {n}")
        for f in ("pooled_FE.bdg", "pooled_treat_pileup.bdg", "pooled_control_lambda.bdg"):
            (tmp / f).unlink(missing_ok=True)
    step(pd_, "fc_bigwig", fc)

    def idr():
        for n in ("pooled", "pr1", "pr2"):
            sh(f"gzip -dc {pd_}/{n}.pval0.01.300K.narrowPeak.gz > {tmp}/{n}.np")
        out = f"{tmp}/idr.out"
        sh(f"{IDR_BIN} --samples {tmp}/pr1.np {tmp}/pr2.np --peak-list {tmp}/pooled.np "
           f"--input-file-type narrowPeak --output-file {out} --rank p.value "
           f"--soft-idr-threshold {IDR_THRESH} --use-best-multisummit-IDR "
           f"--log-output-file {pd_}/idr.log > {pd_}/idr.stdout 2>&1")
        thr = -math.log10(IDR_THRESH)
        sh(f"awk 'BEGIN{{OFS=\"\\t\"}} NR==FNR {{s[$1]=$2; next}} ($1 in s) && $3<=s[$1] && $12>={thr} "
           f"{{if ($2<0) $2=0; print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}}' {cs} {out} "
           f"| sort | uniq | sort -grk8,8 | gzip -nc > {pd_}/pooled.idr0.05.12col.narrowPeak.gz")
        sh(f"gzip -dc {pd_}/pooled.idr0.05.12col.narrowPeak.gz | cut -f 1-10 | gzip -nc > {tmp}/idr.np.gz")
        bfilt(f"{tmp}/idr.np.gz", f"{sd}/{sample}.idr.narrowPeak.gz", bl)
        sh(f"gzip -c {out} > {pd_}/idr.out.gz")
        for f in ("pooled.np", "pr1.np", "pr2.np", "idr.out", "idr.np.gz"):
            (tmp / f).unlink(missing_ok=True)
    step(pd_, "idr", idr)

    def overlap():
        ov = f"{tmp}/overlap.np.gz"
        naive_overlap(f"{pd_}/pr1.pval0.01.300K.narrowPeak.gz", f"{pd_}/pr2.pval0.01.300K.narrowPeak.gz",
                      f"{pd_}/pooled.pval0.01.300K.narrowPeak.gz", ov, tmp)
        bfilt(ov, f"{sd}/{sample}.overlap.narrowPeak.gz", bl)
        bfilt(f"{pd_}/pooled.pval0.01.300K.narrowPeak.gz", f"{sd}/{sample}.pooled.pval0.01.300K.bfilt.narrowPeak.gz", bl)
        Path(ov).unlink(missing_ok=True)
    step(pd_, "overlap", overlap)

    def fr():
        res = {}
        for kind in ("idr", "overlap"):
            pk = f"{sd}/{sample}.{kind}.narrowPeak.gz"
            K = N = 0
            for r, ta in zip(reps, tas):
                k, n = frip_counts(ta, pk)
                res[f"frip_{r}_{kind}"] = k / n
                K += k; N += n
            res[f"frip_pooled_{kind}"] = K / N
        (pd_ / "frip.json").write_text(json.dumps(res, indent=1))
    step(pd_, "frip", fr)

    # cleanup pseudorep tagAligns
    for f in (pr1, pr2):
        Path(f).unlink(missing_ok=True)
    shutil.rmtree(tmp, ignore_errors=True)
    try:
        (sd / "tmp").rmdir()
    except OSError:
        pass

    qc = summarize(sample, reps)
    (sd / "qc.json").write_text(json.dumps(qc, indent=1))
    log(json.dumps(qc["gate"], indent=1))
    return qc


def _ge(x, thr) -> bool:
    """Gate comparison; NaN/None/inf never pass."""
    return x is not None and bool(np.isfinite(x)) and x >= thr


def summarize(sample, reps) -> dict:
    sd = OUT_ROOT / sample
    pd_ = sd / "pooled"
    reps_qc = {r: rep_qc(sample, r) for r in reps}
    prof = sum(np.load(f"{sd}/{r}/{r}.tss_profile.npy") for r in reps)
    pooled = atac.tss_enrichment(prof)
    fr = json.loads((pd_ / "frip.json").read_text())
    n_idr = atac.count_lines(f"{sd}/{sample}.idr.narrowPeak.gz")
    n_ov = atac.count_lines(f"{sd}/{sample}.overlap.narrowPeak.gz")
    n_pool = atac.count_lines(f"{sd}/{sample}.pooled.pval0.01.300K.bfilt.narrowPeak.gz")
    merged = sum(q["nondup_nonmito_fragments"] for q in reps_qc.values())
    for r, q in reps_qc.items():
        q["frip_idr"] = fr[f"frip_{r}_idr"]; q["frip_overlap"] = fr[f"frip_{r}_overlap"]
    gate = {
        "tss_enrichment": {r: (q["tss_enrichment"], _ge(q["tss_enrichment"], GATE["tss_enrichment_min"])) for r, q in reps_qc.items()},
        "tss_enrichment_pooled": (pooled["tss_enrichment"], _ge(pooled["tss_enrichment"], GATE["tss_enrichment_min"])),
        "NRF": {r: (q["NRF"], _ge(q["NRF"], GATE["nrf_min"])) for r, q in reps_qc.items()},
        "PBC1": {r: (q["PBC1"], _ge(q["PBC1"], GATE["pbc1_min"])) for r, q in reps_qc.items()},
        "nondup_nonmito_fragments_merged": (merged, _ge(merged, GATE["nondup_nonmito_fragments_merged_min"])),
        "idr_peaks": (n_idr, _ge(n_idr, GATE["idr_peaks_min"])),
    }
    ok = all(v[1] for k, v in gate.items() if isinstance(v, tuple)) and \
        all(x[1] for k, v in gate.items() if isinstance(v, dict) for x in v.values())
    gate["PASS"] = bool(ok)
    return {
        "sample": sample, "replicates": reps_qc,
        "pooled": {"tss_enrichment": pooled["tss_enrichment"], "tss_enrichment_bin10": pooled["tss_enrichment_bin10"],
                   "nondup_nonmito_fragments_merged": merged,
                   "n_idr_peaks": n_idr, "n_overlap_peaks": n_ov, "n_pooled_peaks_bfilt": n_pool,
                   "frip_idr": fr["frip_pooled_idr"], "frip_overlap": fr["frip_pooled_overlap"]},
        "gate_thresholds": GATE, "gate": gate,
        "outputs": {k: f"data/atac_local/{sample}/{sample}.{k}" for k in
                    ("fc.bigWig", "idr.narrowPeak.gz", "overlap.narrowPeak.gz", "pooled.pval0.01.300K.bfilt.narrowPeak.gz")},
        "generated": time.strftime("%Y-%m-%d %H:%M"),
    }


def finalize_storage(sample, reps, limit_gb=10.0):
    """Keep the final filtered BAMs only if they total < limit_gb; tagAligns are always kept."""
    sd = OUT_ROOT / sample
    bams = [sd / r / f"{r}.nodup.no_chrM.bam" for r in reps]
    tot = sum(b.stat().st_size for b in bams if b.exists()) / 1e9
    if tot >= limit_gb:
        log(f"final BAMs total {tot:.1f} GB >= {limit_gb} GB: deleting BAMs, keeping tagAligns")
        for b in bams:
            b.unlink(missing_ok=True); Path(str(b) + ".bai").unlink(missing_ok=True)
    else:
        log(f"final BAMs total {tot:.1f} GB < {limit_gb} GB: kept")
    return tot


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("rep"); a.add_argument("--sample", required=True); a.add_argument("--rep", required=True)
    a.add_argument("--r1", required=True); a.add_argument("--r2", required=True); a.add_argument("--threads", type=int, default=10)
    b = sub.add_parser("pool"); b.add_argument("--sample", required=True); b.add_argument("--reps", nargs="+", required=True)
    b.add_argument("--threads", type=int, default=10)
    c = sub.add_parser("sheet"); c.add_argument("--sheet", required=True); c.add_argument("--sample", required=True)
    c.add_argument("--threads", type=int, default=10)
    args = ap.parse_args()
    if args.cmd == "rep":
        run_rep(args.sample, args.rep, args.r1.split(","), args.r2.split(","), args.threads)
    elif args.cmd == "pool":
        run_pool(args.sample, args.reps, args.threads)
        finalize_storage(args.sample, args.reps)
    else:
        sh_ = pd.read_csv(args.sheet, sep="\t", comment="#")
        sh_ = sh_[sh_["sample"] == args.sample]
        if sh_.empty:
            raise SystemExit(f"no rows for {args.sample}")
        for _, row in sh_.iterrows():
            r1 = ",".join(str(ROOT / x) for x in row["r1"].split(","))
            r2 = ",".join(str(ROOT / x) for x in row["r2"].split(","))
            run_rep(args.sample, row["rep"], r1.split(","), r2.split(","), args.threads)
        reps = list(sh_["rep"])
        run_pool(args.sample, reps, args.threads)
        finalize_storage(args.sample, reps)


if __name__ == "__main__":
    main()
