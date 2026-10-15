# Local ENCODE-equivalent ATAC-seq pipeline

Code: `scripts/atac_pipeline.py` (driver), `src/cdtest/atac.py` (TSS enrichment, bedGraph→bigWig,
bigWig means; tests in `tests/test_atac.py`), `scripts/atac_validate_encode.py` (comparison with
ENCODE), `scripts/atac_samples.tsv` (sample sheet). Outputs go to `data/atac_local/<sample>/`
(git-ignored). Status: **provisional** stand-in for the ENCODE ATAC pipeline v1.9.1 run on OSC
(a deviation from the analysis plan); the OSC rerun stays confirmatory.

## Environment

| Tool | Version | Where |
|---|---|---|
| fastp | 1.4.0 | Homebrew |
| bowtie2 | 2.5.5 | Homebrew |
| samtools | 1.24 | Homebrew |
| bedtools | 2.31.1 | Homebrew |
| MACS3 | 3.0.5 | `env/venv` |
| IDR (Kundaje lab, nboley/idr) | 2.0.3 | `env/idr-venv/bin/idr` (Python 3.11). Two `numpy.int` → `int` edits in `idr/idr.py` (lines 298–299) were needed for numpy 1.26; numerics unchanged. |
| pyBigWig, numpy, pandas, scipy | — | `env/venv` |

ENCODE v1.9.1 pins (scripts/requirements.txt of that tag): cutadapt 2.5, bowtie2 2.3.4.3,
samtools 1.9, picard 2.20.7, macs2 2.2.4, idr 2.0.4.2.

Reference: GRCh38 no-alt analysis set (ENCODE ENCFF110MCL bowtie2 index, the same reference as
ENCODE ATAC v1.9.1). Chrom sizes taken from the index (`bowtie2-inspect -s`) →
`data/atac_local/GRCh38_no_alt.chrom.sizes`. Blacklist: ENCODE hg38 blacklist v2
(`data/annot/blacklist/hg38-blacklist.v2.bed.gz`). TSSs: Ensembl-canonical TSSs of
protein-coding genes in `data/annot/genes_gencode_v50.tsv` (19,9xx unique TSSs; ~19.6k after window/chrom checks).

## How to run

```bash
cd ~/chromatin-divergence; unset PYTHONPATH; source env/venv/bin/activate
# all replicates + pooled steps of one sample from the sheet
python scripts/atac_pipeline.py sheet --sheet scripts/atac_samples.tsv --sample jurkat_GSE273644
# or step by step
python scripts/atac_pipeline.py rep  --sample S --rep rep1 --r1 R1.fq.gz[,run2_R1.fq.gz] --r2 R2.fq.gz[,...]
python scripts/atac_pipeline.py pool --sample S --reps rep1 rep2 rep3
# validation against ENCODE outputs of the same experiment
python scripts/atac_validate_encode.py --sample k562_ENCSR868FGK \
    --encode-bw data/atac/K562__ENCFF102ARJ.bigWig --encode-peaks data/atac/K562__ENCFF948AFM.bed.gz
```

Each step writes a `.<step>.done` marker, so an interrupted run resumes. Threads: bowtie2 `-p 10`,
samtools/sort `-@ 5`; replicates run one after another.

**RPE1 (not yet run — FASTQs still downloading).** Once all four files are in
`data/manifest_done.txt`:

```bash
python scripts/atac_pipeline.py sheet --sheet scripts/atac_samples.tsv --sample rpe1_4DNES48BA2C8
```

(rep1 = 4DNFIDNZZ2A9 R1 + 4DNFI14DT5NE R2; rep2 = 4DNFIWGX7K9Y R1 + 4DNFIY8ZSQC1 R2.)

## Steps and correspondence with ENCODE ATAC v1.9.1

| # | Step | Local command / rule | ENCODE v1.9.1 | Same? |
|---|---|---|---|---|
| 1 | Adapter trimming | `fastp --adapter_sequence CTGTCTCTTATACACATCT --adapter_sequence_r2 CTGTCTCTTATACACATCT --disable_quality_filtering --disable_trim_poly_g --length_required 5` (JSON report kept). Output streamed (interleaved) to bowtie2; no trimmed FASTQ written. | `cutadapt -e 0.1 -m 5 -a ADAPT -A ADAPT` (adapter auto-detected, Nextera for these libraries) | Equivalent. fastp also uses PE overlap analysis to find adapters; quality filtering and poly-G trimming are switched off to match cutadapt. |
| 2 | Alignment | `bowtie2 -X 2000 --mm -p 10 -x GRCh38_no_alt --interleaved -` | `bowtie2 -k 5 -X2000 --mm` (atac.wdl default `multimapping = 4`) | Differs: no `-k`; see filtering. |
| 3 | Filtering | `samtools view -F 1804 -f 2 -q 30 \| samtools fixmate -r -m \| samtools view -F 1804 -f 2` (bowtie2 output is mate-collated, so ENCODE's name sort is not needed; orphan check asserts R1 = R2 counts) | Same commands when `multimapping = 0`; with the wdl default `multimapping = 4`, ENCODE instead keeps reads with ≤4 alignments via `assign_multimappers.py` and skips the MAPQ filter | Uses the MAPQ ≥ 30 unique-read rule as specified for this project. This drops multi-mapping reads that ENCODE (default config) would assign at random; expected effect is lower signal in repeats/segmental duplications, small at promoters. |
| 4 | Library complexity | `bedtools bamtobed -bedpe` on the filtered, **pre-dedup**, name-collated BAM; key `$1,$2,$4,$6,$9,$10`; chrM excluded; `sort \| uniq -c`; NRF = distinct/total, PBC1 = one-read/distinct, PBC2 = one-read/two-read | `pbc_qc_pe` identical (awk copied) | Same |
| 5 | Duplicates | `samtools sort` → `samtools markdup -r -s` (needs fixmate `-m` tags) | Picard MarkDuplicates, then `samtools view -F 1804 -f 2` | Equivalent criteria (both mates' 5′ ends + orientation); samtools markdup used instead of Picard. |
| 6 | chrM removal | `samtools view -e 'rname!="chrM" && rname!="MT"'` on the deduplicated BAM → `<rep>.nodup.no_chrM.bam` | `filter_chrs = ['chrM','MT']` on nodup BAM | Same |
| 7 | Blacklist | fragments with **either** mate touching blacklist v2 removed (`bedtools pairtobed -type neither`) before tagAlign | Not applied to reads; peaks only (`.bfilt`) | Extra step (project spec). Also applied to peaks as ENCODE does. |
| 8 | tagAlign + Tn5 shift | `bamtobed -bedpe -mate1` → 2 lines/fragment (`N 1000 strand`) → `awk '$6=="+"{$2+=4} $6=="-"{$3-=5}'` | `bam2ta_pe` + `tn5_shift_ta` (same awk) | Same |
| 9 | Pooling | the replicate tagAligns are passed together to MACS3 (`-t rep1 rep2 rep3`), not concatenated on disk | `pool_ta` (cat) | Same reads |
| 10 | Pseudoreplicates | each replicate's fragments (line pairs) shuffled with a seed derived from master seed 20261011 (sha256 of sample/rep) and split into exact halves; rep halves pooled → pooled-pr1 / pooled-pr2 | `spr_pe` per replicate (shuf with openssl seed), then pooled | Same design, different RNG |
| 11 | Peak calling | `macs3 callpeak -f BED -g hs -p 0.01 --shift -75 --extsize 150 --nomodel -B --SPMR --keep-dup all --call-summits` (pr1/pr2 without `-B --SPMR`, which do not change peaks) | `macs2 callpeak` same flags (MACS2 2.2.4) | MACS3 3.0.5 vs MACS2 2.2.4 |
| 12 | Peak post-processing | `sort -k 8gr,8gr`, rename `Peak_N`, fill summit, first 300,000 (`awk 'NR<=300000'`; `head` raises SIGPIPE under pipefail), drop peaks beyond chrom ends | identical awk; `bedClip` | Same |
| 13 | FC signal | `macs3 bdgcmp -t treat_pileup -c control_lambda -m FE` → clip to chrom sizes, sort → bigWig with pyBigWig (`<sample>.fc.bigWig`) | `macs2 bdgcmp -m FE` → `slop -b 0 \| bedClip` → `bedGraphToBigWig` | Same values (float32) |
| 14 | IDR | `idr --samples pr1 pr2 --peak-list pooled --input-file-type narrowPeak --rank p.value --soft-idr-threshold 0.05 --use-best-multisummit-IDR`; keep col12 ≥ −log10(0.05); `sort \| uniq \| sort -grk8`; blacklist + `chr[0-9XY]+` filter → `<sample>.idr.narrowPeak.gz` | Same (`encode_task_idr.py`), output `pooled-pr1_vs_pooled-pr2.idr0.05.bfilt.narrowPeak` | Same |
| 15 | Naive overlap | pooled peaks overlapping a pr1 peak **and** a pr2 peak by ≥ 50 % of either peak; bfilt → `<sample>.overlap.narrowPeak.gz` | `naive_overlap` (same awk) | Same |
| 16 | FRiP | fraction of tagAlign reads overlapping peaks (`intersect -wa -u`), per replicate and pooled, for IDR and overlap peaks | `frip` | Same definition (ENCODE reports per-replicate FRiP on its own peaks) |
| 17 | TSS enrichment | Tn5 insertion sites (shifted read 5′ ends) summed per bp over ±2 kb of every protein-coding canonical TSS (strand-oriented); normalised by the mean of the outer 100 bp at both ends **of the aggregate profile** (so zero flanks at single TSSs cannot occur); 20-bp centred boxcar; max within ±250 bp. Also reported: 10-bp-bin variant (`tss_enrichment_bin10`). A non-finite value raises an error, and the QC gate treats NaN as fail. | ATAqC `make_tss_plot`: metaseq read coverage (reads shifted by −read_len/2) in 400 bins over ±2 kb, normalised by the 100 bp at each end, max of the binned mean; ENCODE's own TSS file | Same normalisation; insertion-based rather than read-coverage-based, different TSS set. Values are calibrated against ENCODE below. |
| 18 | mito fraction | chrM reads / primary mapped reads in the raw BAM | idxstats of raw BAM | Same |

Not implemented (ENCODE extras not needed here): p-value signal track, per-replicate peaks and
true-replicate IDR, cross-correlation, fragment-length/GC/annotation-enrichment QC, IDR plots.

Files kept per sample (`data/atac_local/<sample>/`): `<sample>.fc.bigWig`,
`<sample>.idr.narrowPeak.gz`, `<sample>.overlap.narrowPeak.gz`,
`<sample>.pooled.pval0.01.300K.bfilt.narrowPeak.gz`, `qc.json`; per replicate
`<rep>/<rep>.tn5.tagAlign.gz`, `<rep>.nodup.no_chrM.bam` (+ .bai; deleted automatically if all
final BAMs total ≥ 10 GB), logs, fastp JSON, flagstats, `lib_complexity.qc`, TSS profile;
`pooled/` keeps unfiltered pooled/pr1/pr2 peaks, IDR raw output and logs. Raw/unfiltered BAMs,
name-sorted BAMs, MACS bedGraphs and pseudoreplicate tagAligns are deleted.

## Validation against ENCODE (K562, ENCSR868FGK) — run 2026-10-11/12

Input: all three ENCODE biological replicates, each **subsampled to its first 40,000,000 read
pairs** (2×100 bp, NovaSeq) because the full set is ~45 GB / 400M pairs
(`data/atac_fastq/k562_ENCSR868FGK/get.sh`: whole runs downloaded and md5-checked, `head`
of the concatenated runs; rep1 = run ENCFF962SRP/ENCFF420XJE, rep2 = ENCFF098UCE/ENCFF703BGR
+ ENCFF903FIK/ENCFF924YGE, rep3 = ENCFF354EXH/ENCFF575ZTZ + ENCFF827WXQ/ENCFF742FCG).
ENCODE processed all reads (≈ 1.2–1.7× more usable fragments per replicate) with
multimapping = 4. Reference outputs: ENCFF102ARJ (pooled FC bigWig), ENCFF948AFM (pooled
pseudoreplicate IDR peaks, 181,340). Results: `results/atac_local/k562_ENCSR868FGK.*.json`.

| Comparison (log2(1 + mean FC)) | n | Pearson | Spearman |
|---|---|---|---|
| Union of both IDR peak sets, summit ±250 bp (merged, blacklist removed) | 95,342 | 0.840 | 0.901 |
| TSS ±1 kb, all 53,476 genes (blacklisted TSS windows removed) | 53,476 | 0.972 | **0.967** |
| TSS ±1 kb, protein-coding | 19,685 | 0.980 | 0.986 |

Peak overlap: local IDR peaks 150,128 vs ENCODE 181,340; bp Jaccard (merged) 0.691; 82.6 % of
ENCODE IDR peaks overlap a local IDR peak; 95.5 % of local IDR peaks overlap an ENCODE IDR peak;
97.6 % of ENCODE IDR peaks overlap the local top-300k pooled peaks.

**Pre-specified pass criterion (Spearman ≥ 0.9 over TSS windows): PASS (0.967).**

QC side by side (local on 40M pairs / ENCODE on all reads, ENCODE values from the portal QC
objects of the v1.9.1 analysis ENCAN824UKX):

| Metric | rep1 local / ENCODE | rep2 | rep3 |
|---|---|---|---|
| Alignment rate | 98.3 % / 98.3 % | 98.5 % / 98.5 % | 98.3 % / 98.3 % |
| Fraction chrM (raw) | 0.118 / 0.119 | 0.152 / 0.148 | 0.113 / 0.114 |
| NRF | 0.947 / 0.889 | 0.952 / 0.907 | 0.947 / 0.906 |
| PBC1 | 0.948 / 0.902 | 0.953 / 0.921 | 0.947 / 0.917 |
| TSS enrichment | 15.4 / 13.6 | 15.7 / 14.1 | 15.8 / 14.0 |
| Non-dup non-mito fragments | 28.3M / 114.5M | 27.3M / 76.4M | 28.3M / 82.5M |
| FRiP (pooled IDR peaks) | 0.329 | 0.341 | 0.339 (ENCODE pooled 0.356) |

NRF/PBC are higher locally as expected at lower depth. The local TSS enrichment (insertion
based) runs ~1.12× ENCODE's ATAqC value on the same library, so thresholds such as ≥ 5 are
not materially easier to pass. Lower peak-union Pearson (0.84) mainly reflects weak peaks
present in only one set (depth and multimapping differences); promoter signal agrees closely.

## Jurkat GSE273644 (primary) — QC gate

3 replicates, 2×42 bp, all reads (40.1M / 31.4M / 34.8M pairs). Results:
`results/atac_local/jurkat_GSE273644.qc.json`.

| Metric | rep1 (SRR30058502) | rep2 (SRR30058500) | rep3 (SRR30058498) | Gate | Pass |
|---|---|---|---|---|---|
| Alignment rate | 92.3 % | 93.4 % | 94.1 % | — | |
| Fraction chrM (raw) | 0.229 | 0.176 | 0.164 | — | |
| Duplicate fraction (filtered) | 0.361 | 0.159 | 0.156 | — | |
| TSS enrichment | 36.3 | 29.6 | 28.3 (pooled 31.3) | ≥ 5 | yes |
| NRF | **0.789** | 0.933 | 0.934 | ≥ 0.8 | **rep1 no** (misses by 0.011); reps 2–3 yes |
| PBC1 | 0.827 | 0.939 | 0.940 | ≥ 0.7 | yes |
| PBC2 | 6.39 | 17.40 | 17.53 | — | |
| Non-dup non-mito fragments | 18.46M | 18.61M | 21.34M | merged 58.41M ≥ 25M | yes |
| FRiP (pooled IDR / overlap peaks) | 0.540 / 0.595 | 0.455 / 0.493 | 0.432 / 0.471 | — | |
| Pooled IDR peaks (bfilt) | | | 110,488 | ≥ 70k | yes |
| Pooled naive-overlap peaks | | | 183,045 | — | |

Strict reading of the gate: FAIL on one criterion (rep1 NRF 0.789 < 0.8); every other
criterion passes with a wide margin. All three replicates were kept; this is reported as a
deviation from the analysis plan.

Outputs: `data/atac_local/jurkat_GSE273644/jurkat_GSE273644.{fc.bigWig,idr.narrowPeak.gz,overlap.narrowPeak.gz,pooled.pval0.01.300K.bfilt.narrowPeak.gz}`,
`qc.json`; filtered BAMs kept (5.0 GB total).

## Run notes

* Wall time on the shared Mac (other jobs running): Jurkat ~45–70 min per replicate,
  K562 40M×100 bp ~45–130 min per replicate, pooled steps 10–15 min.
* K562 subsampled FASTQs (11 GB) are kept in `data/atac_fastq/k562_ENCSR868FGK/`; full runs were
  deleted after subsampling. K562 final BAMs kept (8.4 GB < 10 GB).
* MACS3 pooled/pr1/pr2 run concurrently only for ≤ 250M pooled tags (memory), else sequentially.
