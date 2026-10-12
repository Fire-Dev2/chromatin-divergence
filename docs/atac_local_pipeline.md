# Local ENCODE-equivalent ATAC-seq pipeline

Code: `scripts/atac_pipeline.py` (driver), `src/cdtest/atac.py` (TSS enrichment, bedGraph→bigWig,
bigWig means; tests in `tests/test_atac.py`), `scripts/atac_validate_encode.py` (comparison with
ENCODE), `scripts/atac_samples.tsv` (sample sheet). Outputs go to `data/atac_local/<sample>/`
(git-ignored). Status: **provisional** stand-in for the ENCODE ATAC pipeline v1.9.1 run on OSC
(log/RESEARCH_LOG.md, "Compute route change"); the OSC rerun stays confirmatory.

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
protein-coding genes in `data/annot/genes_gencode_v50.tsv` (19,9xx unique TSSs).

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
| 12 | Peak post-processing | `sort -k 8gr,8gr`, rename `Peak_N`, fill summit, `head -n 300000`, drop peaks beyond chrom ends | identical awk; `bedClip` | Same |
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
