# Accessibility features (ANALYSIS_PLAN.md §3, signal per plan amendment v1.2)

Code: `src/cdtest/features.py`. Sanity checks and figures: `scripts/features_sanity.py`.
Tests: `tests/test_features.py`. One-command rebuild of both versions: `scripts/build_features.sh`
(~1.5 min on 10 processes).

## What this does

For every gene we want a number that says how open its chromatin is in each cell line.
We make three of them:

* **Promoter**: how open the DNA is right around the gene's start (TSS ± 1 kb).
* **Distal**: how much open DNA (likely enhancers) there is 2.5–100 kb away. Nearer and
  more open regions count more.
* **Combined**: the average of the two after putting both on the same scale.

The test (Aim 2) uses the *difference* of a feature between two lines,
`pair_feature(A, B, kind)` = feature in A − feature in B.

All lines are processed together, because two steps (quantile normalisation and the
pooled z-score) compare lines with each other. **Adding a line changes every line's
values a little**, so the whole build is re-run whenever a line is added or replaced.

## Two signal versions

| version | `signal=` | raw signal per region / promoter window | outputs |
|---|---|---|---|
| **primary** (plan v1.2) | `"cpm"` (default) | Tn5 insertions in the window per million insertions of the line (chr1–22,X), replicates pooled | `gene_features.parquet`, `region_signal.parquet`, `region_signal_raw.parquet`, `build_info.json`, `promoter_cpm_replicates.parquet` |
| sensitivity (v1.0 track) | `"fc"` | mean MACS2 fold-change bigWig value | `gene_features_fc.parquet`, `region_signal_fc.parquet`, `region_signal_raw_fc.parquet`, `build_info_fc.json` |

Why: the fold-change track's local background compresses strong promoters, so promoter FC
was *negatively* rank-correlated with expression among expressed genes (see
`log/RESEARCH_LOG.md`, "Signal-track decision"). Insertion CPM fixes the sign.

Everything after the raw signal (log, quantile normalisation, CN, distal weighting, z) is
identical between the versions, and both use the same union regions.

## Which lines and files

`config/feature_lines.tsv` lists, per line, a glob of replicate insertion files
(`insertions`), a fold-change bigWig (`fc_bigwig`), a peak file and the name of its
copy-number track. A line can have several candidate rows: the build uses the row with the
lowest `priority` whose files *for the requested signal* exist (cpm: ≥ 1 insertion file +
peaks; fc: bigWig + peaks). Lines without usable files are skipped and listed in build_info.

* Insertions: `data/atac_ins/<Line>__rep<k>.npz`, one sorted 0-based int32 array of Tn5
  insertion positions per chromosome (chr1–22,X), `.json` side-car with `n_insertions`.
  ENCODE lines: ENCODE v1.9.1 filtered BAMs re-filtered `-f 2 -F 3852 -q 30`, Tn5-shifted.
  Jurkat / RPE1: local pipeline tagAligns (`docs/atac_local_pipeline.md`; fragments touching
  the blacklist were already dropped there — harmless, blacklisted regions are removed from
  the features anyway).
* Peaks: ENCODE IDR thresholded peaks for the 9 ENCODE lines; local pipeline IDR peaks for
  Jurkat (`data/atac_local/jurkat_GSE273644/…idr.narrowPeak.gz`) and RPE1
  (`data/atac_local/rpe1_4DNES48BA2C8/…idr.narrowPeak.gz`). The 4DN RPE1 bigWig/bigBed row
  remains only as a priority-2 fallback.
* FC bigWigs: ENCODE pooled FC; local pipeline `.fc.bigWig` for Jurkat and RPE1.

Replicates: K562, HepG2, GM12878, A549, Jurkat 3; HCT116, MCF-7, IMR-90, PC-3, Panc1, RPE1 2.

## Steps and every parameter

All coordinates are 0-based half-open (BED). Only chr1–22 and chrX are used.

1. **Summits.** narrowPeak column 10 (summit offset) added to the start; if it is −1, the
   peak midpoint. bigBed narrowPeak files are read with pyBigWig `entries()`.
2. **Union region set.** Window `[summit − 250, summit + 250)`, clipped to the chromosome.
   Windows of all 11 lines are pooled and merged when they overlap or touch (`bedtools merge
   -d 0`). Merged regions overlapping any ENCODE blacklist v2 region are dropped entirely.
3. **Raw region signal.**
   * cpm: for each replicate, insertions with `start ≤ p < end` are counted with
     `np.searchsorted` on the sorted per-chromosome arrays (npz keys are read lazily, one
     process per replicate). Counts and chr1–22,X totals are summed over the line's
     replicates; CPM = count / total × 10⁶. Not length-normalised (region lengths are the same
     in every line, so this cancels in between-line differences).
   * fc: exact mean FC (`pyBigWig stats(type="sum", exact=True)` / length; missing = 0).
4. **Normalisation** (both versions): `log2(x + 1)`, quantile normalisation across lines
   (reference = mean of sorted columns; ties get the mean of their ranks' reference values),
   then **minus the line's log2 relative copy number** over the window (bp-weighted mean of
   segments; gaps take the nearest segment; `annot.CNLookup`) → `s_r`.
5. **Promoter feature.** Steps 3–4 on `[tss0 − 1000, tss0 + 1000)` read directly from the
   insertions / bigWig. Quantile normalisation uses the genes × lines promoter matrix itself.
   TSS = Ensembl-canonical transcript start, GENCODE v50 (54,710 protein-coding + lncRNA genes).
6. **Distal feature.** `log2(1 + Σ_r 2^{s_r} · exp(−d_r / 25,000))` over union regions whose
   centre lies at 2,500 ≤ |centre − tss0| ≤ 100,000 bp (both sides). No region → 0.
7. **Combined.** `0.5·z(promoter) + 0.5·z(distal)`; each z uses one mean and SD (ddof = 0)
   pooled over all genes × all 11 lines (constants stored in build_info).

### Copy number

| feature line | CN track (`cn_line`) | source |
|---|---|---|
| K562, HepG2, HCT116, A549 | same name | DepMap 24Q4 segments |
| Jurkat | Jurkat | DepMap 24Q4 (WES) |
| MCF-7, PC-3, Panc1 | MCF7, PC3, PANC1 | DepMap 24Q4 segments (names differ) |
| RPE1 | RPE1 | karyotype track: +log2(3/2) on 10q21.2–qter, 0 elsewhere |
| GM12878, IMR-90 | GM12878, IMR90 | **not in DepMap 24Q4 → 0 everywhere (no correction; flagged)** |

## Outputs (`data/features/`, git-ignored)

| file | content |
|---|---|
| `region_signal[_fc].parquet` | region_id, chrom, start, end, one column per line = `s_r` (float32) |
| `region_signal_raw[_fc].parquet` | same, raw CPM (or raw mean FC) |
| `gene_features[_fc].parquet` | long: gene_id, gene_name, line, promoter, distal, combined, promoter_raw (CPM or FC), promoter_cn_log2, n_distal_regions |
| `promoter_cpm_replicates.parquet` | gene_id + one column per replicate (`<Line>__rep<k>`): TSS ± 1 kb CPM using that replicate's own total |
| `build_info[_fc].json` | signal, chosen files + md5 (insertion npz files for cpm), insertion totals per line / replicate, skipped lines, counts, parameters, pooled z constants |

Python access:

```python
from cdtest.features import pair_feature, feature_matrix
x = pair_feature("RPE1", "K562", "combined")               # primary (cpm)
x_fc = pair_feature("RPE1", "K562", "combined", signal="fc")  # sensitivity
P = feature_matrix("promoter")                             # genes x lines
```

CLI: `python -m cdtest.features --signal {cpm,fc,both} [--procs 10] [--no-md5]`.

## Sanity checks (results/features/sanity.md, build of 2026-10-12, 11 lines)

* 324,842 union regions (1.96 M summit windows, 327,145 merged, 2,303 blacklisted).
* Replicate promoter CPM Spearman 0.945–0.995 over all genes (Jurkat lowest, 0.945–0.953).
* **CPM (primary): 11/14 checks pass.** Promoter vs baseline expression over the 5,262 Aim 1
  genes Spearman +0.111 (RPE1), +0.113 (K562) (raw CPM before log/QN/CN: +0.115, +0.155);
  GAPDH, ACTB, EEF1A1 ≥ 94th percentile in every line, RPL13A 73rd–94th (fail); HBG1, GATA1,
  ALB, APOA1, CD3E, LCK pass; KRT8 / KRT18 fail only because K562 ≥ RPE1 (both lines express
  them; RPE1 and HepG2 do exceed Jurkat). Between-line promoter r 0.80–0.93.
* **FC (sensitivity): 7/14.** Expression Spearman −0.147 / −0.090, all four housekeeping
  genes fail, between-line promoter r 0.64–0.86.
* CPM vs FC agreement of x (all genes, Pearson): RPE1−K562 0.72–0.74, K562−Jurkat 0.71–0.74,
  RPE1−HepG2 0.68–0.70 over the three kinds; lower within the Aim 1 genes (promoter 0.33–0.45),
  i.e. the two versions are not interchangeable for the expressed genes the test uses.
