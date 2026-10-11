# Accessibility features (ANALYSIS_PLAN.md §3)

Code: `src/cdtest/features.py`. Sanity checks and figure: `scripts/features_sanity.py`.
Tests: `tests/test_features.py`. One-command rebuild: `scripts/build_features.sh`.

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

## Which lines and files

`config/feature_lines.tsv` lists, per line, a fold-change (FC) bigWig, a peak file and the
name of its copy-number track. A line can have several candidate rows: the build uses the
row with the lowest `priority` whose two files both exist. Lines without existing files are
skipped and listed in `build_info.json`.

* RPE1: priority 1 = local pipeline output `data/atac_local/rpe1_4DNES48BA2C8/…`
  (not yet present); priority 2 = provisional 4DN FC bigWig + conservative-peak bigBed.
* Jurkat: `data/atac_local/jurkat_GSE273644/…` (not yet present → skipped).
* K562, HepG2 and 7 decoys: ENCODE pooled FC bigWig + IDR thresholded peaks.

When Jurkat / local RPE1 appear, just run `scripts/build_features.sh` (no code edits).
Make sure the pipeline has finished writing both files first.

## Steps and every parameter

All coordinates are 0-based half-open (BED). Only chr1–22 and chrX are used (chrY and
unplaced contigs are dropped: lines differ in sex and CN tracks cover chr1–22,X).

1. **Summits.** narrowPeak column 10 (summit offset) added to the start; if it is −1, the
   peak midpoint. bigBed narrowPeak files are read with pyBigWig `entries()`.
2. **Union region set.** Window `[summit − 250, summit + 250)` (`SUMMIT_HALF = 250`),
   clipped to the chromosome. Windows of all lines are pooled and merged when they overlap
   or touch (like `bedtools merge -d 0`). Merged regions that overlap any ENCODE blacklist v2
   region are dropped (the whole region, not just the overlapping part).
3. **Region signal.** Exact mean FC per region and line (`pyBigWig stats(type="sum",
   exact=True)` / region length; bases without data count as 0).
   Then `log2(x + 1)` (`PSEUDO = 1`), then quantile normalisation across lines
   (reference = mean of the sorted columns; ties receive the average of their ranks'
   reference values), then **minus the line's log2 relative copy number** over the region
   (bp-weighted mean of segments; gaps take the nearest segment; `annot.CNLookup`).
   Result `s_r` → `data/features/region_signal.parquet`. Raw mean FC →
   `region_signal_raw.parquet`.
4. **Promoter feature.** Same as step 3 on the window `[tss0 − 1000, tss0 + 1000)`
   (`PROMOTER_HALF = 1000`) read straight from the bigWig (not from the union regions).
   Quantile normalisation uses the genes × lines promoter matrix itself. TSS = Ensembl-
   canonical transcript start, GENCODE v50 (`data/annot/genes_gencode_v50.tsv`, 54,710
   protein-coding + lncRNA genes).
5. **Distal feature.** `log2(1 + Σ_r 2^{s_r} · exp(−d_r / 25,000))` over union regions whose
   centre `(start + end)/2` lies at distance `d_r = |centre − tss0|` with
   2,500 ≤ d_r ≤ 100,000 bp (`DISTAL_MIN`, `DISTAL_MAX`, `DISTAL_SCALE`). Both sides of the
   TSS count. A gene with no such region gets 0.
6. **Combined.** `0.5·z(promoter) + 0.5·z(distal)`; each z uses a single mean and SD
   (ddof = 0) pooled over all genes × all configured lines, so between-line differences stay
   on one scale. The pooled means / SDs are stored in `build_info.json`.

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
| `region_signal.parquet` | region_id, chrom, start, end, one column per line = `s_r` (float32) |
| `region_signal_raw.parquet` | same, raw mean FC |
| `gene_features.parquet` | long: gene_id, gene_name, line, promoter, distal, combined, promoter_raw_fc, promoter_cn_log2, n_distal_regions |
| `build_info.json` | config, chosen files + md5, skipped lines, region / gene counts, parameters, pooled z constants |

Python access:

```python
from cdtest.features import pair_feature, feature_matrix
x = pair_feature("RPE1", "K562", "combined")   # pd.Series indexed by gene_id
P = feature_matrix("promoter")                  # genes x lines
```

## Sanity checks (results/features/sanity.md)

Current build (10 lines, Jurkat pending, RPE1 provisional): 338,483 union regions.
Promoter correlation between lines 0.66–0.87. Lineage markers behave as expected (HBG1 and
GATA1 higher in K562 than RPE1; ALB highest in HepG2, APOA1 second to MCF-7 but above K562
and RPE1). **Two pre-stated checks fail:** promoter feature vs baseline expression within
the 5,262 Aim 1 genes has Spearman −0.09 (K562) and −0.14 (RPE1), not > 0.3; GAPDH/ACTB are
above the median in every line but not ≥ 90th percentile. The same happens with the raw FC
track at any window size, so it is not a bug in the steps above. Promoter accessibility does
separate expressed from non-expressed genes (AUC 0.85–0.88). The likely cause is that MACS2
fold-change, with its local background, does not rank very strong promoters higher. Changing
the signal track would be a plan deviation.
