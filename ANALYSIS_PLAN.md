# Frozen analysis plan, v1.0

Fixed on 2026-10-11, before any response difference is computed. Every choice
below is made without looking at the Perturb-seq response data. A change after
this point needs a dated entry in log/RESEARCH_LOG.md with the reason, and both
the original and the changed analysis are reported.

Git tag: `plan-v1.0`.

---

## 0. Pairs and orientation

The response difference M is always line A minus line B.

| Pair | A | B | Role |
|---|---|---|---|
| K562–RPE1 | RPE1 | K562 | Pipeline development (Aim 1). Expected p53 contrast |
| K562–Jurkat | K562 | Jurkat | **Matched, primary (H1)** |
| RPE1–HepG2 | RPE1 | HepG2 | **Matched, primary (H1)** |
| K562–HepG2, Jurkat–RPE1, Jurkat–HepG2 | first | second | Mixed, pooled contrast |
| K562 essential – K562 genome-wide | essential | genome-wide | Same-line baseline |

Matched pairs are not touched until the pipeline is frozen on K562–RPE1
(tag `pipeline-v1.0`). Before that, the Jurkat and HepG2 files are opened only
to count cells, controls and knockdowns.

## 1. Cells, knockdowns and genes

- **Cells.** Use the provided raw UMI counts. Keep cells with a single assigned
  guide target (`gene` field) and with mitochondrial fraction below the line's
  99th percentile among control cells.
- **Controls.** Cells labelled non-targeting.
- **Normalization.** Scale each cell to the median total UMI count of that
  line's control cells (CP-median). No log transform before averaging.
- **Knockdowns kept.** A knockdown is kept for a pair if, in both lines:
  1. it has at least 20 cells;
  2. its target gene passes the expression threshold below; and
  3. on-target repression is at least 30%, meaning the target's pseudobulk
     log2 fold change is at most log2(0.7).
  The analysis uses only knockdowns kept in both lines of the pair.
- **Genes.** A gene is kept for a pair if, in both lines' control cells, mean
  normalized expression is at least 0.2 UMI per cell and it is detected in at
  least 5% of cells. Genes must also have a GENCODE v50 protein-coding or
  lncRNA annotation with a TSS on chr1–22 or chrX.
- **On-target entries.** The entry for knockdown p and its own target gene is
  excluded from every average that defines M.

## 2. Outcome (Aim 1)

1. **Effect.** For knockdown p, gene g and line k:
   β̂ = log2(mean_p + ε) − log2(mean_ctrl + ε), using pseudobulk means of
   normalized counts and ε = 0.1.
2. **Noise null.** Build K = 50 fake screens per line. In each fake screen,
   draw n_p control cells without replacement as a pseudo-knockdown for every
   real knockdown p, matching its cell count; the remaining controls form the
   pseudo-control. Compute the fake β with the same formula.
3. **Generic response.** For each gene, subtract its mean β over the kept
   knockdowns, giving β̃. The same subtraction is applied within each fake
   screen. s̃²(p,g,k) is the mean over the K fake screens of β̃_fake².
4. **M_g.** M_g = (1/P) Σ_p [(β̃²_A − s̃²_A) − (β̃²_B − s̃²_B)].
   - Unsigned secondary outcome: (1/P) Σ_p [(β̃_A − β̃_B)² − s̃²_A − s̃²_B].
   - Generic share = (mean_p β)² divided by mean_p β², per gene, summarized as
     the median over genes. It is reported for every pair and the same-line
     baseline.
5. **Reliability.** Split cells at random into halves, separately within each
   knockdown and within controls, and recompute M on each half. Compute the
   Pearson r between the halves; reliability = 2r/(1+r). Repeat over 20 splits
   and report the median.
6. **Group-mean variant (sensitivity).** Subtract the mean within preset
   knockdown groups instead of the global mean. Groups are the cytosolic
   ribosome, mitochondrial ribosome, proteasome, spliceosome, RNA polymerase II,
   Mediator, nuclear pore, and "other". Membership comes from HGNC gene groups,
   with CORUM 5.3 as tiebreaker.

**p53 gate.** Fischer 2017 116-gene core targets, restricted to kept genes.
Statistic: the mean of M over these targets for K562–RPE1 (A = RPE1). It is
compared against 1,000 random gene sets of the same size, matched on mean
expression deciles in both lines. The gate passes if the empirical one-sided p
is below 0.05. If it fails, the pipeline is debugged before anything else runs.
The two matched pairs serve as the negative arm of this check, but only in Aim 4.

## 3. Accessibility features

- **Region set.** Take the union of IDR peaks over all 11 lines: K562, RPE1,
  Jurkat, HepG2 and 7 decoys. Use summit ±250 bp, merge overlaps, and drop
  ENCODE blacklist v2 regions.
- **Region signal.** Mean fold-change bigWig signal over each region, per line.
  Then log2(x + 1) and quantile normalization across the 11 lines.
- **Copy-number correction.** Subtract the line's relative log2 copy number at
  each region. Sources: DepMap 24Q4 segments for K562, HepG2 and Jurkat; the
  karyotype-based track for RPE1 (+log2(3/2) on 10q21.2–qter); DepMap where
  available for decoys, otherwise none, which is noted.
- **TSS.** Ensembl-canonical transcript TSS (GENCODE v50).
- **Promoter feature.** Mean normalized, CN-corrected signal over TSS ±1 kb,
  computed directly from the bigWig on the same log2 / quantile / CN scheme.
  Quantile normalization uses the promoter-window matrix itself.
- **Distal feature.** log2(1 + Σ_r 2^{s_r} · exp(−d_r / 25 kb)) over union
  regions whose centre lies 2.5–100 kb from the TSS. Here s_r is the
  normalized, CN-corrected region signal and d_r is the distance to the TSS.
- **Combined feature (H1).** 0.5·z(promoter) + 0.5·z(distal). The z-scores use
  the mean and SD pooled over all genes and all 11 lines, so differences
  between lines stay on a common scale.
- **Feature difference.** x_g = feature_A − feature_B.

## 4. Covariates C (entered before x)

In both lines: log2 mean control expression, and log2 control variance of
normalized expression. Also the expression difference (log2 A − log2 B), the
copy-number difference at the TSS, and promoter CpG observed/expected ratio
over TSS ±500 bp (GRCh38 primary assembly).

f(C): a natural cubic spline with 4 df on each line's log2 mean expression;
all other covariates enter linearly.

## 5. Test (Aim 2)

- **Statistic T.** The OLS coefficient on x in M ~ f(C) + x. H1 is one-sided:
  θ > 0.
- **Main null: block permutation.**
  1. Residualize x on C, then order genes by genomic position.
  2. Cut each chromosome arm into consecutive blocks of k genes; a leftover
     shorter than k joins the previous block.
  3. Permute the residuals at block level among blocks of equal size, keeping
     order within each block, and add back the fitted part.
  4. k is chosen from the feature only: the smallest gene lag at which the
     autocorrelation of the x residuals falls below 0.05, bounded to [5, 50].
  5. Blocks of unequal size are permuted only within their size class.
- **Sensitivity 1: stratified permutation.** Shuffle x within cells of a grid:
  10 bins of expression difference × 3 copy-number-difference states
  (loss < −0.2, neutral, gain > 0.2 in relative log2).
- **Sensitivity 2: decoys.** Each decoy x comes from a pair of decoy lines
  (21 pairs). Similarity is the Pearson correlation between the decoy x and the
  pair's x, both residualized on C.
  - Negative control: decoys with |similarity| < 0.05.
  - H4: all decoys plus the other pairs' x.
- **Secondary block definition.** Co-accessibility clusters: genes are
  clustered on promoter + distal signal across the 11 lines (hierarchical,
  average linkage, 1 − Pearson, cut at 0.5). This is reported, not primary.
- **Permutation counts.** B = 9,999 for primary p-values and B = 1,999 per grid
  point for confidence intervals.
- **p-value.** p = (1 + #{T* ≥ T}) / (1 + B).
- **Confidence interval.** Test inversion over a θ grid of 201 points spanning
  ±4 × |naive OLS θ| (at least ±0.5 in M units per SD of x). This gives the
  95% two-sided interval.
- **Partial R².** For x, as a share of reproducible variance (divided by
  reliability). Its 95% interval comes from 500 bootstrap resamples over
  knockdowns, recomputing reliability in each.
- **Comparators.** Four standard tests are run on the same data:
  1. OLS treating genes as independent;
  2. unrestricted permutation of x;
  3. permutation of x within expression-difference bins only;
  4. a CAMERA-style variance-inflation correction, using inter-gene
     correlation estimated from control cells (VIF = 1 + (n − 1)ρ̄).
- **Multiple testing.** H1 on the combined feature in the two matched pairs,
  with Holm correction at family-wise α = 0.05. Everything else is labelled
  secondary.

## 6. Benchmarks (Aim 3)

- **Synthetic grid.** The settings are crossed as follows:
  - response module correlation ρ_M ∈ {0, 0.1, 0.3, 0.5};
  - feature sharing of modules ρ_x ∈ {0, 0.3, 0.6};
  - feature–expression confounding ∈ {0, 0.3};
  - noise model ∈ {homoscedastic, expression- and cell-count-dependent};
  - effect θ ∈ {0, small, moderate, large}, where small, moderate and large
    give partial R² of 0.5%, 2% and 5%;
  - wrong-block settings: block size k/4, 4k, and random gene order.

  There are 1,000 replicates per setting.
- **Planted signal.** Use real M from each pair. Add θ·x_decoy with θ chosen to
  give partial R² ∈ {0, 0.25%, 0.5%, 1%, 2%}, drawing the decoy from the
  negative-control set. Run 500 replicates per level, each with a different
  decoy and random sign.
- **Scores.** False-positive rate at 5% and 1%, power at matched false-positive
  rate, bias, and 95% CI coverage, each with Monte Carlo SE.
- **Detection limit.** The smallest planted partial R² recovered with 80%
  power by the main test, for each matched pair.

## 7. Aim 4 extras

- **Mixed pairs.** The mean of the four T statistics, with one shared block
  permutation applied to all four at once.
- **H4.** Regress each pair's θ̂ for every available feature on that feature's
  similarity to the pair's own x. Report the slope, the intercept and the plot.
- **H5.** For TFs in the kept knockdown set with ENCODE IDR ChIP in K562 or
  HepG2, or a JASPAR 2024 CORE motif (FIMO, p < 1e-5, inside union regions):
  1. Count accessible binding sites within 50 kb of each gene, where accessible
     means the region signal is above the line's median.
  2. Take the difference between lines as the feature.
  3. Compute the Spearman correlation with β̃²_A − β̃²_B for that knockdown.
  4. Swap control: repeat with 20 random other TFs' sites.

## 8. Data choices fixed here

- **Jurkat ATAC.** GSE273644 is primary and GSE302221 is the sensitivity set.
  QC gate and fallbacks are as logged on 2026-10-11.
- **RPE1 ATAC.** The 4DN FC bigWig and conservative peaks are used
  provisionally. They will be replaced by ENCODE v1.9.1 reprocessing of the
  same FASTQs before any matched-pair analysis.
- **Random seeds.** Master seed 20261011; each job derives its seed from it.
