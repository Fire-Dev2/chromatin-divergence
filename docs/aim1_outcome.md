# Aim 1: the response-difference outcome M, step by step

Plain-language description of what the Aim 1 code does (`src/cdtest/screens.py`,
`src/cdtest/aim1.py`), with the numbers from the development pair RPE1 − K562
(A = RPE1, B = K562). The specification is ANALYSIS_PLAN.md §1–2.

## Step 0. One streaming pass over each screen (`screens.run`)

The Perturb-seq files are dense matrices of raw UMI counts (one row per cell,
one column per gene). They are too large to load twice, so each is read once in
blocks of 4,000 cells. For every cell:

- **Cell filter.** Keep the cell if its mitochondrial fraction is below the
  99th percentile of that line's control (non-targeting) cells.
- **Normalization.** Scale the cell so that its total UMI count equals the
  median total of the line's kept control cells (CP-median). No log.
- **Aggregation.** Add the normalized counts into one running sum per knockdown
  (all cells carrying guides against the same gene) and one for the controls.
  The normalized control cells are also written out individually
  (`ctrl_norm.npy`), because the noise null in step 2 needs single cells.
- **Random halves.** Each cell gets 20 independent coin flips. For each flip the
  pass also keeps the sum over "heads" cells per knockdown; the "tails" sum is
  the total minus that. These give the reliability in step 5 without rereading.

RPE1: 242,296 kept cells, 11,370 controls, 2,393 knockdowns, 8,749 genes.
K562 essential: 302,187 kept cells, 10,584 controls, 2,057 knockdowns, 8,563 genes.

For the K562 genome-wide screen (same-line baseline) the pass aggregates only
the knockdowns that also exist in the K562 essential screen, so the stored
sums stay small (`screens.RESTRICT_TO`).

## Step 1. Which genes and knockdowns are used (`aim1.select`)

- **Genes.** Measured in both screens, mean normalized expression ≥ 0.2 UMI per
  control cell and detected in ≥ 5% of control cells in both lines, and a
  GENCODE v50 protein-coding or lncRNA gene. RPE1 − K562: **5,262 genes**.
- **Knockdowns.** Present in both screens with ≥ 20 cells in each, target gene
  among the kept genes, and the target repressed by ≥ 30% (log2 fold change ≤
  log2 0.7) in both lines. RPE1 − K562: **1,314 of 2,055** shared knockdowns
  kept (210 too few cells, 510 target not expressed, 360 weak repression; a
  knockdown can fail more than one).
- A knockdown's own target gene is masked out of every average below.

## Step 2. Effect of each knockdown on each gene

β = log2(knockdown mean + 0.1) − log2(control mean + 0.1), from the pseudobulk
means of normalized counts.

## Step 3. Remove the generic response

Many knockdowns of essential genes move the same genes the same way (growth
arrest, stress). That shared part is removed: for each gene, subtract its mean
β over all kept knockdowns, giving β̃. What remains is the part of the response
specific to each knockdown.

**Generic share** = (mean β)² / mean β² per gene, i.e. how much of a gene's
response is the shared part; the median over genes is reported. RPE1 − K562:
13% in RPE1, 3% in K562.

## Step 4. Noise null from fake screens

Even with no biology, β̃² is positive because of sampling noise, and noise
differs between lines (cell numbers, depth). To estimate it, the code builds 50
fake screens per line: for every real knockdown with n cells it draws n control
cells at random as a "fake knockdown" (the remaining controls are the fake
control), computes β with the same formula and removes the generic mean the same
way. s̃² is the average of the fake β̃² over the 50 fake screens, per knockdown
and gene.

## Step 5. The outcome M

For each gene,

M = mean over knockdowns of [(β̃²_RPE1 − s̃²_RPE1) − (β̃²_K562 − s̃²_K562)].

Positive M: the gene's knockdown-specific responses are larger in RPE1 than in
K562 after removing noise. RPE1 − K562: mean M 0.048, SD 0.100
(`results/aim1/RPE1_vs_K562/gene_outcome.tsv`; figure
`figures/aim1_M_distribution.png`). The unsigned secondary outcome
(β̃_A − β̃_B)² − s̃²_A − s̃²_B is stored in the same table.

## Step 6. Reliability

Repeat steps 2–5 separately on the two random halves of the cells (within each
knockdown and within the controls; 10 fake screens per half), correlate the two
half-M's over genes, and apply Spearman–Brown 2r/(1+r). Median over 20 splits.
RPE1 − K562: split-half r 0.991, **reliability 0.996**. A quadratic in both
lines' baseline expression explains only 2.9% of M, so the agreement is not an
expression artifact.

## Step 7. p53 gate

The 116 Fischer 2017 core p53 targets (25 pass the gene filter) should respond
more in p53-wild-type RPE1 than in p53-null K562. Their mean M (0.131) is
compared with 1,000 random gene sets matched on expression deciles in both lines
(0.050 ± 0.020): z = 4.1, p = 0.001. **Gate passed.** Top targets: MDM2,
GADD45A, ZMAT3, RPS27L, TRIAP1.

## Step 8. Group-mean variant (sensitivity; plan §2 item 6)

Step 3 is replaced by subtracting the mean within preset knockdown groups
(cytosolic ribosome, mitochondrial ribosome, proteasome, spliceosome, RNA
polymerase II, Mediator, nuclear pore, other; `data/annot/knockdown_groups.tsv`,
HGNC groups with CORUM tiebreak), in the real screen and in every fake screen.
Run: `aim1.run_pair("RPE1", "K562", generic="group")` →
`results/aim1/RPE1_vs_K562_groupmean/`.

Kept knockdowns per group: other 1,038, spliceosome 84, cytosolic ribosome 72,
mitochondrial ribosome 61, proteasome 32, nuclear pore 16, Mediator 7, RNA pol II 4.
M correlates 0.996 with the global-mean M; mean 0.041, SD 0.085; reliability
0.994; p53 gate still passes (mean M of targets 0.100, z = 3.4, p = 0.002).
Small groups (RNA pol II, 4 knockdowns) have their generic mean estimated from
very few knockdowns; that is a property of the preset groups.

## Step 9. Naive analysis (plan Aim 1 step 7)

OLS of M on the expression difference and the accessibility difference x,
treating genes as independent (`aim2.naive_analysis`,
`results/aim1/naive_RPE1_vs_K562.json`). It is the "what a standard analysis
would say" reference. RPE1 − K562 (CPM features): combined −0.019 (partial R²
0.31%, two-sided p 6e−5), promoter −0.013 (0.61%, t = −5.7, p ≈ 1e−8), distal
−0.002 (0.006%, p 0.58); combined FC −0.001 (p 0.81). The sign is negative; for a mixed pair this
is descriptive only.

## Step 10. Same-line baseline

The same pipeline on two screens of one line, K562 essential (A) vs K562
genome-wide (B), restricted to shared knockdowns and genes
(`scripts/run_same_line.py`, `results/aim1/same_line_baseline.json`). This
shows how large M, its reliability and the generic share are when the cell line
does not change and only the experiment does.

The genome-wide file was checked before use: same obs fields as the essential
screen (`gene`, `gene_id`, `gem_group`, `UMI_count`, `mitopercent`, `sgID_AB`, …),
dense float32 `X` of 1,989,578 × 8,248, 75,328 non-targeting cells, 9,866
targets. The pass kept the 2,053 targets shared with the essential screen
(491,053 kept cells incl. 74,574 controls).

| | same line (K562 ess − K562 gw) | cross line (RPE1 − K562) |
|---|---|---|
| genes | 6,207 | 5,262 |
| knockdowns kept / shared | 1,508 / 2,053 | 1,314 / 2,055 |
| mean M | 0.0024 | 0.048 |
| SD of M | 0.0095 | 0.100 |
| split-half r (median of 20) | 0.567 | 0.991 |
| reliability | 0.72 | 0.996 |
| generic share, median (A / B) | 2.5% / 2.9% | 13.0% / 2.7% |

- The same-line M is about 10 times smaller in SD; its reproducible variance
  (reliability × variance) is 0.65% of the cross-line one.
- It is still reproducible (reliability 0.72): the two K562 screens differ
  systematically (separate experiments; median control UMI 13,920 vs 11,603), so
  not all of a cross-line M is "line" — but at this pair the batch part is ~1% of it.
- Over the 5,044 genes in both analyses the same-line and cross-line M are
  uncorrelated (r = −0.03).
