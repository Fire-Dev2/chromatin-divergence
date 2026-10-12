# Methods: the Aim 2 test and the Aim 3 synthetic benchmark

Code: `src/cdtest/crt.py` (test), `src/cdtest/simulate.py` (benchmark),
`src/cdtest/plot_aim3.py` (central figure). Tests: `tests/`. The specification
is ANALYSIS_PLAN.md sections 5 and 6. Where the code had to make a choice the
plan leaves open, it is marked **[choice]**.

---

## 1. What is being tested

For each gene g there are three things:

- **M_g**, the response difference between the two lines (Aim 1 output);
- **x_g**, the accessibility difference (combined feature, line A minus line B);
- **C_g**, covariates: log2 mean expression in each line, log2 variance in each
  line, expression difference, copy-number difference, promoter CpG o/e.

The effect model is

  M_g = θ·x_g + f(C_g) + ε_g,

where f is a natural cubic spline with 4 degrees of freedom on each line's log2
mean expression, plus linear terms for everything else. θ is the change in
response difference per unit of accessibility difference. H1 says θ > 0.

The null hypothesis is: **given C, x is independent of M.**

## 2. Step by step

### 2.1 The statistic

T is the ordinary least-squares coefficient on x when M is regressed on f(C)
and x together. It is computed with the Frisch–Waugh–Lovell identity: remove
from x and from M everything the covariates explain (r_x and r_M), then

  T = ⟨r_x, M⟩ / ⟨r_x, r_x⟩.

`crt.build_design` makes the covariate matrix. The spline basis
(`natural_spline_basis`) has the same span as R's `ns(x, df = 4)`: knots at the
25/50/75% quantiles, boundary knots at the range, linear outside it. The
expression difference is an exact linear combination of the two lines'
expression, so the design is rank-deficient on purpose. The code works with an
orthonormal basis Q of the design's column space (SVD), which handles that.

### 2.2 The main null: block permutation

We need many "replacement" features x* that look like x would if it had no link
to M, then we recompute T for each.

1. **Residualize x on C.** r_x = x − (part of x that C predicts).
2. **Order genes along the genome** by (chromosome, TSS position).
3. **Choose the block size k from the feature alone.** Compute the
   autocorrelation of r_x along the genome (pairs within the same arm only)
   at lags 1, 2, …, 50. Let k_v1.0 be the first lag where it drops below 0.05,
   kept within [5, 50]. M is never looked at.
   - **Plan v1.1 (primary): k = 4 × k_v1.0, kept within [5, 100]**
     (`choose_k`, default `rule="x4"`).
   - The v1.0 rule, k = k_v1.0 (`rule="v1.0"`), is reported as a sensitivity
     analysis.
   - Why the rule changed: §5.4 and `results/aim3/k_rule_selection.md`.
4. **Cut each chromosome arm into blocks** of k consecutive genes. If the arm
   leaves a shorter remainder, it joins the previous block, so blocks have
   between k and 2k−1 genes. An arm with fewer than k genes forms one block
   (`make_blocks`).
5. **Shuffle whole blocks.** Blocks of the same size swap places at random. A
   block keeps its internal order, so the local pattern of accessibility
   (neighbouring genes going up and down together) stays intact. A block whose
   size is unique cannot move.
6. **Add back the fitted part**: x* = (C-predicted part of x) + shuffled r_x.
   Because the fitted part lies in the span of the covariates, it does not
   change T*. The code keeps it in `block_permuted_features` for clarity but
   skips it in the fast path.
7. Recompute T* for each x*. Repeat B = 9,999 times.

**Speed.** With Q and r_M computed once, T* for a replacement feature v is
⟨v, r_M⟩ / (|v|² − |Q′v|²). A chunk of permutations is one matrix product. On
10,000 genes, 9,999 block permutations take about 0.6 s on the M5 Pro.

### 2.3 p-value

  p = (1 + #{b : T*_b ≥ T}) / (1 + B), one-sided for H1 (θ > 0).

There is also a two-sided version, (1 + #{|T*_b| ≥ |T|}) / (1 + B). The "+1"
counts the observed feature as one of the possible draws. That is what makes the
p-value valid rather than only approximately so (§3).

### 2.4 Confidence interval by test inversion

To test "θ = θ₀", subtract θ₀·x from M and run the same test on what is left.
The 95% interval is every θ₀ that is not rejected. The rule is equal-tailed:
each one-sided p-value must be > 0.025. The candidate θ₀ values are the plan's
grid of 201 points spanning ±max(4·|θ̂|, 0.5/SD(x)).

The statistic is linear in M, so the test at θ₀ needs no new permutations:
T(θ₀) = T − θ₀ and T*_b(θ₀) = T*_b − θ₀·c*_b, where
c*_b = ⟨v_b, r_x⟩ / (|v_b|² − |Q′v_b|²). The whole grid is evaluated in
closed form from the same B draws (`TestResult.p_at`, `.ci`, `.covers`). The
plan's B = 1,999 for CIs can therefore be raised to B = 9,999 at no extra cost.
The grid endpoints round inward by at most one grid step. The simulation scores
coverage exactly: it asks whether the true θ is rejected.

### 2.5 Partial R²

Partial R² = (RSS_C − RSS_{C,x}) / RSS_C = T²·|r_x|² / |r_M|². It is divided by
the Aim 1 reliability so it reads as a share of reproducible variance. Its
bootstrap over knockdowns is in `cdtest.aim2` (§6.3).

### 2.6 Sensitivity nulls

- **Stratified permutation.** Strata are 10 quantile bins of expression
  difference × 3 copy-number states (loss < −0.2, neutral, gain > 0.2). x is
  shuffled freely within each stratum.
- **Decoys.** Each decoy feature (from another pair of lines) is used as x* in
  turn, with p = (1 + #{T_decoy ≥ T}) / (1 + n_decoys). `decoy_similarity`
  computes the plan's similarity: the Pearson correlation of residualized decoy
  and residualized x. `max_abs_similarity=0.05` keeps the negative-control set.
  With 21 decoys the smallest possible p is 1/22 ≈ 0.045.

### 2.7 Comparators

1. **Naive OLS**: a t-test that treats genes as independent.
2. **Unrestricted permutation** of x across all genes.
3. **Expression-bin permutation**: x shuffled within 10 bins of expression
   difference only.
4. **CAMERA-style**: the OLS standard error multiplied by √VIF, with
   VIF = 1 + (n − 1)ρ̄. ρ̄ is the mean inter-gene correlation and n defaults to
   the number of genes **[choice]**. In real data ρ̄ comes from control cells; in
   the simulation it is the true mean off-diagonal correlation of the response
   noise.

For comparators 2–3, "x" is the raw feature, as the plan words it **[choice]**.

---

## 3. Why the randomization p-value is valid

### 3.1 Setting

Hold M and C fixed. Let X₀ = x be the real feature and X₁, …, X_B the
replacement features. Write T_b = T(M, C, X_b). The statistic is any fixed
function of (M, C, feature).

**Assumption (exchangeability under the null).** Under H0, conditional on
(M, C), the sequence (X₀, X₁, …, X_B) is exchangeable. That means its joint
distribution does not change when the B + 1 entries are reordered.

### 3.2 Lemma

Let T₀, …, T_B be exchangeable and N_i = #{j ≠ i : T_j ≥ T_i}. Then for every
α ∈ (0, 1),

  P( (1 + N₀)/(1 + B) ≤ α ) ≤ α.

*Proof.* Let m = ⌊α(1 + B)⌋, and let S be the set of indices i with
1 + N_i ≤ m. Take any outcome and assume S is not empty. Let i* be the index in
S with the smallest T. Every other member of S has T ≥ T_{i*}, so
N_{i*} ≥ |S| − 1. Then |S| ≤ 1 + N_{i*} ≤ m. So no outcome puts more than m
indices in S. Ties only make the N_i larger.

The vector (T₀, …, T_B) is exchangeable, and S is defined the same way for
every index, so each index is equally likely to be in S:

  P(0 ∈ S) = E|S| / (1 + B) ≤ m / (1 + B) ≤ α.  ∎

Applying the lemma conditionally on (M, C) and then averaging over (M, C) gives
P(p ≤ α) ≤ α. This holds at any number of genes, and whatever the joint
distribution of M: correlated, heteroscedastic, non-Gaussian. M is a constant
in the argument. If the T_b have no ties, the bound is attained exactly when
α(1 + B) is an integer.

### 3.3 When the assumption holds

**(a) Conditional randomization (Candès et al. 2018).** Draw X₁, …, X_B
independently from the true law of x given C, independently of M. Under H0,
x given (M, C) has that same law. So X₀, …, X_B are i.i.d. given (M, C), and
hence exchangeable.

**(b) Permutations from a group.** Let 𝒢 be a group of gene relabellings. For
block permutation, 𝒢 is all ways of exchanging blocks of equal size. Suppose
that given C, the law of x is invariant under 𝒢: πx has the same law as x for
every π ∈ 𝒢. Under H0 this invariance also holds given (M, C). Draw
π₁, …, π_B i.i.d. uniform on 𝒢, independently of the data, and set X_b = π_b x.

To see that (x, π₁x, …, π_Bx) is exchangeable, take an extra π₀ uniform on 𝒢
and set x′ = π₀x. Then x′ has the law of x and is independent of π₀. Write the
sequence as (σ₀x′, σ₁x′, …, σ_Bx′), with σ₀ = π₀⁻¹ and σ_b = π_bπ₀⁻¹. For a
group, the σ_b are i.i.d. uniform on 𝒢 and independent of x′ (Hemerik &
Goeman 2018). A sequence that is i.i.d. given x′ is exchangeable.

For block permutation, the assumption is that the residual feature is
exchangeable at block level within each size class. Dependence between genes may
sit inside blocks, but blocks of the same size are interchangeable. The sampled
permutations need not be distinct, and the identity may come up. Both are
allowed by the argument.

**(c) Decoys.** If the decoy features and x are exchangeable given (M, C) under
H0, the lemma applies with B = number of decoys. A decoy that shares x's
structure but not its tie to these two lines is the intended case.

**(d) Confidence interval.** Under the effect model with ε independent of x
given C, M − θ_true·x = f(C) + ε is independent of x given C. The null
therefore holds for the shifted response at the true θ, and that θ is rejected
with probability ≤ α. The interval covers θ_true with probability ≥ 1 − α.

### 3.4 What is exact and what is approximate

- With an intercept-only C, residualizing just centres x. Centring commutes with
  any permutation of genes, so block permutation is exact whenever x is
  block-exchangeable. `tests/test_crt_validity.py` checks this: 2,000
  replicates, false-positive rate within Monte Carlo error of 5% and 1%, while
  OLS on the same data inflates.
- With real covariates, the residual r_x = (I − H)x is only approximately
  block-exchangeable, because the projection mixes genes. The approximation
  error is of order rank(C)/G and was not detectable in the tests (G = 300,
  rank 6).
- The main assumption is that blocks are right: dependence in x does not cross
  block boundaries. Blocks are estimated, so this is only approximately true,
  and Aim 3 measures how much that costs (§5).

## 4. Why the feature is randomized and the response is held fixed

- **Validity then depends only on the feature's distribution.** The argument in
  §3 never uses the distribution of M. Correlation between genes in the
  response is the hard part to model. It comes from co-expression modules,
  shared pathways, knockdown-induced programs, and noise that depends on
  expression and cell count. Holding M fixed makes all of that irrelevant to
  validity.
- **Permuting M would need a model of M's dependence.** That is exactly what we
  do not trust. Treating genes as independent (OLS) or as exchangeable
  (unrestricted permutation) fails when response and feature share modules.
  Aim 3 shows this.
- **The feature's dependence is easier to see and to defend.** It is mostly
  genomic and local: neighbouring genes share regulatory elements and
  accessibility domains. It can be estimated from the feature alone, using
  ATAC from many lines. Choosing k without looking at M also means block choice
  cannot be tuned to the outcome.
- **The frameworks must not be mixed.** Grouping genes by co-expression of M
  and then permuting x (the earlier draft) preserves the wrong dependence. When
  x is randomized, it is x's dependence that must be kept.
- **Cost.** The test is only as good as the model for x given C. That is why
  stratified and decoy nulls run as checks, and why Aim 3 probes wrong blocks.

---

## 5. The Aim 3 synthetic benchmark

### 5.1 Generator (`simulate.simulate_replicate`)

- **Genome.** G = 2,000 genes on 44 arms (22 chromosomes × p/q), with random arm
  sizes **[choice]**. A real pair keeps about 5–10k genes. G = 2,000 keeps the
  1,000-replicate grid within about an hour on 14 cores.
- **Modules.** Contiguous runs of genes within an arm, size 1 + Geometric with
  mean 10 **[choice]**. They stand in for co-regulated neighbourhoods that the
  response and the feature can share. They are redrawn in every replicate.
- **Covariates.** e_A ~ N(2, 1.2²); e_B correlated 0.85 with e_A; log-variances
  ≈ 1.1·e plus noise; d = e_A − e_B; arm-level CN states (loss −1, neutral,
  gain +0.585, with probabilities 0.1/0.8/0.1).
- **Feature.** x = κ·z(d) + √(1−κ²)·[√ρ_x·u_module + √(1−ρ_x)·a_gene].
  - ρ_x is the feature's share of module structure.
  - κ ∈ {0, 0.3} is confounding with expression.
- **Response.** M = f(C) + θx + σ_g·[√ρ_M·v_module + √(1−ρ_M)·e_gene], with
  f(C) = 0.4·tanh(e_A−2) − 0.4·tanh(e_B−2) + 0.5·z(d) + 0.3·cn.
  - u and v are independent. The response and the feature share modules, not
    values, so the null is true at θ = 0.
- **Noise.** Homoscedastic (σ_g = 1), or σ_g² ∝ 0.3 + 100/(n_eff,g · 2^{ē_g}).
  Here n_eff is the per-gene effective cell count, LogNormal(log 100, 0.5), and
  ē_g is the gene's mean expression in the two lines. σ_g² is rescaled to
  mean 1.
- **Effect.** θ = √(R²/(1−R²)/(1−κ²)), giving partial R² = 0, 0.5%, 2%, 5%.
- **Decoys.** 99 independent redraws of x from the same model, so they are
  exchangeable with x under the null **[choice: the real study has ~21]**.

The grid is 4 ρ_M × 3 ρ_x × 2 κ × 2 noise × 4 θ = 192 settings. Each replicate
runs every test on the same data:

- naive OLS and CAMERA;
- unrestricted, expression-bin and stratified permutation;
- decoys;
- block permutation with the plan's k (main);
- three wrong-block variants: k/4, 4k, and blocks of k on a random gene order;
- **oracle** blocks equal to the true modules **[added]**, to separate "block
  permutation" from "the plan's rule for k".

Permutation tests use B = 999 **[choice]**. The seed is derived from the master
seed 20261011, the setting index and the replicate index.

### 5.2 Scores (`simulate.summarize`)

Each score is computed for every setting × test, with a binomial Monte Carlo SE:

- false-positive rate at 5% and 1% (θ = 0);
- rejection rate at 5% (θ > 0);
- power at matched false-positive rate: the threshold is the 5% quantile of that
  test's null p-values in the same setting;
- bias of θ̂ (the same estimator for all tests);
- exact coverage of the 95% interval: is θ_true not rejected by the
  equal-tailed test?

For the main test, the plan's 201-point grid interval is also stored.

### 5.3 Results (200 replicates per setting; see `results/aim3/README.md`)

Central figure: `figures/aim3_fpr_vs_structure.png`. It plots false-positive
rate at 5% against response module correlation, with columns for how much the
feature shares the modules and rows for confounding, one line per test.

- **Naive tests fail exactly where predicted.** OLS, unrestricted, expression-bin
  and stratified permutation hold 5% unless response and feature share modules.
  When both are structured they reach up to 25%.
- **The main block test is much closer to nominal.** It reaches 7–10% under
  strong shared structure, and its CI coverage averages 93%.
- **The excess is in the k rule, not in block permutation.** Oracle blocks and
  4k blocks hold 5%. The plan's autocorrelation < 0.05 rule makes blocks only
  about twice the dependence scale.
- **Decoys hold 5%.** By construction this is their best case.
- **CAMERA is conservative.** Its FPR is about 2%.

### 5.4 Choosing the block-size rule (plan v1.1; `cdtest.k_rule_sim`)

The grid in §5.3 showed that the v1.0 rule for k gives blocks that are too
small when the response and the feature share modules. Eight feature-only rules
were therefore compared in a targeted simulation. The protocol and the
selection criterion were committed before the run
(`results/aim3/k_rule_selection_protocol.md`).

**Candidate rules.** All are bounded to [5, 100]; v1.0 keeps its frozen [5, 50].

- v1.0;
- an autocorrelation threshold of 0.01 instead of 0.05;
- 2×, 3× and 4× k_v1.0;
- ⌈c·τ⌉ for c = 1, 2, where τ = 1 + 2·Σ ACF up to the first lag below 0.05;
- the Politis–White automatic block length (2004; Patton et al. 2009
  correction).

**Grid.**

- The shared-structure settings only: ρ_M {0.3, 0.5} × ρ_x {0.3, 0.6} ×
  confounding {0, 0.3} × 2 noise models.
- Four module designs **[added to `simulate.make_genome`]**: Geometric with
  mean 5, 10 and 25, and a mixed LogNormal design (median 8, σ = 1, clipped to
  [1, 150]). Modules are cut at arm ends.
- 1,000 null replicates per cell; 300 per cell at partial R² 0.5% and 2%.

**Pre-set criterion.**

1. Worst null cell ≤ 0.065 at α = 0.05 and ≤ 0.02 at α = 0.01.
2. Then the highest mean matched power at R² = 0.5%.
3. Ties within 0.02 go to the simpler rule.

**Outcome: no rule met the criterion.**

- FPR at α = 0.05, mean / worst cell:
  - v1.0: 7.9% / 11.5%;
  - 4×: 5.4% / 7.0%;
  - 3×: 5.6% / 7.3%;
  - Politis–White: 5.7% / 7.7%;
  - oracle (true modules): 5.0% / 6.5%.
- The τ rules were worse than v1.0.
- At n = 1,000 even an exact test passes the worst-of-64 limit only about 38%
  of the time.
- Under the protocol's fallback, the rule with the lowest worst-case FPR (4×)
  was adopted as plan v1.1.
- It costs power. Matched power at R² = 0.5% is 0.26 against 0.45 for v1.0.
  Part of this is due to the 45-gene simulated arms: few blocks can be
  exchanged once k approaches the arm length.
- Details, MC SEs and the k chosen per design: `results/aim3/k_rule_selection.md`.


---

## 6. The real-data pipeline (`cdtest.covariates`, `cdtest.aim2`, `cdtest.planted`)

Run on the development pair RPE1 − K562 only, before tag `pipeline-v1.0`.
Scripts: `scripts/run_aim1_dev.sh`, `scripts/run_aim2_dev.py`,
`scripts/run_planted_dev.py`, `scripts/report_dev.py`; figures
`cdtest.plot_dev`. Results: `results/aim2/RPE1_vs_K562_test.{json,md}`,
`results/aim3/planted_RPE1_vs_K562.{csv,md}`.

### 6.1 Covariates (plan §4; `cdtest.covariates`)

For the pair's kept genes, in the order of `gene_outcome.tsv`:

| column | definition |
|---|---|
| log2_meanA / log2_meanB | log2 mean normalized (CP-median) expression over the line's kept control cells (`ctrl_norm.npy`) |
| log2_varA / log2_varB | log2 variance (ddof = 1) of the same values **[choice: ddof 1]** |
| expr_diff | log2_meanA − log2_meanB |
| cn_diff | relative log2 copy number at the TSS, A − B, from `data/annot/gene_cn_log2.tsv`; RPE1 = column `RPE1` (karyotype track, +log2(3/2) on 10q21.2–qter), not `RPE1_depmap` |
| cpg_oe | `cpg_oe_tss500` (GENCODE v50 table) |

f(C) = intercept + `ns(log2_meanA, df 4)` + `ns(log2_meanB, df 4)` + the
other five linear: 14 columns, rank 13 (expr_diff is collinear with the
splines' linear parts; the orthonormal basis handles it). Kept genes all pass
the expression filter, so every variance is positive; non-finite covariates
raise an error rather than being imputed.

Genome order and blocks use `chrom`, `tss1` and `arm` (blocks never cross a
chromosome arm, labels `chrom + p/q`).

### 6.2 Test and CI

- Main test: `crt.block_permutation_test` with k = `choose_k(rule="x4")` on the
  x residuals in genome order within arms; B = 9,999. Sensitivity: the same with
  `rule="v1.0"`.
- **CI resolution [choice].** The plan's 201-point grid spans
  ±max(4|θ̂|, 0.5/SD(x)). On the development pair SD(x) ≈ 0.3 and θ̂ ≈ −0.02, so
  the grid step (≈ 0.017) is as large as θ̂ and the accepted set is one or two
  grid points. The plan-grid interval is reported as specified
  (`ci95_plan_grid`). Alongside it, `ci95` applies the identical acceptance rule
  (each one-sided p > 0.025, the same B = 9,999 draws) on 2,001 points between one
  plan step below and above the plan-grid interval (inside the plan's span).
  No parameter of the plan is changed; only the grid resolution.
- Plan B = 1,999 per CI grid point: the closed-form inversion uses the same
  9,999 draws for every θ₀ (§2.4).
- **Exchangeability diagnostic.** `aim2.block_exchangeability` reports how many
  blocks share their size with at least one other block (only those can move,
  plan §5 step 5) and the share of genes in them. On the development pair at
  k = 100 only 21 of 49 blocks (41% of genes) are movable; the rest keep their
  real x in every permutation, so the null distribution is centred near the
  observed association instead of 0. The test stays valid but loses power. This
  is a consequence of the frozen rules (x4, cap 100, equal-size classes) on real
  gene density; it is reported, not changed.

### 6.3 Partial R² and its bootstrap (approximation)

Partial R² = T²|r_x|²/|r_M|², divided by the Aim 1 reliability. The plan asks
for 500 bootstrap resamples over knockdowns, recomputing reliability in each.
The full version would rerun the generic-mean subtraction and the 50 fake
screens per resample (≈ 4 min each). **Approximation used:**

1. `aim1.run_pair(..., save_half_terms=...)` stores, for each of the 20 random
   half-splits and both halves, the per-(knockdown, gene) term
   D = (β̃²_A − s̃²_A) − (β̃²_B − s̃²_B) (`data/aim1_cache/<pair>_halfD.npy`).
   The full-data terms are in `terms.npz`.
2. Each resample draws multinomial weights over the kept knockdowns. M_b is the
   weighted mean of the full-data D (on-target entries excluded); each half's M is
   the weighted mean of its stored D; reliability_b = Spearman–Brown of the median
   split r over the 20 splits.
3. partial R²_b / reliability_b gives the percentile 95% interval.

β̃ and s̃² are held at their full-data values: the generic mean and the
fake-screen noise are not recomputed per resample. Unweighted, step 2
reproduces the stored split r's exactly (checked). The interval reflects
knockdown sampling only, not gene-level dependence. Bootstrap reliabilities are
lower than the full-data value (duplicated knockdowns), which biases the
share slightly upward.

### 6.4 Comparators on real data

- Naive OLS, unrestricted, expression-bin (10 bins of expr_diff), stratified
  (10 bins × CN loss/neutral/gain at ±0.2): as in §2.6–2.7, on raw x.
- CAMERA: ρ̄ = mean off-diagonal Pearson correlation between the pair's genes
  over control cells, computed per line on the normalized counts (no log, the
  scale on which pseudobulk means are formed) and averaged over the two lines
  **[choice]**; the log2(x+1) version is reported too. n = number of genes.
- Decoy negative control: decoy-line pairs (7 decoys, 21 pairs, first minus second
  in `aim2.DECOY_LINES` order) with |similarity| < 0.05.
- Co-accessibility clusters (secondary blocks): genes clustered on
  [z(promoter), z(distal)] across the 11 lines (pooled z as in the combined
  feature) **[choice]**, average linkage, 1 − Pearson, cut at 0.5.

### 6.5 Naive analysis (Aim 1 step 7)

OLS of M on 1 + expr_diff + x, genes independent: coefficient, SE, partial R²
and t-test p, for promoter, distal, combined (CPM) and combined (FC).
`results/aim1/naive_<pair>.json`.

### 6.6 Planted signal on real M (plan §6; `cdtest.planted`)

- Decoys: the 21 decoy-pair combined-CPM differences; negative-control set =
  |similarity| < 0.05 (similarity per plan §5). If fewer than 5 pass, the 5 with
  the lowest |similarity| are used and this is stated (not needed on the
  development pair: 6 pass).
- Replicate: random decoy from the set, random sign s; v = s·x_decoy is the
  tested feature. θ ≥ 0 is set so that
  θ²|r_v|² / (θ²|r_v|² + RSS(M ~ f(C) + v)) equals the level (raw partial R²,
  not divided by reliability; reliability 0.996 makes the difference
  negligible) **[choice]**. M' = M + θ·v. k is chosen from v's residuals.
- Every test is applied to v with B = 999: block x4, block v1.0, naive OLS,
  CAMERA, unrestricted, expression-bin, stratified. Coverage is exact inversion
  at the true θ (`TestResult.covers`).
- The truth assumed is θ_true = θ planted. The decoy may itself be associated with
  the real M (θ₀ ≠ 0); the random sign cancels θ₀ on average but not per
  replicate, so level 0 is not a clean null. θ₀ is stored per replicate.
- Only 2 × (number of decoys) distinct datasets exist; replicates within one
  differ only in permutation draws. MC SEs are reported binomial (as the plan)
  and cluster-robust over (decoy, sign) configurations.
- Detection limit: the smallest level with rejection rate ≥ 0.8 at nominal 5%
  for the main test (also reported at matched FPR and for v1.0).
