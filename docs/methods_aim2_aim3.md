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
   at lags 1, 2, …, 50. k is the first lag where it drops below 0.05, kept
   within [5, 50] (`choose_k`). M is never looked at.
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
bootstrap over knockdowns belongs to the pipeline, not this module.

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
