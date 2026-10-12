# Aim 3 planted-signal benchmark on real M, RPE1 − K562

Code: `cdtest.planted` (`scripts/run_planted_dev.py`). Plan §6: M' = M + θ·v, v = (random sign) × a negative-control decoy feature, θ set so the planted partial R² (raw, given f(C)) is the level; 500 replicates per level, B = 999, every test applied to v.

- Decoys: 21 decoy-line pairs (combined CPM feature). Rule: |similarity| < 0.05. 6 pass: GM12878-MCF-7, GM12878-A549, GM12878-Panc1, MCF-7-A549, MCF-7-Panc1, A549-Panc1.
- k per decoy (x4 / v1.0): GM12878-MCF-7 100/50, GM12878-A549 100/50, GM12878-Panc1 100/50, MCF-7-A549 100/50, MCF-7-Panc1 100/50, A549-Panc1 100/50.
- Only 12 distinct (decoy, sign) datasets exist; replicates within one differ only in permutation draws. MC SEs are therefore given two ways: binomial over 500 (`_se`) and cluster-robust over configurations (`_se_cluster`, the honest one).

## Rejection rate at nominal 5% (level 0 = false-positive rate); ± cluster SE

| test | R²=0.00% | R²=0.25% | R²=0.50% | R²=1.00% | R²=2.00% |
|---|---|---|---|---|---|
| block_x4 | 0.096 ± 0.069 | 0.574 ± 0.153 | 0.602 ± 0.146 | 0.928 ± 0.073 | 1.000 ± 0.000 |
| block_v1.0 | 0.084 ± 0.082 | 0.766 ± 0.127 | 1.000 ± 0.000 | 1.000 ± 0.000 | 1.000 ± 0.000 |
| naive_ols | 0.272 ± 0.140 | 0.930 ± 0.071 | 1.000 ± 0.000 | 1.000 ± 0.000 | 1.000 ± 0.000 |
| camera | 0.000 ± 0.000 | 0.000 ± 0.000 | 0.236 ± 0.126 | 0.560 ± 0.152 | 1.000 ± 0.000 |
| unrestricted | 0.256 ± 0.131 | 0.930 ± 0.071 | 1.000 ± 0.000 | 1.000 ± 0.000 | 1.000 ± 0.000 |
| expr_bin | 0.172 ± 0.108 | 0.930 ± 0.071 | 1.000 ± 0.000 | 1.000 ± 0.000 | 1.000 ± 0.000 |
| stratified | 0.182 ± 0.121 | 0.930 ± 0.071 | 1.000 ± 0.000 | 1.000 ± 0.000 | 1.000 ± 0.000 |

## FPR at 1% / rejection at nominal 1%

| test | R²=0.00% | R²=0.25% | R²=0.50% | R²=1.00% | R²=2.00% |
|---|---|---|---|---|---|
| block_x4 | 0.000 | 0.306 | 0.590 | 0.774 | 0.994 |
| block_v1.0 | 0.000 | 0.492 | 0.802 | 1.000 | 1.000 |
| naive_ols | 0.076 | 0.750 | 1.000 | 1.000 | 1.000 |
| camera | 0.000 | 0.000 | 0.000 | 0.062 | 0.588 |
| unrestricted | 0.076 | 0.746 | 0.990 | 1.000 | 1.000 |
| expr_bin | 0.076 | 0.746 | 0.974 | 1.000 | 1.000 |
| stratified | 0.076 | 0.746 | 0.940 | 1.000 | 1.000 |

## Power at matched 5% FPR

| test | R²=0.00% | R²=0.25% | R²=0.50% | R²=1.00% | R²=2.00% |
|---|---|---|---|---|---|
| block_x4 | 0.050 | 0.574 | 0.596 | 0.928 | 1.000 |
| block_v1.0 | 0.052 | 0.744 | 0.998 | 1.000 | 1.000 |
| naive_ols | 0.076 | 0.608 | 0.924 | 1.000 | 1.000 |
| camera | 0.076 | 0.608 | 0.924 | 1.000 | 1.000 |
| unrestricted | 0.054 | 0.686 | 0.942 | 1.000 | 1.000 |
| expr_bin | 0.058 | 0.706 | 0.940 | 1.000 | 1.000 |
| stratified | 0.052 | 0.638 | 0.924 | 1.000 | 1.000 |

## 95% CI coverage of θ

| test | R²=0.00% | R²=0.25% | R²=0.50% | R²=1.00% | R²=2.00% |
|---|---|---|---|---|---|
| block_x4 | 1.000 | 1.000 | 0.998 | 1.000 | 1.000 |
| block_v1.0 | 0.976 | 0.990 | 0.988 | 0.990 | 0.984 |
| naive_ols | 0.856 | 0.860 | 0.854 | 0.858 | 0.840 |
| camera | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| unrestricted | 0.856 | 0.858 | 0.850 | 0.858 | 0.840 |
| expr_bin | 0.856 | 0.860 | 0.854 | 0.858 | 0.840 |
| stratified | 0.854 | 0.858 | 0.852 | 0.858 | 0.838 |

## Bias of θ̂ (same estimator for every test)

| level | θ planted | bias | SE (cluster) |
|---|---|---|---|
| 0.00% | 0 | 2.64e-05 | 0.00205 |
| 0.25% | 0.01445 | -6.7e-05 | 0.00185 |
| 0.50% | 0.02068 | 0.000168 | 0.00193 |
| 1.00% | 0.02933 | -0.000458 | 0.00196 |
| 2.00% | 0.04137 | -0.000189 | 0.00198 |

## Level 0 by configuration (θ₀ = the decoy's own coefficient in the real M)

| decoy | sign | n | θ₀ | rej block_x4 | rej block_v1.0 | rej naive |
|---|---|---|---|---|---|---|
| A549-Panc1 | -1 | 31.0 | -0.0050 | 0.00 | 0.00 | 0.00 |
| A549-Panc1 | +1 | 38.0 | 0.0050 | 0.00 | 0.00 | 0.00 |
| GM12878-A549 | -1 | 42.0 | -0.0041 | 0.00 | 0.00 | 0.00 |
| GM12878-A549 | +1 | 33.0 | 0.0041 | 0.00 | 0.00 | 0.00 |
| GM12878-MCF-7 | -1 | 37.0 | 0.0010 | 0.00 | 0.00 | 0.00 |
| GM12878-MCF-7 | +1 | 47.0 | -0.0010 | 0.00 | 0.00 | 0.00 |
| GM12878-Panc1 | -1 | 45.0 | -0.0059 | 0.00 | 0.00 | 0.00 |
| GM12878-Panc1 | +1 | 41.0 | 0.0059 | 0.76 | 1.00 | 1.00 |
| MCF-7-A549 | -1 | 57.0 | -0.0072 | 0.00 | 0.00 | 0.00 |
| MCF-7-A549 | +1 | 57.0 | 0.0072 | 0.00 | 0.00 | 1.00 |
| MCF-7-Panc1 | -1 | 34.0 | -0.0122 | 0.00 | 0.00 | 0.00 |
| MCF-7-Panc1 | +1 | 38.0 | 0.0122 | 0.45 | 0.03 | 1.00 |

## Detection limit

- Main test (x4 rule): smallest planted partial R² with ≥ 80% power at nominal 5% = **1.00%**
  (at matched 5% FPR: 1.00%).
- v1.0 rule (sensitivity): 0.50%.

## Reading

- Level 0 is not a clean null: the decoy features are not guaranteed to be unrelated to the real M. Rejections at level 0 come from two configurations with θ₀ > 0 (GM12878−Panc1 +, MCF-7−Panc1 +); for the other ten the block tests never reject. Naive OLS and plain permutations reject whenever |θ₀| is a few thousandths; the block tests need a larger θ₀. With only 12 configurations the FPR estimates are very uncertain (cluster SE ≈ 0.07).
- The x4 main test has coverage 1.00 at every level (conservative CIs) and loses power relative to v1.0, as the synthetic benchmark predicted. On real M part of that loss comes from the immobile blocks (see results/aim2/RPE1_vs_K562_test.md: only ~41% of genes are in movable blocks at k = 100).
- CAMERA (VIF ≈ 17 from ρ̄ ≈ 0.003) is far too conservative here.
- Bias of θ̂ is negligible relative to θ at every level (random signs cancel θ₀).
