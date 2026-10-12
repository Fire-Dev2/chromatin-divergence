# Aim 3 planted-signal benchmark on real M, RPE1 − K562 — plan v1.3 (arm fixed effects in C)

Code: `cdtest.planted.run(arm_effects=True)` (`scripts/run_planted_dev_v13.py`). Protocol as v1.2 (plan §6): M' = M + θ·v, v = random sign × a negative-control decoy, θ set so the raw planted partial R² given f(C) is the level; 500 replicates/level, B = 999; every test applied to v. Under v1.3, f(C) includes the arm indicators everywhere: decoy similarities (negative-control set), k per decoy, θ calibration and all tests. v1.2 results: `planted_RPE1_vs_K562.{csv,md}` (unchanged).

- Decoys (|similarity| < 0.05 recomputed under v1.3 C): 7 of 21 — GM12878-HCT116, GM12878-MCF-7, HCT116-MCF-7, MCF-7-A549, MCF-7-Panc1, A549-Panc1, IMR-90-PC-3 (v1.2: 6). 14 (decoy, sign) configurations.
- k per decoy (x4 / v1.0): GM12878-HCT116 56/14, GM12878-MCF-7 100/27, HCT116-MCF-7 100/29, MCF-7-A549 100/29, MCF-7-Panc1 100/29, A549-Panc1 72/18, IMR-90-PC-3 44/11.
- MC SEs: binomial over 500 and cluster-robust over (decoy, sign) configurations (± below, the honest one).

## Each decoy's own θ₀ in the real M (v1.3 f(C); two-sided)

| decoy | θ₀ | k x4 | p block x4 (B = 9,999) | p naive OLS |
|---|---|---|---|---|
| GM12878-HCT116 | 0.0002 | 56 | 0.9624 | 0.961 |
| GM12878-MCF-7 | -0.0031 | 100 | 0.7091 | 0.402 |
| HCT116-MCF-7 | -0.0071 | 100 | 0.921 | 0.195 |
| MCF-7-A549 | 0.0054 | 100 | 0.5041 | 0.272 |
| MCF-7-Panc1 | 0.0164 | 100 | 0.0976 | 0.000734 |
| A549-Panc1 | 0.0161 | 72 | 0.018 | 0.006 |
| IMR-90-PC-3 | -0.0182 | 44 | 0.0186 | 0.000613 |

Significantly nonzero θ₀ by the main test (two-sided p ≤ 0.05): A549-Panc1, IMR-90-PC-3; by naive OLS: A549-Panc1, IMR-90-PC-3, MCF-7-Panc1. Exclusion drops both signs of such a decoy.

## Rejection at nominal 5% (level 0 = FPR) ± cluster SE — all configurations

| test | R²=0.00% | R²=0.25% | R²=0.50% | R²=1.00% | R²=2.00% |
|---|---|---|---|---|---|
| block_x4 | 0.122 ± 0.085 | 0.652 ± 0.132 | 0.782 ± 0.113 | 1.000 ± 0.000 | 1.000 ± 0.000 |
| block_v1.0 | 0.226 ± 0.122 | 0.792 ± 0.113 | 1.000 ± 0.000 | 1.000 ± 0.000 | 1.000 ± 0.000 |
| naive_ols | 0.226 ± 0.122 | 0.792 ± 0.113 | 1.000 ± 0.000 | 1.000 ± 0.000 | 1.000 ± 0.000 |
| camera | 0.000 ± 0.000 | 0.126 ± 0.088 | 0.210 ± 0.113 | 0.556 ± 0.140 | 1.000 ± 0.000 |
| unrestricted | 0.226 ± 0.122 | 0.792 ± 0.113 | 1.000 ± 0.000 | 1.000 ± 0.000 | 1.000 ± 0.000 |
| expr_bin | 0.228 ± 0.122 | 0.792 ± 0.113 | 1.000 ± 0.000 | 1.000 ± 0.000 | 1.000 ± 0.000 |
| stratified | 0.238 ± 0.121 | 0.792 ± 0.113 | 0.994 ± 0.006 | 1.000 ± 0.000 | 1.000 ± 0.000 |

## Rejection at nominal 5% ± cluster SE — excluding decoys with significant θ₀ (main test: A549-Panc1, IMR-90-PC-3)

| test | R²=0.00% | R²=0.25% | R²=0.50% | R²=1.00% | R²=2.00% |
|---|---|---|---|---|---|
| block_x4 | 0.000 ± 0.000 | 0.720 ± 0.147 | 0.794 ± 0.137 | 1.000 ± 0.000 | 1.000 ± 0.000 |
| block_v1.0 | 0.136 ± 0.131 | 0.914 ± 0.087 | 1.000 ± 0.000 | 1.000 ± 0.000 | 1.000 ± 0.000 |
| naive_ols | 0.136 ± 0.131 | 0.914 ± 0.087 | 1.000 ± 0.000 | 1.000 ± 0.000 | 1.000 ± 0.000 |
| camera | 0.000 ± 0.000 | 0.086 ± 0.087 | 0.083 ± 0.085 | 0.567 ± 0.168 | 1.000 ± 0.000 |
| unrestricted | 0.136 ± 0.131 | 0.914 ± 0.087 | 1.000 ± 0.000 | 1.000 ± 0.000 | 1.000 ± 0.000 |
| expr_bin | 0.138 ± 0.130 | 0.914 ± 0.087 | 1.000 ± 0.000 | 1.000 ± 0.000 | 1.000 ± 0.000 |
| stratified | 0.151 ± 0.129 | 0.914 ± 0.087 | 0.991 ± 0.009 | 1.000 ± 0.000 | 1.000 ± 0.000 |

## Same, excluding by the naive-OLS θ₀ test (A549-Panc1, IMR-90-PC-3, MCF-7-Panc1)

| test | R²=0.00% | R²=0.25% | R²=0.50% | R²=1.00% | R²=2.00% |
|---|---|---|---|---|---|
| block_x4 | 0.000 ± 0.000 | 0.766 ± 0.156 | 0.753 ± 0.164 | 1.000 ± 0.000 | 1.000 ± 0.000 |
| block_v1.0 | 0.000 ± 0.000 | 1.000 ± 0.000 | 1.000 ± 0.000 | 1.000 ± 0.000 | 1.000 ± 0.000 |
| naive_ols | 0.000 ± 0.000 | 1.000 ± 0.000 | 1.000 ± 0.000 | 1.000 ± 0.000 | 1.000 ± 0.000 |
| camera | 0.000 ± 0.000 | 0.000 ± 0.000 | 0.000 ± 0.000 | 0.609 ± 0.187 | 1.000 ± 0.000 |
| unrestricted | 0.000 ± 0.000 | 1.000 ± 0.000 | 1.000 ± 0.000 | 1.000 ± 0.000 | 1.000 ± 0.000 |
| expr_bin | 0.003 ± 0.004 | 1.000 ± 0.000 | 1.000 ± 0.000 | 1.000 ± 0.000 | 1.000 ± 0.000 |
| stratified | 0.021 ± 0.021 | 1.000 ± 0.000 | 1.000 ± 0.000 | 1.000 ± 0.000 | 1.000 ± 0.000 |

## Main test vs v1.2 (rejection at nominal 5%, all configurations)

| test | plan | R²=0.00% | R²=0.25% | R²=0.50% | R²=1.00% | R²=2.00% |
|---|---|---|---|---|---|---|
| block_x4 | v1.3 | 0.122 | 0.652 | 0.782 | 1.000 | 1.000 |
| block_x4 | v1.2 | 0.096 | 0.574 | 0.602 | 0.928 | 1.000 |
| block_v1.0 | v1.3 | 0.226 | 0.792 | 1.000 | 1.000 | 1.000 |
| block_v1.0 | v1.2 | 0.084 | 0.766 | 1.000 | 1.000 | 1.000 |

## Rejection at nominal 1% (all configurations)

| test | R²=0.00% | R²=0.25% | R²=0.50% | R²=1.00% | R²=2.00% |
|---|---|---|---|---|---|
| block_x4 | 0.002 | 0.392 | 0.656 | 0.874 | 1.000 |
| block_v1.0 | 0.080 | 0.666 | 0.862 | 1.000 | 1.000 |
| naive_ols | 0.226 | 0.740 | 0.860 | 1.000 | 1.000 |
| camera | 0.000 | 0.000 | 0.000 | 0.194 | 0.586 |
| unrestricted | 0.226 | 0.790 | 0.860 | 1.000 | 1.000 |
| expr_bin | 0.226 | 0.766 | 0.860 | 1.000 | 1.000 |
| stratified | 0.226 | 0.768 | 0.860 | 1.000 | 1.000 |

## Power at matched 5% FPR (all configurations)

| test | R²=0.00% | R²=0.25% | R²=0.50% | R²=1.00% | R²=2.00% |
|---|---|---|---|---|---|
| block_x4 | 0.066 | 0.400 | 0.704 | 0.880 | 1.000 |
| block_v1.0 | 0.050 | 0.582 | 0.806 | 1.000 | 1.000 |
| naive_ols | 0.058 | 0.580 | 0.792 | 1.000 | 1.000 |
| camera | 0.058 | 0.580 | 0.792 | 1.000 | 1.000 |
| unrestricted | 0.196 | 0.634 | 0.800 | 0.996 | 1.000 |
| expr_bin | 0.166 | 0.590 | 0.790 | 1.000 | 1.000 |
| stratified | 0.172 | 0.576 | 0.796 | 0.996 | 1.000 |

## 95% CI coverage of the planted θ (all configurations)

| test | R²=0.00% | R²=0.25% | R²=0.50% | R²=1.00% | R²=2.00% |
|---|---|---|---|---|---|
| block_x4 | 0.770 | 0.728 | 0.712 | 0.736 | 0.726 |
| block_v1.0 | 0.580 | 0.602 | 0.584 | 0.600 | 0.602 |
| naive_ols | 0.574 | 0.598 | 0.582 | 0.594 | 0.598 |
| camera | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| unrestricted | 0.574 | 0.598 | 0.582 | 0.594 | 0.598 |
| expr_bin | 0.574 | 0.598 | 0.582 | 0.594 | 0.598 |
| stratified | 0.574 | 0.598 | 0.582 | 0.594 | 0.598 |

## 95% CI coverage of the planted θ, excluding decoys with significant θ₀ (main test)

| test | R²=0.00% | R²=0.25% | R²=0.50% | R²=1.00% | R²=2.00% |
|---|---|---|---|---|---|
| block_x4 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| block_v1.0 | 0.749 | 0.828 | 0.831 | 0.814 | 0.833 |
| naive_ols | 0.749 | 0.828 | 0.831 | 0.814 | 0.833 |
| camera | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| unrestricted | 0.749 | 0.828 | 0.831 | 0.814 | 0.833 |
| expr_bin | 0.749 | 0.828 | 0.831 | 0.814 | 0.833 |
| stratified | 0.749 | 0.828 | 0.831 | 0.814 | 0.833 |

## Bias of θ̂ (same estimator for every test)

| level | θ planted | bias | SE (cluster) |
|---|---|---|---|
| 0.00% | 0 | 0.000328 | 0.00328 |
| 0.25% | 0.01743 | -0.000219 | 0.00307 |
| 0.50% | 0.025 | -0.000238 | 0.00322 |
| 1.00% | 0.0352 | -0.000388 | 0.00315 |
| 2.00% | 0.05008 | 0.000177 | 0.0031 |

## Level 0 by configuration (θ₀ = the signed decoy's own coefficient in the real M)

| decoy | sign | n | θ₀ | rej block_x4 | rej block_v1.0 | rej naive |
|---|---|---|---|---|---|---|
| A549-Panc1 | -1 | 27.0 | -0.0161 | 0.00 | 0.00 | 0.00 |
| A549-Panc1 | +1 | 32.0 | 0.0161 | 1.00 | 1.00 | 1.00 |
| GM12878-HCT116 | -1 | 32.0 | -0.0002 | 0.00 | 0.00 | 0.00 |
| GM12878-HCT116 | +1 | 42.0 | 0.0002 | 0.00 | 0.00 | 0.00 |
| GM12878-MCF-7 | -1 | 36.0 | 0.0031 | 0.00 | 0.00 | 0.00 |
| GM12878-MCF-7 | +1 | 28.0 | -0.0031 | 0.00 | 0.00 | 0.00 |
| HCT116-MCF-7 | -1 | 34.0 | 0.0071 | 0.00 | 0.00 | 0.00 |
| HCT116-MCF-7 | +1 | 37.0 | -0.0071 | 0.00 | 0.00 | 0.00 |
| IMR-90-PC-3 | -1 | 29.0 | 0.0182 | 1.00 | 1.00 | 1.00 |
| IMR-90-PC-3 | +1 | 29.0 | -0.0182 | 0.00 | 0.00 | 0.00 |
| MCF-7-A549 | -1 | 44.0 | -0.0054 | 0.00 | 0.00 | 0.00 |
| MCF-7-A549 | +1 | 34.0 | 0.0054 | 0.00 | 0.00 | 0.00 |
| MCF-7-Panc1 | -1 | 44.0 | -0.0164 | 0.00 | 0.00 | 0.00 |
| MCF-7-Panc1 | +1 | 52.0 | 0.0164 | 0.00 | 1.00 | 1.00 |

## Headline (main test = block permutation, x4 rule)

- **Detection limit (80% power at nominal 5%): 1.00%** with all configurations; 1.00% excluding decoys with significant θ₀ (main-test criterion); 1.00% with the naive-OLS criterion. At matched 5% FPR: 1.00%. v1.2: 1.00%.
  Power at 0.5%: 0.782 (all) / 0.794 (excluded); v1.2 0.602.
- **Level-0 FPR at nominal 5%: 0.122 ± 0.085** (all 14 configurations); **0.000 ± 0.000** excluding A549-Panc1, IMR-90-PC-3 (10 configurations, n = 383); 0.000 with the naive criterion. v1.2: 0.096.
- v1.0 rule (sensitivity): detection limit 0.50% (all) / 0.25% (excluded); level-0 FPR 0.226 / 0.136.

## Reading

- All level-0 rejections of the main test come from configurations whose decoy is itself associated with the real M with θ₀ > 0 (A549−Panc1 +, IMR-90−PC-3 −). Both decoys have two-sided main-test p ≈ 0.02 for θ₀ by themselves; without them the main test never rejects at level 0. Level 0 is therefore not a clean null; the decoys are not guaranteed signal-free.
- Under v1.3 the negative-control set changed (7 decoys: 3 new, 2 of the v1.2 six dropped) and three of them carry |θ₀| ≈ 0.016–0.018, comparable to θ at the 0.25–0.5% levels. Because every test estimates θ + θ₀, CI coverage of the planted θ falls (main test 0.71–0.77 over all configurations); excluding the significant-θ₀ decoys the main test covers 1.00 again. Coverage here measures decoy contamination as much as calibration.
- The main test gains power over v1.2 at low levels (0.25%: 0.65 vs 0.57; 0.5%: 0.78 vs 0.60) but just misses 80% at 0.5%, so the detection limit stays 1%. Decoy residuals keep long-range dependence even with arm terms (k x4 = 44–100), unlike the real features (k 20–24).
- Only 14 distinct datasets: cluster SEs (≈ 0.1 at level 0) are the honest uncertainty.
- CAMERA (VIF from ρ̄ ≈ 0.003) remains far too conservative.
