# Aim 2 test, development pair RPE1 − K562 (A = RPE1, B = K562)

**Status: descriptive.** RPE1 and K562 differ in p53 status, tissue and growth, so this is a
mixed pair; θ̂ here is not a test of H1. It is the pipeline's development run, made before
tag `pipeline-v1.0`. Code: `cdtest.aim2` (`scripts/run_aim2_dev.py`).

- Genes 5,262; knockdowns 1,314; reliability of M 0.9955 (knockdown-bootstrap 95% [0.9906, 0.9927]).
- Covariates C: log2_meanA, log2_meanB, log2_varA, log2_varB, expr_diff, cn_diff, cpg_oe; f(C) has 14 columns (intercept, ns(df 4) on each line's log2 mean, 5 linear; rank one less because the expression difference is collinear).
- Permutations B = 9,999 per test; partial-R² bootstrap 500 resamples of knockdowns.
- CAMERA ρ̄ (control cells, normalized counts, mean of the two lines) = 0.00304 (RPE1 0.00331, K562 0.00277; on log2(x+1): 0.00929 / 0.00676).

## Main test (block permutation, k rule x4 = plan v1.1)

| feature | θ̂ | k (x4 / v1.0) | one-sided p | two-sided p | 95% CI (refined) | 95% CI (plan grid) | partial R² / reliability [95% boot] |
|---|---|---|---|---|---|---|---|
| combined_cpm | -0.02239 | 100 / 50 | 0.9976 | 0.0025 | [-0.03184, -0.006324] | [-0.01709, -0.01709] | 0.674% [0.6233, 0.7212]% |
| promoter_cpm | -0.01276 | 100 / 50 | 0.9975 | 0.0026 | [-0.01618, -0.005037] | [-0.01603, -0.008014] | 1.084% [1.027, 1.131]% |
| distal_cpm | -0.003407 | 100 / 32 | 0.6885 | 0.3116 | [-0.0141, 0.006844] | [-0.01021, 0] | 0.043% [0.03108, 0.05651]% |
| combined_fc | -0.00254 | 92 / 23 | 0.8745 | 0.2142 | [-0.03647, 0.008111] | [-0.02458, 0] | 0.005% [0.001384, 0.009719]% |

Block exchangeability (main x4 rule; v1.0 in brackets): share of genes in blocks that have at
least one other block of the same size, i.e. that the permutation can move.

| feature | blocks | movable blocks | genes movable | null mean of T* | null SD of T* |
|---|---|---|---|---|---|
| combined_cpm | 49 [92] | 21 [68] | 41% [67%] | -0.01836 | 0.0015 |
| promoter_cpm | 49 [92] | 21 [68] | 41% [67%] | -0.01049 | 0.000829 |
| distal_cpm | 49 [143] | 21 [126] | 41% [86%] | -0.002771 | 0.00122 |
| combined_fc | 49 [207] | 17 [197] | 30% [93%] | -0.000187 | 0.00201 |

## Sensitivity analyses and comparators (one-sided p for θ > 0; two-sided p; 95% CI)

| feature | block_v1.0 | naive_ols | unrestricted | expr_bin | stratified | camera | block_coaccess_clusters | decoy_neg_control |
|---|---|---|---|---|---|---|---|---|
| combined_cpm | 1 / 0.0002 [-0.03044, -0.01002] | 1 / 2.78e-09 [-0.02976, -0.01502] | 1 / 0.0001 [-0.02933, -0.01557] | 1 / 0.0001 [-0.03037, -0.01554] | 1 / 0.0001 [-0.03141, -0.01632] | 0.926 / 0.149 [-0.05277, 0.007995] | 1 / 0.0002 [-0.04276, -0.01346] | p=1 / 0.1429 (n=6) |
| promoter_cpm | 1 / 0.0001 [-0.0159, -0.007164] | 1 / 4.53e-14 [-0.01607, -0.009456] | 1 / 0.0001 [-0.01594, -0.009526] | 1 / 0.0001 [-0.0162, -0.009464] | 1 / 0.0001 [-0.01672, -0.009817] | 0.967 / 0.0665 [-0.0264, 0.0008689] | 1 / 0.0002 [-0.02042, -0.006708] | p=1 / 0.1429 (n=6) |
| distal_cpm | 0.789 / 0.225 [-0.008396, 0.0038] | 0.934 / 0.133 [-0.00785, 0.001035] | 0.948 / 0.102 [-0.007538, 0.0007355] | 0.962 / 0.122 [-0.00811, 0.0005005] | 0.965 / 0.147 [-0.00857, 0.0003167] | 0.642 / 0.715 [-0.02172, 0.0149] | 0.983 / 0.0176 [-0.01926, -0.0009091] | p=1 / 0.3636 (n=10) |
| combined_fc | 0.647 / 0.615 [-0.01326, 0.009321] | 0.687 / 0.625 [-0.01273, 0.007649] | 0.847 / 0.309 [-0.007484, 0.002323] | 0.709 / 0.376 [-0.006563, 0.003515] | 0.745 / 0.434 [-0.008443, 0.004314] | 0.547 / 0.906 [-0.04454, 0.03946] | 0.914 / 0.193 [-0.04237, 0.007957] | p=0.7778 / 0.4444 (n=8) |

## Naive analysis (Aim 1 step 7: M ~ 1 + expression difference + x, genes independent)

| feature | coefficient | SE | partial R² | one-sided p | two-sided p |
|---|---|---|---|---|---|
| combined_cpm | -0.01944 | 0.00484 | 0.305% | 1 | 6.1e-05 |
| promoter_cpm | -0.01275 | 0.00225 | 0.609% | 1 | 1.45e-08 |
| distal_cpm | -0.001577 | 0.00288 | 0.006% | 0.708 | 0.583 |
| combined_fc | -0.0008153 | 0.00347 | 0.001% | 0.593 | 0.814 |

## Reading

- The primary feature (combined, CPM) has a **negative** coefficient, θ̂ = -0.0224, so the one-sided test of θ > 0 gives p = 0.9976. Genes more accessible in RPE1 than in K562 have a *smaller* RPE1-minus-K562 response difference. Two-sided, the block null puts it at p = 0.0025; every comparator agrees on the sign. The promoter feature carries it; distal and the FC version are near zero.
- Size: partial R² 0.67% of reproducible variance for combined, 1.08% for promoter.
- This is the development pair, which differs in p53, tissue and growth; it is not evidence for or against H1.
- **Pipeline issue found here (for the author, before the matched pairs).** At G ≈ 5,300 genes on 41 arms the x residuals stay autocorrelated past lag 50, so k_v1.0 hits its cap (50) and the x4 rule its cap (100). Blocks are then 100–199 genes, and because the plan permutes blocks only within equal-size classes, only ~41% of genes (21 of 49 blocks) can move. The fixed blocks keep their real x, so the block null is centred near the observed association (null mean of T* ≈ −0.017 for combined), not at 0. The test stays valid (exchangeable units), but it tests only the movable part and loses power; see the planted benchmark (results/aim3/planted_RPE1_vs_K562.md). No parameter was changed.
- CI: the plan's 201-point grid spans ±max(4|θ̂|, 0.5/SD(x)); with SD(x) ≈ 0.3 its step (≈0.017) is about the size of θ̂, so the plan-grid CI collapses to one or two grid points. The table also gives the same inversion on a 2,001-point grid inside the plan span (docs/methods_aim2_aim3.md §6.2).
- Partial-R² bootstrap: approximate (β̃ and s̃² held at full-data values; see methods §6.3). It reflects knockdown sampling only, not gene-level dependence.
