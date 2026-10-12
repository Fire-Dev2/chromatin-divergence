# Aim 2 test, development pair RPE1 − K562 — plan v1.3 (arm fixed effects in C)

**Status: descriptive** (mixed pair: p53, tissue, growth; not an H1 test). Plan v1.3 adds one indicator per chromosome arm to f(C) (reference arm dropped); the arm terms enter the model M ~ f(C) + x, the residualisation of x that sets k and is block-permuted, every comparator and the decoy similarities. Everything else as plan v1.1/v1.2. Code: `cdtest.aim2.run(arm_effects=True)` (`scripts/run_aim2_dev_v13.py`). The v1.2 results (`RPE1_vs_K562_test.{json,md}`) are kept unchanged and shown in brackets.

- Genes 5,262 on 41 arms; knockdowns 1,314; reliability 0.9955.
- f(C): 54 columns (14 of v1.2 + 40 arm indicators), rank 53.
- B = 9,999 per test; partial-R² bootstrap 500 knockdown resamples; CAMERA ρ̄ 0.00304.
- v1.2 reproduction: `arm_effects=False` reruns reproduce `RPE1_vs_K562_test.json` exactly (max |Δ| = 0, null draws identical) and the planted benchmark replicates exactly (`results/aim2/v12_reproduction.json`).

## Main test (block permutation, x4 rule) — v1.3 [v1.2]

| feature | θ̂ | k x4 | k v1.0 | one-sided p | two-sided p | 95% CI (refined) | partial R² / reliability [95% boot] |
|---|---|---|---|---|---|---|---|
| combined_cpm | -0.02235 [-0.02239] | 20 [100] | 5 [50] | 1 [0.9976] | 0.0001 [0.0025] | [-0.03152, -0.01345] [[-0.03184, -0.006324]] | 0.545% [0.4986, 0.5944]% [0.674%] |
| promoter_cpm | -0.01301 [-0.01276] | 20 [100] | 5 [50] | 1 [0.9975] | 0.0001 [0.0026] | [-0.01688, -0.008983] [[-0.01618, -0.005037]] | 0.936% [0.8823, 0.9894]% [1.084%] |
| distal_cpm | -0.002215 [-0.003407] | 20 [100] | 5 [32] | 0.7918 [0.6885] | 0.3431 [0.3116] | [-0.007365, 0.003095] [[-0.0141, 0.006844]] | 0.017% [0.009546, 0.0252]% [0.043%] |
| combined_fc | -0.002605 [-0.00254] | 24 [92] | 6 [23] | 0.7602 [0.8745] | 0.617 [0.2142] | [-0.01534, 0.00757] [[-0.03647, 0.008111]] | 0.005% [0.001382, 0.009591]% [0.005%] |

Plan-grid CIs (201 points, as specified): combined_cpm [-0.01709, -0.01709]; promoter_cpm [-0.01603, -0.01603]; distal_cpm [0, 0]; combined_fc [-0.01229, 0].

## Block diagnostics — v1.3 [v1.2]

Exchangeability = share of genes in blocks with at least one other block of the same size (only those move, plan §5 step 5). Null mean/SD = mean/SD of the B block-permutation T*.

| feature | blocks (x4) | genes movable, x4 | genes movable, v1.0 | null mean T* | null SD T* | observed T | (T − null mean)/null SD |
|---|---|---|---|---|---|---|---|
| combined_cpm | 241 [49] | 96.1% [41.1%] | 100.0% [67.2%] | -0.002589 [-0.01836] | 0.00408 [0.0015] | -0.02235 [-0.02239] | -4.8 [-2.7] |
| promoter_cpm | 241 [49] | 96.1% [41.1%] | 100.0% [67.2%] | -0.001502 [-0.01049] | 0.00181 [0.000829] | -0.01301 [-0.01276] | -6.4 [-2.7] |
| distal_cpm | 241 [49] | 96.1% [41.1%] | 100.0% [85.7%] | -0.0003396 [-0.002771] | 0.00233 [0.00122] | -0.002215 [-0.003407] | -0.8 [-0.5] |
| combined_fc | 201 [49] | 94.2% [30.2%] | 100.0% [93.1%] | 0.0009855 [-0.000187] | 0.00509 [0.00201] | -0.002605 [-0.00254] | -0.7 [-1.2] |

Within-arm ACF of the x residuals (after f(C); v1.3 includes the arm terms) — v1.3 [v1.2]:

| series | lag 1 | lag 5 | lag 10 | lag 25 | lag 50 | lag 100 |
|---|---|---|---|---|---|---|
| x: combined_cpm | 0.229 [0.377] | 0.049 [0.227] | 0.014 [0.194] | 0.011 [0.169] | -0.047 [0.055] | 0.025 [0.031] |
| x: promoter_cpm | 0.115 [0.269] | 0.015 [0.177] | 0.010 [0.168] | 0.021 [0.160] | -0.038 [0.061] | 0.017 [0.030] |
| x: distal_cpm | 0.260 [0.328] | 0.046 [0.137] | 0.008 [0.100] | 0.013 [0.093] | -0.040 [0.015] | -0.006 [0.013] |
| x: combined_fc | 0.191 [0.237] | 0.053 [0.105] | 0.021 [0.075] | 0.006 [0.063] | 0.005 [0.063] | 0.024 [0.056] |
| M (diagnostic only) | 0.045 [0.057] | -0.011 [0.002] | -0.033 [-0.019] | 0.014 [0.025] | -0.006 [0.001] | 0.002 [0.011] |

## Sensitivity analyses and comparators, v1.3 (one-sided p / two-sided p [95% CI])

| feature | block_v1.0 | naive_ols | unrestricted | expr_bin | stratified | camera | block_coaccess_clusters | decoy_neg_control |
|---|---|---|---|---|---|---|---|---|
| combined_cpm | 1 / 0.0001 [-0.03082, -0.01407] | 1 / 1.01e-07 [-0.03057, -0.01413] | 1 / 0.0001 [-0.02917, -0.01545] | 1 / 0.0001 [-0.03032, -0.01562] | 1 / 0.0001 [-0.03035, -0.01542] | 0.902 / 0.196 [-0.05623, 0.01153] | 1 / 0.0003 [-0.0484, -0.01323] | p=1 / 0.125 (n=7) |
| promoter_cpm | 1 / 0.0001 [-0.01681, -0.009312] | 1 / 2.88e-12 [-0.01665, -0.009364] | 1 / 0.0001 [-0.01616, -0.009721] | 1 / 0.0001 [-0.01646, -0.009681] | 1 / 0.0001 [-0.0166, -0.009737] | 0.955 / 0.0895 [-0.02802, 0.002008] | 0.999 / 0.001 [-0.02202, -0.005542] | p=1 / 0.25 (n=3) |
| distal_cpm | 0.818 / 0.357 [-0.006926, 0.002625] | 0.824 / 0.352 [-0.006878, 0.002448] | 0.858 / 0.29 [-0.006343, 0.001951] | 0.886 / 0.315 [-0.006752, 0.001696] | 0.879 / 0.344 [-0.007242, 0.001869] | 0.589 / 0.821 [-0.02144, 0.017] | 0.985 / 0.0415 [-0.02068, -0.001073] | p=1 / 0.625 (n=7) |
| combined_fc | 0.689 / 0.632 [-0.01348, 0.008289] | 0.687 / 0.626 [-0.01308, 0.00787] | 0.852 / 0.303 [-0.007607, 0.002237] | 0.734 / 0.346 [-0.006698, 0.00349] | 0.675 / 0.485 [-0.008553, 0.005506] | 0.547 / 0.906 [-0.04578, 0.04057] | 0.858 / 0.221 [-0.0458, 0.01285] | p=0.8 / 0.9 (n=9) |

Two-sided p under v1.2 for the same comparators:

| feature | block_v1.0 | naive_ols | unrestricted | expr_bin | stratified | camera | block_coaccess_clusters | decoy_neg_control |
|---|---|---|---|---|---|---|---|---|
| combined_cpm | 0.0002 | 2.78e-09 | 0.0001 | 0.0001 | 0.0001 | 0.149 | 0.0002 | 0.143 |
| promoter_cpm | 0.0001 | 4.53e-14 | 0.0001 | 0.0001 | 0.0001 | 0.0665 | 0.0002 | 0.143 |
| distal_cpm | 0.225 | 0.133 | 0.102 | 0.122 | 0.147 | 0.715 | 0.0176 | 0.364 |
| combined_fc | 0.615 | 0.625 | 0.309 | 0.376 | 0.434 | 0.906 | 0.193 | 0.444 |

## Reading

- Arm terms fix the block problem. k for combined CPM falls from 100 to 20 (x4) and from 50 to 5 (v1.0); genes in movable blocks rise from 41% to 96%; the block null is now centred near 0 (mean -0.0026, v1.2 -0.0184) and wider (SD 0.0041 vs 0.0015).
- The estimate barely moves: θ̂ = -0.0224 (v1.2 -0.0224); the association is within arms, not an arm-scale artefact. One-sided p (θ > 0) = 1; two-sided 0.0001 (the minimum at B = 9,999). Partial R² / reliability 0.55% (v1.2 0.67%); promoter 0.94%.
- Sign and pattern as in v1.2: genes more accessible in RPE1 than K562 have a smaller RPE1-minus-K562 response difference, carried by the promoter feature; distal and FC near zero. All comparators agree on the sign.
- Not evidence for or against H1 (mixed pair).
