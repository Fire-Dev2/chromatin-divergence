# H1 sensitivity analyses (labelled sensitivity; not part of the confirmatory test)

Combined feature unless stated; one-sided p is for θ > 0.

## K562_vs_Jurkat

| analysis | θ̂ | one-sided p | two-sided p | 95% CI |
|---|---|---|---|---|
| main (v1.3, x4, CPM) | 0.0001979 | 0.6447 | 0.8381 | [-0.00212, 0.00137] |
| v1.2 (no arm terms) | -8.936e-05 | 0.1977 | 0.8871 | [-0.00299, 0.00523] |
| FC features (combined_fc) | -0.0002331 | 0.8124 | 0.8497 | [-0.00327, 0.0012] |
| v1.0 k rule | 0.0001979 | 0.4478 | 0.808 | [-0.00153, 0.00173] |
| stratified null | 0.0001979 | 0.4708 | 0.8071 | [-0.00153, 0.00169] |
| decoy null (n = 13) | 0.0001979 | 0.2857 | 1 | — |
| co-accessibility-cluster blocks | 0.0001979 | 0.0943 | 0.7854 | [-0.00131, 0.00455] |
| group-mean M (v1.3) | -5.267e-05 | 0.7735 | 0.9559 | [-0.00207, 0.000899] |
| group-mean M (v1.2) | -0.0002374 | 0.2403 | 0.766 | [-0.00273, 0.00416] |

Target-promoter accessibility, knockdown level (exploratory, not pre-specified in detail): Spearman between each knockdown's response divergence (mean over genes of (β̃²−s̃²)_A − (β̃²−s̃²)_B) and its target's promoter CPM difference (A − B), both residualised on the target's expression difference and on-target log2FC difference: ρ = -0.017 over 1,414 knockdowns, one-sided p 0.742, two-sided p 0.526 (9,999 knockdown permutations; knockdowns treated as exchangeable).

## RPE1_vs_HepG2

| analysis | θ̂ | one-sided p | two-sided p | 95% CI |
|---|---|---|---|---|
| main (v1.3, x4, CPM) | -0.01734 | 1 | 0.0001 | [-0.0263, -0.014] |
| v1.2 (no arm terms) | -0.01394 | 1 | 0.0001 | [-0.0407, -0.0189] |
| FC features (combined_fc) | -0.00429 | 0.7619 | 0.2672 | [-0.00972, 0.00474] |
| v1.0 k rule | -0.01734 | 1 | 0.0001 | [-0.0229, -0.0122] |
| stratified null | -0.01734 | 1 | 0.0001 | [-0.0227, -0.013] |
| decoy null (n = 6) | -0.01734 | 1 | 0.1429 | — |
| co-accessibility-cluster blocks | -0.01734 | 0.9588 | 0.0413 | [-0.0271, 0.002] |
| group-mean M (v1.3) | -0.01862 | 1 | 0.0001 | [-0.0267, -0.016] |
| group-mean M (v1.2) | -0.01462 | 1 | 0.0001 | [-0.0414, -0.0212] |

Target-promoter accessibility, knockdown level (exploratory, not pre-specified in detail): Spearman between each knockdown's response divergence (mean over genes of (β̃²−s̃²)_A − (β̃²−s̃²)_B) and its target's promoter CPM difference (A − B), both residualised on the target's expression difference and on-target log2FC difference: ρ = -0.052 over 1,424 knockdowns, one-sided p 0.977, two-sided p 0.0501 (9,999 knockdown permutations; knockdowns treated as exchangeable).

## Target-promoter check, all pairs (exploratory)

| pair | n kd | ρ (resid.) | p₁ | p₂ |
|---|---|---|---|---|
| K562_vs_Jurkat | 1414 | -0.017 | 0.742 | 0.526 |
| RPE1_vs_HepG2 | 1424 | -0.052 | 0.977 | 0.0501 |
| RPE1_vs_K562 | 1314 | 0.001 | 0.493 | 0.982 |
| K562_vs_HepG2 | 1389 | -0.018 | 0.753 | 0.497 |
| Jurkat_vs_RPE1 | 1413 | -0.052 | 0.978 | 0.0452 |
| Jurkat_vs_HepG2 | 1482 | -0.111 | 1 | 0.0002 |

Reading: no sensitivity analysis turns either matched pair positive. K562 − Jurkat stays ≈ 0 under every variant; RPE1 − HepG2 stays negative (two-sided p ≤ 0.01 under every variant that permutes x; the decoy null, with few decoys, has a coarse p).
