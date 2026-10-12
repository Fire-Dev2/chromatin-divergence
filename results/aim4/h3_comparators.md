# H3 (primary in the proposal) — naive tests vs the main test on real data, all six pairs

Combined CPM feature. Comparators as in plan §5 (all with the v1.3 f(C)); 'naive step 7' is the Aim 1 naive analysis M ~ 1 + expr_diff + x (OLS, genes independent, no f(C)). The naive-vs-main gap = θ̂(naive step 7) − θ̂(main) with a 500-resample knockdown bootstrap (β̃, s̃² held fixed; same approximation as the partial-R² bootstrap). The OLS-with-f(C) comparator has the same θ̂ as the main test; its inflation is in the p-value (SE ratio = naive SE / SD of the block-permutation null).

| pair | role | θ̂ main | p₂ main | θ̂ naive step 7 | p₂ naive step 7 | gap [boot 95%] p | naive SE / block-null SD |
|---|---|---|---|---|---|---|---|
| K562_vs_Jurkat | matched | 0.0001979 | 0.838 | 5.872e-05 | 0.942 | -0.0001392 [-0.00042, 0.000154] p=0.36 | 1.03 |
| RPE1_vs_HepG2 | matched | -0.01734 | 0.0001 | -0.009736 | 0.000318 | 0.007604 [0.00603, 0.00901] p=<0.002 | 1.06 |
| RPE1_vs_K562 | mixed | -0.02235 | 0.0001 | -0.01944 | 6.1e-05 | 0.002919 [0.00233, 0.00345] p=<0.002 | 1.03 |
| K562_vs_HepG2 | mixed | 0.007052 | 0.0551 | -0.002592 | 0.411 | -0.009644 [-0.0106, -0.00863] p=<0.002 | 1.06 |
| Jurkat_vs_RPE1 | mixed | -0.0219 | 0.0002 | -0.01291 | 0.00359 | 0.008984 [0.0083, 0.00977] p=<0.002 | 1.07 |
| Jurkat_vs_HepG2 | mixed | 0.002405 | 0.244 | 0.001321 | 0.592 | -0.001084 [-0.00153, -0.000673] p=<0.002 | 1.18 |

Two-sided p by comparator (combined CPM):

| pair | main (block x4) | naive OLS | unrestricted | expr-bin | stratified | CAMERA | block v1.0 | co-access. clusters | decoy (n) |
|---|---|---|---|---|---|---|---|---|---|
| K562_vs_Jurkat | 0.838 | 0.809 | 0.777 | 0.811 | 0.807 | 0.95 | 0.808 | 0.785 | 1 (13) |
| RPE1_vs_HepG2 | 0.0001 | 8.25e-11 | 0.0001 | 0.0001 | 0.0001 | 0.154 | 0.0001 | 0.0413 | 0.143 (6) |
| RPE1_vs_K562 | 0.0001 | 1.01e-07 | 0.0001 | 0.0001 | 0.0001 | 0.196 | 0.0001 | 0.0003 | 0.125 (7) |
| K562_vs_HepG2 | 0.0551 | 0.0215 | 0.008 | 0.0102 | 0.0157 | 0.611 | 0.0221 | 0.032 | 0.571 (6) |
| Jurkat_vs_RPE1 | 0.0002 | 2.87e-07 | 0.0001 | 0.0001 | 0.0001 | 0.194 | 0.0001 | 0.068 | 0.2 (4) |
| Jurkat_vs_HepG2 | 0.244 | 0.32 | 0.241 | 0.26 | 0.281 | 0.819 | 0.33 | 0.107 | 0.75 (7) |

Reading:
- Simulation half of H3 (Aim 3, unchanged): naive tests reject a true null well above 5% under shared structure; on the matched-pair planted benchmarks the naive OLS level-0 FPR is K562_vs_Jurkat 0.32; RPE1_vs_HepG2 0.43 vs main test 0.08; 0.32.
- Real-data half: the naive step-7 estimate is *not* larger than the main-test estimate. In absolute value it is smaller in 4 of 6 pairs (RPE1 − HepG2 −0.0097 vs −0.0173; Jurkat − RPE1 −0.0129 vs −0.0219; RPE1 − K562 −0.0194 vs −0.0224; Jurkat − HepG2 0.0013 vs 0.0024), has the opposite sign in K562 − HepG2 (−0.0026 vs +0.0071) and is equal within noise in K562 − Jurkat. With the full f(C) the naive OLS SE is within 3–18% of the SD of the block-permutation null, so naive OLS p-values are smaller mainly because a t-test has no 1/(B+1) floor; the only pair where naive-type tests reject at 5% and the main test does not is K562 − HepG2 (mixed; naive OLS p₂ 0.022, unrestricted 0.008, expr-bin 0.010, stratified 0.016 vs main 0.055). CAMERA is non-significant everywhere (ρ̄ 0.0025–0.0035, VIF ≈ 15–20, SE × ~4).
- H3 verdict: the simulation half holds (Aim 3; and on these pairs' planted benchmarks naive OLS level-0 FPR is 0.27–0.31 even after excluding decoys with significant θ₀, vs 0.00 for the main test); the real-data half ('a larger real-data effect than the main test') is not supported.
