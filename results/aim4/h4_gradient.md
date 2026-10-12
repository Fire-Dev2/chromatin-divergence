# H4 (secondary) — similarity gradient (plan §7)

For each pair: θ̂ (M ~ f(C) + v, v3 f(C) with arm terms) for every available combined-CPM feature difference v — the pair's own x, the other five pairs' x (each in its plan orientation) and the 21 decoy-line differences — against similarity = Pearson correlation of v and the pair's own x, both residualised on C. OLS of θ̂ on similarity per pair (n = 27) and pooled over pairs (n = 162). H4 'supported if' effect increases with similarity across other pairs and decoys and is near zero at low similarity: slope > 0 and intercept ≈ 0. 'θ_std' = θ̂ × SD(r_v)/SD(r_x), which equals θ_own × similarity exactly under the effect model (diagnostic).

| pair | slope [95% CI] p | intercept [95% CI] p | R² | slope excl. own [p] | slope θ_std [p] |
|---|---|---|---|---|---|
| K562_vs_Jurkat | -0.0001183 [-0.00268, 0.00245] 0.925 | -0.0006994 [-0.00139, -1.23e-05] 0.0463 | 0.00 | -0.001467 [0.453] | 6.506e-05 [0.952] |
| RPE1_vs_HepG2 | -0.03102 [-0.0454, -0.0167] 0.000156 | 0.007049 [0.00313, 0.011] 0.00105 | 0.44 | -0.03802 [0.000866] | -0.0302 [0.000373] |
| RPE1_vs_K562 | -0.0328 [-0.0548, -0.0108] 0.00515 | 0.005356 [-0.000913, 0.0116] 0.0907 | 0.27 | -0.03817 [0.0232] | -0.03094 [0.00253] |
| K562_vs_HepG2 | 0.01313 [-0.00399, 0.0302] 0.127 | -0.004502 [-0.00889, -0.000109] 0.045 | 0.09 | 0.01515 [0.24] | 0.01279 [0.143] |
| Jurkat_vs_RPE1 | -0.03253 [-0.0554, -0.00963] 0.00722 | -0.002061 [-0.00863, 0.00451] 0.524 | 0.26 | -0.04296 [0.00736] | -0.0316 [0.00245] |
| Jurkat_vs_HepG2 | 0.007204 [-0.00859, 0.023] 0.357 | -0.002611 [-0.00708, 0.00186] 0.24 | 0.03 | 0.009131 [0.392] | 0.00584 [0.36] |
| pooled_all_pairs | -0.01587 [-0.0234, -0.0084] 4.55e-05 | 0.0008883 [-0.00117, 0.00295] 0.396 | 0.10 | -0.02197 [4.54e-05] | -0.01495 [2.28e-05] |

Pooled OLS treats the 162 points as independent (they are not: features and Ms are shared across pairs); descriptive.

## Decoy negative control (|similarity| < 0.05; each decoy's θ̂ vs its own block-permutation null, x4, B = 9,999, two-sided)

| pair | n decoys | n with p ≤ 0.05 | decoys significant (θ̂) |
|---|---|---|---|
| K562_vs_Jurkat | 13 | 2 | HCT116-IMR-90 (-0.0033), A549-IMR-90 (-0.0038) |
| RPE1_vs_HepG2 | 6 | 4 | GM12878-HCT116 (0.0125), GM12878-MCF-7 (0.0063), GM12878-A549 (0.0112), HCT116-MCF-7 (-0.0125) |
| RPE1_vs_K562 | 7 | 2 | A549-Panc1 (0.0161), IMR-90-PC-3 (-0.0182) |
| K562_vs_HepG2 | 6 | 4 | GM12878-IMR-90 (-0.0107), MCF-7-PC-3 (-0.0126), A549-PC-3 (-0.0060), IMR-90-Panc1 (0.0169) |
| Jurkat_vs_RPE1 | 4 | 1 | GM12878-HCT116 (0.0108) |
| Jurkat_vs_HepG2 | 7 | 4 | HCT116-IMR-90 (-0.0259), HCT116-Panc1 (-0.0143), MCF-7-IMR-90 (-0.0209), IMR-90-Panc1 (0.0145) |

Overall 17 of 43 low-similarity decoy tests reject at 5% (expected ≈ 2.1 if decoys carried no signal and the null were exact). Low similarity to the pair's x does not make a feature signal-free in M.

Reading: the pooled slope is negative (θ̂ falls as similarity rises) and the intercept is ≈ 0. A negative slope is what the effect model predicts when the pair's own θ is negative: the gradient follows the own effect, whose sign is opposite to H1. H4 as stated (effect increases with similarity) is not supported; the near-zero intercept says the features do not carry a generic, similarity-independent association.
