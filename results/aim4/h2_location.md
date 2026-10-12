# H2 (secondary) — promoter vs distal CPM, matched pairs

Main test per feature (v1.3, B = 9,999). θ is in M units per unit of the feature difference; the two features are on different scales, so the scale-free comparison is θ per SD of x and partial R² / reliability (500-resample knockdown bootstrap). Supported only if distal is larger with non-overlapping intervals.

| pair | feature | θ̂ | 95% CI | θ̂ × SD(x) [CI × SD(x)] | one-/two-sided p | partial R²/rel [boot 95%] |
|---|---|---|---|---|---|---|
| K562_vs_Jurkat | promoter_cpm | -0.0004862 | [-0.00154, -0.000114] | -0.0003008 [-0.000951, -7.08e-05] | 0.988 / 0.231 | 0.041% [0.008%, 0.108%] |
| K562_vs_Jurkat | distal_cpm | 0.001018 | [5.64e-05, 0.00202] | 0.0004764 [2.64e-05, 0.000946] | 0.0212 / 0.0283 | 0.092% [0.036%, 0.169%] |
| RPE1_vs_HepG2 | promoter_cpm | -0.005902 | [-0.00769, -0.00327] | -0.00431 [-0.00562, -0.00239] | 1 / 0.0001 | 0.585% [0.443%, 0.691%] |
| RPE1_vs_HepG2 | distal_cpm | -0.005561 | [-0.011, -0.00342] | -0.002825 [-0.00557, -0.00174] | 1 / 0.0001 | 0.230% [0.169%, 0.286%] |

Mixed pairs (descriptive):

| pair | promoter θ̂ [CI] p₂ | distal θ̂ [CI] p₂ | pR²/rel promoter / distal |
|---|---|---|---|
| RPE1_vs_K562 | -0.01301 [-0.0169, -0.00898] 0.0001 | -0.002215 [-0.00736, 0.0031] 0.343 | 0.936% / 0.017% |
| K562_vs_HepG2 | 0.00198 [-0.0032, 0.00335] 0.483 | 0.003743 [6.89e-05, 0.00728] 0.0287 | 0.035% / 0.086% |
| Jurkat_vs_RPE1 | -0.01575 [-0.0193, -0.0112] 0.0001 | 0.0009169 [-0.00401, 0.00655] 0.663 | 1.285% / 0.003% |
| Jurkat_vs_HepG2 | 1.781e-05 [-0.00223, 0.00194] 0.983 | 0.002273 [-0.00102, 0.00467] 0.119 | 0.000% / 0.048% |

**Verdict (matched pairs):** K562_vs_Jurkat: not supported; RPE1_vs_HepG2: not supported.

- K562 − Jurkat: distal θ̂ > 0 (one-sided p 0.021) and promoter slightly negative, but the partial-R² bootstrap intervals overlap (distal 0.092% [0.036, 0.169] vs promoter 0.041% [0.008, 0.108]) and both shares are below 0.1%.
- RPE1 − HepG2: both features negative (opposite to H1); promoter carries the larger share (0.59% vs 0.23%, non-overlapping intervals) — the reverse of H2.
- Across all six pairs: the promoter feature carries a negative association in every RPE1-containing pair; in the three pairs without RPE1 the distal feature is weakly positive (K562 − Jurkat p₂ 0.028, K562 − HepG2 0.029, Jurkat − HepG2 0.12) and the promoter ≈ 0. Descriptive (secondary).
