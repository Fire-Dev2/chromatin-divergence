# Calibration of the main test on the real outcomes (v0.2.0; not part of the frozen analysis)

## 1. Decoy-free null: circular shifts of the pair's own feature within chromosome arms

Main test: block permutation, x4 rule, arm terms, B = 999; 2,000 shifted features per pair. Rejection at two-sided 5% (Wilson 95% interval).

The last column uses the shifted features as a decoy-free reference distribution for the observed θ̂ (two-sided: share of shifted |θ̂| at least the observed |θ̂|).

| pair | main test | naive OLS | observed θ̂ | circular-shift p |
|---|---|---|---|---|
| K562_vs_Jurkat | 5.0% [4.1, 6.0] | 5.5% [4.5, 6.5] | 0.0002 | 0.7996 |
| RPE1_vs_HepG2 | 5.2% [4.4, 6.3] | 5.2% [4.4, 6.3] | -0.0173 | 0.0004998 |
| RPE1_vs_K562 | 3.2% [2.5, 4.1] | 3.8% [3.0, 4.7] | -0.0224 | 0.0004998 |
| K562_vs_HepG2 | 4.7% [3.9, 5.7] | 5.2% [4.4, 6.3] | 0.0071 | 0.01499 |
| Jurkat_vs_RPE1 | 6.9% [5.9, 8.1] | 7.3% [6.2, 8.5] | -0.0219 | 0.0004998 |
| Jurkat_vs_HepG2 | 4.2% [3.4, 5.2] | 4.4% [3.6, 5.4] | 0.0024 | 0.2924 |
| K562_vs_K562gw | 4.8% [3.9, 5.8] | 4.9% [4.0, 5.9] |  |  |
| six pairs pooled | 4.9% [4.5, 5.3] | 5.2% [4.8, 5.6] | | |

## 2. Planted signal on the same-line M (K562 essential − K562 genome-wide)

Rejection rate of the main test (and naive OLS) at two-sided 5% by planted partial R²; level 0 is the false-positive rate on this M.

| null features | planted partial R² | main test | naive OLS |
|---|---|---|---|
| circular | 0% | 4.6% [3.1, 6.8] | 5.4% |
| circular | 0.25% | 96.2% [94.1, 97.6] | 97.0% |
| circular | 0.5% | 100.0% [99.2, 100.0] | 100.0% |
| circular | 1% | 100.0% [99.2, 100.0] | 100.0% |
| circular | 2% | 100.0% [99.2, 100.0] | 100.0% |
| decoy | 0% | 25.6% [22.0, 29.6] | 53.8% |
| decoy | 0.25% | 61.8% [57.5, 66.0] | 70.2% |
| decoy | 0.5% | 78.6% [74.8, 82.0] | 93.6% |
| decoy | 1% | 95.0% [92.7, 96.6] | 100.0% |
| decoy | 2% | 100.0% [99.2, 100.0] | 100.0% |

Decoy-line differences tested against the same-line M (B = 9,999): 6 of 21 reject with the main test (11 of 21 with naive OLS); the six pair features: 3 of 6 (Figure 4A).
