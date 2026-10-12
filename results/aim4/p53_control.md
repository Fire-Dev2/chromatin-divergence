# p53 control — all six pairs

Aim 1 p53 gate statistic: mean M over the Fischer 2017 core targets kept for the pair vs 1,000 random gene sets matched on expression deciles in both lines (`aim1.p53_gate` algorithm and seed; RPE1 − K562 reproduces the stored gate exactly). Oriented so that positive = more response in the p53-WT line (K562 and Jurkat p53-mutant/null; RPE1 and HepG2 WT). Matched pairs keep A − B and are read two-sided. Expectation (proposal): divergence in the four mixed pairs, none in the two matched pairs.

| pair | p53 differs | n targets | stat (oriented) | null mean | z | one-sided p (oriented) | two-sided p | as expected? |
|---|---|---|---|---|---|---|---|---|
| K562_vs_Jurkat | False | 30 | -0.0080 | -0.0058 | -0.72 | 0.802 | 0.397 | yes |
| RPE1_vs_HepG2 | False | 40 | 0.0536 | 0.0220 | 3.34 | 0.000999 | 0.000999 | NO |
| RPE1_vs_K562 | True | 25 | 0.1312 | 0.0497 | 4.05 | 0.000999 | 0.000999 | yes |
| K562_vs_HepG2 | True | 36 | 0.1030 | 0.0297 | 6.95 | 0.000999 | 0.000999 | yes |
| Jurkat_vs_RPE1 | True | 25 | 0.1239 | 0.0383 | 5.30 | 0.002 | 0.002 | yes |
| Jurkat_vs_HepG2 | True | 32 | 0.0768 | 0.0186 | 6.14 | 0.000999 | 0.000999 | yes |

Reading: all four mixed pairs diverge at p53 targets in the expected direction (z 4–7). K562 − Jurkat (both p53-deficient) shows none, as expected. **RPE1 − HepG2 (both p53 WT) diverges strongly (z = 3.3, RPE1 > HepG2)** — the matched pair is not matched in p53 pathway activity, so the RPE1 − HepG2 comparison does not fully hold the p53 axis constant.
