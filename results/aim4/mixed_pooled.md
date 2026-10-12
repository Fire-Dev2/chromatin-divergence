# Mixed pairs — pooled contrast (plan §7) and matched-vs-mixed gap

Statistic: mean of the four mixed-pair θ̂ (combined CPM; RPE1 − K562, K562 − HepG2, Jurkat − RPE1, Jurkat − HepG2; each pair's own signed x = A − B and M = A − B, so under the effect model every θ is positive regardless of orientation).

**Shared permutation (documented choice).** One block permutation is applied identically to all four pairs. Exact sharing requires the same genes, so the statistic is computed on the genes kept in all four pairs (4,688 genes = genes kept in all four lines; the same set for the two matched pairs and for all six). Genes in genome order, blocks cut within arms with k = the largest of the four pairs' x4 k on these genes (k = 28; 90.7% of genes in movable blocks), one permutation of equal-size blocks per draw applied to every pair's x residuals (each pair keeps its own f(C) with arm terms), B = 9,999. CI by inverting the pooled statistic with one common θ. The union-order alternative requested (each pair on its own genes, blocks on the union order) cannot move the same gene to the same gene in every pair; an implementation that couples the pairs through a shared ordering of union blocks gave null correlations ≈ 0 even for a 90% subset of one pair against itself, i.e. it behaves like independent permutations. It is reported as secondary.

| contrast | genes | θ̂ pooled | one-sided p | two-sided p | 95% CI | null mean / SD |
|---|---|---|---|---|---|---|
| mixed (4 pairs), primary | 4,688 | -0.009614 | 1 | 0.0001 | [-0.0187, -0.00702] | 0.000462 / 0.00239 |
| matched (2 pairs) | 4,688 | -0.005461 | 1 | 0.0001 | [-0.0104, -0.00414] | 0.000161 / 0.00126 |
| all six | 4,688 | -0.00823 | 1 | 0.0001 | [-0.0152, -0.00671] | 0.000245 / 0.00177 |
| mixed, union-rank coupling (secondary; full gene sets) | 6,546 union | -0.008699 | 1 | 0.0001 | [-0.013, -0.00491] | -0.00147 / 0.00168 |

Per-pair θ̂ on the common genes (shared draws) vs on the pair's full gene set (own main test):

| pair | θ̂ common genes | p₂ shared draws | θ̂ full set | p₂ own main test |
|---|---|---|---|---|
| K562_vs_Jurkat | 0.001562 | 0.346 | 0.0001979 | 0.838 |
| RPE1_vs_HepG2 | -0.01248 | 0.0001 | -0.01734 | 0.0001 |
| RPE1_vs_K562 | -0.02304 | 0.0001 | -0.02235 | 0.0001 |
| K562_vs_HepG2 | 0.007954 | 0.0397 | 0.007052 | 0.0551 |
| Jurkat_vs_RPE1 | -0.02551 | 0.0001 | -0.0219 | 0.0002 |
| Jurkat_vs_HepG2 | 0.002137 | 0.403 | 0.002405 | 0.244 |

Null correlation of the four mixed-pair T* under the shared permutation: [[1.0, 0.225, 0.477, -0.236], [0.225, 1.0, -0.015, 0.362], [0.477, -0.015, 1.0, 0.343], [-0.236, 0.362, 0.343, 1.0]].

**Gap (mixed − matched)** on the common genes, one permutation shared by all six pairs: -0.004153 (mixed -0.009614, matched -0.005461); null SD 0.00228; approximate 95% interval [-0.00863, 0.000325]; two-sided p (joint null) 0.0564.

Reading: the mixed pooled θ̂ is negative (opposite to H1) and clearly nonzero; the matched pooled θ̂ is also negative, driven by RPE1 − HepG2. The proposal expected a larger positive effect in the mixed pairs; instead both are negative and the gap is modest and not clearly different from 0. Within the mixed pairs the sign splits by line: the two RPE1-containing pairs are strongly negative, K562 − HepG2 is positive (p₂ ≈ 0.05), Jurkat − HepG2 ≈ 0.
