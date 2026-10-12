# Aim 4 step 1 — Aim 1 outcome M for all six pairs (frozen pipeline pipeline-v1.1)

Streaming pseudobulk pass (`screens.run`) on Jurkat and HepG2, then `aim1.run_pair(A, B)` with plan §0 orientation (M = A − B). RPE1 − K562 is the development-pair run (reused). Logs: `results/aim1/screens_*.log`, `results/aim4/logs/aim1_*.log`.

- Jurkat: 262,956 cells, 257,157 kept, 11,892 controls, 2,393 knockdowns, 8,882 genes.
- HepG2: 145,473 cells, 142,813 kept, 4,926 controls, 2,393 knockdowns, 9,624 genes.

| pair | role | genes | knockdowns kept / shared | few cells | target not expr. | weak repression | split-half r | reliability | generic share A / B | mean M | SD M |
|---|---|---|---|---|---|---|---|---|---|---|---|
| K562_vs_Jurkat | matched | 5,675 | 1,414 / 2,055 | 153 | 385 | 400 | 0.759 | 0.863 | 2.4% / 3.2% | -0.0054 | 0.0166 |
| RPE1_vs_HepG2 | matched | 5,631 | 1,424 / 2,392 | 386 | 620 | 416 | 0.958 | 0.978 | 11.9% / 4.3% | 0.0214 | 0.0619 |
| RPE1_vs_K562 | mixed | 5,262 | 1,314 / 2,055 | 210 | 510 | 360 | 0.991 | 0.996 | 13.0% / 2.7% | 0.0476 | 0.1004 |
| K562_vs_HepG2 | mixed | 5,785 | 1,389 / 2,055 | 245 | 407 | 321 | 0.964 | 0.982 | 2.7% / 4.8% | -0.0273 | 0.0597 |
| Jurkat_vs_RPE1 | mixed | 5,222 | 1,413 / 2,392 | 297 | 646 | 503 | 0.990 | 0.995 | 3.1% / 12.1% | -0.0405 | 0.0963 |
| Jurkat_vs_HepG2 | mixed | 5,587 | 1,482 / 2,392 | 322 | 542 | 469 | 0.954 | 0.977 | 3.1% / 4.2% | -0.0172 | 0.0535 |

Same-line baseline (K562 essential − K562 genome-wide; `results/aim1/same_line_baseline.json`): SD M 0.0095, reliability 0.72, generic share 2.5% / 2.9%.

Reproducible variance of M (SD² × reliability) relative to the same-line baseline (0.0095² × 0.72 = 6.5e-5):

- K562_vs_Jurkat: 2.37e-04 (3.6× baseline)
- RPE1_vs_HepG2: 3.75e-03 (57.2× baseline)
- RPE1_vs_K562: 1.00e-02 (152.8× baseline)
- K562_vs_HepG2: 3.50e-03 (53.4× baseline)
- Jurkat_vs_RPE1: 9.22e-03 (140.6× baseline)
- Jurkat_vs_HepG2: 2.80e-03 (42.7× baseline)

Notes: K562 − Jurkat has the smallest and least reliable M (SD 0.017, reliability 0.86), about 4× the same-line baseline in reproducible variance; the other pairs are 40–150×. The generic share is below 13% in every line (largest in RPE1, 12–13%).
