# Aim 4 — summary (frozen pipeline `pipeline-v1.1`, plan v1.3; H1 test run once)

Date: 2026-10-12. All analyses: `scripts/run_aim4_aim1.sh`, `scripts/run_aim4.sh` (`scripts/run_aim4.py`), `scripts/report_aim4.py`; new code `src/cdtest/aim4.py` (pooled contrast, H4, p53 control, target-promoter check); every test, rule and parameter as frozen (B = 9,999; x4 block rule; arm terms in f(C); CPM features; 500 bootstrap resamples; planted 500 reps/level, B = 999).

## Pairs (Aim 1; M = A − B)

| pair | role | genes | knockdowns kept | reliability | generic share A / B | SD M |
|---|---|---|---|---|---|---|
| K562_vs_Jurkat | matched | 5,675 | 1,414 / 2,055 | 0.863 | 2.4% / 3.2% | 0.0166 |
| RPE1_vs_HepG2 | matched | 5,631 | 1,424 / 2,392 | 0.978 | 11.9% / 4.3% | 0.0619 |
| RPE1_vs_K562 | mixed | 5,262 | 1,314 / 2,055 | 0.996 | 13.0% / 2.7% | 0.1004 |
| K562_vs_HepG2 | mixed | 5,785 | 1,389 / 2,055 | 0.982 | 2.7% / 4.8% | 0.0597 |
| Jurkat_vs_RPE1 | mixed | 5,222 | 1,413 / 2,392 | 0.995 | 3.1% / 12.1% | 0.0963 |
| Jurkat_vs_HepG2 | mixed | 5,587 | 1,482 / 2,392 | 0.977 | 3.1% / 4.2% | 0.0535 |

## Headline table (combined CPM feature, main test)

| pair | role | θ̂ | 95% CI | p (θ > 0) | p two-sided | Holm p | pR²/rel [boot] | detection limit | p53 control z (oriented) |
|---|---|---|---|---|---|---|---|---|---|
| K562_vs_Jurkat | matched | 0.0002 | [-0.00212, 0.00137] | 0.6447 | 0.8381 | 1 | 0.001% [0.000%, 0.023%] | 0.5% / 0.5% | -0.7 |
| RPE1_vs_HepG2 | matched | -0.0173 | [-0.0263, -0.014] | 1 | 0.0001 | 1 | 0.770% [0.607%, 0.905%] | 2% / 0.5% | 3.3 |
| RPE1_vs_K562 | mixed | -0.0224 | [-0.0315, -0.0135] | 1 | 0.0001 | — | 0.545% [0.499%, 0.594%] | — | 4.1 |
| K562_vs_HepG2 | mixed | 0.0071 | [-0.00129, 0.0123] | 0.0551 | 0.0551 | — | 0.094% [0.069%, 0.118%] | — | 6.9 |
| Jurkat_vs_RPE1 | mixed | -0.0219 | [-0.0303, -0.0104] | 0.9999 | 0.0002 | — | 0.511% [0.465%, 0.565%] | — | 5.3 |
| Jurkat_vs_HepG2 | mixed | 0.0024 | [-0.00248, 0.00815] | 0.1369 | 0.2443 | — | 0.018% [0.009%, 0.033%] | — | 6.1 |
| pooled mixed (4 pairs, 4,688 common genes, one shared block permutation) | mixed | -0.0096 | [-0.0187, -0.00702] | 1 | 0.0001 | — | — | — | — |
| pooled matched (same genes) | matched | -0.0055 | [-0.0104, -0.00414] | 1 | 0.0001 | — | — | — | — |

Detection limit = smallest planted partial R² recovered with 80% power at nominal 5% (all decoys / excluding decoys whose own θ₀ is significant). p53 z oriented to the p53-WT line (matched: A − B).

## Hypotheses

- **H1 (primary): not supported** — no matched pair shows a positive effect at 5% (Holm) under the main test. K562 − Jurkat: θ̂ ≈ 0 with adequate power (detection limit 0.5%). RPE1 − HepG2: θ̂ = -0.0173, opposite sign, two-sided p 1e-4. Holm p = 1 for both.
- **H2 (secondary): not supported.** Distal is not larger with non-overlapping intervals in either matched pair; in RPE1 − HepG2 the promoter carries more (0.59% vs 0.23%).
- **H3: simulation half supported (Aim 3, planted naive FPR 0.27–0.43), real-data half not supported** — the naive step-7 estimate is smaller in magnitude than the main estimate in 4/6 pairs and flips sign in K562 − HepG2; naive and main tests agree except borderline K562 − HepG2.
- **H4 (secondary): not supported as stated.** Pooled slope -0.0159 (p 4.5e-05), intercept 0.0009 (p 0.4); per-pair slopes are negative where the own θ is negative and ≈ 0 where it is ≈ 0 — the gradient tracks the own effect (consistent with the effect model, but with negative θ). Decoy negative control: 17 of 43 low-similarity decoys reject at 5% — low similarity does not imply no association.
- **Mixed vs matched:** mixed pooled -0.0096 vs matched pooled -0.0055; gap -0.0042 (approx. 95% [-0.00863, 0.000325], p 0.056).
- **p53 control:** all four mixed pairs diverge (z 4.1–7.0, oriented as expected); K562 − Jurkat none (z −0.7); **RPE1 − HepG2 diverges (z 3.3)** although both are p53 WT.
- **Sensitivity (H1):** v1.2, FC features, v1.0 k rule, stratified, decoy null, co-accessibility blocks and group-mean M: no variant gives a positive effect in either matched pair (`sensitivity_h1.md`). The RPE1 − HepG2 negative association holds under every CPM variant but not with FC features (θ̂ −0.004, p₂ 0.27), and is weaker with co-accessibility blocks (p₂ 0.04) and the 6-decoy null (p₂ 0.14). Target-promoter knockdown-level check: no positive association (ρ −0.05 to 0.00).

## Outcome classification (proposal 'Interpreting outcomes')

| outcome | supported? | evidence |
|---|---|---|
| Positive effect in both matched pairs, rising with similarity | no | K562 − Jurkat θ̂ ≈ 0; RPE1 − HepG2 θ̂ < 0; H4 slopes ≤ 0 where the own effect is nonzero. |
| Effect in the mixed pairs only | no (as stated) | Mixed pooled θ̂ -0.0096 (CI [-0.0187, -0.00702]) is nonzero but negative, and a matched pair (RPE1 − HepG2) shows the same negative sign; gap mixed − matched -0.0042 (approx. CI [-0.00863, 0.000325], p 0.056). |
| The generic response accounts for most of the raw difference | no | Median generic share 2–13% per line in every pair. |
| Naive tests find an effect; the main test does not | no (one borderline mixed pair) | Naive and main agree in sign and significance in 5 of 6 pairs; K562 − HepG2 (mixed): permutation-free naive tests p₂ 0.008–0.022 vs main 0.055. The naive step-7 estimate is not larger than the main one. |
| No effect under any test, with adequate power | yes, for K562 − Jurkat | θ̂ ≈ 0 under every test and sensitivity; planted detection limit 0.5% partial R²; pR²/rel upper bound 0.023%. |
| Cross-line numbers sit near the same-line baseline | partly, for K562 − Jurkat | K562 − Jurkat reproducible variance of M is 3.6× the same-line baseline (other pairs 43–153×); reliability 0.86. |
| The two matched pairs disagree | YES — lead | K562 − Jurkat null; RPE1 − HepG2 significantly negative (opposite to H1). What differs: RPE1 − HepG2 still diverges at p53 targets (z 3.3) and every RPE1-containing pair shows the same promoter-driven negative association, while the three pairs without RPE1 do not. |
| Signal in distal accessibility only (secondary) | no (weak hint only) | Distal θ̂ weakly positive in the pairs without RPE1 (K562 − Jurkat one-sided p 0.021) but H2's non-overlap criterion fails in both matched pairs. |

**The report should lead with:** the pair-by-pair comparison — *the two matched pairs disagree*. K562 − Jurkat gives a null with adequate power (state its detection limit, 0.5% partial R²); RPE1 − HepG2 gives a significant association of the opposite sign to H1. Then 'what differs': the negative, promoter-carried association appears in every pair containing RPE1 (RPE1 − K562, Jurkat − RPE1, RPE1 − HepG2) and in none without it, and RPE1 − HepG2 is not p53-matched in practice. Candidate explanations to state, not test: RPE1-specific feature/CN processing (local reprocessing of 4DN ATAC, karyotype-only CN track), RPE1 biology (strong p53 response), cross-study batch. H1 as pre-specified is not supported.

## Files

- `aim1_pairs.{json,md}`, `h1_primary.{json,md}`, `planted_<pair>.{json,md}`, `h2_location.{json,md}`, `h3_comparators.{json,md}`, `mixed_pooled.{json,md}`, `h4_gradient.{json,md}`, `p53_control.{json,md}`, `sensitivity_h1.{json,md}`, `outcome_classification.json`; raw stage outputs in `tests/`, `sens_v12/`, `sens_groupmean/`, `planted/`, `extras/`, logs in `logs/`.
- Figures: `figures/aim4_h1_forest.png`, `aim4_h2_promoter_distal.png`, `aim4_h3_naive_vs_main.png`, `aim4_similarity_gradient.png`, `aim4_p53_control.png`.
