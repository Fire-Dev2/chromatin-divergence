# H1 (primary) — matched pairs, combined CPM feature, main test

Run once on the frozen pipeline (tag `pipeline-v1.1`, plan v1.3): block permutation, x4 rule, chromosome-arm terms in f(C), B = 9,999, one-sided θ > 0, Holm across the two matched pairs. 95% CI by test inversion (refined grid inside the plan span; plan-grid interval also given). Partial R² / reliability with 500-resample knockdown bootstrap. Detection limit from the planted-signal benchmark on each pair (`planted_<pair>.md`).

| pair (A − B) | θ̂ | k | one-sided p (θ > 0) | two-sided p [equal-tailed] | Holm p | 95% CI (refined) | plan-grid CI | partial R²/rel [boot 95%] | detection limit (all / excl. sig. θ₀ decoys) |
|---|---|---|---|---|---|---|---|---|---|
| K562_vs_Jurkat | 0.0001979 | 20 | 0.6447 | 0.8381 [0.7108] | 1 | [-0.002123, 0.001367] | [0, 0] | 0.001% [0.000%, 0.023%] | 0.5% / 0.5% |
| RPE1_vs_HepG2 | -0.01734 | 28 | 1 | 0.0001 [0.0002] | 1 | [-0.02625, -0.01404] | [-0.0159, -0.0159] | 0.770% [0.607%, 0.905%] | 2% / 0.5% |

**Verdict: H1 not supported: no matched pair shows a positive effect at 5% (Holm) under the main test.**

Proposal 'Supported if': the effect of the signed accessibility difference on the signed response difference is positive beyond expression, at 5% under the main test. Neither matched pair meets it.

- **K562 − Jurkat:** θ̂ ≈ 0 (2e-4; 95% CI [−0.0021, 0.0014]; partial R²/rel 0.001%, bootstrap upper bound 0.023%). The planted benchmark on this pair detects 0.5% partial R² with ≥ 80% power (level-0 FPR 0.08 all decoys, 0.00 without the two decoys whose own θ₀ is significant), so an effect of the size the study can see is absent. 'If not' column: accessibility adds nothing detectable here once p53, tissue and growth mode are held constant.
- **RPE1 − HepG2:** a clear association of the *opposite* sign (θ̂ = −0.017, CI [−0.026, −0.014], two-sided p 1e-4 at the B = 9,999 floor, partial R²/rel 0.77%): genes more accessible in one line respond *less* in that line. The one-sided H1 test cannot reject (p = 1). The p53 control shows RPE1 − HepG2 still diverges at p53 targets (z = 3.3), so p53 activity is not held constant in this pair. The planted benchmark here is contaminated by decoys with their own θ₀ (4 of 6 significant): detection limit 2% with all decoys, 0.5% with the 2 clean decoys (4 configurations only).
- Two-sided p in the frozen code is (1 + #|T*| ≥ |T|)/(1 + B); the equal-tailed version 2·min(p>, p<) is given in brackets. They differ when the block null is not centred at 0.
