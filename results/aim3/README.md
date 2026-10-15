# Aim 3 synthetic grid: run record

## Run

| item | value |
|---|---|
| Date | 2026-10-11 |
| Command | `PYTHONPATH=src python -m cdtest.simulate --reps 200 --B 999 --G 2000 --n-jobs 14` |
| Replicates per setting | **200, not the plan's 1,000** (see below) |
| Settings | 192 = ρ_M {0, .1, .3, .5} × ρ_x {0, .3, .6} × confounding {0, .3} × noise {homo, hetero} × θ {0, small 0.5%, moderate 2%, large 5% partial R²} |
| Genes per replicate | 2,000 |
| Permutations per test | B = 999 |
| Decoys | 99 |
| Tests per replicate | 11 |
| Seeds | master 20261011 → (setting index, replicate) via `SeedSequence` |
| Wall time | 997 s on 14 workers |

The tests are: naive OLS, CAMERA VIF, unrestricted / expression-bin / stratified
permutation, decoys, the main block test, and blocks at k/4, 4k, random order
and oracle (true modules).

**Why 200.** The 1,000-replicate run was started first, at 15:25. Other jobs
then loaded the machine (load average 90–160), and after 94 minutes only 23% of
the batches were done. It was stopped, and the 200-replicate run was done
instead.

The runner now checkpoints every 25-replicate batch to `results/aim3/_parts/`.
Replicate r of a setting always gets the same seed, so the 200 replicates here
are exactly the first 200 of the full run. To finish the plan's 1,000 (for
example on OSC), run:

```
PYTHONPATH=src python -m cdtest.simulate --reps 1000 --B 999 --G 2000 --n-jobs <cores>
```

It resumes from the existing parts. Uncontended, 14 cores need about 75 min.

## Files

| file | content | in git |
|---|---|---|
| `synthetic_grid.parquet` | one row per replicate: p-values (one- and two-sided), CI coverage per test, θ̂, partial R², k, grid CI of the main test | no (regenerate) |
| `_parts/` | checkpoints | no |
| `synthetic_grid_summary.csv` | per setting × test: FPR@5%/1%, power, power at matched FPR, bias, coverage, all with MC SE | yes |
| `fpr_vs_structure.csv` | data behind `figures/aim3_fpr_vs_structure.*` (θ = 0, noise models pooled, n = 400 per point) | yes |
| `headline_by_test.csv` | the table below | yes |
| `k_rule_selection*.{md,csv}`, `k_rule_selection_protocol.md` | plan v1.1 k-rule selection | yes |
| `k_rule_grid.parquet`, `_parts_krule/`, `k_rule_run.log` | its per-replicate output | no |

Read the summary with `pd.read_csv(..., keep_default_na=False)`. Effect level
"zero" is θ = 0.

## Headline (200 reps; 400 per point when the two noise models are pooled)

Column definitions:

- **"shared"**: cells with ρ_M > 0 and ρ_x > 0.
- **"no shared"**: cells with ρ_M = 0 or ρ_x = 0.
- **FPR**: rejection rate at θ = 0.
- **Power**: mean over the θ > 0 settings.
- **Coverage**: exact coverage of the equal-tailed 95% inversion interval.

Monte Carlo SE of a 5% rate: 1.1 points per cell (n = 400), 1.5 points
(n = 200).

| test | FPR@5% no shared | FPR@5% shared | FPR@5% worst cell | FPR@1% mean | FPR@1% worst | power@5%, R²=0.5% | power at matched FPR, R²=0.5% | same, R²=2% | coverage mean | coverage min |
|---|---|---|---|---|---|---|---|---|---|---|
| **block (main, k rule)** | 0.051 | **0.079** | 0.102 | 0.017 | 0.045 | 0.76 | 0.71 | 0.95 | 0.928 | 0.85 |
| block = true modules (oracle) | 0.053 | 0.051 | 0.075 | 0.013 | 0.030 | 0.71 | 0.69 | 0.97 | 0.952 | 0.91 |
| block 4k | 0.049 | 0.053 | 0.078 | 0.010 | 0.030 | 0.48 | 0.47 | 0.77 | 0.951 | 0.89 |
| block k/4 | 0.052 | 0.119 | 0.170 | 0.030 | 0.090 | 0.85 | 0.77 | 0.98 | 0.896 | 0.73 |
| block, random gene order | 0.049 | 0.149 | 0.250 | 0.042 | 0.175 | 0.88 | 0.78 | 0.99 | 0.866 | 0.57 |
| decoys (99) | 0.053 | 0.050 | 0.070 | 0.009 | 0.030 | 0.77 | 0.78 | 0.98 | 0.961 | 0.93 |
| stratified (expr × CN) | 0.049 | 0.148 | 0.238 | 0.042 | 0.160 | 0.88 | 0.78 | 0.99 | 0.864 | 0.57 |
| expression-bin permutation | 0.051 | 0.150 | 0.250 | 0.042 | 0.160 | 0.88 | 0.78 | 0.99 | 0.861 | 0.56 |
| unrestricted permutation | 0.054 | 0.154 | 0.258 | 0.047 | 0.180 | 0.89 | 0.79 | 0.99 | 0.854 | 0.55 |
| naive OLS | 0.050 | 0.151 | 0.248 | 0.043 | 0.160 | 0.88 | 0.77 | 0.98 | 0.862 | 0.55 |
| CAMERA VIF | 0.027 | 0.018 | 0.068 | 0.003 | 0.025 | 0.53 | 0.77 | 0.98 | 0.983 | 0.92 |

## What it shows

1. **The predicted failure condition holds.** OLS, unrestricted, expression-bin
   and stratified permutation are at nominal level whenever either the response
   or the feature has no module structure. They inflate only when both share
   modules, reaching 20–26% at 5% and up to 18% at 1%. Expression-bin and
   stratified permutation give no protection, because the shared structure lies
   within bins. This is H3's prediction.
2. **The main block test is much better but not exact.** With strong shared
   structure it reaches 7–10% at the 5% level (mean 7.9% over shared cells) and
   up to 4.5% at 1%. Its CI coverage is 92.8% on average, with a minimum of 85%.
3. **Block permutation itself is valid when blocks are right.** Oracle blocks
   (the true modules) and blocks of 4k hold 5% within MC error. The excess in
   the main test therefore comes from the plan's rule for k (autocorrelation
   < 0.05), which gives blocks only about 1.5–2.3× the module scale (mean k:
   5, 15, 23 for ρ_x = 0, .3, .6; module mean 10). Dependence across block
   boundaries is then lost. k/4 and random order fail, as designed.
4. **Decoys are valid** in this simulation, because the decoys are exchangeable
   with x by construction. This is the best case for that null.
5. **CAMERA-style VIF is conservative.** Its FPR is about 2% and its power at
   the nominal level is low, but its coverage is about 98%. With n = all genes
   and the true ρ̄, VIF ≈ 1 + ρ_M·(mean module size − 1), whether or not the
   feature shares the modules.
6. **Power.** At matched FPR, the main test loses little (0.71 vs 0.77 at
   R² = 0.5%; 0.95 vs 0.98 at 2%). 4k pays heavily (0.47 and 0.77). θ̂ is
   unbiased: |mean bias| < 0.007, within MC SE.

## Follow-up: plan v1.1 k rule (2026-10-11)

A targeted simulation compared 8 block-size rules over 4 module designs
(`k_rule_selection.md`, `k_rule_selection.csv`, `k_rule_selection_by_cell.csv`,
`figures/aim3_k_rule_selection.png`). No rule met the pre-set criterion. 4×k
had the best size (mean 5.4%, worst 7.0%) and was adopted as plan v1.1. v1.0
is kept as a sensitivity analysis.

`synthetic_grid*` keeps the v1.0 rule as `block_main`: `simulate.analyze_replicate`
pins `rule="v1.0"`.

`k_rule_grid.parquet` holds the per-replicate output and is not in git.
Regenerate it with `python -m cdtest.k_rule_sim run`.

## Implications (written for v1.0; see the follow-up above)

- The plan's k rule is mildly anti-conservative under strong shared local
  structure. In the worst simulated cell it doubles the 5% FPR. This is a
  pre-specified Aim 3 finding to report, together with the 4k / oracle
  contrast.
- If a change is made, it needs a dated plan amendment in ANALYSIS_PLAN.md, with the
  original and changed analyses both reported. Candidates are a larger k
  multiplier, or a lower autocorrelation threshold such as 0.01.
- How large the inflation is depends on the simulated module sizes and
  strengths, which are modelling choices. The planted-signal benchmark on real
  M is the more realistic check.
