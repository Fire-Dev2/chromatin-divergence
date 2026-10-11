# Block-size rule selection for plan v1.1 (simulation only)

Date 2026-10-11. Protocol (committed before the run, commit 2690811):
`results/aim3/k_rule_selection_protocol.md`. Code: `cdtest.k_rule_sim`
(`python -m cdtest.k_rule_sim run --n-jobs 8`, then `report`). No real-data
test was run.

## Design

- 64 null cells: ρ_M ∈ {0.3, 0.5} × ρ_x ∈ {0.3, 0.6} × confounding ∈ {0, 0.3} ×
  noise ∈ {homo, hetero} × 4 module designs. 1,000 replicates per null cell.
- Power: partial R² 0.5% and 2%, 300 replicates per cell (128 cells).
- 102,400 replicates in total. G = 2,000, B = 999, seeds
  `SeedSequence([20261011, 11, cell, replicate])`. Wall time 73 min on 8 workers.
- Module designs. Modules are contiguous within an arm and are cut at arm ends.
  The simulated arms average 45 genes, so the realized sizes are smaller than
  the nominal ones:

  | design | nominal | realized mean | gene-weighted mean | largest per replicate (mean / max) |
  |---|---|---|---|---|
  | geom5 | Geometric, mean 5 | 4.6 | 8.1 | 27 / 61 |
  | geom10 | Geometric, mean 10 | 8.4 | 15.2 | 45 / 109 |
  | geom25 | Geometric, mean 25 | 16.4 | 28.3 | 71 / 147 |
  | mixed | LogNormal, median 8, σ = 1, clipped to [1, 150] | 9.9 | 18.9 | 60 / 150 |

- Per replicate, every rule's k was computed from the x residuals only. The
  block test was run once per distinct k, so rules that pick the same k share
  the same draw. Oracle blocks (the true modules) are a reference only.
- Matched power: the threshold is the 5% lower quantile of the rule's 1,000 null
  p-values in the same cell. Power is averaged over the 64 cells. The SE shown is
  binomial and ignores error in the threshold.
- CI coverage: exact coverage of the equal-tailed 95% inversion interval,
  averaged over null and power replicates. Min is the worst cell × effect.

## Result per rule

| rule | FPR@5% mean (SE) | FPR@5% worst (SE) | cells > 0.065 | FPR@1% mean (SE) | FPR@1% worst (SE) | CI coverage mean / min | matched power R²=0.5% (SE) | matched power R²=2% (SE) | eligible |
|---|---|---|---|---|---|---|---|---|---|
| v1.0 | 0.0791 (0.0011) | 0.115 (0.010) | 55/64 | 0.0228 (0.0006) | 0.037 (0.006) | 0.909 / 0.843 | 0.449 (0.003) | 0.813 (0.002) | no |
| acf01 | 0.0711 (0.0010) | 0.107 (0.010) | 47/64 | 0.0190 (0.0005) | 0.034 (0.006) | 0.920 / 0.870 | 0.409 (0.003) | 0.757 (0.003) | no |
| x2 | 0.0625 (0.0010) | 0.085 (0.009) | 22/64 | 0.0148 (0.0005) | 0.030 (0.005) | 0.934 / 0.893 | 0.345 (0.003) | 0.680 (0.003) | no |
| x3 | 0.0561 (0.0009) | 0.073 (0.008) | 9/64 | 0.0116 (0.0004) | 0.021 (0.005) | 0.942 / 0.880 | 0.289 (0.003) | 0.604 (0.003) | no |
| x4 | 0.0538 (0.0009) | 0.070 (0.008) | 3/64 | 0.0110 (0.0004) | 0.022 (0.005) | 0.945 / 0.897 | 0.257 (0.003) | 0.560 (0.003) | no |
| tau1 | 0.1063 (0.0012) | 0.148 (0.011) | 64/64 | 0.0369 (0.0007) | 0.059 (0.007) | 0.865 / 0.790 | 0.523 (0.003) | 0.894 (0.002) | no |
| tau2 | 0.0830 (0.0011) | 0.115 (0.010) | 58/64 | 0.0249 (0.0006) | 0.048 (0.007) | 0.903 / 0.843 | 0.470 (0.003) | 0.828 (0.002) | no |
| pw | 0.0572 (0.0009) | 0.077 (0.008) | 11/64 | 0.0120 (0.0004) | 0.023 (0.005) | 0.941 / 0.893 | 0.291 (0.003) | 0.602 (0.003) | no |
| oracle | 0.0499 (0.0009) | 0.065 (0.008) | 0/64 | 0.0104 (0.0004) | 0.017 (0.004) | 0.950 / 0.920 | 0.487 (0.003) | 0.888 (0.002) | ref |

Typical k chosen (median [IQR] over all replicates of a design):

| rule | geom5 | geom10 | geom25 | mixed |
|---|---|---|---|---|
| v1.0 | 9 [8–11] | 17 [14–22] | 36 [27–49] | 22 [17–30] |
| acf01 | 13 [10–16] | 23 [18–30] | 44 [34–59] | 30 [23–41] |
| x2 | 18 [16–22] | 34 [28–44] | 72 [54–98] | 44 [34–60] |
| x3 | 27 [24–33] | 51 [42–66] | 100 [81–100] | 66 [51–90] |
| x4 | 36 [32–44] | 68 [56–88] | 100 [100–100] | 88 [68–100] |
| tau1 | 5 [5–6] | 8 [6–11] | 17 [11–25] | 10 [7–14] |
| tau2 | 9 [6–11] | 15 [11–21] | 34 [21–50] | 19 [13–27] |
| pw | 28 [24–33] | 47 [41–57] | 80 [70–87] | 58 [48–74] |

FPR at α = 0.05 by module design (mean / worst of the 16 null cells):

| rule | geom5 | geom10 | geom25 | mixed |
|---|---|---|---|---|
| v1.0 | 0.088 / 0.115 | 0.082 / 0.106 | 0.067 / 0.087 | 0.080 / 0.100 |
| acf01 | 0.080 / 0.107 | 0.073 / 0.093 | 0.060 / 0.072 | 0.071 / 0.087 |
| x2 | 0.073 / 0.082 | 0.063 / 0.085 | 0.055 / 0.062 | 0.059 / 0.082 |
| x3 | 0.063 / 0.073 | 0.055 / 0.069 | 0.053 / 0.063 | 0.053 / 0.065 |
| x4 | 0.058 / 0.070 | 0.054 / 0.068 | 0.052 / 0.062 | 0.051 / 0.063 |
| tau1 | 0.096 / 0.132 | 0.110 / 0.148 | 0.106 / 0.140 | 0.114 / 0.135 |
| tau2 | 0.087 / 0.105 | 0.089 / 0.115 | 0.071 / 0.092 | 0.086 / 0.108 |
| pw | 0.066 / 0.077 | 0.056 / 0.070 | 0.053 / 0.060 | 0.053 / 0.061 |
| oracle | 0.050 / 0.056 | 0.050 / 0.065 | 0.052 / 0.060 | 0.048 / 0.055 |

Figure: `figures/aim3_k_rule_selection.png`. Per-cell numbers:
`k_rule_selection_by_cell.csv`. Per-rule numbers: `k_rule_selection.csv`.

## Applying the pre-set criterion

**No rule meets it.** Every candidate's worst-case FPR at 5% is above 0.065.
Closest:

- x4: 0.070 (SE 0.008), 3/64 cells above 0.065. Its worst FPR at 1% is 0.022,
  also above the 0.02 limit.
- x3: 0.073.
- pw: 0.077.

By the protocol's no-eligible-rule clause, the rule with the lowest worst-case
FPR at 5% is **x4** (4 × the v1.0 k, bounded to [5, 100]).

## Reading the result

1. **v1.0 is anti-conservative across all module designs.** Its mean FPR at 5%
   is 7.9% and the worst cell is 11.5%. CI coverage averages 0.91. This
   reproduces the main grid's finding (7.9% / 10.2%) and extends it to smaller,
   larger and mixed modules. Small modules (geom5) are the hardest case for
   every ACF-based rule.
2. **Integrated-autocorrelation rules are worse than v1.0.** tau1 and tau2 pick
   smaller blocks: τ is small when ρ_x is moderate. Feature-side τ measures how
   much x is smoothed, not how far the shared structure reaches.
3. **Size improves steadily with k.** Mean FPR goes v1.0 7.9% → acf01 7.1% →
   x2 6.3% → x3 5.6% → x4 5.4% (SE about 0.1 points). x4 is still slightly
   anti-conservative on average (z ≈ 4 against 0.05), mostly in geom5 (mean
   5.8%). Politis–White sits between x2 and x3 in size, with x3's power.
4. **The 0.065 worst-case threshold is close to Monte Carlo noise at n = 1,000.**
   Take an exactly valid test. Each cell exceeds 0.065 with probability 0.015,
   so the test passes all 64 cells only with probability 0.38. The expected
   maximum over 64 cells is 0.066 (90% range 0.062–0.073).
   - The oracle reference lands exactly on the limit (0.065).
   - x4's worst cell (0.070) is within the range expected for an exact test.
     Its mean (5.4%) is not, so it is close to nominal but not exact.
   - This is a weakness of the pre-set criterion. It was not changed after the
     results were seen.
5. **Size is bought with power.** At partial R² 0.5%, matched power is 0.449
   for v1.0, 0.289 for x3 and 0.257 for x4. The oracle gets 0.487. Part of this
   loss comes from G = 2,000:
   - With 44 arms of about 45 genes, k ≥ 68 makes almost every arm a single
     block. Blocks are exchanged only among equal sizes, so only about 20 blocks
     can move. At k = 36 about 30 can move.
   - At a realistic G of about 8,000 (≈180 genes per arm), about 34 blocks can
     move at k = 100 and 73 at k = 68.
   - So the simulated power of x3, x4 and pw is likely pessimistic for the real
     test. This is not verified here.

## Recommendation (adopted as plan v1.1, flagged for author/mentor sign-off)

- **Main test: x4.** k = 4 × the v1.0 k, bounded to [5, 100].
  - It has the best size of the candidates (mean 5.4%, worst 7.0%; at 1%, mean
    1.1% and worst 2.2%) and the best coverage (mean 0.945).
  - It is the protocol's fallback choice.
  - It is simple: one multiplier on the frozen rule.
- **Keep v1.0 as a reported sensitivity analysis.**
- **Before any confirmatory real-data test,** the planted-signal benchmark on
  real M (plan §6) should report the main test under both rules. Power with x4
  is the main cost, and the real G per pair is larger than simulated.
- **Not done here: a rule that meets the criterion.** None of the
  feature-side-only rules tested meets it. Exact size in these simulations needs
  blocks matched to the response's shared structure (the oracle), which the
  feature alone does not reveal.
