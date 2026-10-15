# API reference

The public interface is the test engine in `cdtest.crt`. Everything else in `src/cdtest/` is
the analysis that applies it to Perturb-seq and ATAC-seq data (see the repository map in the
README). All functions take NumPy arrays (one entry per gene) unless stated otherwise.

## Main test

### `crt.block_permutation_test(M, x, Z=None, chrom=None, pos=None, arm=None, k=None, B=9999, rng=None, k_rule=None, prep=None, order=None, ...) -> TestResult`

Feature-side block permutation. Fits `M ~ f(C) + θ·x` by OLS (`Z` = design matrix of f(C)),
residualizes `x` on f(C), orders genes by `(chrom, pos)`, cuts the residuals into blocks of `k`
consecutive genes within chromosome arms, and permutes blocks of equal size `B` times. `M` is
never permuted, and the block size is chosen from the feature only.

- `k`: block size; if `None`, chosen by `k_rule` (`"x4"`, the rule used in the paper: four times
  the integrated autocorrelation time of the x residuals, bounded below by 5).
- `rng`: seed or `numpy.random.Generator`.
- `prep`, `order`: reuse a `crt.prepare(...)` result and a `crt.genome_order(...)` ordering
  across calls (faster in loops).

### `TestResult`

| attribute / method | meaning |
|---|---|
| `T` | effect estimate θ̂ |
| `p` | one-sided randomization p for θ > 0, (1 + #{T* ≥ T}) / (1 + B) |
| `p_two_sided` | two-sided p |
| `Tstar` | the B null statistics |
| `ci(grid=None, alpha=0.05, sd_x=None)` | test-inversion confidence interval `(lo, hi, info)` over a θ grid (default: 201 points over ± max(4·abs(θ̂), 0.5 per SD of x)) |
| `p_at(theta0, alternative="greater")` | p-value for H0: θ = θ0 (vectorized) |
| `covers(theta0, alpha=0.05)` | whether θ0 is inside the equal-tailed interval |
| `info` | dict with `k`, block counts and other diagnostics |

## Design and helpers

| function | purpose |
|---|---|
| `crt.build_design(C, spline_cols=(), df=4, intercept=True)` | f(C): natural cubic spline (df 4) on `spline_cols`, other columns linear; `C` is a DataFrame |
| `covariates.arm_indicators(arm)` | one 0/1 column per chromosome arm (first arm dropped); append to the design to absorb arm-scale copy-number error, as in the paper |
| `crt.prepare(M, x, Z=None, Q=None)` | residualized M and x, θ̂ and the orthonormal basis Q of the design |
| `crt.choose_k(x_resid_in_genome_order, groups=arm, rule="x4")` | block size from the feature residuals only |
| `crt.genome_order(chrom, pos)`, `crt.make_blocks(arm_sorted, k)` | ordering and block ids |
| `crt.partial_r2(M, x, Z=None, reliability=1.0)` | partial R² of x given f(C), optionally divided by the outcome's split-half reliability |
| `crt.fit_effect_model(M, x, Z)` | OLS fit with naive SE, p and partial R² |
| `crt.decoy_similarity(x, decoys, Z=None)` | correlation of decoy features with x after residualization on f(C) |

## Comparators (used in the benchmarks)

| function | null |
|---|---|
| `crt.naive_ols_test` | OLS t-test, genes independent |
| `crt.unrestricted_permutation_test` | free permutation of x |
| `crt.expression_bin_permutation_test` | permutation within expression-difference bins |
| `crt.stratified_permutation_test` | permutation within expression × copy-number strata |
| `crt.camera_vif_test` | CAMERA-style variance-inflated OLS |
| `crt.decoy_test` | decoy features as the null distribution |

## Benchmarks

- `cdtest.simulate`: synthetic genomes with shared module structure (`python -m cdtest.simulate`).
- `cdtest.k_rule_sim`: block-size rule selection (`python -m cdtest.k_rule_sim`).
- `cdtest.planted`: planted signal in a real outcome (`planted.run(a, b, reps, B)`).

## Minimal example

```python
import numpy as np
from cdtest import crt
from cdtest.covariates import arm_indicators

# M, x: arrays per gene; C: covariate DataFrame; genes: DataFrame with chrom, pos, arm
Z = np.column_stack([crt.build_design(C, spline_cols=("expr_A", "expr_B"), df=4),
                     arm_indicators(genes["arm"].to_numpy())[0]])
res = crt.block_permutation_test(M, x, Z=Z, chrom=genes["chrom"].to_numpy(), pos=genes["pos"].to_numpy(),
                                 arm=genes["arm"].to_numpy(), B=9999, k_rule="x4", rng=20261011)
lo, hi, _ = res.ci(sd_x=float(np.std(x)))
print(res.T, res.p_two_sided, (lo, hi), res.info["k"])
```
