# cdtest: dependence-aware feature-randomization test

`cdtest` tests whether a gene-level feature (for example, chromatin accessibility) explains a
gene-level response (for example, how differently the same CRISPRi knockdown changes a gene's
expression in two cell lines) while staying valid when genes are not independent.

Neighbouring genes share regulatory elements and copy-number segments, so gene-level tests that
treat genes as independent reject far too often. `cdtest` keeps the response fixed and permutes
the **feature** in genome-ordered blocks (a conditional-randomization-style test), with
covariate adjustment, a feature-only block-size rule, test-inversion confidence intervals and
planted-signal power benchmarks.

This repository also contains the complete, pre-registered analysis that applies `cdtest` to
CRISPRi Perturb-seq screens in K562, RPE1, Jurkat and HepG2 with bulk ATAC-seq features.

## Install

Requires Python ≥ 3.10.

```bash
git clone https://github.com/Fire-Dev2/chromatin-divergence.git
cd chromatin-divergence
python -m venv env/venv && source env/venv/bin/activate
pip install -e ".[test]"          # core test + simulations
pip install -e ".[full,test]"     # also h5py, pyBigWig, pyfaidx for the real-data pipeline
```

## Quick start (about 5 seconds, no downloads)

```bash
python examples/quickstart.py
```

It simulates 2,000 genes whose response and feature share correlated gene modules with **no
true effect**, compares a naive gene-level OLS test with the `cdtest` block-permutation test,
then plants a real effect. Expected output (seed 20261011):

```
Null simulation (100 replicates, 2,000 genes, true effect = 0)
  false-positive rate at 5%   naive OLS: 39.0%   cdtest block permutation: 7.0%

Planted effect (100 replicates, true theta = 0.229)
  cdtest power at 5%: 62.0%   mean estimate: 0.225
  last replicate: estimate 0.189, 95% CI [0.091, 0.423], block size k = 100
```

With 100 replicates the Monte Carlo SE of a 5% rate is about 2.2 points. The full calibration
(192 settings × 200 replicates, and 64 block-size cells × 1,000 replicates) is reproduced with
`python -m cdtest.simulate` and `python -m cdtest.k_rule_sim`.

## Run the tests

```bash
pytest
```

69 tests pass and 1 is skipped on a fresh clone; the skipped test checks the frozen
development-pair numbers and needs the regenerated Aim 1 intermediates.

## Use it on your own data

```python
import numpy as np
from cdtest import crt

# M: response per gene; x: feature per gene; C: covariate DataFrame (one row per gene)
# genes: DataFrame with chrom, pos (TSS) and arm (chromosome arm) per gene
Z = crt.build_design(C, spline_cols=("expr_A", "expr_B"), df=4)   # f(C)
res = crt.block_permutation_test(M, x, Z=Z, chrom=genes["chrom"], pos=genes["pos"],
                                 arm=genes["arm"], B=9999, k_rule="x4", rng=20261011)
lo, hi, _ = res.ci(sd_x=float(np.std(x)))
print(res.T, res.p_two_sided, (lo, hi), res.info["k"])
```

`res.T` is the effect estimate (θ̂), `res.p` the one-sided p for θ > 0, `res.p_two_sided` the
two-sided p, and `res.ci()` a test-inversion 95% interval. As in the paper (plan v1.3), add one
indicator column per chromosome arm to `C` to absorb arm-scale copy-number error. Comparators
(`naive_ols_test`, `unrestricted_permutation_test`, `stratified_permutation_test`,
`decoy_test`, `camera_vif_test`) live in the same module. Method details:
[`docs/methods_aim2_aim3.md`](docs/methods_aim2_aim3.md).

## Reproduce the paper

All input data are public. Large files are not in git; they are downloaded into `data/`
(about 100 GB in total). The shell scripts expect the environment at `env/venv`.

| Step | Command / document |
|---|---|
| 1. Download screens, ATAC-seq, annotations | `python data/fetch_all.py` (manifests and accessions in `recon/`) |
| 2. Gene annotations, copy number, p53 sets | `python scripts/build_annotations.py` ([`docs/annotations.md`](docs/annotations.md)) |
| 3. Local ATAC pipeline (Jurkat, RPE1; K562 validation) | `scripts/atac_pipeline.py`, `scripts/run_rpe1_sub60M.sh` ([`docs/atac_local_pipeline.md`](docs/atac_local_pipeline.md)) |
| 4. Accessibility features | `scripts/build_features.sh` ([`docs/features.md`](docs/features.md)) |
| 5. Outcome M | `scripts/run_aim4_aim1.sh` ([`docs/aim1_outcome.md`](docs/aim1_outcome.md)) |
| 6. Main test, controls, sensitivities | `scripts/run_aim4.sh`, `scripts/run_same_line_control.py`, `scripts/aim4x_rpe1_diagnostic.py` |
| 7. Figures | `python scripts/make_figures.py` (reads `results/`; needs step 2 for the gene table) |

Result tables and figures from the published run are already in `results/` and `figures/final/`.

## Pre-registration and record

- [`ANALYSIS_PLAN.md`](ANALYSIS_PLAN.md): the frozen analysis plan (git tag `plan-v1.0`;
  amendments tagged `plan-v1.1` and `plan-v1.3`).
- `pipeline-v1.0`, `pipeline-v1.1`: the frozen pipeline tags used for the matched pairs.
- Every deviation from the plan is listed with its date and reason in the paper;
  commit timestamps in this repository record when each step was made.

## Repository layout

| Folder | Contents | In git? |
|---|---|---|
| `src/cdtest/` | The package (`crt.py` is the test; `simulate.py`, `k_rule_sim.py`, `planted.py` the benchmarks) | yes |
| `examples/` | Quick-start | yes |
| `tests/` | Unit tests, including the randomization-validity check | yes |
| `scripts/` | Pipeline, analysis and figure scripts | yes |
| `results/`, `figures/` | Outputs (small files) | yes |
| `docs/`, `recon/` | Method notes; data inventories and accessions | yes |
| `data/` | Downloaded inputs (~100 GB) | no |
| `env/` | Python environment | no |

## License

MIT; see [`LICENSE`](LICENSE).
