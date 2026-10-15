# cdtest: dependence-aware feature-randomization test

[![tests](https://github.com/Fire-Dev2/chromatin-divergence/actions/workflows/tests.yml/badge.svg)](https://github.com/Fire-Dev2/chromatin-divergence/actions/workflows/tests.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

`cdtest` tests whether a gene-level feature (for example, chromatin accessibility) explains a
gene-level response (for example, how differently the same CRISPRi knockdown changes a gene's
expression in two cell lines) while staying valid when genes are not independent.

Neighbouring genes share regulatory elements and copy-number segments, so gene-level tests that
treat genes as independent reject far too often. `cdtest` keeps the response fixed and permutes
the **feature** in genome-ordered blocks (a conditional-randomization-style test), with
covariate adjustment, a feature-only block-size rule, test-inversion confidence intervals and
planted-signal power benchmarks.

This repository also contains the complete analysis, run under a frozen analysis plan, that
applies `cdtest` to CRISPRi Perturb-seq screens in K562, RPE1, Jurkat and HepG2 with bulk
ATAC-seq features.

## Install

Requires Python ≥ 3.10.

```bash
pip install "git+https://github.com/Fire-Dev2/chromatin-divergence.git@v0.2.0"   # the package only
```

or, to run the tests, the quick-start and the paper's analyses:

```bash
git clone https://github.com/Fire-Dev2/chromatin-divergence.git
cd chromatin-divergence
python -m venv env/venv && source env/venv/bin/activate
pip install -e ".[full,test]"     # full = h5py, pyBigWig, pyfaidx for the real-data pipeline
```

The command-line tools of the local ATAC-seq pipeline (fastp, Bowtie2, SAMtools, BEDTools,
MACS3, IDR) are listed in [`environment-atac.yml`](environment-atac.yml)
(`conda env create -f environment-atac.yml`); the [`Dockerfile`](Dockerfile) builds a container
with them and cdtest. They are needed only to rebuild the accessibility features from reads.

## Quick start (no downloads)

```bash
python examples/quickstart.py
```

It simulates 2,000 genes whose response and feature share correlated gene modules with **no
true effect** (response module correlation 0.5, feature module sharing 0.6: the most strongly
shared setting of the paper's benchmark), compares a naive gene-level OLS test with the `cdtest`
block-permutation test (x4 block-size rule, B = 999), then plants a real effect. About 20–40 s.
Expected output (seed 20261011):

```
Null simulation (500 replicates, 2,000 genes, true effect = 0)
  false-positive rate at 5%   naive OLS: 37.8%   cdtest block permutation: 4.4%
  Monte Carlo SE of a 5% rate: 1.0%; block size k: median 84 (range 44-100)

Planted effect (500 replicates, true theta = 0.229)
  cdtest power at 5%: 70.6%   mean estimate: 0.230
  last replicate: estimate 0.225, 95% CI [-0.270, 0.379], block size k = 76
```

Across the full benchmark (192 settings × 200 replicates) the x4 rule held 4.9–5.3% on average
and 7.8% in the worst cell; it is reproduced with `python -m cdtest.simulate` and
`python -m cdtest.k_rule_sim`. The simulated genome has no arm-scale artefact (arm-level copy
number is a covariate), so the demo uses no chromosome-arm terms; `--arm-terms` adds them. Without
arm-scale structure to absorb, arm terms shorten the blocks the x4 rule chooses, and in this
setting the false-positive rate rises to about 9%. Arm terms are meant for real data whose
residual feature carries arm-scale structure, as in the paper; their calibration on the real
outcomes is checked in `results/calibration/`.

## Run the tests

```bash
pytest
```

On a fresh clone all tests run except the annotation sanity checks, which need the full
annotation build (step 2 below). `tests/test_intermediates.py` recomputes the primary (H1)
numbers of both matched pairs from the deposited intermediates and checks them against the
committed results, including the B = 9,999 permutation p-values.

## Use it on your own data

```python
import numpy as np
from cdtest import crt
from cdtest.covariates import arm_indicators

# M: response per gene; x: feature per gene; C: covariate DataFrame (one row per gene)
# genes: DataFrame with chrom, pos (TSS) and arm (chromosome arm) per gene
Z = np.column_stack([crt.build_design(C, spline_cols=("expr_A", "expr_B"), df=4),   # f(C)
                     arm_indicators(genes["arm"].to_numpy())[0]])                    # arm terms
res = crt.block_permutation_test(M, x, Z=Z, chrom=genes["chrom"].to_numpy(), pos=genes["pos"].to_numpy(),
                                 arm=genes["arm"].to_numpy(), B=9999, k_rule="x4", rng=20261011)
lo, hi, _ = res.ci(sd_x=float(np.std(x)))
print(res.T, res.p_two_sided, (lo, hi), res.info["k"])
```

`res.T` is the effect estimate (θ̂), `res.p` the one-sided p for θ > 0, `res.p_two_sided` the
two-sided p, and `res.ci()` a test-inversion 95% interval. Comparators (`naive_ols_test`,
`unrestricted_permutation_test`, `stratified_permutation_test`, `decoy_test`,
`camera_vif_test`) live in the same module. API reference: [`docs/api.md`](docs/api.md);
method details: [`docs/methods_aim2_aim3.md`](docs/methods_aim2_aim3.md).

## Reproduce the paper

### From the deposited intermediates (under a minute, no large downloads)

```bash
python scripts/reproduce_from_intermediates.py
```

[`intermediates/`](intermediates/) holds the processed inputs of every confirmatory analysis
(~10 MB): gene annotation, copy number, gene-level accessibility features for the four screened
and seven decoy lines, and control-cell expression moments for the covariates; the outcome
vectors M are in `results/aim1/<pair>/`. The script recomputes the main test in all six pairs,
Holm p, the pooled contrasts and the same-line control with the frozen code and seeds, writes
Table 2 (`results_repro/table2.md`) and Figures 3–4 (`results_repro/figures/`), and compares
every number with the committed results (`results_repro/COMPARISON.md`). The same files are
deposited in GigaDB together with this release.

The decoy-free calibration analyses added in v0.2.0 run from the same inputs:
`python scripts/run_calibration_nulls.py circular`, `... sameline`, `... report` (a few minutes;
results in [`results/calibration/calibration.md`](results/calibration/calibration.md)). Circular
shifts of each pair's own feature within chromosome arms (2,000 per pair) were rejected at
two-sided 5% in 4.9% of cases over the six pairs and 4.8% on the same-line outcome.

### From the raw data (about 100 GB of downloads)

All input data are public. Large files are not in git; they are downloaded into `data/`.
The shell scripts expect the environment at `env/venv`.

| Step | Command / document |
|---|---|
| 1. Download screens, ATAC-seq, annotations | `python data/fetch_all.py` (manifests and accessions in `recon/`) |
| 2. Gene annotations, copy number, p53 sets | `python scripts/build_annotations.py` ([`docs/annotations.md`](docs/annotations.md)) |
| 3. Local ATAC pipeline (Jurkat, RPE1; K562 validation) | `scripts/atac_pipeline.py`, `scripts/run_rpe1_sub60M.sh` ([`docs/atac_local_pipeline.md`](docs/atac_local_pipeline.md)) |
| 4. Accessibility features | `scripts/build_features.sh` ([`docs/features.md`](docs/features.md)) |
| 5. Outcome M | `scripts/run_aim4_aim1.sh` ([`docs/aim1_outcome.md`](docs/aim1_outcome.md)) |
| 6. Main test, controls, sensitivities | `scripts/run_aim4.sh`, `scripts/run_same_line_control.py`, `scripts/aim4x_rpe1_diagnostic.py` |
| 7. Figures | `python scripts/make_figures.py` (reads `results/`) |

Steps 1–5 can be skipped by starting from the intermediates. Only the covariates need the
screens, and `scripts/fetch_control_moments.py` reads just their control cells (~4 GB).
Result tables and figures from the published run are already in `results/` and `figures/final/`.

## Frozen analysis plan and record

- [`ANALYSIS_PLAN.md`](ANALYSIS_PLAN.md): the analysis plan, frozen before the confirmatory
  analysis (git tag `plan-v1.0`; amendments tagged `plan-v1.1` and `plan-v1.3`). The tags live in
  this repository; the plan was not registered with an external registry.
- `pipeline-v1.0`, `pipeline-v1.1`: the frozen pipeline tags used for the matched pairs.
- Every deviation from the plan is listed with its date and reason in the paper; commit
  timestamps record when each step was made. Changes after v0.1.0 are in
  [`CHANGELOG.md`](CHANGELOG.md).

## Repository map

The analysis was organised in four stages, and file names keep the stage labels of the plan:
**Aim 1** = outcome M from the Perturb-seq screens, **Aim 2** = the test on the development
pair (RPE1 − K562), **Aim 3** = calibration benchmarks (synthetic and planted signal),
**Aim 4** = the confirmatory run on the matched and mixed pairs. `aim4x_*` and
`results/aim4_exploratory/` are the post hoc (exploratory) RPE1 diagnostic.

| Path | Contents | In git? |
|---|---|---|
| `src/cdtest/` | The package: `crt.py` (the test), `simulate.py`, `k_rule_sim.py`, `planted.py` (benchmarks), `aim1.py`/`screens.py` (outcome), `features.py`/`atac.py` (accessibility), `covariates.py`/`annot.py`, `intermediates.py` | yes |
| `examples/` | Quick-start | yes |
| `tests/` | Unit tests, randomization-validity check, reproduction from intermediates | yes |
| `scripts/` | Pipeline, analysis, reproduction and figure scripts | yes |
| `intermediates/` | Deposited processed inputs (see above) | yes |
| `results/`, `figures/` | Outputs of the published run (small files); `results/calibration/` added in v0.2.0 | yes |
| `docs/` | Method notes and the API reference | yes |
| `recon/` | Data reconnaissance done before the plan was frozen: dataset inventories, accessions, candidate-dataset searches | yes |
| `proposal/` | The study proposal (11 October 2026) that preceded the analysis plan | yes |
| `status.sh` | Developer helper that printed download and job status during the run; not needed to use or reproduce cdtest | yes |
| `data/` | Downloaded inputs (~100 GB) | no |
| `env/` | Python environment | no |

## Citation

See [`CITATION.cff`](CITATION.cff) (GitHub shows a "Cite this repository" button).

## License

MIT; see [`LICENSE`](LICENSE).
