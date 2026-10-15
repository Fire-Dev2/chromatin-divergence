# Changelog

All notable changes to cdtest. Versions follow [semantic versioning](https://semver.org/).

## [0.2.0] - 2026-10-15

### Added
- **Reproduce from intermediates.** `intermediates/` holds the processed inputs of every
  confirmatory analysis (gene annotation, copy number, gene-level accessibility features for the
  4 screened and 7 decoy lines, control-cell expression moments; ~10 MB, gene subset), with
  `MANIFEST.tsv` md5 sums. `scripts/reproduce_from_intermediates.py` recomputes Table 2 and
  Figures 3–4 (except Figure 4C and the 4DN rows of 4B, drawn from committed results) from them in
  under a minute and compares every number with the committed results
  (`results_repro/COMPARISON.md`). `cdtest.intermediates.install()` puts the files in place.
- `scripts/fetch_control_moments.py`: control-cell expression moments read from the public h5ad
  files over HTTP range requests (~4 GB instead of ~95 GB).
- `scripts/run_calibration_nulls.py`: decoy-free real-data null (circular shifts of each pair's
  feature within chromosome arms) and planted-signal calibration on the same-line outcome
  (`results/calibration/`).
- `tests/test_intermediates.py`: the H1 numbers of both matched pairs (θ̂, block size, partial R²,
  B = 9,999 p-values) are recomputed from the intermediates on a fresh clone; the development-pair
  plan v1.2 test no longer skips.
- `CITATION.cff`, this changelog, `environment-atac.yml` and `Dockerfile` for the ATAC-seq
  tools, GitHub Actions CI, `docs/api.md`.

### Changed
- Quick-start: 500 replicates per scenario (was 100); see the README for the settings shown.
- `covariates.control_moments` reads the deposited moments when the pseudobulk matrix is absent
  (identical numbers; frozen analysis code otherwise unchanged).
- `make_figures.py`: optional overlay directory of recomputed results (`CDTEST_RESULTS_OVERLAY`)
  and output directory (`CDTEST_FIG_OUT`).
- README: "frozen analysis plan" wording; repository map explaining the project-stage names.

## [0.1.0] - 2026-10-14

First public release, accompanying the GigaByte Technical Release submission: the `cdtest`
package (feature-side block permutation, x4 block-size rule, test-inversion intervals, comparators,
synthetic and planted-signal benchmarks), the full analysis of four CRISPRi Perturb-seq screens with
bulk ATAC-seq features, frozen analysis plan tags (`plan-v1.0`, `plan-v1.1`, `plan-v1.3`) and
pipeline tags (`pipeline-v1.0`, `pipeline-v1.1`), results and figures.
