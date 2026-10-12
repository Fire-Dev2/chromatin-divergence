# H5 mechanism case studies — K562 − HepG2

Run 2026-10-12 (`scripts/run_h5_cases.py`). Factor by factor, **no pooled test** (plan §7 H5). K562 − HepG2 is a mixed pair, so these cases illustrate a mechanism and do not establish it.

- TFs: sequence-specific, kept knockdowns of K562_vs_HepG2, ENCODE IDR thresholded peaks (GRCh38, no audit ERROR; ENCODE4 > ENCODE3 > ENCODE2, then fewer NOT_COMPLIANT, more reps, newest).
- Site: IDR peak overlapping a union region accessible in that line (signal > line median); counted if its summit is within 50 kb of the TSS. One-line ChIP: that peak set used for both lines, gated by each line's accessibility.
- Feature: count_K562 − count_HepG2. Outcome: d_g = (β̃A² − s²A) − (β̃B² − s²B) for that knockdown; Spearman ρ over genes, own gene excluded.
- Swap: 20 random other TFs from the same list (seed 20261012); rank 1 = own ρ highest of 21; p_emp = (1 + #swap ρ ≥ own)/21. Expected pattern = own ρ > 0 and rank 1.

## Main cases (Lambert 2018 sequence-specific TFs)

| TF | ChIP source | genes with sites | genes Δ≠0 | n genes | ρ own | p(ρ) | ρ swap median [min, max] | rank /21 | p_emp |
|---|---|---|---|---|---|---|---|---|---|
| MYBL2 | both | 5269 | 5217 | 5784 | +0.036 | 0.0056 | +0.020 [-0.016, +0.058] | 4 | 0.190 |
| ATF4 | both | 3866 | 3353 | 5784 | +0.031 | 0.019 | +0.011 [-0.023, +0.049] | 4 | 0.190 |
| PREB ✓ | both | 3134 | 2776 | 5784 | +0.028 | 0.033 | +0.010 [-0.026, +0.026] | 1 | 0.048 |
| THAP1 | K562 only | 3126 | 20 | 5784 | +0.026 | 0.044 | +0.013 [-0.026, +0.030] | 2 | 0.095 |
| CTCF | both | 5667 | 4440 | 5784 | +0.022 | 0.099 | +0.004 [-0.030, +0.052] | 4 | 0.190 |
| SNAPC5 | both | 3155 | 3107 | 5784 | +0.021 | 0.11 | +0.001 [-0.041, +0.043] | 6 | 0.286 |
| USF2 | both | 4533 | 3712 | 5784 | +0.020 | 0.14 | +0.008 [-0.023, +0.046] | 9 | 0.429 |
| TBP | HepG2 only | 1868 | 23 | 5784 | +0.018 | 0.17 | +0.010 [-0.018, +0.025] | 7 | 0.333 |
| MAX | both | 5751 | 5165 | 5784 | +0.009 | 0.49 | +0.014 [-0.013, +0.028] | 14 | 0.667 |
| ARID3A | both | 5624 | 5357 | 5784 | +0.007 | 0.57 | +0.009 [-0.024, +0.034] | 13 | 0.619 |
| CEBPZ | both | 1169 | 862 | 5784 | +0.002 | 0.88 | -0.007 [-0.024, +0.025] | 7 | 0.333 |
| NFYB | HepG2 only | 4939 | 322 | 5784 | +0.001 | 0.97 | +0.002 [-0.022, +0.040] | 11 | 0.524 |
| YBX1 | both | 262 | 250 | 5784 | -0.000 | 0.98 | +0.018 [-0.034, +0.048] | 16 | 0.762 |
| TFDP1 | both | 5657 | 4495 | 5784 | -0.001 | 0.97 | +0.008 [-0.014, +0.023] | 16 | 0.762 |
| DRAP1 | HepG2 only | 5499 | 696 | 5784 | -0.001 | 0.94 | +0.006 [-0.028, +0.058] | 14 | 0.667 |
| ATF5 | HepG2 only | 2761 | 66 | 5784 | -0.001 | 0.91 | +0.012 [-0.033, +0.028] | 17 | 0.810 |
| GABPA | both | 5630 | 4702 | 5784 | -0.003 | 0.84 | +0.004 [-0.015, +0.031] | 16 | 0.762 |
| CDC5L | K562 only | 2023 | 118 | 5784 | -0.006 | 0.62 | +0.006 [-0.015, +0.038] | 17 | 0.810 |
| THAP11 | HepG2 only | 5659 | 1090 | 5784 | -0.026 | 0.049 | +0.012 [-0.032, +0.023] | 20 | 0.952 |
| E2F6 | K562 only | 5572 | 552 | 5784 | -0.027 | 0.037 | -0.002 [-0.031, +0.037] | 20 | 0.952 |

## Supplementary: ENCODE-flagged factors that are not Lambert TFs (cofactors/RBPs; swap-pool members)

| TF | ChIP source | genes with sites | genes Δ≠0 | n genes | ρ own | p(ρ) | ρ swap median [min, max] | rank /21 | p_emp |
|---|---|---|---|---|---|---|---|---|---|
| XRCC5 | both | 5419 | 4646 | 5784 | +0.010 | 0.46 | +0.004 [-0.012, +0.022] | 9 | 0.429 |
| GPN1 | HepG2 only | 1774 | 10 | 5784 | +0.007 | 0.57 | +0.013 [-0.007, +0.031] | 15 | 0.714 |
| PHB2 | both | 2141 | 2135 | 5784 | +0.007 | 0.59 | +0.009 [-0.020, +0.031] | 12 | 0.571 |
| ZNF720 | HepG2 only | 1294 | 24 | 5784 | +0.006 | 0.66 | +0.009 [-0.021, +0.032] | 13 | 0.619 |
| HNRNPK | both | 4493 | 3598 | 5784 | +0.006 | 0.67 | +0.013 [-0.024, +0.048] | 15 | 0.714 |
| CBX1 | both | 5088 | 4172 | 5784 | -0.004 | 0.77 | -0.001 [-0.024, +0.040] | 13 | 0.619 |
| SFPQ | both | 4821 | 4807 | 5784 | -0.025 | 0.056 | +0.008 [-0.038, +0.040] | 20 | 0.952 |
| BCLAF1 | both | 944 | 944 | 5784 | -0.029 | 0.026 | +0.011 [-0.017, +0.054] | 21 | 1.000 |

Not evaluable (no accessible sites near any tested gene): TFAM.

**Summary.** 20 TFs (12 with ChIP in both lines). Own ρ > 0 for 12/20; expected pattern (own ρ > 0 and own sites best of 21) for 1: PREB. No pooled claim is made.

Figure: `figures/aim4_h5_tf_cases.png`. ChIP accessions and md5s: `h5_tf_cases.json`.

**Reading.** Only PREB shows the expected pattern (own ρ = +0.028, best of 21, p_emp 0.048 — one of 20 factor-level looks, so not beyond chance). MYBL2, ATF4, CTCF (rank 4) and THAP1 (rank 2) have positive own ρ but are not uniquely best. All |ρ| ≤ 0.04: site differences explain < 0.2% of rank variance. Swap ρ medians are mostly slightly positive, i.e. any TF's K562-minus-HepG2 accessible-site count tracks a generic accessibility difference that weakly co-varies with d_g; own-vs-swap rank is the informative comparison. THAP11 and E2F6 (one-line ChIP) and the cofactors SFPQ, BCLAF1 sit at the bottom (negative ρ).

**Caveats.** (1) Mixed pair: line identity, not chromatin alone, differs. (2) Two-line ChIP mixes binding differences with ChIP depth/quality differences (e.g. MYBL2 758 vs 24,426 peaks; BCLAF1 1,987 vs 50). (3) One-line ChIP features vary only through accessibility gating, so few genes have Δ ≠ 0. (4) The recon table lists only 'IDR thresholded peaks' (all GRCh38); conservative vs optimal is not distinguished there. (5) The ENCODE sequence-specific flag includes cofactors/RBPs; main cases are restricted to Lambert 2018 TFs (`data/annot/lambert2018_humantfs_v1.01.csv`).
