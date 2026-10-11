# Accessibility feature sanity checks

Built 2026-10-11T17:38:57 from `config/feature_lines.tsv`. Lines (10): K562 (target), RPE1 (target_provisional), HepG2 (target), GM12878 (decoy), HCT116 (decoy), MCF-7 (decoy), A549 (decoy), IMR-90 (decoy), PC-3 (decoy), Panc1 (decoy).
Skipped (inputs not yet available): Jurkat.

## Union regions

- Summit windows pooled: 1,932,956; merged: 340,993; blacklist-overlapping dropped: 2,510.
- **Union regions: 338,483.**
- Region length: median 725 bp, 90th pct 1482, max 6,195 bp.
- Genes: 54,710; gene-region distal pairs 1,582,423; genes with no union region 2.5-100 kb away: 848.
- Pooled z: promoter mean 1.201 SD 0.844; distal mean 3.867 SD 1.157.

Per-line summit counts: K562 181,336, RPE1 277,348, HepG2 208,885, GM12878 207,873, HCT116 169,178, MCF-7 173,570, A549 208,895, IMR-90 174,980, PC-3 193,312, Panc1 137,579.

## Promoter feature correlation across lines (Pearson, all genes)

| | K562 | RPE1 | HepG2 | GM12878 | HCT116 | MCF-7 | A549 | IMR-90 | PC-3 | Panc1 |
|---|---|---|---|---|---|---|---|---|---|---|
| K562 | 1.00 | 0.68 | 0.69 | 0.73 | 0.74 | 0.66 | 0.71 | 0.72 | 0.66 | 0.72 |
| RPE1 | 0.68 | 1.00 | 0.72 | 0.80 | 0.79 | 0.66 | 0.77 | 0.87 | 0.71 | 0.73 |
| HepG2 | 0.69 | 0.72 | 1.00 | 0.75 | 0.80 | 0.70 | 0.79 | 0.77 | 0.66 | 0.73 |
| GM12878 | 0.73 | 0.80 | 0.75 | 1.00 | 0.79 | 0.67 | 0.77 | 0.84 | 0.68 | 0.73 |
| HCT116 | 0.74 | 0.79 | 0.80 | 0.79 | 1.00 | 0.74 | 0.82 | 0.81 | 0.74 | 0.76 |
| MCF-7 | 0.66 | 0.66 | 0.70 | 0.67 | 0.74 | 1.00 | 0.77 | 0.70 | 0.71 | 0.71 |
| A549 | 0.71 | 0.77 | 0.79 | 0.77 | 0.82 | 0.77 | 1.00 | 0.82 | 0.74 | 0.79 |
| IMR-90 | 0.72 | 0.87 | 0.77 | 0.84 | 0.81 | 0.70 | 0.82 | 1.00 | 0.73 | 0.77 |
| PC-3 | 0.66 | 0.71 | 0.66 | 0.68 | 0.74 | 0.71 | 0.74 | 0.73 | 1.00 | 0.73 |
| Panc1 | 0.72 | 0.73 | 0.73 | 0.73 | 0.76 | 0.71 | 0.79 | 0.77 | 0.73 | 1.00 |

Off-diagonal range 0.66-0.87 (median 0.73).
Distal feature off-diagonal correlation range 0.85-0.95. Promoter vs distal within line (median over lines): 0.37.

## Promoter feature vs baseline expression (expect Spearman > 0.3)

Source: `results/aim1/RPE1_vs_K562/gene_outcome.tsv`, columns gene_id, mean_ctrlA (RPE1), mean_ctrlB (K562) only.

| line | n genes | rho promoter | rho distal | rho combined | pass (promoter > 0.3) |
|---|---|---|---|---|---|
| RPE1 | 5,262 | -0.143 | 0.097 | 0.005 | FAIL |
| K562 | 5,262 | -0.091 | 0.005 | -0.040 | FAIL |

Supplementary (uses only gene_id membership): the 5,262 genes in gene_outcome passed the Aim 1 expression filter in both lines. AUC of the promoter feature, expressed set vs all other genes:
- RPE1: AUC 0.881 (median 2.27 vs 0.86).
- K562: AUC 0.851 (median 2.19 vs 0.93).

Between-line: Spearman(promoter K562-RPE1, log2 expr K562/RPE1) = -0.242 (descriptive only).

## Housekeeping promoters (percentile of promoter feature within each line; expect high)

| gene | K562 | RPE1 | HepG2 | GM12878 | HCT116 | MCF-7 | A549 | IMR-90 | PC-3 | Panc1 | pass (all >= 90th) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| GAPDH | 93 | 68 | 76 | 76 | 66 | 85 | 67 | 66 | 83 | 85 | FAIL |
| ACTB | 72 | 70 | 65 | 63 | 70 | 63 | 59 | 64 | 54 | 53 | FAIL |

## Lineage markers (promoter feature; value [within-line percentile])

| gene | K562 | RPE1 | HepG2 | GM12878 | HCT116 | MCF-7 | A549 | IMR-90 | PC-3 | Panc1 | check | pass |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| HBG1 | 1.20 [55] | 0.19 [5] | 0.52 [25] | 0.62 [31] | 1.10 [52] | 1.23 [55] | 0.38 [19] | 0.80 [41] | 0.18 [12] | -0.42 [1] | K562 > RPE1 (all-line rank 2/10) | PASS |
| GATA1 | 1.80 [71] | 0.85 [44] | 1.85 [73] | 0.87 [45] | 2.41 [94] | 1.49 [63] | 1.23 [57] | 0.60 [30] | 1.31 [59] | 1.35 [58] | K562 > RPE1 (all-line rank 3/10) | PASS |
| ALB | 0.23 [13] | 0.33 [12] | 1.25 [57] | 0.46 [20] | 0.40 [15] | 0.54 [26] | 0.62 [32] | 0.31 [10] | 0.16 [12] | 0.72 [35] | HepG2 > K562, RPE1 (all-line rank 1/10) | PASS |
| APOA1 | 0.52 [25] | 0.70 [37] | 1.48 [63] | 0.49 [23] | 0.44 [18] | 1.69 [68] | 1.12 [54] | 0.66 [34] | 0.24 [14] | 0.25 [14] | HepG2 > K562, RPE1 (all-line rank 2/10) | PASS |

## Summary

**4/8 checks pass.**

Diagnosis of the expression / housekeeping failures (2026-10-11 build, 10 lines): the same Spearman computed from the raw fold-change bigWig over TSS +/-250, +/-500 and +/-1000 bp, and the max over +/-1 kb, is -0.05 to -0.07 (K562) and -0.13 to -0.14 (RPE1), so the failure is not caused by the window, quantile normalisation or CN step. Promoter accessibility separates expressed from non-expressed genes strongly (AUC above), but within the expression-filtered set the MACS2 fold-change track (local-lambda background) does not rank strong promoters higher: very highly expressed genes (RPs, EEF1A1, GAPDH, ACTB) sit at FC ~1.5-3.5 over +/-1 kb. Any change of the signal track would be a plan deviation and is left to the PI.

Caveats: RPE1 uses the provisional 4DN fold-change bigWig + conservative peaks until the local pipeline output exists; GM12878 and IMR-90 have no DepMap CN (CN correction = 0). Jurkat is absent until its local pipeline finishes; rerun `scripts/build_features.sh`.
