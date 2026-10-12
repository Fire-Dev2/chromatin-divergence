# RPE1 diagnostic — EXPLORATORY (post hoc)

**EXPLORATORY (post hoc).** H1 was tested once on the frozen pipeline (`pipeline-v1.1`, plan v1.3; `results/aim4/`). Everything below was decided after seeing that result, is not confirmatory, and changes no primary file or tag. Question: is the negative accessibility-response association that appears only in pairs containing RPE1 an RPE1 processing artefact or biology?

Method: for each RPE1 pair the frozen main test (block permutation, x4 k rule, chromosome-arm terms in f(C), B = 9,999, the aim2 seed of the primary feature for every variant, refined test-inversion CI) is rerun with one input changed. Feature variants are rebuilt from the saved raw signal with the frozen feature steps (log2 → quantile normalisation across the 11 lines → − log2 CN → distal sum → combined z); only the RPE1 column (or RPE1 CN track) changes. The rebuild reproduces the primary feature to ≤ 3e-8 and the primary θ̂, CI and p exactly. pR²/rel = partial R² / the pair's full-data reliability (no bootstrap). Gene-subset variants refit f(C) (splines, arm terms) and k on the subset. Code: `scripts/aim4x_rpe1_diagnostic.py` (stages features, fc4dn, tests, extra, joint), `scripts/aim4x_rpe1_report.py`; logs in `results/aim4_exploratory/logs/`.

## RPE1 − K562 (x = A − B)

| group | variant | genes | θ̂ | 95% CI | p two-sided | pR²/rel |
|---|---|---|---|---|---|---|
| primary | primary (frozen main test; reproduced) | 5,262 | -0.0224 | [-0.0315, -0.0135] | 0.0001 | 0.545% |
| copy number | CN: RPE1 = DepMap RPE1-clone consensus (features + cn_diff) | 5,262 | -0.0224 | [-0.0315, -0.0134] | 0.0001 | 0.546% |
| copy number | … features only | 5,262 | -0.0224 | [-0.0314, -0.0134] | 0.0001 | 0.544% |
| copy number | … cn_diff covariate only | 5,262 | -0.0224 | [-0.0316, -0.0135] | 0.0001 | 0.547% |
| copy number | CN: no CN correction for RPE1 (features + cn_diff) | 5,262 | -0.0224 | [-0.0315, -0.0135] | 0.0001 | 0.547% |
| accessibility source | RPE1 CPM from rep1 only | 5,262 | -0.0224 | [-0.0314, -0.0135] | 0.0001 | 0.552% |
| accessibility source | RPE1 CPM from rep2 only | 5,262 | -0.0219 | [-0.0311, -0.0130] | 0.0001 | 0.527% |
| accessibility source | RPE1 = 4DN-processed FC bigWig (others CPM) | 5,262 | 0.0044 | [-0.0084, 0.0050] | 0.6682 | 0.075% |
| accessibility source | RPE1 = local-pipeline FC bigWig (others CPM; control) | 5,262 | 0.0036 | [-0.0082, 0.0054] | 0.6752 | 0.053% |
| accessibility source | RPE1 CPM thinned to Jurkat depth (115.8M ins.), mean of 3 seeds (seed range -0.0226 … -0.0222) | 5,262 | -0.0224 | [-0.0316, -0.0136] | 0.0001 | 0.549% |
| accessibility source | RPE1 CPM thinned to 58M insertions, mean of 3 seeds (seed range -0.0224 … -0.0221) | 5,262 | -0.0222 | [-0.0312, -0.0132] | 0.0001 | 0.541% |
| accessibility source | all lines FC (frozen FC sensitivity) | 5,262 | -0.0026 | [-0.0153, 0.0076] | 0.617 | 0.005% |
| accessibility source | all lines FC, RPE1 = 4DN FC + 4DN peaks | 5,262 | 3.77e-04 | [-0.0168, 0.0091] | 0.9545 | 0.000% |
| accessibility source | partner = its FC, RPE1 kept CPM (control) | 5,262 | -0.0172 | [-0.0238, -0.0129] | 0.0001 | 0.860% |
| accessibility source | partner and RPE1 = FC (others CPM) | 5,262 | -0.0072 | [-0.0159, -1.15e-04] | 0.0291 | 0.079% |
| accessibility source | RPE1 CPM part shared with RPE1 FC (RPE1 side only) | 5,262 | 0.0074 | [-0.0117, 0.0251] | 0.3877 | 0.014% |
| accessibility source | RPE1 CPM part NOT shared with RPE1 FC (RPE1 side only) | 5,262 | -0.0154 | [-0.0211, -0.0017] | 0.0121 | 0.541% |
| promoter vs distal | promoter feature | 5,262 | -0.0130 | [-0.0169, -0.0090] | 0.0001 | 0.936% |
| promoter vs distal | distal feature | 5,262 | -0.0022 | [-0.0074, 0.0031] | 0.3431 | 0.017% |
| chromosomes | exclude chr10 | 5,053 | -0.0227 | [-0.0305, -0.0125] | 0.0001 | 0.563% |
| chromosomes | exclude chrX | 5,118 | -0.0224 | [-0.0316, -0.0135] | 0.0001 | 0.536% |
| chromosomes | exclude chr10 and chrX | 4,909 | -0.0228 | [-0.0312, -0.0125] | 0.0001 | 0.554% |
| p53 / expression | exclude Fischer p53 core targets | 5,237 | -0.0228 | [-0.0298, -0.0119] | 0.0001 | 0.572% |
| p53 / expression | exclude HALLMARK_P53_PATHWAY | 5,193 | -0.0214 | [-0.0305, -0.0128] | 0.0001 | 0.514% |
| p53 / expression | exclude both p53 sets | 5,181 | -0.0215 | [-0.0322, -0.0145] | 0.0001 | 0.518% |
| p53 / expression | bottom tercile of RPE1 expression | 1,754 | -0.0090 | [-0.0176, 0.0043] | 0.1187 | 0.227% |
| p53 / expression | middle tercile of RPE1 expression | 1,754 | -0.0162 | [-0.0394, -0.0069] | 0.0025 | 0.343% |
| p53 / expression | top tercile of RPE1 expression | 1,754 | -0.0340 | [-0.0643, -0.0170] | 0.0011 | 0.829% |
| which side | RPE1 feature alone (oriented, resid. on C) | 5,262 | -0.0128 | [-0.0183, -0.0057] | 0.0001 | 0.425% |
| which side | partner feature alone (oriented, resid. on C) | 5,262 | 0.0030 | [-0.0089, 0.0068] | 0.596 | 0.026% |

Leave-one-chromosome-out (θ̂ only): range -0.0245 (drop chr19) to -0.0203 (drop chr2); jackknife SE 0.0046; 0/23 drops give θ̂ > 0.

## Jurkat − RPE1 (x = A − B)

| group | variant | genes | θ̂ | 95% CI | p two-sided | pR²/rel |
|---|---|---|---|---|---|---|
| primary | primary (frozen main test; reproduced) | 5,222 | -0.0219 | [-0.0303, -0.0104] | 0.0002 | 0.511% |
| copy number | CN: RPE1 = DepMap RPE1-clone consensus (features + cn_diff) | 5,222 | -0.0219 | [-0.0302, -0.0102] | 0.0002 | 0.509% |
| copy number | … features only | 5,222 | -0.0219 | [-0.0302, -0.0101] | 0.0002 | 0.507% |
| copy number | … cn_diff covariate only | 5,222 | -0.0218 | [-0.0303, -0.0104] | 0.0002 | 0.507% |
| copy number | CN: no CN correction for RPE1 (features + cn_diff) | 5,222 | -0.0220 | [-0.0305, -0.0105] | 0.0002 | 0.518% |
| accessibility source | RPE1 CPM from rep1 only | 5,222 | -0.0221 | [-0.0304, -0.0104] | 0.0002 | 0.523% |
| accessibility source | RPE1 CPM from rep2 only | 5,222 | -0.0213 | [-0.0299, -0.0101] | 0.0002 | 0.489% |
| accessibility source | RPE1 = 4DN-processed FC bigWig (others CPM) | 5,222 | 0.0023 | [-0.0082, 0.0137] | 0.2377 | 0.023% |
| accessibility source | RPE1 = local-pipeline FC bigWig (others CPM; control) | 5,222 | 0.0020 | [-0.0110, 0.0117] | 0.3616 | 0.016% |
| accessibility source | RPE1 CPM thinned to Jurkat depth (115.8M ins.), mean of 3 seeds (seed range -0.0221 … -0.0217) | 5,222 | -0.0219 | [-0.0304, -0.0104] | 0.0002 | 0.514% |
| accessibility source | RPE1 CPM thinned to 58M insertions, mean of 3 seeds (seed range -0.0220 … -0.0215) | 5,222 | -0.0218 | [-0.0301, -0.0102] | 0.0002 | 0.509% |
| accessibility source | all lines FC (frozen FC sensitivity) | 5,222 | -0.0046 | [-0.0168, 0.0097] | 0.4643 | 0.011% |
| accessibility source | all lines FC, RPE1 = 4DN FC + 4DN peaks | 5,222 | -0.0014 | [-0.0112, 0.0107] | 0.7938 | 0.001% |
| accessibility source | partner = its FC, RPE1 kept CPM (control) | 5,222 | -0.0146 | [-0.0213, -0.0033] | 0.0038 | 0.559% |
| accessibility source | partner and RPE1 = FC (others CPM) | 5,222 | -0.0059 | [-0.0131, 0.0030] | 0.1476 | 0.045% |
| accessibility source | RPE1 CPM part shared with RPE1 FC (RPE1 side only) | 5,222 | 0.0082 | [-0.0068, 0.0236] | 0.24 | 0.024% |
| accessibility source | RPE1 CPM part NOT shared with RPE1 FC (RPE1 side only) | 5,222 | -0.0108 | [-0.0218, -0.0028] | 0.0075 | 0.289% |
| promoter vs distal | promoter feature | 5,222 | -0.0158 | [-0.0193, -0.0112] | 0.0001 | 1.285% |
| promoter vs distal | distal feature | 5,222 | 9.17e-04 | [-0.0040, 0.0066] | 0.6634 | 0.003% |
| chromosomes | exclude chr10 | 5,002 | -0.0233 | [-0.0326, -0.0124] | 0.0001 | 0.584% |
| chromosomes | exclude chrX | 5,071 | -0.0223 | [-0.0310, -0.0107] | 0.0001 | 0.525% |
| chromosomes | exclude chr10 and chrX | 4,851 | -0.0239 | [-0.0329, -0.0125] | 0.0002 | 0.602% |
| p53 / expression | exclude Fischer p53 core targets | 5,197 | -0.0220 | [-0.0321, -0.0131] | 0.0001 | 0.519% |
| p53 / expression | exclude HALLMARK_P53_PATHWAY | 5,164 | -0.0225 | [-0.0329, -0.0122] | 0.0001 | 0.545% |
| p53 / expression | exclude both p53 sets | 5,150 | -0.0225 | [-0.0322, -0.0129] | 0.0001 | 0.548% |
| p53 / expression | bottom tercile of RPE1 expression | 1,741 | -0.0087 | [-0.0223, -6.74e-05] | 0.0309 | 0.228% |
| p53 / expression | middle tercile of RPE1 expression | 1,740 | -0.0204 | [-0.0456, -0.0112] | 0.0009 | 0.505% |
| p53 / expression | top tercile of RPE1 expression | 1,741 | -0.0344 | [-0.0626, -0.0151] | 0.0005 | 0.752% |
| which side | RPE1 feature alone (oriented, resid. on C) | 5,222 | -0.0083 | [-0.0166, -0.0027] | 0.0027 | 0.196% |
| which side | partner feature alone (oriented, resid. on C) | 5,222 | 1.06e-04 | [-0.0133, 0.0076] | 0.9583 | 0.000% |

Leave-one-chromosome-out (θ̂ only): range -0.0239 (drop chr14) to -0.0197 (drop chr15); jackknife SE 0.0053; 0/23 drops give θ̂ > 0.

## RPE1 − HepG2 (x = A − B)

| group | variant | genes | θ̂ | 95% CI | p two-sided | pR²/rel |
|---|---|---|---|---|---|---|
| primary | primary (frozen main test; reproduced) | 5,631 | -0.0173 | [-0.0263, -0.0140] | 0.0001 | 0.770% |
| copy number | CN: RPE1 = DepMap RPE1-clone consensus (features + cn_diff) | 5,631 | -0.0174 | [-0.0255, -0.0135] | 0.0001 | 0.777% |
| copy number | … features only | 5,631 | -0.0174 | [-0.0255, -0.0135] | 0.0001 | 0.780% |
| copy number | … cn_diff covariate only | 5,631 | -0.0173 | [-0.0263, -0.0140] | 0.0001 | 0.768% |
| copy number | CN: no CN correction for RPE1 (features + cn_diff) | 5,631 | -0.0170 | [-0.0251, -0.0131] | 0.0001 | 0.750% |
| accessibility source | RPE1 CPM from rep1 only | 5,631 | -0.0174 | [-0.0262, -0.0141] | 0.0001 | 0.782% |
| accessibility source | RPE1 CPM from rep2 only | 5,631 | -0.0170 | [-0.0258, -0.0137] | 0.0001 | 0.745% |
| accessibility source | RPE1 = 4DN-processed FC bigWig (others CPM) | 5,631 | 0.0039 | [-0.0025, 0.0070] | 0.1608 | 0.127% |
| accessibility source | RPE1 = local-pipeline FC bigWig (others CPM; control) | 5,631 | 0.0031 | [-0.0031, 0.0060] | 0.2568 | 0.082% |
| accessibility source | RPE1 CPM thinned to Jurkat depth (115.8M ins.), mean of 3 seeds (seed range -0.0174 … -0.0173) | 5,631 | -0.0173 | [-0.0263, -0.0141] | 0.0001 | 0.770% |
| accessibility source | RPE1 CPM thinned to 58M insertions, mean of 3 seeds (seed range -0.0179 … -0.0171) | 5,631 | -0.0175 | [-0.0263, -0.0141] | 0.0001 | 0.788% |
| accessibility source | all lines FC (frozen FC sensitivity) | 5,631 | -0.0043 | [-0.0097, 0.0047] | 0.2672 | 0.029% |
| accessibility source | all lines FC, RPE1 = 4DN FC + 4DN peaks | 5,631 | -7.12e-04 | [-0.0071, 0.0095] | 0.8661 | 0.001% |
| accessibility source | partner = its FC, RPE1 kept CPM (control) | 5,631 | -0.0096 | [-0.0145, -0.0041] | 0.0009 | 0.761% |
| accessibility source | partner and RPE1 = FC (others CPM) | 5,631 | -0.0027 | [-0.0068, 0.0037] | 0.3179 | 0.022% |
| accessibility source | RPE1 CPM part shared with RPE1 FC (RPE1 side only) | 5,631 | 0.0081 | [-0.0034, 0.0214] | 0.1528 | 0.036% |
| accessibility source | RPE1 CPM part NOT shared with RPE1 FC (RPE1 side only) | 5,631 | -0.0113 | [-0.0154, -0.0043] | 0.0002 | 0.707% |
| promoter vs distal | promoter feature | 5,631 | -0.0059 | [-0.0077, -0.0033] | 0.0001 | 0.585% |
| promoter vs distal | distal feature | 5,631 | -0.0056 | [-0.0110, -0.0034] | 0.0001 | 0.230% |
| chromosomes | exclude chr10 | 5,378 | -0.0179 | [-0.0272, -0.0151] | 0.0001 | 0.834% |
| chromosomes | exclude chrX | 5,471 | -0.0183 | [-0.0270, -0.0144] | 0.0001 | 0.845% |
| chromosomes | exclude chr10 and chrX | 5,218 | -0.0189 | [-0.0279, -0.0151] | 0.0001 | 0.920% |
| p53 / expression | exclude Fischer p53 core targets | 5,591 | -0.0176 | [-0.0252, -0.0129] | 0.0001 | 0.812% |
| p53 / expression | exclude HALLMARK_P53_PATHWAY | 5,542 | -0.0178 | [-0.0242, -0.0122] | 0.0001 | 0.835% |
| p53 / expression | exclude both p53 sets | 5,520 | -0.0179 | [-0.0248, -0.0131] | 0.0001 | 0.851% |
| p53 / expression | bottom tercile of RPE1 expression | 1,877 | -0.0040 | [-0.0105, 0.0036] | 0.1884 | 0.086% |
| p53 / expression | middle tercile of RPE1 expression | 1,877 | -0.0155 | [-0.0317, -0.0092] | 0.0002 | 0.632% |
| p53 / expression | top tercile of RPE1 expression | 1,877 | -0.0284 | [-0.0396, -0.0118] | 0.0001 | 1.458% |
| which side | RPE1 feature alone (oriented, resid. on C) | 5,631 | -0.0094 | [-0.0128, -0.0036] | 0.0003 | 0.546% |
| which side | partner feature alone (oriented, resid. on C) | 5,631 | 0.0021 | [-0.0045, 0.0045] | 0.4998 | 0.028% |

Leave-one-chromosome-out (θ̂ only): range -0.0187 (drop chr19) to -0.0154 (drop chr1); jackknife SE 0.0032; 0/23 drops give θ̂ > 0.

## (5) Leave-one-chromosome-out jackknife, RPE1 − HepG2

| dropped | chr1 | chr2 | chr3 | chr4 | chr5 | chr6 | chr7 | chr8 | chr9 | chrX | chr10 | chr11 | chr12 | chr13 | chr14 | chr15 | chr16 | chr17 | chr18 | chr19 | chr20 | chr21 | chr22 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| θ̂ | -0.0154 | -0.0181 | -0.0177 | -0.0169 | -0.0167 | -0.0174 | -0.0174 | -0.0174 | -0.0169 | -0.0183 | -0.0179 | -0.0166 | -0.0168 | -0.0175 | -0.0173 | -0.0170 | -0.0179 | -0.0171 | -0.0177 | -0.0187 | -0.0179 | -0.0174 | -0.0167 |

## (6) p53 gene sets

| pair | Fischer core in pair | Hallmark in pair | union | share of θ̂ numerator from union | share of genes |
|---|---|---|---|---|---|
| RPE1_vs_K562 | 25 | 69 | 81 | 4.8% | 1.5% |
| Jurkat_vs_RPE1 | 25 | 58 | 72 | -1.6% | 1.4% |
| RPE1_vs_HepG2 | 40 | 89 | 111 | -0.4% | 2.0% |

Fischer: 116 core targets listed (`fischer2017_p53_core_mapped.tsv`); Hallmark: 200 symbols (`HALLMARK_P53_PATHWAY.v2025.1`), matched by gene name.

## (7) Which side — joint two-line model and pairs without RPE1

Joint OLS of M (residualised on f(C)) on both lines' residualised features, each oriented as in x (A feature, −B feature); descriptive, no permutation. The two residualised features correlate −0.79 to −0.85 (shared accessibility), so the marginal one-line tests above are the cleaner read.

| pair | coef A | coef B | r (resid., oriented) |
|---|---|---|---|
| RPE1_vs_K562 | -0.0293 | -0.0191 | -0.81 |
| Jurkat_vs_RPE1 | -0.0205 | -0.0271 | -0.84 |
| RPE1_vs_HepG2 | -0.0210 | -0.0143 | -0.79 |
| K562_vs_Jurkat | -0.0001 | 0.0004 | -0.82 |
| K562_vs_HepG2 | 0.0066 | 0.0081 | -0.85 |
| Jurkat_vs_HepG2 | 0.0014 | 0.0047 | -0.81 |

One-line-alone main tests in the pairs without RPE1 (context):

| pair | variant | θ̂ | p two-sided |
|---|---|---|---|
| K562_vs_Jurkat | primary | 1.98e-04 | 0.8381 |
| K562_vs_Jurkat | side_K562_alone | -4.80e-04 | 0.496 |
| K562_vs_Jurkat | side_Jurkat_alone | 5.20e-04 | 0.7598 |
| K562_vs_HepG2 | primary | 0.0071 | 0.0551 |
| K562_vs_HepG2 | side_K562_alone | 6.12e-05 | 0.9553 |
| K562_vs_HepG2 | side_HepG2_alone | 0.0022 | 0.2577 |
| Jurkat_vs_HepG2 | primary | 0.0024 | 0.2443 |
| Jurkat_vs_HepG2 | side_Jurkat_alone | -0.0021 | 0.2868 |
| Jurkat_vs_HepG2 | side_HepG2_alone | 0.0035 | 0.2178 |

## Data-quality context

Insertion TSS enrichment (TSS ± 50 bp vs flanks 1.9–2.0 kb, protein-coding, replicates pooled): RPE1 10.6, K562 13.7, HepG2 14.3, A549 15.5, GM12878 15.7, Panc1 15.8, PC-3 17.9, IMR-90 23.4, HCT116 25.1, Jurkat 25.2, MCF-7 28.3. CPM vs FC combined feature (the two frozen builds), per line over the RPE1 − HepG2 genes: HepG2 0.39, PC-3 0.41, A549 0.44, Panc1 0.46, IMR-90 0.51, HCT116 0.53, RPE1 0.55, MCF-7 0.57, K562 0.57, GM12878 0.65, Jurkat 0.70 — RPE1 is not an outlier; FC is a substantially different summary for every line. Within the hybrid (RPE1 local FC quantile-normalised against the CPM lines) the RPE1 CPM vs FC correlation over each pair's genes is RPE1_vs_K562 0.29, Jurkat_vs_RPE1 0.33, RPE1_vs_HepG2 0.28.

## Link to the planned same-line control (`results/aim4/same_line_control.md`)

M = K562 essential − K562 genome-wide (two screens of the same line: technical divergence only). Frozen main test, seed keyed on the row label (so p differs slightly from `same_line_control.json`; θ̂ reproduces). partial r = θ̂·SD(resid. x)/SD(resid. M). For comparison, the cross-line partial r of the pair x in the RPE1 pairs is −0.074 / −0.071 / −0.087, and of the RPE1 CPM-specific part −0.073 / −0.054 / −0.083 (oriented RPE1 side).

| feature tested against same-line M | θ̂ | 95% CI | p two-sided | partial r |
|---|---|---|---|---|
| x_RPE1_vs_K562 (CPM; reproduces same_line_control) | -0.0013 | [-0.0021, -4.96e-04] | 0.0012 | -0.046 |
| x_RPE1_vs_HepG2 (CPM; reproduces same_line_control) | -9.32e-04 | [-0.0017, -7.60e-05] | 0.02 | -0.031 |
| x_RPE1_vs_K562, RPE1 = local FC hybrid | -0.0018 | [-0.0027, -0.0015] | 0.0001 | -0.103 |
| x_RPE1_vs_K562, RPE1 = 4DN FC hybrid | -0.0019 | [-0.0028, -0.0016] | 0.0001 | -0.102 |
| x_RPE1_vs_HepG2, RPE1 = 4DN FC hybrid | -0.0018 | [-0.0025, -0.0013] | 0.0001 | -0.096 |
| RPE1 CPM alone | 9.51e-04 | [8.63e-05, 0.0015] | 0.0143 | 0.046 |
| RPE1 local FC (hybrid) alone | -0.0011 | [-0.0017, -3.64e-04] | 0.0025 | -0.046 |
| RPE1 CPM part shared with FC | -0.0018 | [-0.0030, -6.06e-04] | 0.0026 | -0.046 |
| RPE1 CPM part NOT shared with FC | 0.0019 | [8.20e-04, 0.0027] | 0.0005 | 0.080 |
| K562 CPM alone | 0.0016 | [0.0012, 0.0025] | 0.0001 | 0.079 |
| HepG2 CPM alone | 0.0013 | [5.68e-04, 0.0019] | 0.0002 | 0.067 |
| Jurkat CPM alone | 0.0011 | [3.99e-04, 0.0016] | 0.0009 | 0.056 |

## Reading

- **Copy number (1): not the cause.** DepMap RPE1-clone consensus CN or no RPE1 CN correction changes θ̂ by < 0.0004 in every pair (features, covariate, or both). Arm-level CN is already absorbed by the arm terms in f(C).
- **Replicates and depth (2a, 2c): not the cause.** Each RPE1 replicate alone, and the pooled insertions thinned to Jurkat's depth (115.8M insertions = ~58M fragments) or to 58M insertions (3 seeds each), give θ̂ within 0.0007 of the primary, all p₂ ≤ 0.0002. (RPE1 at 152M insertions is not deep relative to the other lines anyway.)
- **Signal type (2b): the effect disappears when RPE1 accessibility is measured as fold-change.** Swapping only the RPE1 column for the 4DN-processed FC bigWig gives θ̂ ≈ +0.002 to +0.004 (p₂ 0.16–0.67) in all three pairs; the local-pipeline FC bigWig from the same reads does the same — so it is not 4DN-vs-local processing but CPM-vs-FC for RPE1. All-FC features (RPE1 local or 4DN FC/peaks) are also n.s. Swapping only the partner to FC keeps the effect (θ̂ −0.010 to −0.017, p₂ ≤ 0.004). With both RPE1 and the partner as FC, a smaller residual remains in RPE1 − K562 (θ̂ −0.007, p₂ 0.03) but not in the other two pairs (−0.006, −0.003, n.s.). Splitting the RPE1 CPM feature into the part shared with RPE1 FC and the remainder: the shared part gives θ̂ ≈ +0.007 to +0.008 (n.s.); the CPM-specific remainder carries all of the negative association (θ̂ −0.011 to −0.015, p₂ 0.0002–0.012).
- **Location (3):** promoter-driven in all three pairs; distal is n.s. in RPE1 − K562 and Jurkat − RPE1 and negative (smaller pR²) in RPE1 − HepG2.
- **Chromosomes (4, 5): not driven by chr10 (10q gain) or chrX, nor by any single chromosome.** Exclusions slightly strengthen θ̂; every leave-one-chromosome-out θ̂ stays negative (RPE1 − HepG2: −0.0154 to −0.0187).
- **p53 (6): not driven by p53 genes.** Excluding Fischer core targets, Hallmark p53 pathway genes or both changes θ̂ by < 0.001 (all p₂ ≤ 0.0001); in two of three pairs the p53 genes contribute less than their gene share to the θ̂ numerator.
- **Expression (6):** the association grows with RPE1 expression (bottom tercile −0.004 to −0.009, n.s. or marginal; top tercile −0.028 to −0.034, p₂ ≤ 0.0011).
- **Side (7): it is RPE1's own accessibility.** RPE1 feature alone: θ̂ −0.008 to −0.013, p₂ ≤ 0.003 in every pair; partner alone: ≈ 0 (p₂ ≥ 0.5). The same one-line tests in the three pairs without RPE1 are all ≈ 0, and the joint two-line coefficients are ≈ 0 or positive there but negative for both terms in the RPE1 pairs (the partner term picks up RPE1's deviation from the shared accessibility profile).
- **Data quality:** RPE1 has the lowest insertion TSS enrichment of the 11 lines (10.6 vs 13.7–28.3).
- **Same-line control (planned, added after this diagnostic was drafted):** accessibility does associate with purely technical between-screen divergence. Every line's own CPM feature correlates positively with K562 essential − genome-wide M (partial r +0.046 to +0.079; K562 strongest). The RPE1_vs_K562 x reaches a partial r of −0.046, ~60% of its cross-line size (−0.074), mostly through its −K562 term. But this technical channel does not have the structure of the RPE1 effect. Replacing RPE1 CPM with FC removes the cross-line RPE1 association, yet makes the same-line association stronger (partial r about −0.10, p₂ 1e-4). The RPE1 CPM-specific part, which carries the whole cross-line effect, also relates to same-line M (+0.080), and the part shared with FC has the opposite sign (−0.046); the pattern does not match the cross-line one. So the same-line result is a general warning that non-biological between-screen divergence tracks accessibility at a magnitude comparable to the cross-line effects. It is not a mechanistic explanation of the RPE1-specific pattern.

## Conclusion (plain language)

EXPLORATORY (post hoc). The evidence moderately favours a non-biological explanation of the negative association in the RPE1 pairs over RPE1 biology. The planned same-line control adds weight against biology but makes the exact artefact less certain. Within the RPE1 pairs, the association is not explained by copy-number handling, either replicate, sequencing depth, the 10q gain, chrX, any single chromosome, or p53 target / p53-pathway genes. It belongs to RPE1's own accessibility, is strongest in highly expressed genes, and exists only in the insertion-count (CPM) summary of the RPE1 reads. When the same reads are summarised as fold-change over local background (by the local pipeline or by 4DN's), it vanishes in all three pairs. It sits entirely in the part of the RPE1 CPM signal that FC does not share, and RPE1 has the lowest TSS enrichment of the 11 libraries. A real accessibility-response relation should survive both summaries of the same reads, so this points to a CPM-specific background/normalisation component of the RPE1 track. The same-line control (two K562 screens, no biological difference) shows that accessibility differences involving RPE1 also associate with purely technical between-screen divergence. In standardised terms that is about 60% of the cross-line size, so cross-line θ̂ values of this magnitude cannot be read as biology by default. That technical association is not the same signal, though: it is mostly carried by the partner/K562 accessibility and by the FC-shared part of RPE1, and it gets stronger, not weaker, when RPE1 is switched to FC. The two non-biological channels therefore look distinct, and these data cannot assign the RPE1-pair effect wholly to one of them. Biology is not excluded: FC and CPM agree only modestly for every line (r ≈ 0.4–0.7; RPE1 0.55), and no independent RPE1 ATAC dataset was available. The RPE1-pair negatives and the matched-pair disagreement in H1 should be reported as most likely technical. They are not evidence of a biological negative accessibility effect.
