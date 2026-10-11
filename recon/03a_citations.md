# 03a — Verified citations (recon)

Verified 2026-10-11 against the Crossref API (`api.crossref.org/works/<doi>`), the arXiv API, the bioRxiv API, and full text where it could be fetched. BibTeX keys refer to `references.bib`.

Project question: how much of the difference in CRISPRi Perturb-seq knockdown responses between K562, RPE1, Jurkat and HepG2 can bulk ATAC-seq chromatin accessibility explain?

---

## 1. Nadig et al. 2025 — TRADE  `[nadig2025trade]`

**Citation.** Nadig A, Replogle JM, Pogson AN, Murthy M, McCarroll SA, Weissman JS, Robinson EB, O'Connor LJ. Transcriptome-wide analysis of differential expression in perturbation atlases. *Nature Genetics* 57(5):1228–1237 (2025; online 21 Apr 2025). DOI: **10.1038/s41588-025-02169-3**.
Preprint: "Transcriptome-wide characterization of genetic perturbations", bioRxiv 2024, DOI 10.1101/2024.07.03.601903.

**DOI corrections.**
- 10.1038/s41588-025-02177-3 is a different item: Peidli S., "Quantifying reproducible multivariate differential expression in genetic perturbation screens", *Nat Genet* 57:1066–1067 (2025). It is a single-author News & Views on TRADE. Cite it only if we want the commentary.
- 10.1101/2024.09.16.613342 is **not** the TRADE preprint. It is the Ahlmann-Eltze/Huber/Anders preprint (work 9).

**Summary.** TRADE (TRanscriptome-wide Analysis of Differential Expression) uses adaptive shrinkage to estimate the transcriptome-wide distribution of true DE effects from noisy gene-level Perturb-seq estimates. It defines a "transcriptome-wide impact" metric that stays stable across sampling depth. That stability allows comparisons across contexts that differ in power. The authors reanalysed K562-GenomeWide (9,866 expressed genes), K562-Essential (**2,057** knockdowns) and RPE1-Essential (**2,393**). They also generated new Jurkat-Essential and HepG2-Essential screens (**2,393** knockdowns each). Cross-line correlations of transcriptome-wide effects were estimated for **2,053 shared essential genes**, restricted to the **1,660** perturbations with nonzero impact in all four lines. Median correlations were **K562/Jurkat 0.74**, **HepG2/RPE1 0.75**, K562/HepG2 0.64, Jurkat/HepG2 0.69, RPE1/Jurkat 0.56 and **K562/RPE1 0.40**. The authors write that "K562 and Jurkat are hematopoietic cell lines that are p53 mutant and grow in suspension, while HepG2 and RPE1 are epithelial cell lines that are p53 wild-type and are adherent." They identify 241 perturbations with cell-type-specific excess effects (K562 47, RPE1 118, Jurkat 10, HepG2 66), and 56% of perturbations are highly correlated across all four lines.

*Provenance note:* these numbers come from the bioRxiv v1 full text. The published Nature Genetics version is paywalled and could not be checked, so figures may have shifted slightly in revision.

**What it leaves open.** The paper explains the correlation structure only through coarse line labels (lineage, p53 status, suspension vs. adherent). It does not test any gene-level or chromatin feature as an explanation, and it does not model which target genes or which downstream genes drive the divergence. The p53, tissue and growth-mode labels are confounded across only four lines, so they cannot be separated. Both gaps leave room for an accessibility-based explanation.

---

## 2. Molina & Zhang 2026 — Perturbation response decomposition  `[molina2026decomposition]`

**Citation.** Molina A, Zhang X. Perturbation response decomposition enables biologically aligned generalization to unseen perturbations and cellular contexts. *bioRxiv* (posted 27 Jul 2026; "new results", bioinformatics, CC-BY). DOI: **10.64898/2026.07.24.740459**. Corresponding author: Xinyi Zhang (AITHYRA, Austrian Academy of Sciences, Vienna). Not peer-reviewed.

**Summary.** The authors decompose pseudobulk responses as δ(c,p,g) = μ_g + α_{c,g} + β_{p,g} + γ_{c,p,g}. The terms are a global mean, a cell-line template, a perturbation effect conserved across lines, and a **cell-line × perturbation interaction**. The data are the four Replogle/Nadig CRISPRi screens (K562 and RPE1 from Replogle 2022; HepG2 and Jurkat from Nadig 2025), built into a balanced tensor of **4 lines × 638 perturbations × 2,000 genes**. After noise correction, response variance splits into cell-line template ~27.8–27.9%, conserved perturbation 29.4%, **interaction 23.5%** and noise ~19.2–19.3%. The two figures differ slightly between the main text and Supplementary Table 1. Once the template is removed, the interaction accounts for up to ~27.7% (32.6% under ANOVA accounting). Its per-line magnitude is largest in RPE1. The gene-level features used are **DepMap CRISPRGeneEffect coessentiality profiles** (1,208 lines; PCA50) as perturbation descriptors, and **control (baseline) expression**. Neither feature recovers γ in zero-shot cross-line transfer: correlations with γ are ≈0 or negative for Ridge, MLP, STATE and MORPH. γ becomes partly learnable only when ~30% of target-line perturbations are added, and this also works with disjoint perturbation sets, which indicates learnable, generalizable interaction structure. **No chromatin or accessibility features are used**: the full text has no hits for "chromatin" or "ATAC".

**What it leaves open.** The paper establishes that roughly a quarter of reproducible response variance is a line × perturbation interaction. It also shows that this interaction cannot be predicted from baseline expression or DepMap. It does not test whether line-specific regulatory or chromatin state, such as bulk ATAC at target or response-gene loci, explains γ. That is exactly the gap this project targets, and its γ definition and balanced 638-perturbation tensor are a direct template for us.

---

## 3. Ruan et al. 2026 — ChromaPert  `[ruan2026chromapert]`

**Citation.** Ruan J, He C, Quan R, Xue F, Yang Z, Yang Y. Beyond the Transcriptome: Chromatin-Informed Prediction of Cell-State-Dependent Perturbation Responses. *arXiv* 2610.09054 [q-bio.GN], v1 6 Oct 2026. DOI: 10.48550/arXiv.2610.09054 (arXiv DOI; we did not resolve it via DataCite). Not peer-reviewed.

**Summary.** ChromaPert predicts perturbation responses by combining molecular and DNA-informed target-locus priors with paired control RNA+ATAC features. A "Chromatin-Guided Response Router" retrieves control-derived expression offsets by perturbation/locus or control-ATAC similarity, and conditional transport flow models the remaining response. It is evaluated on two paired RNA–ATAC perturbation datasets:
- **Perturb-Multiome** (Martin-Rufino et al. 2025, GSE274113): 19 TF targets in *primary human hematopoietic cells undergoing erythroid differentiation*, 7 held-out states. **This is not K562.**
- **CAT-ATAC** (Shevade et al. 2025, GSE288996): 18 targets in DMSO-treated **K562**, roughly 12k cells.

ChromaPert reports the best all-gene response correlation among the benchmarked methods. On held-out states it improves correlation by 23.2% across all genes and 46.2% for the top-20 responsive genes. Correct RNA–ATAC pairing beats shuffled pairing in recovering chromatin-associated response slopes. The authors state that the main limitation is the need for paired scRNA+scATAC controls with a promoter-accessibility mapping. They also state: "Future work will extend ChromaPert to unpaired or single-modality inputs and incorporate distal regulatory information."

*Correction to prior notes:* the "within-line cell states" are erythroid differentiation states in primary cells (Perturb-Multiome). The K562 data are CAT-ATAC.

**What it leaves open.** ChromaPert works at small scale (about 18–19 TF targets), within a single lineage or line, and needs paired single-cell multiome data. It does not address cross-cell-line divergence across K562, RPE1, Jurkat and HepG2, genome-scale essential-gene panels, bulk ATAC, or distal (enhancer) accessibility. It is a prediction method, not an inferential test of how much variance accessibility explains.

---

## 4. Metzner, Southard & Norman 2025 — Multiome Perturb-seq  `[metzner2025multiome]`

**Citation.** Metzner E, Southard KM, Norman TM. Multiome Perturb-seq unlocks scalable discovery of integrated perturbation effects on the transcriptome and epigenome. *Cell Systems* 16(1):101161 (2025). DOI: **10.1016/j.cels.2024.12.002**. Preprint: bioRxiv 10.1101/2024.07.26.605307.

**Summary.** Multiome Perturb-seq extends single-cell CRISPR screens so that perturbation-induced changes in gene expression and chromatin accessibility are measured in the same nuclei (10x Multiome). An improved capture of sgRNA barcode transcripts from nuclear RNA makes guide assignment efficient. The authors ran a CRISPRi screen of **13 chromatin remodelers in human RPE-1 cells**. They organized expression and accessibility changes into joint programs and found that ARID1A and SUZ12 knockdowns induce programs enriched for developmental features. Modeling perturbation-induced heterogeneity linked accessibility changes to expression changes. (Details from the Cell Systems abstract/highlights; the full text was not fetched.)

**What it leaves open.** The study measures how a perturbation *changes* accessibility, in one line (RPE1) and for 13 targets. It does not ask whether *baseline* accessibility differences between lines explain differences in transcriptional responses to the same perturbation. Its RPE1 multiome data are a possible within-line validation resource for our project.

---

## 5. Replogle et al. 2022 — Genome-scale Perturb-seq  `[replogle2022mapping]`

**Citation.** Replogle JM, Saunders RA, Pogson AN, Hussmann JA, Lenail A, Guna A, Mascibroda L, Wagner EJ, Adelman K, Lithwick-Yanai G, Iremadze N, Oberstrass F, Lipson D, Bonnar JL, Jost M, Norman TM, Weissman JS. Mapping information-rich genotype-phenotype landscapes with genome-scale Perturb-seq. *Cell* 185(14):2559–2575.e28 (2022). DOI: **10.1016/j.cell.2022.05.013**.

**Summary.** This is the first genome-scale CRISPRi Perturb-seq. It used multiplexed (dual-sgRNA) CRISPRi in K562, targeting all expressed genes (9,866 per Nadig et al.), and >2.5 million cells were profiled. It also includes essential-gene screens in K562 (2,057) and RPE1 (2,393). The paper uses the transcriptional phenotypes to assign functions to poorly characterized genes and to define complexes and pathways. It also studies phenotypes such as transcriptional and chromosomal instability and RNA-processing defects (e.g., 3′ processing and Integrator). The data are the K562 and RPE1 inputs for both TRADE and the decomposition preprint. (Scale and focus are from our knowledge of the paper; Crossref has no abstract. The knockdown counts are confirmed through the Nadig preprint.)

**What it leaves open.** The paper compares K562 and RPE1 only descriptively. It does not explain line-specific response differences mechanistically and does not connect them to chromatin context.

---

## 6. Wu & Smyth 2012 — CAMERA  `[wu2012camera]`

**Citation.** Wu D, Smyth GK. Camera: a competitive gene set test accounting for inter-gene correlation. *Nucleic Acids Research* 40(17):e133 (2012). DOI: **10.1093/nar/gks461**.

**Summary.** CAMERA is a competitive gene-set test: are genes in the set more differentially expressed than genes outside it? It estimates the average inter-gene correlation within the set and inflates the variance of the set statistic by the variance inflation factor 1 + (m−1)ρ̄. This corrects the severe false-positive inflation that standard competitive tests show when set members are correlated. It works within limma's linear-model framework, is fast, and remains valid for small sample sizes. It is implemented as `camera()` (and later `cameraPR()`) in limma.

**What it leaves open.** It assumes a single exchangeable average correlation and tests sets defined in advance. It does not give a variance-explained estimate for a continuous gene-level covariate such as accessibility. Our use would be as a correlation-aware competitive test, for example whether high-accessibility gene sets are enriched for line-divergent responses.

---

## 7. Efron & Tibshirani 2007 — Testing significance of gene sets  `[efron2007testing]`

**Citation.** Efron B, Tibshirani R. On testing the significance of sets of genes. *The Annals of Applied Statistics* 1(1):107–129 (2007). DOI: **10.1214/07-AOAS101**.

**Summary.** The paper starts from GSEA and analyses gene-set testing in generality. It distinguishes the null hypotheses implied by permuting samples from those implied by randomizing genes. It proposes the **maxmean** statistic for summarizing a gene set and **restandardization**, which combines gene randomization with sample permutation to get more accurate inference. It also covers extensions, including using gene-set scores for class prediction, and introduces the R package GSA.

**What it leaves open.** This is the statistical basis for the choice between competitive and self-contained nulls in our tests. It does not handle inter-gene correlation in the way CAMERA does, and it does not handle continuous covariates.

---

## 8. Candès, Fan, Janson & Lv 2018 — Model-X knockoffs / CRT  `[candes2018panning]`

**Citation.** Candès E, Fan Y, Janson L, Lv J. Panning for gold: 'model-X' knockoffs for high dimensional controlled variable selection. *Journal of the Royal Statistical Society Series B: Statistical Methodology* 80(3):551–577 (2018). DOI: **10.1111/rssb.12265**.

**Summary.** The paper introduces model-X knockoffs, which give finite-sample FDR control for variable selection when the conditional distribution of Y|X is arbitrary and unknown. The condition is that the joint distribution of the covariates X is known (or well estimated), and the result holds regardless of the number of covariates. Knockoffs are built probabilistically rather than geometrically. The paper also introduces the **conditional randomization test (CRT)**: resample X_j from its conditional distribution given X_{−j} to get an exact p-value for conditional independence. The authors give preliminary evidence of robustness when the X distribution is estimated.

**What it leaves open.** Validity depends on modeling the covariate distribution, and the CRT is expensive to compute. For our project the question is how to model accessibility given confounders (e.g., expression level, gene length, essentiality) so that a CRT for "accessibility explains divergence beyond confounders" is valid.

---

## 9. Ahlmann-Eltze, Huber & Anders 2025 — DL vs linear baselines  `[ahlmanneltze2025deep]`

**Citation.** Ahlmann-Eltze C, Huber W, Anders S. Deep-learning-based gene perturbation effect prediction does not yet outperform simple linear baselines. *Nature Methods* 22(8):1657–1661 (2025). DOI: **10.1038/s41592-025-02772-6**. Preprint: bioRxiv 10.1101/2024.09.16.613342 (titled "Deep learning-based predictions of gene perturbation effects do not yet outperform simple linear baselines").

**Summary.** The authors benchmarked five single-cell foundation models and two other deep-learning models (e.g., GEARS, scGPT, scFoundation) against deliberately simple baselines on predicting transcriptome changes after single or double genetic perturbations. The baselines were the "additive" model for double perturbations, the mean, and a linear model with gene and perturbation embeddings. No deep model beat the baselines. A linear model using embeddings pretrained on perturbation data did as well as or better than the foundation models. The paper argues for critical benchmarking in method development.

**What it leaves open.** It concerns prediction accuracy within a cell line, not cross-line divergence. It also does not ask whether any biological prior, chromatin included, carries information beyond expression. It motivates our use of simple linear or variance-partitioning models as the primary analysis.

---

## Verification table

| # | Work | DOI | Status |
|---|---|---|---|
| 1 | Nadig et al. 2025 Nat Genet 57:1228–1237 | 10.1038/s41588-025-02169-3 | Verified (Crossref). Numbers from preprint 10.1101/2024.07.03.601903. 02177-3 = Peidli News & Views |
| 2 | Molina & Zhang 2026 bioRxiv | 10.64898/2026.07.24.740459 | Verified (Crossref + bioRxiv API + full-text PDF) |
| 3 | Ruan et al. 2026 arXiv 2610.09054 | 10.48550/arXiv.2610.09054 | Verified (arXiv API + PDF). DOI is the arXiv convention, not resolved |
| 4 | Metzner, Southard & Norman 2025 Cell Syst 16:101161 | 10.1016/j.cels.2024.12.002 | Verified (Crossref) |
| 5 | Replogle et al. 2022 Cell 185:2559–2575.e28 | 10.1016/j.cell.2022.05.013 | Verified (Crossref) |
| 6 | Wu & Smyth 2012 NAR 40:e133 | 10.1093/nar/gks461 | Verified (Crossref) |
| 7 | Efron & Tibshirani 2007 AOAS 1:107–129 | 10.1214/07-AOAS101 | Verified (Crossref) |
| 8 | Candès et al. 2018 JRSS-B 80:551–577 | 10.1111/rssb.12265 | Verified (Crossref) |
| 9 | Ahlmann-Eltze et al. 2025 Nat Methods 22:1657–1661 | 10.1038/s41592-025-02772-6 | Verified (Crossref) |
