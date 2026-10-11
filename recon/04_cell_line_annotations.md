# 04 — Cell-line annotations: TP53 status, copy-number sources, p53 targets, gene/CpG/functional-group annotations

Compiled 2026-10-11. Lines: K562 (CVCL_0004, DepMap ACH-000551), hTERT-RPE1 (CVCL_4388; not in DepMap as the parental line), Jurkat E6-1 (CVCL_0367; DepMap ACH-000995 is the parental "JURKAT", CVCL_0065), HepG2 (CVCL_0027, DepMap ACH-000739).

Downloaded files live in `data/annot/` (git-ignored via `.gitignore: data/`). MD5 + size + date for every file is in `data/annot/MANIFEST_md5.tsv`; key entries are repeated below.

---

## 1. TP53 status per line

| Line | TP53 status | Variant (HGVS, Cellosaurus) | Functional consequence | Primary citation(s) | Secondary / confirmatory |
|---|---|---|---|---|---|
| **K562** | **Null** (homozygous frameshift) | c.406_407insC, p.Gln136fs*13 (exon 5); homozygous | No full-length p53; p53 pathway inactive | Ikediobi ON et al. *Mol Cancer Ther* 2006;5:2606 (NCI-60 mutation survey) doi:10.1158/1535-7163.MCT-06-0433 (PMID 17088437) | Berglind H et al. *Cancer Biol Ther* 2008;7:699 (PMID 18277095); Zhou B et al. *Genome Res* 2019;29:472 (K562 WGS; PMID 30737237); Cellosaurus CVCL_0004 |
| **hTERT-RPE1** | **Wild-type** (functional) | none reported | Intact p53-p21 axis: p53-dependent arrest after centrosome loss/prolonged mitosis and after ribosome/essential-gene knockdown | Lambrus BG et al. *J Cell Biol* 2016;214:143 doi:10.1083/jcb.201604054 (PMID 27432896) — p53/p21-dependent arrest in RPE1 | Bodnar AG et al. *Science* 1998;279:349 (origin; PMID 9454332); Replogle JM et al. *Cell* 2022;185:2559 (RPE1 Perturb-seq; PMID 35688146); Nadig A et al. *Nat Genet* 2025;57:1228 (TRADE; PMID 40259084) — uses the Replogle K562/RPE1 atlases and contrasts RPE1 (p53-WT, non-transformed) with K562 (p53-null) — **see blocker B5: wording not verified from full text** |
| **Jurkat E6-1** | **Mutant** (heterozygous nonsense; functionally p53-deficient) | c.586C>T, p.Arg196Ter (R196*), exon 6; rs397516435; allele fraction ≈0.25 on a hypotetraploid background | Truncated DNA-binding-domain product; no canonical p53 transcriptional response reported; possible p53ψ-like isoform | Cheng J & Haas M. *Mol Cell Biol* 1990;10:5502 doi:10.1128/MCB.10.10.5502 (first report of R196*) | Gioia L et al. *BMC Genomics* 2018;19:334 doi:10.1186/s12864-018-4718-6 (PMID 29739316; WGS confirms heterozygous R196*; note: this is **BMC Genomics, not Sci Data**); Cellosaurus CVCL_0367 |
| **HepG2** | **Wild-type** | none; NRAS p.Gln61Leu and TERT promoter C228T are the drivers | Functional p53 (classic wt-p53 control line vs Hep3B null) | Bressac B et al. *PNAS* 1990;87:1973 doi:10.1073/pnas.87.5.1973 (PMID 2155427) | Zhou B et al. *Nucleic Acids Res* 2019;47:3846 (HepG2 WGS; PMID 30864654); Cellosaurus CVCL_0027; DepMap ACH-000739 (no TP53 call) |

Caveats to carry into the analysis:
- Jurkat is not p53-*null*: WGS shows R196* at ~1/4 allelic fraction with the remaining alleles apparently WT at that position (Gioia 2018), yet Jurkat is consistently reported p53-deficient functionally. Treat as "p53-mutant/deficient" and, if a p53-target score is used as a covariate, verify empirically (e.g. CDKN1A/MDM2 response to ribosomal-protein knockdown in the Jurkat Perturb-seq data).
- K562 TP53 variant coordinates in Cellosaurus are protein-level; the DepMap OmicsSomaticMutations table (339 MB, not downloaded) carries the GRCh38 genomic coordinates if needed.

---

## 2. Copy-number sources per line (for ATAC / expression CN-correction, GRCh38)

### Chosen primary source for the three cancer lines: DepMap 24Q4 WGS/WES relative CN segments (GRCh38)
- Release: **DepMap Public 24Q4** (Broad), Figshare+ DOI 10.25452/figshare.plus.27993248.v1, published 2024-12-10. **Licence: CC BY 4.0.** Citation: "DepMap, Broad (2024). DepMap 24Q4 Public. Figshare+."
- Assembly: GRCh38/hg38 (README: GATK CN pipeline aligned to hg38).
- Files (all retrieved 2026-10-11 from `https://ndownloader.figshare.com/files/<id>`; md5 matches Figshare `computed_md5`):

| File | Figshare id | bytes | md5 | content |
|---|---|---|---|---|
| `depmap_24q4/OmicsCNSegmentsProfile.csv` | 51065333 | 41,098,979 | 79fc5cf5c46f215ef085dddaa0e56860 | segment-level *relative* copy ratio (`SegmentMean`, 1.0 = ploidy-normalised diploid), all profiles |
| `depmap_24q4/OmicsAbsoluteCNSegmentsProfile.csv` | 51065321 | 33,904,573 | e5d3a4425a3ce6d77c42ca75cd13b71b | PureCN absolute total + minor-allele CN, LOH status (high-ploidy samples omitted) |
| `depmap_24q4/Model.csv` | 51065297 | 645,696 | 675210d17675f3517b0ce39a3c274f16 | model metadata (ModelID ↔ RRID) |
| `depmap_24q4/OmicsProfiles.csv` | 51065723 | 254,733 | 7195088a5fc41eddb3f896c6e984a992 | ProfileID ↔ ModelID ↔ datatype (wgs/wes) |
| `depmap_24q4/OmicsDefaultModelProfiles.csv` | 51065339 | 90,080 | 4d342a03701c0a4c0742fedb56c1de40 | default DNA profile per model |
| `depmap_24q4/README.txt` | 51065795 | 43,103 | 54628c15a10b9ff6536db245cfed5231 | column definitions |

- Subsets written for this project (same licence):
  - `depmap24q4_CNsegments_K562_HepG2_Jurkat_RPE1clones.csv` — md5 9959eb6398299aadd75085c941acdffd (176,899 B). Labels: `K562_wgs` (PR-aheaZL, 817 seg; default), `K562_wes_broad`, `K562_wes_sanger`, `HEPG2_wgs` (PR-ozp1kx, 165 seg; default), `JURKAT_wes_broad` (PR-Ds9QdK, 69 seg; default), `JURKAT_wes_sanger`, `RPE1-ss6_wgs`, `RPE1-ss77_wgs`, `RPE1-ss119_wgs`, `RPE1-ss48_wes`, `RPE1-ss51_wes`, `RPE1-ss111_wes`.
  - `depmap24q4_AbsoluteCNsegments_K562_HepG2_Jurkat_RPE1clones.csv` — md5 60163c03f766c5b3891bec1676dde20a (80,045 B). **No absolute CN for K562_wgs or HepG2_wgs** (PureCN omitted as likely high-ploidy); Jurkat absolute CN was fitted as diploid (chr1 = 2.0), which contradicts the hypotetraploid karyotype — use relative CN only.
- Sanity checks on the relative segments (this run): K562 shows 9p loss (CDKN2A region mean 0.74) and 22q/BCR-ABL-region gain; HepG2 ~20 % of genome gained, 7 % lost; Jurkat nearly flat (as expected for a tetraploid genome measured as relative CN); RPE1-ss6 shows the known 10q gain (10q21.2→qter segment means 1.5–1.6 ≈ 3 copies).
- Gene-level matrix (`OmicsCNGene.csv`, 1.39 GB; `PortalOmicsCNGeneLog2.csv`, 1.40 GB; Figshare ids 51065324 / 51065765) was **not** downloaded (size). For gene-level CN, intersect the segment table with GENCODE gene bodies/TSS (weighted mean, as DepMap does) — trivial with bedtools/pyranges.

### Per-line recommendation

| Line | Best source (chosen) | Format / assembly / licence | Fallback |
|---|---|---|---|
| **K562** | DepMap 24Q4 WGS profile PR-aheaZL (817 segments; relative CN) — in `depmap24q4_CNsegments_*.csv` | CSV segments; GRCh38; CC BY 4.0 | (1) Zhou B et al. *Genome Res* 2019 (PMID 30737237) — ENCODE K562 WGS/10x phased genome; CN by 10-kb bins in paper supplement (hg19; needs liftOver). ENCODE deposits only VCFs for these experiments: ENCSR025GPQ, ENCSR711UNY (Illumina WGS), ENCSR053AXS (10x linked reads) — all hg19. (2) ENCODE HAIB genotyping-array CNV bed **downloaded**: `encode/ENCFF486MJU_K562_CNV_hg19.bed.gz` (ENCSR000BBO; 537 segments; hg19; md5 b4a78fa77789291cb5e3d528b32aa55a; ENCODE data: free reuse, CC0-equivalent "ENCODE data release policy"). Requires UCSC liftOver hg19→hg38 (`hg19ToHg38.over.chain.gz`). (3) 4DN K562 Hi-C (4DNSRP5XBV7W / 4DNSRWQFKJYO) for SV-aware checks. |
| **HepG2** | DepMap 24Q4 WGS profile PR-ozp1kx (165 segments) | CSV segments; GRCh38; CC BY 4.0 | (1) Zhou B et al. *NAR* 2019;47:3846 (PMID 30864654) — haplotype-resolved HepG2 genome; CN of chromosomal segments in supplement (hg19). ENCODE WGS VCFs ENCSR276ECO, ENCSR350RUO, ENCSR319QHO (hg19). (2) ENCODE HAIB array CNV bed **downloaded**: `encode/ENCFF074XLG_HepG2_CNV_hg19.bed.gz` (ENCSR000BBN; 428 segments; md5 8fa76109277f2f2e970c92dab43a5d7a; hg19 → liftOver). |
| **Jurkat E6-1** | DepMap 24Q4 WES profile PR-Ds9QdK (ACH-000995; relative CN; 69 segments). Caveat: WES-derived and relative, on a tetraploid baseline; profile is of parental Jurkat, not the E6-1 subclone (near-identical genomes) | CSV segments; GRCh38; CC BY 4.0 | Gioia L et al. *BMC Genomics* 2018 (PMID 29739316): WGS (SRA SRP101994) with CNVnator CN calls (modal CN 4); variant/CNV calls on Zenodo 10.5281/zenodo.400615 — `jurkat_final_variant_calls.tar.gz` 338 MB, CC BY 4.0, **hg19** (not downloaded: >50 MB). Second fallback: DepMap Sanger WES profile PR-Yz6AIC (95 segments) in the same subset file. |
| **hTERT-RPE1** | **No parental-line CN table exists in DepMap/CCLE.** Chosen: treat as near-diploid 46,XX with a single large gain — **10q21.2→qter at 3 copies** via der(X)t(X;10)(Xq28;10q21.2), plus corresponding Xq28-ter loss on one X — per Volpe M et al. *Nat Commun* 2025 "The reference genome of the human diploid cell line RPE-1" (doi:10.1038/s41467-025-62428-z; karyotype 46,X,der(X)(Xpter→Xq28::10q21.2→10qter), confirmed by G/R-banding and HiFi coverage). Operationally use the DepMap RPE1 single-cell-clone WGS profiles (RPE1-ss6 PR-g5QF67, RPE1-ss77 PR-Yf6hXD, RPE1-ss119 PR-l4NGGv; from Cohen-Sharir Y et al. *Nature* 2021;590:486, PMID 33505028; Cellosaurus CVCL_C8DQ etc.) in `depmap24q4_CNsegments_*.csv` as an empirical GRCh38 profile — **take the consensus of the three clones** (segments shared by all; clone-private aneuploidies are artefacts of single-cell cloning). | CSV segments; GRCh38; CC BY 4.0 (DepMap) / karyotype from CC BY paper | (1) RPE1v1.1 diploid assembly + Illumina/HiFi coverage (Volpe 2025; data in ENA/GEO per paper — accession not captured here, see blocker B3). (2) 4DN hTERT-RPE1 biosource 4DNSRHGVFSRJ (Hi-C) for SV sanity. (3) Simplest fallback: synthetic CN track = 2 everywhere, 3 on chr10:≈60.8 Mb (GRCh38 10q21.2 boundary observed in RPE1-ss6 at 60,799,376)–qter, 1 on chrX:q28-ter. |

Recommended processing: convert relative `SegmentMean` to log2 ratio (log2(SegmentMean)); map to GENCODE TSS ±1 kb (ATAC promoter peaks) and gene bodies (expression) by weighted overlap; use as an offset/covariate. For Jurkat and K562, do not interpret relative CN as absolute ploidy.

---

## 3. p53 target gene sets

| Set | n genes | Source / URL | File downloaded | md5 | Licence |
|---|---|---|---|---|---|
| **Fischer 2017 "116 core" p53-activated targets (Table 1: identified in ≥6 of 16 genome-wide datasets)** — *recommended primary set* | 116 | Fischer M. *Oncogene* 2017;36:3943 doi:10.1038/onc.2016.502; table https://www.nature.com/articles/onc2016502/tables/1 | `fischer2017_oncogene_table1_116_p53_targets.tsv` (symbol + no. of datasets; parsed from the HTML table 2026-10-11) | 02ed713780a61919f485e36c79b07108 | Article text © Springer Nature (not OA); gene-list facts are reusable; cite the paper. The full 343-gene list is in Supplementary Table S3 of the same paper (not downloaded; Springer ESM behind the article page). |
| Fischer meta-analysis direct p53 targets (MSigDB C2:CGP `FISCHER_DIRECT_P53_TARGETS_META_ANALYSIS`, M61) — **note: derived from Fischer 2016 *NAR* Table S3 (PMID 27280975), not the 2017 Oncogene paper** | 311 | https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/FISCHER_DIRECT_P53_TARGETS_META_ANALYSIS | `msigdb/FISCHER_DIRECT_P53_TARGETS_META_ANALYSIS.v2025.1.Hs.json` | f655ea02cb1587632447ff4094844d25 | MSigDB CC BY 4.0 (v2022.1+) |
| MSigDB `HALLMARK_P53_PATHWAY` (M5939; Liberzon 2015) — broad pathway, use only as a secondary check | 200 | https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_P53_PATHWAY | `msigdb/HALLMARK_P53_PATHWAY.v2025.1.Hs.json` | 8ae500379223ea1aaa6c922e357468cd | MSigDB CC BY 4.0 |
| Andrysik 2017 core TP53 program ("∼100" direct targets activated across HCT116/MCF7/SJSA via Nutlin; GRO-seq+ChIP-seq) — alternative with mechanistic (direct-binding) support | ~100 | Andrysik Z et al. *Genome Res* 2017;27:1645 doi:10.1101/gr.220533.117 (PMID 28904012); Supplemental Table S? at genome.cshlp.org | not downloaded (supplement behind CSHL site; see B6) | — | Genome Research OA (CC BY-NC after 6 mo) |

Download URL pattern for any single MSigDB set without login: `https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?geneSetName=<NAME>&fileType=json` (or `gmt`/`grp`).

---

## 4. Gene annotation, TSS, promoter CpG

| Item | Choice | URL | Size / status | Licence |
|---|---|---|---|---|
| Gene annotation / TSS | **GENCODE Human Release 50 (GRCh38.p14, June 2026; current)**. Use `gencode.v50.basic.annotation.gtf.gz` (80.0 MB) or the comprehensive `gencode.v50.annotation.gtf.gz` (124.5 MB). TSS = 5′ end of each transcript; for a per-gene TSS take Ensembl canonical (`tag "Ensembl_canonical"`) transcript. Match to the Perturb-seq reference: Replogle 2022 used GENCODE v32-era Cell Ranger refs — map by Ensembl gene ID (version-stripped). | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.basic.annotation.gtf.gz ; https://www.gencodegenes.org/human/ | HEAD verified 2026-10-11 (Last-Modified 2026-06-09); **not downloaded (>50 MB)** | GENCODE: freely available, no restrictions (Ensembl/EMBL-EBI terms; effectively public domain/CC0-like) |
| CpG islands | UCSC hg38 `cpgIslandExt` (Gardiner-Garden & Frommer criteria; 32,038 islands) | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz | `ucsc/hg38_cpgIslandExt.txt.gz` 717,984 B, md5 17da03fd8be6e33d2c09228df677a0f8 | UCSC annotation tracks: free for all uses (CpG island track has no restrictions) |
| Promoter CpG content (continuous) | Compute directly: CpG observed/expected and GC% in TSS ±500 bp (Saxonov 2006 HCP/ICP/LCP classes: HCP if CpG o/e >0.75 & GC >0.55 in a 500-bp window) from the GRCh38 primary assembly FASTA (`GRCh38.primary_assembly.genome.fa.gz`, same GENCODE release dir, ~800 MB) | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/GRCh38.primary_assembly.genome.fa.gz | not downloaded | as above |

---

## 5. Knockdown functional groups (ribosome, proteasome, spliceosome, …)

| Resource | Use | URL | File downloaded | md5 | Licence |
|---|---|---|---|---|---|
| **CORUM 5.3** (2026-04-14) human complexes | Complex-level groups: 26S proteasome, 40S/60S ribosome, spliceosome (U1/U2/tri-snRNP, NTC), Mediator, Pol II, exosome, cohesin, condensin, NPC, etc. Subunits by gene symbol (`subunits_gene_name`) | API (no login): https://mips.helmholtz-muenchen.de/fastapi-corum/public/file/download_current_file?file_id=human&file_format=txt ; site https://mips.helmholtz-muenchen.de/corum/ | `corum/humanComplexes_corum5.3.txt` (6,268,940 B; 7,799 complexes) | 5cc8f6a98de30907f99f80adafce704a | **CC BY-NC 4.0** (non-commercial) |
| **Reactome** (release current as of 2026-06-21) human pathway GMT | Pathway-level groups: Eukaryotic Translation Initiation/Elongation/Termination (R-HSA-72613/156842/72764), mRNA Splicing – Major Pathway (R-HSA-72163), Proteasome assembly (R-HSA-9907900), Mitochondrial translation (R-HSA-5368287), Respiratory electron transport (R-HSA-611105), NMD (R-HSA-927802), RNA Pol II Transcription Initiation (R-HSA-75953), Golgi/ER transport … | https://reactome.org/download/current/ReactomePathways.gmt.zip | `reactome/ReactomePathways.gmt` (1,032,186 B; 2,868 sets) + `.zip` | 26b6c62fa4ac0904738987b64d0575ca | **CC0 1.0** (Reactome data) |
| MSigDB 2025.1.Hs C2:CP:REACTOME GMT | Same as above but frozen with MSigDB versioning | https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/c2.cp.reactome.v2025.1.Hs.symbols.gmt | `msigdb/c2.cp.reactome.v2025.1.Hs.symbols.gmt` (1,787 sets) | 25813c4bed9d944b00c51b7532cf14f2 | CC BY 4.0 |
| MSigDB 2025.1.Hs C5:GO:CC GMT | Cellular-component groups (GOCC_CYTOSOLIC_RIBOSOME, GOCC_PROTEASOME_COMPLEX, GOCC_SPLICEOSOMAL_COMPLEX, GOCC_MITOCHONDRIAL_RIBOSOME…) | https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2025.1.Hs/c5.go.cc.v2025.1.Hs.symbols.gmt | `msigdb/c5.go.cc.v2025.1.Hs.symbols.gmt` (1,042 sets) | f2dd8896a671d015be4a247f32567198 | CC BY 4.0 |
| MSigDB 2025.1.Hs C2:CP:KEGG_MEDICUS GMT | Optional cross-check | …/c2.cp.kegg_medicus.v2025.1.Hs.symbols.gmt | `msigdb/c2.cp.kegg_medicus.v2025.1.Hs.symbols.gmt` (658 sets) | 6896df38436f9f4f6a22ee99c2762876 | KEGG MEDICUS sets carry **additional KEGG terms** (non-commercial; see MSigDB licence page) |
| **HGNC gene groups** | Clean, curated membership lists: "L ribosomal proteins" (id 729), "S ribosomal proteins" (728), "Ribosomal proteins" (1054), "Mitochondrial ribosomal proteins" (646), "Proteasome" (690), "Sm spliceosomal proteins" (1255), "Mediator complex" (1061), "Nucleoporins" (1051), "Exosome complex" (817), "Cohesin complex" (1060), "Condensin subunits" (1454), "SAGA complex" (1059), "Anaphase promoting complex" (402) | https://storage.googleapis.com/public-download-files/hgnc/csv/csv/genefamily_db_tables/family.csv and `gene_has_family.csv` (join on family id → HGNC id; map HGNC id → symbol via `hgnc_complete_set`) | `hgnc/family.csv` (1,028,349 B), `hgnc/gene_has_family.csv` (519,501 B) | 1f1528d0228ecd99f5c515d7beadaf54 / c7861b43a465a5a1d5d549defb1f0b20 | **CC0** (HGNC data) |
| Replogle 2022 perturbation clusters | The Perturb-seq papers' own functional cluster labels (Cell 2022 Table S?) — use for direct comparability of "knockdown groups" with the published atlas | https://doi.org/10.1016/j.cell.2022.05.013 (supplementary tables; also gwps.wi.mit.edu) | not downloaded (Cell ESM, size/licence unverified) | Cell Press, CC BY 4.0 (OA article) |

Suggested canonical group definitions (to fix in `src/`): ribosome = HGNC 729 ∪ 728 (cytosolic RPL/RPS, excluding RPL22L1-type paralogues as desired); mito-ribosome = HGNC 646; proteasome = HGNC 690 (PSMA/B/C/D/E/F/G); spliceosome = CORUM "spliceosome"-named complexes ∪ Reactome R-HSA-72163; translation initiation = Reactome R-HSA-72613 minus ribosomal proteins; Pol II/Mediator = HGNC 1061 + POLR2*; OXPHOS = Reactome R-HSA-611105.

---

## 6. Blockers / open items

- **B1 — DepMap portal is Cloudflare-gated for scripts.** `depmap.org/portal/api/download/files` and `data_page` return a Turnstile "Verification" page to curl; the newest releases (25Q2–26Q1, current = 26Q1 dated 2026-04-01) are *only* there (Figshare+ stops at 24Q4). We used 24Q4 from Figshare+ (open API, CC BY 4.0). CN pipelines are stable across releases for these old lines; if 26Q1 is wanted, download `OmicsCNSegmentsProfile.csv` manually in a browser and drop it in `data/annot/depmap_26q1/`.
- **B2 — No absolute copy number for K562/HepG2 WGS in DepMap** (PureCN omitted high-ploidy samples) and Jurkat absolute CN mis-fitted as diploid. Relative CN is sufficient for CN-correction; absolute ploidy must come from karyotype literature (K562 near-triploid; HepG2 hyperdiploid modal 55; Jurkat hypotetraploid; RPE1 diploid).
- **B3 — RPE1 parental CN track is not in any portal.** We rely on (i) the Volpe 2025 *Nat Commun* karyotype/assembly paper for the single structural event (10q gain via t(X;10)) and (ii) DepMap RPE1 single-cell clones (Cohen-Sharir 2021) as an empirical GRCh38 CN profile. Clone-private aneuploidies must be removed by taking the 3-clone consensus. The raw-data accession for the RPE1v1.1 assembly was not captured (article text truncated in extraction) — check the paper's Data Availability if coverage-based CN from the parental line is required.
- **B4 — ENCODE K562/HepG2 CN tracks are hg19 only** (HAIB array CNV beds ENCSR000BBO/BBN; Zhou 2019 WGS supplements). They are downloaded but need UCSC liftOver to GRCh38 before use; they serve as fallbacks only. ENCODE did not deposit CN segment files from the Zhou 2019 WGS experiments (only VCFs).
- **B5 — Nadig 2025 (TRADE) full text not retrieved** (Europe PMC XML returned empty; Nature page paywalled for text extraction). The paper is cited as describing the K562/RPE1 atlases; the exact wording about p53 status should be checked in the PDF before quoting it. The K562-null / RPE1-WT facts themselves are supported by the primary citations in §1.
- **B6 — Supplementary tables not auto-downloadable:** Fischer 2017 Table S3 (343-gene list) and Andrysik 2017 supplement are behind publisher ESM pages; only Fischer Table 1 (116 genes) was captured from the HTML. Gioia 2018 Jurkat CNV calls (Zenodo, 338 MB, hg19) exceed the 50 MB cap.
- **B7 — Licence constraint:** CORUM is CC BY-NC 4.0 and MSigDB KEGG_MEDICUS carries KEGG non-commercial terms. Fine for an academic project; avoid redistributing derived complex tables commercially. Everything else used here is CC BY 4.0 / CC0 / unrestricted (ENCODE, GENCODE, UCSC CpG, HGNC, Reactome, DepMap).
- **B8 — Jurkat E6-1 vs parental Jurkat:** DepMap ACH-000995 is parental JURKAT (CVCL_0065); the Perturb-seq line is E6-1 (CVCL_0367). Gioia 2018 sequenced E6-1-derived cells (Cellosaurus assigns their variants to CVCL_0367). Differences are expected to be negligible for CN-correction purposes.
