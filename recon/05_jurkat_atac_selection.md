# 05 — Jurkat bulk ATAC-seq dataset selection (for ENCODE ATAC pipeline v1.9.1, GRCh38)

Date: 2026-10-11

## Search performed
- NCBI eutils `esearch db=sra term='(Jurkat[All]) AND (ATAC-seq[Strategy]) AND "Homo sapiens"[Organism]'` → **1003 SRA experiments** (all pages fetched via history server, retstart 0/500/1000) → **1139 runs in 50 SRA studies** (`recon/jurkat_atac_search/all_runinfo.json`; supersedes the 299-run `jurkat_atac_sra_runinfo.json`).
- GEO `acc.cgi?...&form=text` pulled for all 489 GSMs (excluding three scATAC/low-depth studies: SRP126342, SRP401017, SRP141347) and 24 key GSEs (`recon/jurkat_atac_search/gsm/`, `gse/`, `gsm_info.json`).
- Excluded: single-cell/multiplexed (GSE330650, GSE331391, GSE291229, GSE107816, GSE214812, GSE113414, GSE311565 ORF-screen), single-end (GSE132209, GSE199691, GSE276690), HIV-infected-only clones (2D10, J-Lat, Cm6 HIV clones), samples where every Jurkat arm is stimulated / KO / knock-in (GSE90718 has untreated arm → kept; GSE207607, GSE206376, GSE224204 (no Jurkat), GSE290905 GFP-TAL1 knock-in, GSE267359 FKBP-MYB knock-in, etc.), <25M or single-replicate controls (GSE162690, GSE206985, GSE302716, GSE136686, GSE115438, GSE129086, GSE145783, GSE135320, GSE130140, GSE294051 DMSO 16–19M).

## Top 5 candidates (ranked)

| Rank | GSE / BioProject (SRP) | Paper (DOI) | Condition used | Cell line | Reps | Control runs (GSM → SRR) | Read pairs / run | Read length | Protocol / platform |
|---|---|---|---|---|---|---|---|---|---|
| **1** | **GSE273644 / PRJNA1142562 (SRP523435)** | Gottardi Zamperla et al. 2024 *Clin Epigenetics* — 10.1186/s13148-024-01725-8 (PMID 39155390) | Vehicle (DMSO, 16 h) | **Jurkat clone E6-1 (ATCC TIB-152)** stated in GEO | 3 | GSM8434460→SRR30058502; GSM8434462→SRR30058500; GSM8434464→SRR30058498 | 40.1M; 31.4M; 34.8M | PE 2×42 (NextSeq 500) | Buenrostro-style Nextera Tn5 ATAC (Active Motif service; cryopreserved cells, DNase pre-treatment to remove dead-cell DNA), 10 PCR cycles |
| 2 | GSE155555 / PRJNA650287 (SRP275565) | Nasser et al. 2021 *Nature* (ABC enhancer maps) — 10.1038/s41586-021-03446-x | Unstimulated ("ctrl") | Jurkat (clone not stated) | 2 | GSM4706084→SRR12368302+SRR12368303; GSM4706085→SRR12368304+SRR12368305 | 27.4M (13.7+13.7); 35.8M (17.9+17.9) | PE 2×50 (HiSeq 2500) | Buenrostro 2013 w/ Nuclei EZ lysis, 50k cells |
| 3 | GSE90718 / PRJNA355577 (SRP094125) | Brignall et al. 2017 *J Immunol* — 10.4049/jimmunol.1602033 | Untreated | Jurkat **clone A3** (not E6-1) | 3 | GSM2411156→SRR5063984; GSM2411157→SRR5063985; GSM2411158→SRR5063986 | 33.3M; 29.9M; 47.6M | PE 2×101 (HiSeq 2500) | Buenrostro-style, Nuclei EZ lysis, 50k cells |
| 4 | GSE184345 / PRJNA764194 (SRP337615) | Shah, Jung, Corces lab — no PMID linked in GEO ("Multiomics profiling of HIV-infected cells") | Untreated parental Jurkat | Jurkat (clone not stated) | 2 | GSM5584797→SRR15931092; GSM5584798→SRR15931093 | 93.7M; 122.7M | PE 2×50 (HiSeq 2500) | Omni-ATAC (Corces 2017) |
| 5 | GSE173868 / PRJNA727329 (SRP318468) | Zhou et al. 2022 *Mol Cell* — 10.1016/j.molcel.2022.01.015 | Cas9 + control (non-targeting) sgRNA, sorted d6 | **Jurkat E6-1** | 3 | GSM5282145→SRR14423624; GSM5282146→SRR14423625; GSM5282147→SRR14423626 | 50.7M; 41.4M; 36.2M | PE 2×38 (NextSeq 500) | Buenrostro 2013 |

Honourable mentions (all PE, ≥2 control reps, ≥25M):
- GSE327360 / PRJNA1450595 (Mullighan lab, 2026, no PMID yet): Jurkat BCL11B-WT + DMSO 6 h, GSM9624240→SRR37814640 (87.8M), GSM9624241→SRR37814639 (153.1M); Omni-ATAC, 2×~34 NovaSeq.
- GSE284218 / PRJNA1198009 (Chang et al. 2025 *Cell*, 10.1016/j.cell.2025.02.002): "Pre-exposure" Jurkat, GSM8679858→SRR31708097 (64.0M), GSM8679859→SRR31708096 (77.6M); 2×150; TruePrep Tn5 kit.
- GSE206220 / PRJNA849647 (Lenardo lab, EBF4): empty-vector Jurkat ×3 (69–79M, 2×150) — lentiviral vector control.
- GSE229927 (Collora & Ho 2023 *Genome Res*): uninfected "Jurkat Cm6" ×3 (GSM7181451-3; 24.2/25.8/20.4M) — marginal depth.

## Bonus: Jurkat dCas9-KRAB / CRISPRi non-targeting controls
| GSE / BioProject | Paper | Samples | Read pairs | Notes |
|---|---|---|---|---|
| **GSE302221 / PRJNA1289670** | Weissman et al. 2026 *mBio* — 10.1128/mbio.02755-25 | Jurkat_Jur_**Uninfected**_CRISPRi_nontargeting rep1–3: GSM9100014→SRR34472975, GSM9100018→SRR34473029, GSM9100022→SRR34473025 | 56.5M; 46.6M; 59.0M | dCas9-KRAB + non-targeting gRNA, HIV-uninfected Jurkat, 3 reps, 2×151 NovaSeq. Best dCas9-KRAB control match. |
| GSE261695 (SuperSeries GSE220976) / PRJNA1088439 | Pacalin et al. 2025 *Nat Biotechnol* — 10.1038/s41587-024-02213-3 | ATAC NTC Rep1/Rep2: GSM8148870→SRR28353591, GSM8148871→SRR28353590 | 39.4M; 28.8M | "Jurkat CRISPRi", dCas9 + NTC gRNA, dox; 2×75 NextSeq 550. Effector (KRAB vs ZIM3-KRAB) not stated in GEO — check paper methods. |
| GSE203235 / PRJNA839152 | Ho lab | GSM6165071 dCas9-KRAB_NT (137.3M) | — | single rep, HIV-infected 1D7 clone → not usable as control set |

ENA FASTQ links for the bonus runs were verified (all PAIRED, two files each); see `recon/jurkat_atac_search/` if needed.

## Recommendation
**GSE273644 (PRJNA1142562), DMSO-vehicle arm, 3 biological replicates: SRR30058502, SRR30058500, SRR30058498.**

Justification against the priority list:
1. Control condition: DMSO vehicle only (16 h) — allowed.
2. Cell line: the only Jurkat ATAC dataset whose GEO metadata explicitly states **Jurkat clone E6-1 (ATCC TIB-152)**, unmodified (the CRISPR HDAC6-KO arm in this series is B16F10 mouse, not Jurkat).
3. 3 biological replicates (more than the ENCODE minimum of 2).
4. Paired-end Illumina (NextSeq 500), bulk Nextera-Tn5 ATAC (Buenrostro-style).
5. 31–40M read pairs per replicate (all ≥25M).
6. Published, peer-reviewed (Clin Epigenetics 2024).

Caveats: short reads (PE42) — fine for bowtie2 in the ENCODE pipeline but slightly lower mappability than 2×50+; Active Motif service protocol on cryopreserved cells with DNase pre-treatment (not Omni-ATAC); expect TSS enrichment to be checked after the run. If QC fails (e.g. TSS enrichment <5 or low library complexity), fall back to **GSE155555** (Nasser 2021, 2 reps, 2×50) or **GSE90718** (3 reps, 2×101, but clone A3). For a matched CRISPRi-background control, add **GSE302221** uninfected CRISPRi-nontargeting (3 reps, 47–59M).

## Selected runs
`recon/jurkat_atac_selected_runs.tsv` (from ENA filereport):

| replicate | GSM | SRR | spots | bases | layout | FASTQ (ENA) |
|---|---|---|---|---|---|---|
| rep1 | GSM8434460 (Jurkat-S13) | SRR30058502 | 40,111,050 | 3,275,490,692 | PAIRED | ftp.sra.ebi.ac.uk/vol1/fastq/SRR300/002/SRR30058502/SRR30058502_{1,2}.fastq.gz |
| rep2 | GSM8434462 (Jurkat-S14) | SRR30058500 | 31,419,225 | 2,606,213,687 | PAIRED | ftp.sra.ebi.ac.uk/vol1/fastq/SRR300/000/SRR30058500/SRR30058500_{1,2}.fastq.gz |
| rep3 | GSM8434464 (Jurkat-S15) | SRR30058498 | 34,816,240 | 2,886,651,267 | PAIRED | ftp.sra.ebi.ac.uk/vol1/fastq/SRR300/098/SRR30058498/SRR30058498_{1,2}.fastq.gz |

Total download ≈ 7.35 GB (6 fastq.gz). md5s are in the TSV.

## ENCODE Jurkat DNase-seq fallback
Only one Jurkat accessibility experiment exists on ENCODE (searchTerm=Jurkat; no ATAC):
- **ENCSR000EOS** — DNase-seq, Jurkat clone E6-1, UW (Stamatoyannopoulos), released 2011-05-04, 2 biological replicates, **single-end 36 bp** (FASTQ ENCFF001DPF, ENCFF001DPG).
- GRCh38 processed files (latest 2025-08 reprocessing): peaks bed3+ ENCFF768GGN / ENCFF262WSM (rep1), ENCFF932KNJ / ENCFF504GVU (rep2); hotspots ENCFF356OEJ (rep1), ENCFF773COO (rep2); signal bigWig ENCFF799GDL (rep1), ENCFF150LUL (rep2).
- 2020 pipeline: narrowPeak ENCFF613ANR / ENCFF340IYE (rep1), ENCFF957LSS / ENCFF035NVP (rep2); filtered BAM ENCFF361WBN (rep1), ENCFF455EQZ (rep2); bigWig ENCFF794WMH (rep1), ENCFF999ECC (rep2).
- 2017/18: BAM ENCFF415QST/ENCFF825XSS; narrowPeak ENCFF785KPP/ENCFF552VOA; hotspots ENCFF591NSJ/ENCFF361FEQ.

## Files
- `recon/05_jurkat_atac_selection.md` (this file)
- `recon/jurkat_atac_selected_runs.tsv`
- `recon/jurkat_atac_search/` — full SRA runinfo (1139 runs), GSM/GSE text dumps, parsed `gsm_info.json`, `study_gse.json`
