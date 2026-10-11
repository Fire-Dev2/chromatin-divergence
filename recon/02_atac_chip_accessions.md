# 02 — ATAC-seq and TF ChIP-seq accession inventory (ENCODE / 4DN)

Generated 2026-10-11 from live ENCODE REST API (`www.encodeproject.org`) and 4DN data portal API (`data.4dnucleome.org`) metadata; no data files downloaded. Scripts: `fetch_encode_atac.py`, `fetch_encode_files.py`, `fetch_encode_tf_chip.py`, `build_inventory.py` (this directory). Raw JSON: `encode_atac_candidates.json`, `atac_probe.json`, `fourdn_atac_sets.json`, `fourdn_4DNES48BA2C8.json`, `fourdn_file_*.json`.

All ENCODE experiments below are ENCODE4 Snyder-lab bulk ATAC-seq, processed by the ENCODE uniform ATAC-seq pipeline (`ENCPL344QWT` + `ENCPL787FUN`, lab `encode-processing-pipeline`) on GRCh38; the *released* default analysis is pipeline **v1.9.1** for all four target/decoy core experiments and most extras; exceptions: NCI-H929 ENCSR382LBS = **v1.8.0**, WTC11 ENCSR541KFY = **v1.10.0** (older v1.8.0 analyses on the others are `revoked`/`deleted`). File selection rule used throughout: files from the experiment's `default_analysis` only; bigWigs = MACS2 `signal p-value` and `fold change over control` (pooled and per replicate); peaks = `IDR thresholded peaks` from the pooled-pseudoreplicate IDR step (`pooled-pr1_vs_pooled-pr2.idr0.05.bfilt.narrowPeak`; this is the ENCODE ATAC 'IDR' set) plus the `pseudoreplicated peaks` pooled overlap set (`preferred_default=True` on the portal) and, where present, the `conservative IDR thresholded peaks` (true-replicate IDR).

## Summary table

| Role | Cell line | Source | Experiment | Lab | Bio reps | Assembly | Pipeline | Pooled signal p-value bigWig | Pooled fold-change bigWig | IDR peaks (pooled-pr) narrowPeak | Pseudoreplicated peaks narrowPeak | Notes |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| target | K562 | ENCODE | ENCSR868FGK | Michael Snyder, Stanford | 3 | GRCh38 | ENCODE ATAC v1.9.1 (ENCAN824UKX) | ENCFF357GNC (2.12 GB) | ENCFF102ARJ (4.51 GB) | ENCFF948AFM (4.2 MB) | ENCFF333TAT (7.1 MB) |  |
| target | hTERT-RPE1 | 4DN | 4DNES48BA2C8 (exps 4DNEX7A3MJJV, 4DNEX5R8DYKC) | Bing Ren, UCSD | 2 | GRCh38 | 4DN run of ENCODE ATAC pipeline (`encode-atacseq-aln` 1.1.1 + `encode-atacseq-postaln` 1.1.1; 4DNWF2BBZ1EE/4DNWF2VTG8II) | — (not produced) | 4DNFIS4IDZNE (4.41 GB) | 4DNFI4H1MT8U 'conservative peaks' bigBed (15.2 MB) | 4DNFIMNHSGU1 'peaks' bigBed (15.2 MB) | Not ENCODE-uniform; pipeline version differs (4DN CWL port of ENCODE atac-seq-pipeline v1.1.x vs ENCODE v1.9.1); peaks are bigBed not bed; no p-value track; **contact Ren lab before publication** (see policy below) |
| target | Jurkat (clone E6-1) | — | **none** | — | — | — | — | — | — | — | — | No ATAC-seq on ENCODE (any assembly/status) or 4DN; see blockers |
| target | HepG2 | ENCODE | ENCSR291GJU | Michael Snyder, Stanford | 3 | GRCh38 | ENCODE ATAC v1.9.1 (ENCAN877DUG) | ENCFF262URW (1.86 GB) | ENCFF645OHB (4.49 GB) | ENCFF935GLR (4.8 MB) | ENCFF439EIO (7.2 MB) |  |
| decoy | GM12878 | ENCODE | ENCSR637XSC | Michael Snyder, Stanford | 3 | GRCh38 | ENCODE ATAC v1.9.1 (ENCAN038LHH) | ENCFF667MDI (1.43 GB) | ENCFF603BJO (4.21 GB) | ENCFF429WTW (4.7 MB) | ENCFF748UZH (7.2 MB) |  |
| decoy | HCT116 | ENCODE | ENCSR872WGW | Michael Snyder, Stanford | 2 | GRCh38 | ENCODE ATAC v1.9.1 (ENCAN867NRM) | ENCFF624HRW (728.0 MB) | ENCFF259PSA (1.06 GB) | ENCFF054UBO (3.8 MB) | ENCFF296ZZB (5.3 MB) | only untreated/unmodified HCT116 ATAC (16 others are degron-tagged ± auxin/5-Ph-IAA) |
| decoy | MCF-7 | ENCODE | ENCSR422SUG | Michael Snyder, Stanford | 2 | GRCh38 | ENCODE ATAC v1.9.1 (ENCAN088EGP) | ENCFF976UNK (725.4 MB) | ENCFF782BVX (1.05 GB) | ENCFF882OVP (3.8 MB) | ENCFF821OEF (5.9 MB) |  |
| decoy | A549 | ENCODE | ENCSR032RGS | Michael Snyder, Stanford | 3 | GRCh38 | ENCODE ATAC v1.9.1 (ENCAN065GGH) | ENCFF399KCR (1.97 GB) | ENCFF872SDF (5.00 GB) | ENCFF964AUZ (4.8 MB) | ENCFF899OMR (7.2 MB) |  |
| decoy | IMR-90 | ENCODE | ENCSR200OML | Michael Snyder, Stanford | 2 | GRCh38 | ENCODE ATAC v1.9.1 (ENCAN754KOE) | ENCFF770EAV (763.9 MB) | ENCFF282RNO (1.12 GB) | ENCFF114GDS (3.8 MB) | ENCFF243NTP (6.5 MB) |  |
| decoy | PC-3 | ENCODE | ENCSR499ASS | Michael Snyder, Stanford | 2 | GRCh38 | ENCODE ATAC v1.9.1 (ENCAN146VQB) | ENCFF145UAD (672.7 MB) | ENCFF124ZQJ (1.12 GB) | ENCFF369VMF (4.2 MB) | ENCFF811MOZ (6.7 MB) |  |
| decoy | Panc1 | ENCODE | ENCSR591PIX | Michael Snyder, Stanford | 2 | GRCh38 | ENCODE ATAC v1.9.1 (ENCAN507JXY) | ENCFF794CNJ (781.7 MB) | ENCFF757YRZ (1.28 GB) | ENCFF481ELM (3.1 MB) | ENCFF182SSP (5.2 MB) |  |
| decoy | SK-N-SH | ENCODE | ENCSR587TRP | Michael Snyder, Stanford | 3 | GRCh38 | ENCODE ATAC v1.9.1 (ENCAN752KCW) | ENCFF370SJJ (1.69 GB) | ENCFF855TQX (4.00 GB) | ENCFF855CQR (4.9 MB) | ENCFF716JUM (7.1 MB) | **treated: all-trans-retinoic acid 6 µM 48 h** (only SK-N-SH ATAC on ENCODE); use with caution or drop |
| decoy | GM23338 | ENCODE | ENCSR485TLP | Michael Snyder, Stanford | 3 | GRCh38 | ENCODE ATAC v1.9.1 (ENCAN295SCS) | ENCFF234AYB (1.53 GB) | ENCFF378ZGI (4.49 GB) | ENCFF738KYF (4.8 MB) | ENCFF567ZCX (7.4 MB) | iPSC line (male, from GM23248 fibroblasts) |
| decoy | NCI-H929 | ENCODE | ENCSR382LBS | Michael Snyder, Stanford | 3 | GRCh38 | ENCODE ATAC v1.8.0 (ENCAN245MNX) | ENCFF571JLT (1.64 GB) | ENCFF573SLV (1.91 GB) | ENCFF369SMB (3.9 MB) | ENCFF060HHC (6.5 MB) |  |
| decoy | DND-41 | ENCODE | ENCSR660WSB | Michael Snyder, Stanford | 3 | GRCh38 | ENCODE ATAC v1.9.1 (ENCAN022AWY) | ENCFF631KUG (588.2 MB) | ENCFF679DXT (980.4 MB) | ENCFF225RUL (2.5 MB) | ENCFF095ATB (3.9 MB) | T-ALL line — same lineage as Jurkat; candidate Jurkat stand-in/control |
| decoy | RWPE2 | ENCODE | ENCSR080SNF | Michael Snyder, Stanford | 3 | GRCh38 | ENCODE ATAC v1.9.1 (ENCAN367YCE) | ENCFF881UWW (997.3 MB) | ENCFF047EBV (1.52 GB) | ENCFF646WSM (3.8 MB) | ENCFF729MMJ (6.3 MB) |  |
| decoy | WTC11 | ENCODE | ENCSR541KFY | Michael Snyder, Stanford | 2 | GRCh38 | ENCODE ATAC v1.10.0 (ENCAN873DEH) | ENCFF123YPY (1.12 GB) | ENCFF056MYW (1.46 GB) | ENCFF492GYV (3.8 MB) | ENCFF321VDH (6.6 MB) | iPSC; a second WTC11 ATAC (ENCSR089KIJ, Yin Shen) exists |
| alt | K562 | ENCODE | ENCSR483RKN | Michael Snyder, Stanford | 2 | GRCh38 | ENCODE ATAC v1.9.1 | see detail | | | | second untreated Snyder K562 ATAC; 62 other K562 ATAC are Greenleaf-lab drug/perturbation series (e.g. ENCSR165JXS = JQ1 5 µM 4 h) |
| alt | HepG2 | ENCODE | ENCSR042AWH | Michael Snyder, Stanford | 2 | GRCh38 | ENCODE ATAC v1.9.1 | see detail | | | | |
| alt | GM12878 | ENCODE | ENCSR095QNB | Michael Snyder, Stanford | 2 | GRCh38 | ENCODE ATAC v1.9.1 | see detail | | | | |

Cell lines requested as decoys but **absent** from ENCODE ATAC-seq (any status, any assembly, checked by `biosample_ontology.term_name` and free-text `searchTerm`): HeLa-S3, HEK293/HEK293T. 4DN has HEK293 ATAC only as CRISPR-perturbation sets from Taekjip Ha (JHU) (untreated control = 4DNES694VKB2, 2 exps); not ENCODE-uniform.

## TF ChIP-seq availability (ENCODE, released, GRCh38, IDR thresholded narrowPeak present)

Query: `type=Experiment&assay_title=TF ChIP-seq&status=released&biosample_ontology.term_name=<cell>&assembly=GRCh38&files.output_type=IDR thresholded peaks&files.file_type=bed narrowPeak`, then each experiment re-checked to have ≥1 *released* GRCh38 `IDR thresholded peaks` bed narrowPeak. 'Sequence-specific TF' = ENCODE target `investigated_as` contains `transcription factor` (excludes cofactors, chromatin remodelers, RBPs, Pol II, cohesin, DNA repair/replication).

| Cell line | Experiments (all TF-ChIP targets) | Distinct targets (all) | Experiments with TF targets | **Distinct sequence-specific TF targets** | CSV |
|---|---|---|---|---|---|
| K562 | 769 | 525 | 531 | **352** | `recon/encode_tf_chip_K562.csv` |
| HepG2 | 823 | 720 | 632 | **552** | `recon/encode_tf_chip_HepG2.csv` |

TF targets shared by K562 and HepG2: **205**; K562-only 147; HepG2-only 347. CSV columns: target, investigated_as, is_sequence_specific_TF, experiment, lab, award_rfa, date_released, n_bio_reps, chosen IDR peak file accession/replicates/size/URL, number of GRCh38 IDR files, audit ERROR / NOT_COMPLIANT categories (e.g. 'extremely low read depth' — filter on these before use). Where a target has several experiments, every experiment is a row; dedupe on `target` for counts.

## Per-experiment file detail (ENCODE)

### Targets

#### K562 — ENCSR868FGK

- Portal: https://www.encodeproject.org/experiments/ENCSR868FGK/ · lab Michael Snyder, Stanford · ENCODE4 · released 2020-08-05 · biosample summary: Homo sapiens K562 · biosamples ENCBS030HJI, ENCBS049FLD, ENCBS254OQS · bio reps [1, 2, 3]
- Treatments: none · genetic modifications: none
- Default analysis ENCAN824UKX (/analyses/ENCAN824UKX/): GRCh38, ENCODE ATAC-seq pipeline v1.9.1, status released; other analyses: ENCAN623SSG v1.8.0 revoked
- Audits: WARNING: moderate library complexity

| File | Output type | Replicates | Size | Pipeline file name | URL |
|---|---|---|---|---|---|
| ENCFF357GNC | signal p-value | rep1+2+3 | 2.12 GB | `rep.pooled.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF357GNC/@@download/ENCFF357GNC.bigWig |
| ENCFF102ARJ | fold change over control | rep1+2+3 | 4.51 GB | `rep.pooled.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF102ARJ/@@download/ENCFF102ARJ.bigWig |
| ENCFF670QXU | signal p-value | rep1 | 1.58 GB | `ENCFF962SRP.trim.merged.srt.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF670QXU/@@download/ENCFF670QXU.bigWig |
| ENCFF905MUY | signal p-value | rep2 | 993.5 MB | `ENCFF098UCE.trim.merged.srt.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF905MUY/@@download/ENCFF905MUY.bigWig |
| ENCFF043TYW | signal p-value | rep3 | 1.38 GB | `ENCFF354EXH.trim.merged.srt.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF043TYW/@@download/ENCFF043TYW.bigWig |
| ENCFF093IIW | fold change over control | rep1 | 2.19 GB | `ENCFF962SRP.trim.merged.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF093IIW/@@download/ENCFF093IIW.bigWig |
| ENCFF019IPA | fold change over control | rep2 | 1.69 GB | `ENCFF098UCE.trim.merged.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF019IPA/@@download/ENCFF019IPA.bigWig |
| ENCFF656ZKM | fold change over control | rep3 | 1.77 GB | `ENCFF354EXH.trim.merged.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF656ZKM/@@download/ENCFF656ZKM.bigWig |
| ENCFF948AFM | IDR thresholded peaks | rep1+2+3 | 4.2 MB | `pooled-pr1_vs_pooled-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF948AFM/@@download/ENCFF948AFM.bed.gz |
| ENCFF855PCP | conservative IDR thresholded peaks | rep1+2+3 | 4.0 MB | `rep1_vs_rep2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF855PCP/@@download/ENCFF855PCP.bed.gz |
| ENCFF333TAT | pseudoreplicated peaks | rep1+2+3 | 7.1 MB | `pooled-pr1_vs_pooled-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF333TAT/@@download/ENCFF333TAT.bed.gz |
| ENCFF057UYP | replicated peaks | rep1+2+3 | 6.8 MB | `rep2_vs_rep3.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF057UYP/@@download/ENCFF057UYP.bed.gz |
| ENCFF590CPE | replicated peaks | rep1+2+3 | 6.9 MB | `rep1_vs_rep3.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF590CPE/@@download/ENCFF590CPE.bed.gz |
| ENCFF678GBM | replicated peaks | rep1+2+3 | 6.9 MB | `rep1_vs_rep2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF678GBM/@@download/ENCFF678GBM.bed.gz |
| ENCFF993BAP | IDR thresholded peaks | rep1 | 3.7 MB | `rep1-pr1_vs_rep1-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF993BAP/@@download/ENCFF993BAP.bed.gz |
| ENCFF771HDN | IDR thresholded peaks | rep2 | 3.5 MB | `rep2-pr1_vs_rep2-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF771HDN/@@download/ENCFF771HDN.bed.gz |
| ENCFF433EPT | IDR thresholded peaks | rep3 | 3.4 MB | `rep3-pr1_vs_rep3-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF433EPT/@@download/ENCFF433EPT.bed.gz |
| ENCFF349AKO | pseudoreplicated peaks | rep1 | 6.4 MB | `rep1-pr1_vs_rep1-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF349AKO/@@download/ENCFF349AKO.bed.gz |
| ENCFF581TPM | pseudoreplicated peaks | rep2 | 6.1 MB | `rep2-pr1_vs_rep2-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF581TPM/@@download/ENCFF581TPM.bed.gz |
| ENCFF973QQZ | pseudoreplicated peaks | rep3 | 6.1 MB | `rep3-pr1_vs_rep3-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF973QQZ/@@download/ENCFF973QQZ.bed.gz |

Additional true-replicate pairwise IDR files (bed narrowPeak, all labelled `IDR thresholded peaks`, reps 1+2+3): ENCFF223QDM (`rep1_vs_rep3.idr0.05.bfilt.narrowPeak`, 4.0 MB), ENCFF135AEX (`rep2_vs_rep3.idr0.05.bfilt.narrowPeak`, 4.0 MB)

#### HepG2 — ENCSR291GJU

- Portal: https://www.encodeproject.org/experiments/ENCSR291GJU/ · lab Michael Snyder, Stanford · ENCODE4 · released 2020-08-05 · biosample summary: Homo sapiens HepG2 · biosamples ENCBS353INV, ENCBS710QMG, ENCBS827CSM · bio reps [1, 2, 3]
- Treatments: none · genetic modifications: none
- Default analysis ENCAN877DUG (/analyses/ENCAN877DUG/): GRCh38, ENCODE ATAC-seq pipeline v1.9.1, status released; other analyses: ENCAN223PFH v1.8.0 revoked
- Audits: WARNING: mild to moderate bottlenecking, moderate library complexity

| File | Output type | Replicates | Size | Pipeline file name | URL |
|---|---|---|---|---|---|
| ENCFF262URW | signal p-value | rep1+2+3 | 1.86 GB | `rep.pooled.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF262URW/@@download/ENCFF262URW.bigWig |
| ENCFF645OHB | fold change over control | rep1+2+3 | 4.49 GB | `rep.pooled.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF645OHB/@@download/ENCFF645OHB.bigWig |
| ENCFF254ZJG | signal p-value | rep1 | 1.53 GB | `ENCFF529IMC.trim.merged.srt.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF254ZJG/@@download/ENCFF254ZJG.bigWig |
| ENCFF380DRH | signal p-value | rep2 | 1.53 GB | `ENCFF134LRP.trim.merged.srt.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF380DRH/@@download/ENCFF380DRH.bigWig |
| ENCFF692UVH | signal p-value | rep3 | 1.18 GB | `ENCFF978ARG.trim.merged.srt.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF692UVH/@@download/ENCFF692UVH.bigWig |
| ENCFF622FRD | fold change over control | rep1 | 2.01 GB | `ENCFF529IMC.trim.merged.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF622FRD/@@download/ENCFF622FRD.bigWig |
| ENCFF024GLW | fold change over control | rep2 | 1.94 GB | `ENCFF134LRP.trim.merged.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF024GLW/@@download/ENCFF024GLW.bigWig |
| ENCFF240VVR | fold change over control | rep3 | 1.85 GB | `ENCFF978ARG.trim.merged.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF240VVR/@@download/ENCFF240VVR.bigWig |
| ENCFF935GLR | IDR thresholded peaks | rep1+2+3 | 4.8 MB | `pooled-pr1_vs_pooled-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF935GLR/@@download/ENCFF935GLR.bed.gz |
| ENCFF536RJV | conservative IDR thresholded peaks | rep1+2+3 | 4.7 MB | `rep1_vs_rep3.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF536RJV/@@download/ENCFF536RJV.bed.gz |
| ENCFF439EIO | pseudoreplicated peaks | rep1+2+3 | 7.2 MB | `pooled-pr1_vs_pooled-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF439EIO/@@download/ENCFF439EIO.bed.gz |
| ENCFF963FUC | replicated peaks | rep1+2+3 | 7.0 MB | `rep1_vs_rep2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF963FUC/@@download/ENCFF963FUC.bed.gz |
| ENCFF069GAC | replicated peaks | rep1+2+3 | 7.0 MB | `rep1_vs_rep3.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF069GAC/@@download/ENCFF069GAC.bed.gz |
| ENCFF471EFV | replicated peaks | rep1+2+3 | 7.0 MB | `rep2_vs_rep3.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF471EFV/@@download/ENCFF471EFV.bed.gz |
| ENCFF919RKQ | IDR thresholded peaks | rep1 | 4.2 MB | `rep1-pr1_vs_rep1-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF919RKQ/@@download/ENCFF919RKQ.bed.gz |
| ENCFF576UEM | IDR thresholded peaks | rep2 | 3.8 MB | `rep2-pr1_vs_rep2-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF576UEM/@@download/ENCFF576UEM.bed.gz |
| ENCFF161ZZX | IDR thresholded peaks | rep3 | 4.2 MB | `rep3-pr1_vs_rep3-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF161ZZX/@@download/ENCFF161ZZX.bed.gz |
| ENCFF258UBH | pseudoreplicated peaks | rep1 | 6.6 MB | `rep1-pr1_vs_rep1-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF258UBH/@@download/ENCFF258UBH.bed.gz |
| ENCFF046XTK | pseudoreplicated peaks | rep2 | 6.3 MB | `rep2-pr1_vs_rep2-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF046XTK/@@download/ENCFF046XTK.bed.gz |
| ENCFF474POM | pseudoreplicated peaks | rep3 | 6.7 MB | `rep3-pr1_vs_rep3-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF474POM/@@download/ENCFF474POM.bed.gz |

Additional true-replicate pairwise IDR files (bed narrowPeak, all labelled `IDR thresholded peaks`, reps 1+2+3): ENCFF915FZC (`rep1_vs_rep2.idr0.05.bfilt.narrowPeak`, 4.5 MB), ENCFF609BSU (`rep2_vs_rep3.idr0.05.bfilt.narrowPeak`, 4.6 MB)

### Decoys

#### GM12878 — ENCSR637XSC

- Portal: https://www.encodeproject.org/experiments/ENCSR637XSC/ · lab Michael Snyder, Stanford · ENCODE4 · released 2020-08-05 · biosample summary: Homo sapiens GM12878 · biosamples ENCBS737CLB, ENCBS836FXC, ENCBS911HVB · bio reps [1, 2, 3]
- Treatments: none · genetic modifications: none
- Default analysis ENCAN038LHH (/analyses/ENCAN038LHH/): GRCh38, ENCODE ATAC-seq pipeline v1.9.1, status released; other analyses: ENCAN632TSJ v1.8.0 revoked
- Audits: WARNING: mild to moderate bottlenecking, moderate library complexity

| File | Output type | Replicates | Size | Pipeline file name | URL |
|---|---|---|---|---|---|
| ENCFF667MDI | signal p-value | rep1+2+3 | 1.43 GB | `rep.pooled.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF667MDI/@@download/ENCFF667MDI.bigWig |
| ENCFF603BJO | fold change over control | rep1+2+3 | 4.21 GB | `rep.pooled.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF603BJO/@@download/ENCFF603BJO.bigWig |
| ENCFF695YUX | signal p-value | rep1 | 994.7 MB | `ENCFF453JCV.trim.merged.srt.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF695YUX/@@download/ENCFF695YUX.bigWig |
| ENCFF425HXV | signal p-value | rep2 | 1.11 GB | `ENCFF948NNH.trim.merged.srt.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF425HXV/@@download/ENCFF425HXV.bigWig |
| ENCFF191AAC | signal p-value | rep3 | 1.31 GB | `ENCFF856KEW.trim.merged.srt.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF191AAC/@@download/ENCFF191AAC.bigWig |
| ENCFF165WGA | fold change over control | rep1 | 1.76 GB | `ENCFF453JCV.trim.merged.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF165WGA/@@download/ENCFF165WGA.bigWig |
| ENCFF987JEB | fold change over control | rep2 | 1.57 GB | `ENCFF948NNH.trim.merged.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF987JEB/@@download/ENCFF987JEB.bigWig |
| ENCFF267VKW | fold change over control | rep3 | 1.85 GB | `ENCFF856KEW.trim.merged.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF267VKW/@@download/ENCFF267VKW.bigWig |
| ENCFF429WTW | IDR thresholded peaks | rep1+2+3 | 4.7 MB | `pooled-pr1_vs_pooled-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF429WTW/@@download/ENCFF429WTW.bed.gz |
| ENCFF654BVB | conservative IDR thresholded peaks | rep1+2+3 | 4.6 MB | `rep1_vs_rep3.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF654BVB/@@download/ENCFF654BVB.bed.gz |
| ENCFF748UZH | pseudoreplicated peaks | rep1+2+3 | 7.2 MB | `pooled-pr1_vs_pooled-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF748UZH/@@download/ENCFF748UZH.bed.gz |
| ENCFF687HBT | replicated peaks | rep1+2+3 | 7.0 MB | `rep1_vs_rep2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF687HBT/@@download/ENCFF687HBT.bed.gz |
| ENCFF061QPG | replicated peaks | rep1+2+3 | 7.1 MB | `rep1_vs_rep3.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF061QPG/@@download/ENCFF061QPG.bed.gz |
| ENCFF972TMS | replicated peaks | rep1+2+3 | 7.0 MB | `rep2_vs_rep3.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF972TMS/@@download/ENCFF972TMS.bed.gz |
| ENCFF914YHG | IDR thresholded peaks | rep1 | 4.1 MB | `rep1-pr1_vs_rep1-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF914YHG/@@download/ENCFF914YHG.bed.gz |
| ENCFF452ZBX | IDR thresholded peaks | rep2 | 4.1 MB | `rep2-pr1_vs_rep2-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF452ZBX/@@download/ENCFF452ZBX.bed.gz |
| ENCFF838ISC | IDR thresholded peaks | rep3 | 4.2 MB | `rep3-pr1_vs_rep3-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF838ISC/@@download/ENCFF838ISC.bed.gz |
| ENCFF764XLY | pseudoreplicated peaks | rep1 | 6.5 MB | `rep1-pr1_vs_rep1-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF764XLY/@@download/ENCFF764XLY.bed.gz |
| ENCFF078VWH | pseudoreplicated peaks | rep2 | 6.4 MB | `rep2-pr1_vs_rep2-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF078VWH/@@download/ENCFF078VWH.bed.gz |
| ENCFF798HKX | pseudoreplicated peaks | rep3 | 6.5 MB | `rep3-pr1_vs_rep3-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF798HKX/@@download/ENCFF798HKX.bed.gz |

Additional true-replicate pairwise IDR files (bed narrowPeak, all labelled `IDR thresholded peaks`, reps 1+2+3): ENCFF871TST (`rep1_vs_rep2.idr0.05.bfilt.narrowPeak`, 4.5 MB), ENCFF620WFJ (`rep2_vs_rep3.idr0.05.bfilt.narrowPeak`, 4.6 MB)

#### HCT116 — ENCSR872WGW

- Portal: https://www.encodeproject.org/experiments/ENCSR872WGW/ · lab Michael Snyder, Stanford · ENCODE4 · released 2020-10-08 · biosample summary: Homo sapiens HCT116 · biosamples ENCBS162XQL, ENCBS666SZS · bio reps [1, 2]
- Treatments: none · genetic modifications: none
- Default analysis ENCAN867NRM (/analyses/ENCAN867NRM/): GRCh38, ENCODE ATAC-seq pipeline v1.9.1, status released; other analyses: ENCAN524UIW v1.8.0 revoked
- Audits: WARNING: mild to moderate bottlenecking, moderate library complexity

| File | Output type | Replicates | Size | Pipeline file name | URL |
|---|---|---|---|---|---|
| ENCFF624HRW | signal p-value | rep1+2 | 728.0 MB | `rep.pooled.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF624HRW/@@download/ENCFF624HRW.bigWig |
| ENCFF259PSA | fold change over control | rep1+2 | 1.06 GB | `rep.pooled.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF259PSA/@@download/ENCFF259PSA.bigWig |
| ENCFF807VJD | signal p-value | rep1 | 463.7 MB | `ENCFF121EPT.trim.srt.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF807VJD/@@download/ENCFF807VJD.bigWig |
| ENCFF962GFP | signal p-value | rep2 | 486.9 MB | `ENCFF624DNH.trim.srt.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF962GFP/@@download/ENCFF962GFP.bigWig |
| ENCFF480SPP | fold change over control | rep1 | 579.7 MB | `ENCFF121EPT.trim.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF480SPP/@@download/ENCFF480SPP.bigWig |
| ENCFF725NBM | fold change over control | rep2 | 644.6 MB | `ENCFF624DNH.trim.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF725NBM/@@download/ENCFF725NBM.bigWig |
| ENCFF054UBO | IDR thresholded peaks | rep1+2 | 3.8 MB | `pooled-pr1_vs_pooled-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF054UBO/@@download/ENCFF054UBO.bed.gz |
| ENCFF070LSQ | conservative IDR thresholded peaks | rep1+2 | 3.7 MB | `rep1_vs_rep2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF070LSQ/@@download/ENCFF070LSQ.bed.gz |
| ENCFF296ZZB | pseudoreplicated peaks | rep1+2 | 5.3 MB | `pooled-pr1_vs_pooled-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF296ZZB/@@download/ENCFF296ZZB.bed.gz |
| ENCFF107FCE | replicated peaks | rep1+2 | 5.2 MB | `rep1_vs_rep2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF107FCE/@@download/ENCFF107FCE.bed.gz |
| ENCFF647AGG | IDR thresholded peaks | rep1 | 2.9 MB | `rep1-pr1_vs_rep1-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF647AGG/@@download/ENCFF647AGG.bed.gz |
| ENCFF429LWF | IDR thresholded peaks | rep2 | 3.1 MB | `rep2-pr1_vs_rep2-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF429LWF/@@download/ENCFF429LWF.bed.gz |
| ENCFF515RJQ | pseudoreplicated peaks | rep1 | 4.3 MB | `rep1-pr1_vs_rep1-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF515RJQ/@@download/ENCFF515RJQ.bed.gz |
| ENCFF662SAQ | pseudoreplicated peaks | rep2 | 4.6 MB | `rep2-pr1_vs_rep2-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF662SAQ/@@download/ENCFF662SAQ.bed.gz |

#### MCF-7 — ENCSR422SUG

- Portal: https://www.encodeproject.org/experiments/ENCSR422SUG/ · lab Michael Snyder, Stanford · ENCODE4 · released 2020-08-05 · biosample summary: Homo sapiens MCF-7 · biosamples ENCBS105SMU, ENCBS352FBE · bio reps [1, 2]
- Treatments: none · genetic modifications: none
- Default analysis ENCAN088EGP (/analyses/ENCAN088EGP/): GRCh38, ENCODE ATAC-seq pipeline v1.9.1, status released; other analyses: ENCAN662XCF v1.8.0 revoked
- Audits: WARNING: mild to moderate bottlenecking, moderate library complexity

| File | Output type | Replicates | Size | Pipeline file name | URL |
|---|---|---|---|---|---|
| ENCFF976UNK | signal p-value | rep1+2 | 725.4 MB | `rep.pooled.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF976UNK/@@download/ENCFF976UNK.bigWig |
| ENCFF782BVX | fold change over control | rep1+2 | 1.05 GB | `rep.pooled.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF782BVX/@@download/ENCFF782BVX.bigWig |
| ENCFF864CTV | signal p-value | rep1 | 430.2 MB | `ENCFF343ABN.trim.srt.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF864CTV/@@download/ENCFF864CTV.bigWig |
| ENCFF463PDI | signal p-value | rep2 | 408.7 MB | `ENCFF053NCI.trim.srt.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF463PDI/@@download/ENCFF463PDI.bigWig |
| ENCFF194DID | fold change over control | rep1 | 630.9 MB | `ENCFF343ABN.trim.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF194DID/@@download/ENCFF194DID.bigWig |
| ENCFF854ZKH | fold change over control | rep2 | 584.8 MB | `ENCFF053NCI.trim.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF854ZKH/@@download/ENCFF854ZKH.bigWig |
| ENCFF882OVP | IDR thresholded peaks | rep1+2 | 3.8 MB | `pooled-pr1_vs_pooled-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF882OVP/@@download/ENCFF882OVP.bed.gz |
| ENCFF905GGG | conservative IDR thresholded peaks | rep1+2 | 3.7 MB | `rep1_vs_rep2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF905GGG/@@download/ENCFF905GGG.bed.gz |
| ENCFF821OEF | pseudoreplicated peaks | rep1+2 | 5.9 MB | `pooled-pr1_vs_pooled-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF821OEF/@@download/ENCFF821OEF.bed.gz |
| ENCFF291WTA | replicated peaks | rep1+2 | 5.8 MB | `rep1_vs_rep2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF291WTA/@@download/ENCFF291WTA.bed.gz |
| ENCFF882FLC | IDR thresholded peaks | rep1 | 2.9 MB | `rep1-pr1_vs_rep1-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF882FLC/@@download/ENCFF882FLC.bed.gz |
| ENCFF204NII | IDR thresholded peaks | rep2 | 3.1 MB | `rep2-pr1_vs_rep2-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF204NII/@@download/ENCFF204NII.bed.gz |
| ENCFF536UWD | pseudoreplicated peaks | rep1 | 4.9 MB | `rep1-pr1_vs_rep1-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF536UWD/@@download/ENCFF536UWD.bed.gz |
| ENCFF679PKT | pseudoreplicated peaks | rep2 | 5.2 MB | `rep2-pr1_vs_rep2-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF679PKT/@@download/ENCFF679PKT.bed.gz |

#### A549 — ENCSR032RGS

- Portal: https://www.encodeproject.org/experiments/ENCSR032RGS/ · lab Michael Snyder, Stanford · ENCODE4 · released 2020-08-05 · biosample summary: Homo sapiens A549 · biosamples ENCBS003HWJ, ENCBS024PJA, ENCBS049MEE · bio reps [1, 2, 3]
- Treatments: none · genetic modifications: none
- Default analysis ENCAN065GGH (/analyses/ENCAN065GGH/): GRCh38, ENCODE ATAC-seq pipeline v1.9.1, status released; other analyses: ENCAN170AAI v1.8.0 revoked
- Audits: WARNING: mild to moderate bottlenecking, moderate library complexity

| File | Output type | Replicates | Size | Pipeline file name | URL |
|---|---|---|---|---|---|
| ENCFF399KCR | signal p-value | rep1+2+3 | 1.97 GB | `rep.pooled.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF399KCR/@@download/ENCFF399KCR.bigWig |
| ENCFF872SDF | fold change over control | rep1+2+3 | 5.00 GB | `rep.pooled.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF872SDF/@@download/ENCFF872SDF.bigWig |
| ENCFF561BCI | signal p-value | rep1 | 1.25 GB | `ENCFF047YJC.trim.merged.srt.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF561BCI/@@download/ENCFF561BCI.bigWig |
| ENCFF950CJJ | signal p-value | rep2 | 1.54 GB | `ENCFF837LSF.trim.merged.srt.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF950CJJ/@@download/ENCFF950CJJ.bigWig |
| ENCFF022BCI | signal p-value | rep3 | 1.39 GB | `ENCFF965LPB.trim.merged.srt.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF022BCI/@@download/ENCFF022BCI.bigWig |
| ENCFF536RWL | fold change over control | rep1 | 2.02 GB | `ENCFF047YJC.trim.merged.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF536RWL/@@download/ENCFF536RWL.bigWig |
| ENCFF550BVZ | fold change over control | rep2 | 2.18 GB | `ENCFF837LSF.trim.merged.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF550BVZ/@@download/ENCFF550BVZ.bigWig |
| ENCFF059JAV | fold change over control | rep3 | 2.45 GB | `ENCFF965LPB.trim.merged.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF059JAV/@@download/ENCFF059JAV.bigWig |
| ENCFF964AUZ | IDR thresholded peaks | rep1+2+3 | 4.8 MB | `pooled-pr1_vs_pooled-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF964AUZ/@@download/ENCFF964AUZ.bed.gz |
| ENCFF876UEM | conservative IDR thresholded peaks | rep1+2+3 | 4.8 MB | `rep2_vs_rep3.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF876UEM/@@download/ENCFF876UEM.bed.gz |
| ENCFF899OMR | pseudoreplicated peaks | rep1+2+3 | 7.2 MB | `pooled-pr1_vs_pooled-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF899OMR/@@download/ENCFF899OMR.bed.gz |
| ENCFF819ZZO | replicated peaks | rep1+2+3 | 6.8 MB | `rep1_vs_rep3.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF819ZZO/@@download/ENCFF819ZZO.bed.gz |
| ENCFF648AEN | replicated peaks | rep1+2+3 | 7.0 MB | `rep2_vs_rep3.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF648AEN/@@download/ENCFF648AEN.bed.gz |
| ENCFF325MIK | replicated peaks | rep1+2+3 | 6.9 MB | `rep1_vs_rep2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF325MIK/@@download/ENCFF325MIK.bed.gz |
| ENCFF681TUR | IDR thresholded peaks | rep1 | 3.8 MB | `rep1-pr1_vs_rep1-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF681TUR/@@download/ENCFF681TUR.bed.gz |
| ENCFF948AQW | IDR thresholded peaks | rep2 | 4.2 MB | `rep2-pr1_vs_rep2-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF948AQW/@@download/ENCFF948AQW.bed.gz |
| ENCFF237URW | IDR thresholded peaks | rep3 | 4.4 MB | `rep3-pr1_vs_rep3-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF237URW/@@download/ENCFF237URW.bed.gz |
| ENCFF434KBP | pseudoreplicated peaks | rep1 | 6.1 MB | `rep1-pr1_vs_rep1-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF434KBP/@@download/ENCFF434KBP.bed.gz |
| ENCFF998PAN | pseudoreplicated peaks | rep2 | 6.5 MB | `rep2-pr1_vs_rep2-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF998PAN/@@download/ENCFF998PAN.bed.gz |
| ENCFF316ADH | pseudoreplicated peaks | rep3 | 6.6 MB | `rep3-pr1_vs_rep3-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF316ADH/@@download/ENCFF316ADH.bed.gz |

Additional true-replicate pairwise IDR files (bed narrowPeak, all labelled `IDR thresholded peaks`, reps 1+2+3): ENCFF429FOV (`rep1_vs_rep2.idr0.05.bfilt.narrowPeak`, 4.6 MB), ENCFF674RFR (`rep1_vs_rep3.idr0.05.bfilt.narrowPeak`, 4.6 MB)

#### IMR-90 — ENCSR200OML

- Portal: https://www.encodeproject.org/experiments/ENCSR200OML/ · lab Michael Snyder, Stanford · ENCODE4 · released 2020-08-05 · biosample summary: Homo sapiens IMR-90 nuclear fraction and unspecified fraction · biosamples ENCBS666PTP, ENCBS872CEV · bio reps [1, 2]
- Treatments: none · genetic modifications: none
- Default analysis ENCAN754KOE (/analyses/ENCAN754KOE/): GRCh38, ENCODE ATAC-seq pipeline v1.9.1, status released; other analyses: ENCAN564SBA v1.8.0 revoked
- Audits: WARNING: low read depth, mild to moderate bottlenecking, moderate library complexity

| File | Output type | Replicates | Size | Pipeline file name | URL |
|---|---|---|---|---|---|
| ENCFF770EAV | signal p-value | rep1+2 | 763.9 MB | `rep.pooled.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF770EAV/@@download/ENCFF770EAV.bigWig |
| ENCFF282RNO | fold change over control | rep1+2 | 1.12 GB | `rep.pooled.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF282RNO/@@download/ENCFF282RNO.bigWig |
| ENCFF861RID | signal p-value | rep1 | 580.5 MB | `ENCFF846GEV.trim.srt.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF861RID/@@download/ENCFF861RID.bigWig |
| ENCFF619RWM | signal p-value | rep2 | 531.1 MB | `ENCFF841MYT.trim.srt.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF619RWM/@@download/ENCFF619RWM.bigWig |
| ENCFF310KYH | fold change over control | rep1 | 708.2 MB | `ENCFF846GEV.trim.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF310KYH/@@download/ENCFF310KYH.bigWig |
| ENCFF921WHS | fold change over control | rep2 | 542.3 MB | `ENCFF841MYT.trim.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF921WHS/@@download/ENCFF921WHS.bigWig |
| ENCFF114GDS | IDR thresholded peaks | rep1+2 | 3.8 MB | `pooled-pr1_vs_pooled-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF114GDS/@@download/ENCFF114GDS.bed.gz |
| ENCFF982UNH | conservative IDR thresholded peaks | rep1+2 | 3.7 MB | `rep1_vs_rep2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF982UNH/@@download/ENCFF982UNH.bed.gz |
| ENCFF243NTP | pseudoreplicated peaks | rep1+2 | 6.5 MB | `pooled-pr1_vs_pooled-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF243NTP/@@download/ENCFF243NTP.bed.gz |
| ENCFF830JMV | replicated peaks | rep1+2 | 6.5 MB | `rep1_vs_rep2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF830JMV/@@download/ENCFF830JMV.bed.gz |
| ENCFF235JZB | IDR thresholded peaks | rep1 | 3.3 MB | `rep1-pr1_vs_rep1-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF235JZB/@@download/ENCFF235JZB.bed.gz |
| ENCFF383XOM | IDR thresholded peaks | rep2 | 2.7 MB | `rep2-pr1_vs_rep2-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF383XOM/@@download/ENCFF383XOM.bed.gz |
| ENCFF961VMB | pseudoreplicated peaks | rep1 | 6.2 MB | `rep1-pr1_vs_rep1-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF961VMB/@@download/ENCFF961VMB.bed.gz |
| ENCFF503RFW | pseudoreplicated peaks | rep2 | 5.7 MB | `rep2-pr1_vs_rep2-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF503RFW/@@download/ENCFF503RFW.bed.gz |

#### PC-3 — ENCSR499ASS

- Portal: https://www.encodeproject.org/experiments/ENCSR499ASS/ · lab Michael Snyder, Stanford · ENCODE4 · released 2020-10-08 · biosample summary: Homo sapiens PC-3 · biosamples ENCBS132BFL, ENCBS235HQZ · bio reps [1, 2]
- Treatments: none · genetic modifications: none
- Default analysis ENCAN146VQB (/analyses/ENCAN146VQB/): GRCh38, ENCODE ATAC-seq pipeline v1.9.1, status released; other analyses: ENCAN999QIV v1.8.0 revoked
- Audits: WARNING: mild to moderate bottlenecking, moderate library complexity

| File | Output type | Replicates | Size | Pipeline file name | URL |
|---|---|---|---|---|---|
| ENCFF145UAD | signal p-value | rep1+2 | 672.7 MB | `rep.pooled.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF145UAD/@@download/ENCFF145UAD.bigWig |
| ENCFF124ZQJ | fold change over control | rep1+2 | 1.12 GB | `rep.pooled.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF124ZQJ/@@download/ENCFF124ZQJ.bigWig |
| ENCFF285RAK | signal p-value | rep1 | 477.4 MB | `ENCFF561JSL.trim.srt.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF285RAK/@@download/ENCFF285RAK.bigWig |
| ENCFF746GZT | signal p-value | rep2 | 541.4 MB | `ENCFF284TLV.trim.srt.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF746GZT/@@download/ENCFF746GZT.bigWig |
| ENCFF064FOF | fold change over control | rep1 | 583.7 MB | `ENCFF561JSL.trim.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF064FOF/@@download/ENCFF064FOF.bigWig |
| ENCFF895PMJ | fold change over control | rep2 | 695.8 MB | `ENCFF284TLV.trim.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF895PMJ/@@download/ENCFF895PMJ.bigWig |
| ENCFF369VMF | IDR thresholded peaks | rep1+2 | 4.2 MB | `pooled-pr1_vs_pooled-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF369VMF/@@download/ENCFF369VMF.bed.gz |
| ENCFF466CVI | conservative IDR thresholded peaks | rep1+2 | 4.1 MB | `rep1_vs_rep2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF466CVI/@@download/ENCFF466CVI.bed.gz |
| ENCFF811MOZ | pseudoreplicated peaks | rep1+2 | 6.7 MB | `pooled-pr1_vs_pooled-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF811MOZ/@@download/ENCFF811MOZ.bed.gz |
| ENCFF009DVW | replicated peaks | rep1+2 | 6.6 MB | `rep1_vs_rep2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF009DVW/@@download/ENCFF009DVW.bed.gz |
| ENCFF061YKV | IDR thresholded peaks | rep1 | 3.5 MB | `rep1-pr1_vs_rep1-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF061YKV/@@download/ENCFF061YKV.bed.gz |
| ENCFF742PNL | IDR thresholded peaks | rep2 | 3.6 MB | `rep2-pr1_vs_rep2-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF742PNL/@@download/ENCFF742PNL.bed.gz |
| ENCFF906HPY | pseudoreplicated peaks | rep1 | 6.0 MB | `rep1-pr1_vs_rep1-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF906HPY/@@download/ENCFF906HPY.bed.gz |
| ENCFF984BYE | pseudoreplicated peaks | rep2 | 6.3 MB | `rep2-pr1_vs_rep2-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF984BYE/@@download/ENCFF984BYE.bed.gz |

#### Panc1 — ENCSR591PIX

- Portal: https://www.encodeproject.org/experiments/ENCSR591PIX/ · lab Michael Snyder, Stanford · ENCODE4 · released 2020-10-08 · biosample summary: Homo sapiens Panc1 · biosamples ENCBS298BRW, ENCBS630ELJ · bio reps [1, 2]
- Treatments: none · genetic modifications: none
- Default analysis ENCAN507JXY (/analyses/ENCAN507JXY/): GRCh38, ENCODE ATAC-seq pipeline v1.9.1, status released; other analyses: ENCAN556MAX v1.8.0 revoked
- Audits: WARNING: mild to moderate bottlenecking, moderate library complexity

| File | Output type | Replicates | Size | Pipeline file name | URL |
|---|---|---|---|---|---|
| ENCFF794CNJ | signal p-value | rep1+2 | 781.7 MB | `rep.pooled.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF794CNJ/@@download/ENCFF794CNJ.bigWig |
| ENCFF757YRZ | fold change over control | rep1+2 | 1.28 GB | `rep.pooled.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF757YRZ/@@download/ENCFF757YRZ.bigWig |
| ENCFF346FZQ | signal p-value | rep1 | 595.7 MB | `ENCFF020KID.trim.srt.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF346FZQ/@@download/ENCFF346FZQ.bigWig |
| ENCFF055ZEE | signal p-value | rep2 | 473.2 MB | `ENCFF695SZK.trim.srt.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF055ZEE/@@download/ENCFF055ZEE.bigWig |
| ENCFF185TQC | fold change over control | rep1 | 786.0 MB | `ENCFF020KID.trim.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF185TQC/@@download/ENCFF185TQC.bigWig |
| ENCFF174PXJ | fold change over control | rep2 | 676.7 MB | `ENCFF695SZK.trim.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF174PXJ/@@download/ENCFF174PXJ.bigWig |
| ENCFF481ELM | IDR thresholded peaks | rep1+2 | 3.1 MB | `pooled-pr1_vs_pooled-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF481ELM/@@download/ENCFF481ELM.bed.gz |
| ENCFF544MQW | conservative IDR thresholded peaks | rep1+2 | 3.1 MB | `rep1_vs_rep2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF544MQW/@@download/ENCFF544MQW.bed.gz |
| ENCFF182SSP | pseudoreplicated peaks | rep1+2 | 5.2 MB | `pooled-pr1_vs_pooled-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF182SSP/@@download/ENCFF182SSP.bed.gz |
| ENCFF112PVK | replicated peaks | rep1+2 | 5.1 MB | `rep1_vs_rep2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF112PVK/@@download/ENCFF112PVK.bed.gz |
| ENCFF953NZY | IDR thresholded peaks | rep1 | 2.2 MB | `rep1-pr1_vs_rep1-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF953NZY/@@download/ENCFF953NZY.bed.gz |
| ENCFF510OSH | IDR thresholded peaks | rep2 | 2.4 MB | `rep2-pr1_vs_rep2-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF510OSH/@@download/ENCFF510OSH.bed.gz |
| ENCFF835BPL | pseudoreplicated peaks | rep1 | 4.0 MB | `rep1-pr1_vs_rep1-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF835BPL/@@download/ENCFF835BPL.bed.gz |
| ENCFF255BWR | pseudoreplicated peaks | rep2 | 4.1 MB | `rep2-pr1_vs_rep2-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF255BWR/@@download/ENCFF255BWR.bed.gz |

#### SK-N-SH — ENCSR587TRP

- Portal: https://www.encodeproject.org/experiments/ENCSR587TRP/ · lab Michael Snyder, Stanford · ENCODE4 · released 2020-08-05 · biosample summary: Homo sapiens SK-N-SH treated with 6 μM all-trans-retinoic acid for 48 hours · biosamples ENCBS037ALL, ENCBS460YKE, ENCBS612PAQ · bio reps [1, 2, 3]
- Treatments: ['all-trans-retinoic acid 6μM 48hour'] · genetic modifications: none
- Default analysis ENCAN752KCW (/analyses/ENCAN752KCW/): GRCh38, ENCODE ATAC-seq pipeline v1.9.1, status released; other analyses: ENCAN194AXP v1.8.0 revoked
- Audits: WARNING: mild to moderate bottlenecking, moderate library complexity

| File | Output type | Replicates | Size | Pipeline file name | URL |
|---|---|---|---|---|---|
| ENCFF370SJJ | signal p-value | rep1+2+3 | 1.69 GB | `rep.pooled.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF370SJJ/@@download/ENCFF370SJJ.bigWig |
| ENCFF855TQX | fold change over control | rep1+2+3 | 4.00 GB | `rep.pooled.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF855TQX/@@download/ENCFF855TQX.bigWig |
| ENCFF958VOQ | signal p-value | rep1 | 992.0 MB | `ENCFF992SGW.trim.merged.srt.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF958VOQ/@@download/ENCFF958VOQ.bigWig |
| ENCFF038IHB | signal p-value | rep2 | 1.02 GB | `ENCFF111PVC.trim.merged.srt.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF038IHB/@@download/ENCFF038IHB.bigWig |
| ENCFF895UUA | signal p-value | rep3 | 1.27 GB | `ENCFF339SJV.trim.merged.srt.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF895UUA/@@download/ENCFF895UUA.bigWig |
| ENCFF655LOB | fold change over control | rep1 | 1.57 GB | `ENCFF992SGW.trim.merged.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF655LOB/@@download/ENCFF655LOB.bigWig |
| ENCFF923TGK | fold change over control | rep2 | 1.37 GB | `ENCFF111PVC.trim.merged.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF923TGK/@@download/ENCFF923TGK.bigWig |
| ENCFF367OPD | fold change over control | rep3 | 1.58 GB | `ENCFF339SJV.trim.merged.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF367OPD/@@download/ENCFF367OPD.bigWig |
| ENCFF855CQR | IDR thresholded peaks | rep1+2+3 | 4.9 MB | `pooled-pr1_vs_pooled-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF855CQR/@@download/ENCFF855CQR.bed.gz |
| ENCFF685CFY | conservative IDR thresholded peaks | rep1+2+3 | 4.7 MB | `rep1_vs_rep2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF685CFY/@@download/ENCFF685CFY.bed.gz |
| ENCFF716JUM | pseudoreplicated peaks | rep1+2+3 | 7.1 MB | `pooled-pr1_vs_pooled-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF716JUM/@@download/ENCFF716JUM.bed.gz |
| ENCFF687SBG | replicated peaks | rep1+2+3 | 6.9 MB | `rep2_vs_rep3.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF687SBG/@@download/ENCFF687SBG.bed.gz |
| ENCFF795ASL | replicated peaks | rep1+2+3 | 6.9 MB | `rep1_vs_rep2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF795ASL/@@download/ENCFF795ASL.bed.gz |
| ENCFF823UVH | replicated peaks | rep1+2+3 | 6.8 MB | `rep1_vs_rep3.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF823UVH/@@download/ENCFF823UVH.bed.gz |
| ENCFF371ZSI | IDR thresholded peaks | rep1 | 3.7 MB | `rep1-pr1_vs_rep1-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF371ZSI/@@download/ENCFF371ZSI.bed.gz |
| ENCFF464MPR | IDR thresholded peaks | rep2 | 4.4 MB | `rep2-pr1_vs_rep2-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF464MPR/@@download/ENCFF464MPR.bed.gz |
| ENCFF833RCJ | IDR thresholded peaks | rep3 | 4.3 MB | `rep3-pr1_vs_rep3-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF833RCJ/@@download/ENCFF833RCJ.bed.gz |
| ENCFF575GKL | pseudoreplicated peaks | rep1 | 6.1 MB | `rep1-pr1_vs_rep1-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF575GKL/@@download/ENCFF575GKL.bed.gz |
| ENCFF920LQX | pseudoreplicated peaks | rep2 | 6.6 MB | `rep2-pr1_vs_rep2-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF920LQX/@@download/ENCFF920LQX.bed.gz |
| ENCFF290VVZ | pseudoreplicated peaks | rep3 | 6.6 MB | `rep3-pr1_vs_rep3-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF290VVZ/@@download/ENCFF290VVZ.bed.gz |

Additional true-replicate pairwise IDR files (bed narrowPeak, all labelled `IDR thresholded peaks`, reps 1+2+3): ENCFF720DJU (`rep2_vs_rep3.idr0.05.bfilt.narrowPeak`, 4.6 MB), ENCFF461GIA (`rep1_vs_rep3.idr0.05.bfilt.narrowPeak`, 4.6 MB)

#### GM23338 — ENCSR485TLP

- Portal: https://www.encodeproject.org/experiments/ENCSR485TLP/ · lab Michael Snyder, Stanford · ENCODE4 · released 2020-08-05 · biosample summary: Homo sapiens GM23338 · biosamples ENCBS057LKL, ENCBS247UHY, ENCBS994YIB · bio reps [1, 2, 3]
- Treatments: none · genetic modifications: none
- Default analysis ENCAN295SCS (/analyses/ENCAN295SCS/): GRCh38, ENCODE ATAC-seq pipeline v1.9.1, status released; other analyses: ENCAN462LFE v1.8.0 revoked
- Audits: 

| File | Output type | Replicates | Size | Pipeline file name | URL |
|---|---|---|---|---|---|
| ENCFF234AYB | signal p-value | rep1+2+3 | 1.53 GB | `rep.pooled.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF234AYB/@@download/ENCFF234AYB.bigWig |
| ENCFF378ZGI | fold change over control | rep1+2+3 | 4.49 GB | `rep.pooled.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF378ZGI/@@download/ENCFF378ZGI.bigWig |
| ENCFF774QTX | signal p-value | rep1 | 1.46 GB | `ENCFF546FGZ.trim.merged.srt.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF774QTX/@@download/ENCFF774QTX.bigWig |
| ENCFF480VEI | signal p-value | rep2 | 1.39 GB | `ENCFF453RQX.trim.merged.srt.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF480VEI/@@download/ENCFF480VEI.bigWig |
| ENCFF581SQD | signal p-value | rep3 | 1.33 GB | `ENCFF482WXD.trim.merged.srt.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF581SQD/@@download/ENCFF581SQD.bigWig |
| ENCFF095DJH | fold change over control | rep1 | 1.90 GB | `ENCFF546FGZ.trim.merged.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF095DJH/@@download/ENCFF095DJH.bigWig |
| ENCFF276QZQ | fold change over control | rep2 | 2.01 GB | `ENCFF453RQX.trim.merged.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF276QZQ/@@download/ENCFF276QZQ.bigWig |
| ENCFF957TWY | fold change over control | rep3 | 1.90 GB | `ENCFF482WXD.trim.merged.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF957TWY/@@download/ENCFF957TWY.bigWig |
| ENCFF738KYF | IDR thresholded peaks | rep1+2+3 | 4.8 MB | `pooled-pr1_vs_pooled-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF738KYF/@@download/ENCFF738KYF.bed.gz |
| ENCFF870AYG | conservative IDR thresholded peaks | rep1+2+3 | 4.7 MB | `rep1_vs_rep2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF870AYG/@@download/ENCFF870AYG.bed.gz |
| ENCFF567ZCX | pseudoreplicated peaks | rep1+2+3 | 7.4 MB | `pooled-pr1_vs_pooled-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF567ZCX/@@download/ENCFF567ZCX.bed.gz |
| ENCFF180RLK | replicated peaks | rep1+2+3 | 7.3 MB | `rep1_vs_rep3.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF180RLK/@@download/ENCFF180RLK.bed.gz |
| ENCFF977EIN | replicated peaks | rep1+2+3 | 7.2 MB | `rep2_vs_rep3.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF977EIN/@@download/ENCFF977EIN.bed.gz |
| ENCFF213STV | replicated peaks | rep1+2+3 | 7.2 MB | `rep1_vs_rep2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF213STV/@@download/ENCFF213STV.bed.gz |
| ENCFF699ZTV | IDR thresholded peaks | rep1 | 4.1 MB | `rep1-pr1_vs_rep1-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF699ZTV/@@download/ENCFF699ZTV.bed.gz |
| ENCFF166IDS | IDR thresholded peaks | rep2 | 4.3 MB | `rep2-pr1_vs_rep2-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF166IDS/@@download/ENCFF166IDS.bed.gz |
| ENCFF239SMJ | IDR thresholded peaks | rep3 | 4.1 MB | `rep3-pr1_vs_rep3-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF239SMJ/@@download/ENCFF239SMJ.bed.gz |
| ENCFF441YUP | pseudoreplicated peaks | rep1 | 6.8 MB | `rep1-pr1_vs_rep1-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF441YUP/@@download/ENCFF441YUP.bed.gz |
| ENCFF567SUT | pseudoreplicated peaks | rep2 | 6.9 MB | `rep2-pr1_vs_rep2-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF567SUT/@@download/ENCFF567SUT.bed.gz |
| ENCFF018IPG | pseudoreplicated peaks | rep3 | 6.8 MB | `rep3-pr1_vs_rep3-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF018IPG/@@download/ENCFF018IPG.bed.gz |

Additional true-replicate pairwise IDR files (bed narrowPeak, all labelled `IDR thresholded peaks`, reps 1+2+3): ENCFF942QQL (`rep2_vs_rep3.idr0.05.bfilt.narrowPeak`, 4.6 MB), ENCFF426DRI (`rep1_vs_rep3.idr0.05.bfilt.narrowPeak`, 4.6 MB)

#### NCI-H929 — ENCSR382LBS

- Portal: https://www.encodeproject.org/experiments/ENCSR382LBS/ · lab Michael Snyder, Stanford · ENCODE4 · released 2021-04-06 · biosample summary: Homo sapiens NCI-H929 · biosamples ENCBS231ULS, ENCBS379IVU, ENCBS952JSJ · bio reps [1, 2, 3]
- Treatments: none · genetic modifications: none
- Default analysis ENCAN245MNX (/analyses/ENCAN245MNX/): GRCh38, ENCODE ATAC-seq pipeline v1.8.0, status released; other analyses: ENCAN741UVL v1.8.0 deleted
- Audits: NOT_COMPLIANT: insufficient read depth; WARNING: low read depth

| File | Output type | Replicates | Size | Pipeline file name | URL |
|---|---|---|---|---|---|
| ENCFF571JLT | signal p-value | rep1+2+3 | 1.64 GB | `rep.pooled.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF571JLT/@@download/ENCFF571JLT.bigWig |
| ENCFF573SLV | fold change over control | rep1+2+3 | 1.91 GB | `rep.pooled.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF573SLV/@@download/ENCFF573SLV.bigWig |
| ENCFF739NUO | signal p-value | rep1 | 912.6 MB | `ENCFF171VBP.trim.merged.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF739NUO/@@download/ENCFF739NUO.bigWig |
| ENCFF276XUH | signal p-value | rep2 | 353.0 MB | `ENCFF533VFC.trim.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF276XUH/@@download/ENCFF276XUH.bigWig |
| ENCFF195QTQ | signal p-value | rep3 | 613.9 MB | `ENCFF373MPT.trim.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF195QTQ/@@download/ENCFF195QTQ.bigWig |
| ENCFF905FFW | fold change over control | rep1 | 1.23 GB | `ENCFF171VBP.trim.merged.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF905FFW/@@download/ENCFF905FFW.bigWig |
| ENCFF691YCM | fold change over control | rep2 | 464.5 MB | `ENCFF533VFC.trim.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF691YCM/@@download/ENCFF691YCM.bigWig |
| ENCFF031ABQ | fold change over control | rep3 | 572.5 MB | `ENCFF373MPT.trim.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF031ABQ/@@download/ENCFF031ABQ.bigWig |
| ENCFF369SMB | IDR thresholded peaks | rep1+2+3 | 3.9 MB | `pooled-pr1_vs_pooled-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF369SMB/@@download/ENCFF369SMB.bed.gz |
| ENCFF785EFE | conservative IDR thresholded peaks | rep1+2+3 | 3.5 MB | `rep2_vs_rep3.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF785EFE/@@download/ENCFF785EFE.bed.gz |
| ENCFF060HHC | pseudoreplicated peaks | rep1+2+3 | 6.5 MB | `pooled-pr1_vs_pooled-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF060HHC/@@download/ENCFF060HHC.bed.gz |
| ENCFF418XBY | replicated peaks | rep1+2+3 | 5.7 MB | `rep1_vs_rep2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF418XBY/@@download/ENCFF418XBY.bed.gz |
| ENCFF206NEG | replicated peaks | rep1+2+3 | 5.8 MB | `rep1_vs_rep3.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF206NEG/@@download/ENCFF206NEG.bed.gz |
| ENCFF367CZS | replicated peaks | rep1+2+3 | 5.9 MB | `rep2_vs_rep3.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF367CZS/@@download/ENCFF367CZS.bed.gz |
| ENCFF786KPZ | IDR thresholded peaks | rep1 | 3.3 MB | `rep1-pr1_vs_rep1-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF786KPZ/@@download/ENCFF786KPZ.bed.gz |
| ENCFF089MWN | IDR thresholded peaks | rep2 | 1.9 MB | `rep2-pr1_vs_rep2-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF089MWN/@@download/ENCFF089MWN.bed.gz |
| ENCFF386YHU | IDR thresholded peaks | rep3 | 2.4 MB | `rep3-pr1_vs_rep3-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF386YHU/@@download/ENCFF386YHU.bed.gz |
| ENCFF824YZQ | pseudoreplicated peaks | rep1 | 6.0 MB | `rep1-pr1_vs_rep1-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF824YZQ/@@download/ENCFF824YZQ.bed.gz |
| ENCFF680EJZ | pseudoreplicated peaks | rep2 | 4.6 MB | `rep2-pr1_vs_rep2-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF680EJZ/@@download/ENCFF680EJZ.bed.gz |
| ENCFF113TWQ | pseudoreplicated peaks | rep3 | 5.3 MB | `rep3-pr1_vs_rep3-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF113TWQ/@@download/ENCFF113TWQ.bed.gz |

Additional true-replicate pairwise IDR files (bed narrowPeak, all labelled `IDR thresholded peaks`, reps 1+2+3): ENCFF746BDQ (`rep1_vs_rep3.idr0.05.bfilt.narrowPeak`, 3.3 MB), ENCFF045KOG (`rep1_vs_rep2.idr0.05.bfilt.narrowPeak`, 3.1 MB)

#### DND-41 — ENCSR660WSB

- Portal: https://www.encodeproject.org/experiments/ENCSR660WSB/ · lab Michael Snyder, Stanford · ENCODE4 · released 2021-04-06 · biosample summary: Homo sapiens DND-41 · biosamples ENCBS013RPN, ENCBS167JGE, ENCBS308ABW · bio reps [1, 2, 3]
- Treatments: none · genetic modifications: none
- Default analysis ENCAN022AWY (/analyses/ENCAN022AWY/): GRCh38, ENCODE ATAC-seq pipeline v1.9.1, status released; other analyses: ENCAN969UPI v1.8.0 deleted
- Audits: NOT_COMPLIANT: insufficient read depth; WARNING: low read depth

| File | Output type | Replicates | Size | Pipeline file name | URL |
|---|---|---|---|---|---|
| ENCFF631KUG | signal p-value | rep1+2+3 | 588.2 MB | `rep.pooled.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF631KUG/@@download/ENCFF631KUG.bigWig |
| ENCFF679DXT | fold change over control | rep1+2+3 | 980.4 MB | `rep.pooled.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF679DXT/@@download/ENCFF679DXT.bigWig |
| ENCFF296KBZ | signal p-value | rep1 | 288.4 MB | `ENCFF597TFO.trim.srt.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF296KBZ/@@download/ENCFF296KBZ.bigWig |
| ENCFF512SEH | signal p-value | rep2 | 313.2 MB | `ENCFF812FQU.trim.srt.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF512SEH/@@download/ENCFF512SEH.bigWig |
| ENCFF763EUH | signal p-value | rep3 | 275.7 MB | `ENCFF288JLG.trim.srt.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF763EUH/@@download/ENCFF763EUH.bigWig |
| ENCFF239VIO | fold change over control | rep1 | 384.3 MB | `ENCFF597TFO.trim.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF239VIO/@@download/ENCFF239VIO.bigWig |
| ENCFF299LJP | fold change over control | rep2 | 436.8 MB | `ENCFF812FQU.trim.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF299LJP/@@download/ENCFF299LJP.bigWig |
| ENCFF113SPI | fold change over control | rep3 | 391.6 MB | `ENCFF288JLG.trim.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF113SPI/@@download/ENCFF113SPI.bigWig |
| ENCFF225RUL | IDR thresholded peaks | rep1+2+3 | 2.5 MB | `pooled-pr1_vs_pooled-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF225RUL/@@download/ENCFF225RUL.bed.gz |
| ENCFF788MOK | conservative IDR thresholded peaks | rep1+2+3 | 2.4 MB | `rep2_vs_rep3.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF788MOK/@@download/ENCFF788MOK.bed.gz |
| ENCFF095ATB | pseudoreplicated peaks | rep1+2+3 | 3.9 MB | `pooled-pr1_vs_pooled-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF095ATB/@@download/ENCFF095ATB.bed.gz |
| ENCFF381HCM | replicated peaks | rep1+2+3 | 3.6 MB | `rep1_vs_rep2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF381HCM/@@download/ENCFF381HCM.bed.gz |
| ENCFF222REH | replicated peaks | rep1+2+3 | 3.7 MB | `rep1_vs_rep3.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF222REH/@@download/ENCFF222REH.bed.gz |
| ENCFF959GGH | replicated peaks | rep1+2+3 | 3.7 MB | `rep2_vs_rep3.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF959GGH/@@download/ENCFF959GGH.bed.gz |
| ENCFF126TTT | IDR thresholded peaks | rep1 | 1.7 MB | `rep1-pr1_vs_rep1-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF126TTT/@@download/ENCFF126TTT.bed.gz |
| ENCFF455OSM | IDR thresholded peaks | rep2 | 1.9 MB | `rep2-pr1_vs_rep2-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF455OSM/@@download/ENCFF455OSM.bed.gz |
| ENCFF073XWR | IDR thresholded peaks | rep3 | 1.8 MB | `rep3-pr1_vs_rep3-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF073XWR/@@download/ENCFF073XWR.bed.gz |
| ENCFF110QQR | pseudoreplicated peaks | rep1 | 2.8 MB | `rep1-pr1_vs_rep1-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF110QQR/@@download/ENCFF110QQR.bed.gz |
| ENCFF054ATC | pseudoreplicated peaks | rep2 | 3.1 MB | `rep2-pr1_vs_rep2-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF054ATC/@@download/ENCFF054ATC.bed.gz |
| ENCFF147RPW | pseudoreplicated peaks | rep3 | 2.9 MB | `rep3-pr1_vs_rep3-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF147RPW/@@download/ENCFF147RPW.bed.gz |

Additional true-replicate pairwise IDR files (bed narrowPeak, all labelled `IDR thresholded peaks`, reps 1+2+3): ENCFF478CDQ (`rep1_vs_rep3.idr0.05.bfilt.narrowPeak`, 2.4 MB), ENCFF635DAN (`rep1_vs_rep2.idr0.05.bfilt.narrowPeak`, 2.4 MB)

#### RWPE2 — ENCSR080SNF

- Portal: https://www.encodeproject.org/experiments/ENCSR080SNF/ · lab Michael Snyder, Stanford · ENCODE4 · released 2021-04-06 · biosample summary: Homo sapiens RWPE2 · biosamples ENCBS123OBK, ENCBS775VKW, ENCBS961CAJ · bio reps [1, 2, 3]
- Treatments: none · genetic modifications: none
- Default analysis ENCAN367YCE (/analyses/ENCAN367YCE/): GRCh38, ENCODE ATAC-seq pipeline v1.9.1, status released; other analyses: ENCAN195UMD v1.8.0 deleted
- Audits: WARNING: low read depth

| File | Output type | Replicates | Size | Pipeline file name | URL |
|---|---|---|---|---|---|
| ENCFF881UWW | signal p-value | rep1+2+3 | 997.3 MB | `rep.pooled.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF881UWW/@@download/ENCFF881UWW.bigWig |
| ENCFF047EBV | fold change over control | rep1+2+3 | 1.52 GB | `rep.pooled.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF047EBV/@@download/ENCFF047EBV.bigWig |
| ENCFF407MVF | signal p-value | rep1 | 532.4 MB | `ENCFF374LMV.trim.srt.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF407MVF/@@download/ENCFF407MVF.bigWig |
| ENCFF448JGC | signal p-value | rep2 | 534.7 MB | `ENCFF598ZYD.trim.srt.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF448JGC/@@download/ENCFF448JGC.bigWig |
| ENCFF316DBN | signal p-value | rep3 | 448.7 MB | `ENCFF855HKQ.trim.srt.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF316DBN/@@download/ENCFF316DBN.bigWig |
| ENCFF547BBS | fold change over control | rep1 | 678.3 MB | `ENCFF374LMV.trim.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF547BBS/@@download/ENCFF547BBS.bigWig |
| ENCFF808OVI | fold change over control | rep2 | 603.9 MB | `ENCFF598ZYD.trim.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF808OVI/@@download/ENCFF808OVI.bigWig |
| ENCFF051JJI | fold change over control | rep3 | 551.6 MB | `ENCFF855HKQ.trim.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF051JJI/@@download/ENCFF051JJI.bigWig |
| ENCFF646WSM | IDR thresholded peaks | rep1+2+3 | 3.8 MB | `pooled-pr1_vs_pooled-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF646WSM/@@download/ENCFF646WSM.bed.gz |
| ENCFF846VVX | conservative IDR thresholded peaks | rep1+2+3 | 3.6 MB | `rep1_vs_rep2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF846VVX/@@download/ENCFF846VVX.bed.gz |
| ENCFF729MMJ | pseudoreplicated peaks | rep1+2+3 | 6.3 MB | `pooled-pr1_vs_pooled-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF729MMJ/@@download/ENCFF729MMJ.bed.gz |
| ENCFF934EAT | replicated peaks | rep1+2+3 | 5.8 MB | `rep1_vs_rep3.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF934EAT/@@download/ENCFF934EAT.bed.gz |
| ENCFF813RPV | replicated peaks | rep1+2+3 | 5.9 MB | `rep1_vs_rep2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF813RPV/@@download/ENCFF813RPV.bed.gz |
| ENCFF939WLZ | replicated peaks | rep1+2+3 | 5.7 MB | `rep2_vs_rep3.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF939WLZ/@@download/ENCFF939WLZ.bed.gz |
| ENCFF663XSQ | IDR thresholded peaks | rep1 | 2.8 MB | `rep1-pr1_vs_rep1-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF663XSQ/@@download/ENCFF663XSQ.bed.gz |
| ENCFF104IMP | IDR thresholded peaks | rep2 | 2.5 MB | `rep2-pr1_vs_rep2-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF104IMP/@@download/ENCFF104IMP.bed.gz |
| ENCFF652GKL | IDR thresholded peaks | rep3 | 2.4 MB | `rep3-pr1_vs_rep3-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF652GKL/@@download/ENCFF652GKL.bed.gz |
| ENCFF630PFX | pseudoreplicated peaks | rep1 | 5.2 MB | `rep1-pr1_vs_rep1-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF630PFX/@@download/ENCFF630PFX.bed.gz |
| ENCFF043RPD | pseudoreplicated peaks | rep2 | 4.7 MB | `rep2-pr1_vs_rep2-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF043RPD/@@download/ENCFF043RPD.bed.gz |
| ENCFF230JLB | pseudoreplicated peaks | rep3 | 4.7 MB | `rep3-pr1_vs_rep3-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF230JLB/@@download/ENCFF230JLB.bed.gz |

Additional true-replicate pairwise IDR files (bed narrowPeak, all labelled `IDR thresholded peaks`, reps 1+2+3): ENCFF818HSM (`rep1_vs_rep3.idr0.05.bfilt.narrowPeak`, 3.5 MB), ENCFF033QDT (`rep2_vs_rep3.idr0.05.bfilt.narrowPeak`, 3.4 MB)

#### WTC11 — ENCSR541KFY

- Portal: https://www.encodeproject.org/experiments/ENCSR541KFY/ · lab Michael Snyder, Stanford · ENCODE4 · released 2021-06-24 · biosample summary: Homo sapiens WTC11 · biosamples ENCBS503FBI, ENCBS732JEF · bio reps [1, 2]
- Treatments: none · genetic modifications: none
- Default analysis ENCAN873DEH (/analyses/ENCAN873DEH/): GRCh38, ENCODE ATAC-seq pipeline v1.10.0, status released; other analyses: 
- Audits: WARNING: mild to moderate bottlenecking, moderate library complexity

| File | Output type | Replicates | Size | Pipeline file name | URL |
|---|---|---|---|---|---|
| ENCFF123YPY | signal p-value | rep1+2 | 1.12 GB | `rep.pooled.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF123YPY/@@download/ENCFF123YPY.bigWig |
| ENCFF056MYW | fold change over control | rep1+2 | 1.46 GB | `rep.pooled.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF056MYW/@@download/ENCFF056MYW.bigWig |
| ENCFF058HBV | signal p-value | rep1 | 677.5 MB | `ENCFF016TLO.trim.srt.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF058HBV/@@download/ENCFF058HBV.bigWig |
| ENCFF427BHI | signal p-value | rep2 | 730.3 MB | `ENCFF638XYC.trim.srt.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF427BHI/@@download/ENCFF427BHI.bigWig |
| ENCFF938VZR | fold change over control | rep1 | 826.8 MB | `ENCFF016TLO.trim.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF938VZR/@@download/ENCFF938VZR.bigWig |
| ENCFF456IER | fold change over control | rep2 | 858.9 MB | `ENCFF638XYC.trim.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF456IER/@@download/ENCFF456IER.bigWig |
| ENCFF492GYV | IDR thresholded peaks | rep1+2 | 3.8 MB | `pooled-pr1_vs_pooled-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF492GYV/@@download/ENCFF492GYV.bed.gz |
| ENCFF121CAA | conservative IDR thresholded peaks | rep1+2 | 3.8 MB | `rep1_vs_rep2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF121CAA/@@download/ENCFF121CAA.bed.gz |
| ENCFF321VDH | pseudoreplicated peaks | rep1+2 | 6.6 MB | `pooled-pr1_vs_pooled-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF321VDH/@@download/ENCFF321VDH.bed.gz |
| ENCFF704LOO | replicated peaks | rep1+2 | 6.6 MB | `rep1_vs_rep2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF704LOO/@@download/ENCFF704LOO.bed.gz |
| ENCFF770QET | IDR thresholded peaks | rep1 | 3.1 MB | `rep1-pr1_vs_rep1-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF770QET/@@download/ENCFF770QET.bed.gz |
| ENCFF262EMM | IDR thresholded peaks | rep2 | 3.1 MB | `rep2-pr1_vs_rep2-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF262EMM/@@download/ENCFF262EMM.bed.gz |
| ENCFF694YTL | pseudoreplicated peaks | rep1 | 6.0 MB | `rep1-pr1_vs_rep1-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF694YTL/@@download/ENCFF694YTL.bed.gz |
| ENCFF459OUO | pseudoreplicated peaks | rep2 | 6.1 MB | `rep2-pr1_vs_rep2-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF459OUO/@@download/ENCFF459OUO.bed.gz |

### Alternates

#### K562 — ENCSR483RKN

- Portal: https://www.encodeproject.org/experiments/ENCSR483RKN/ · lab Michael Snyder, Stanford · ENCODE4 · released 2020-07-21 · biosample summary: Homo sapiens K562 · biosamples ENCBS321NQQ, ENCBS644OPH · bio reps [1, 2]
- Treatments: none · genetic modifications: none
- Default analysis ENCAN217QUL (/analyses/ENCAN217QUL/): GRCh38, ENCODE ATAC-seq pipeline v1.9.1, status released; other analyses: ENCAN252OFV v1.8.0 revoked
- Audits: WARNING: low read depth, mild to moderate bottlenecking, moderate library complexity

| File | Output type | Replicates | Size | Pipeline file name | URL |
|---|---|---|---|---|---|
| ENCFF600FDO | signal p-value | rep1+2 | 864.2 MB | `rep.pooled.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF600FDO/@@download/ENCFF600FDO.bigWig |
| ENCFF754EAC | fold change over control | rep1+2 | 1.15 GB | `rep.pooled.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF754EAC/@@download/ENCFF754EAC.bigWig |
| ENCFF459KUR | signal p-value | rep1 | 548.0 MB | `ENCFF391BFJ.trim.srt.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF459KUR/@@download/ENCFF459KUR.bigWig |
| ENCFF867QEW | signal p-value | rep2 | 477.2 MB | `ENCFF440UAD.trim.srt.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF867QEW/@@download/ENCFF867QEW.bigWig |
| ENCFF451IKJ | fold change over control | rep1 | 658.6 MB | `ENCFF391BFJ.trim.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF451IKJ/@@download/ENCFF451IKJ.bigWig |
| ENCFF496GRF | fold change over control | rep2 | 647.7 MB | `ENCFF440UAD.trim.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF496GRF/@@download/ENCFF496GRF.bigWig |
| ENCFF925CYR | IDR thresholded peaks | rep1+2 | 2.8 MB | `pooled-pr1_vs_pooled-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF925CYR/@@download/ENCFF925CYR.bed.gz |
| ENCFF695IGF | conservative IDR thresholded peaks | rep1+2 | 2.7 MB | `rep1_vs_rep2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF695IGF/@@download/ENCFF695IGF.bed.gz |
| ENCFF558BLC | pseudoreplicated peaks | rep1+2 | 5.0 MB | `pooled-pr1_vs_pooled-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF558BLC/@@download/ENCFF558BLC.bed.gz |
| ENCFF091NSB | replicated peaks | rep1+2 | 4.9 MB | `rep1_vs_rep2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF091NSB/@@download/ENCFF091NSB.bed.gz |
| ENCFF117MSK | IDR thresholded peaks | rep1 | 2.0 MB | `rep1-pr1_vs_rep1-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF117MSK/@@download/ENCFF117MSK.bed.gz |
| ENCFF976CEI | IDR thresholded peaks | rep2 | 1.9 MB | `rep2-pr1_vs_rep2-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF976CEI/@@download/ENCFF976CEI.bed.gz |
| ENCFF656YON | pseudoreplicated peaks | rep1 | 3.8 MB | `rep1-pr1_vs_rep1-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF656YON/@@download/ENCFF656YON.bed.gz |
| ENCFF858JFG | pseudoreplicated peaks | rep2 | 3.6 MB | `rep2-pr1_vs_rep2-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF858JFG/@@download/ENCFF858JFG.bed.gz |

#### HepG2 — ENCSR042AWH

- Portal: https://www.encodeproject.org/experiments/ENCSR042AWH/ · lab Michael Snyder, Stanford · ENCODE4 · released 2020-08-05 · biosample summary: Homo sapiens HepG2 · biosamples ENCBS495DEU, ENCBS502XQO · bio reps [1, 2]
- Treatments: none · genetic modifications: none
- Default analysis ENCAN074ZCE (/analyses/ENCAN074ZCE/): GRCh38, ENCODE ATAC-seq pipeline v1.9.1, status released; other analyses: ENCAN015VOH v1.8.0 revoked
- Audits: WARNING: low read depth, mild to moderate bottlenecking, moderate library complexity

| File | Output type | Replicates | Size | Pipeline file name | URL |
|---|---|---|---|---|---|
| ENCFF285FQS | signal p-value | rep1+2 | 773.2 MB | `rep.pooled.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF285FQS/@@download/ENCFF285FQS.bigWig |
| ENCFF664EJT | fold change over control | rep1+2 | 1.11 GB | `rep.pooled.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF664EJT/@@download/ENCFF664EJT.bigWig |
| ENCFF699PHI | signal p-value | rep1 | 416.4 MB | `ENCFF664UPL.trim.srt.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF699PHI/@@download/ENCFF699PHI.bigWig |
| ENCFF167CGQ | signal p-value | rep2 | 445.1 MB | `ENCFF094QZE.trim.srt.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF167CGQ/@@download/ENCFF167CGQ.bigWig |
| ENCFF782GKX | fold change over control | rep1 | 602.1 MB | `ENCFF664UPL.trim.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF782GKX/@@download/ENCFF782GKX.bigWig |
| ENCFF029XKY | fold change over control | rep2 | 622.9 MB | `ENCFF094QZE.trim.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF029XKY/@@download/ENCFF029XKY.bigWig |
| ENCFF438JMM | IDR thresholded peaks | rep1+2 | 2.9 MB | `pooled-pr1_vs_pooled-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF438JMM/@@download/ENCFF438JMM.bed.gz |
| ENCFF314RIW | conservative IDR thresholded peaks | rep1+2 | 2.9 MB | `rep1_vs_rep2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF314RIW/@@download/ENCFF314RIW.bed.gz |
| ENCFF913MQB | pseudoreplicated peaks | rep1+2 | 5.3 MB | `pooled-pr1_vs_pooled-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF913MQB/@@download/ENCFF913MQB.bed.gz |
| ENCFF263JJA | replicated peaks | rep1+2 | 5.2 MB | `rep1_vs_rep2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF263JJA/@@download/ENCFF263JJA.bed.gz |
| ENCFF791RKW | IDR thresholded peaks | rep1 | 2.1 MB | `rep1-pr1_vs_rep1-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF791RKW/@@download/ENCFF791RKW.bed.gz |
| ENCFF851SSE | IDR thresholded peaks | rep2 | 1.9 MB | `rep2-pr1_vs_rep2-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF851SSE/@@download/ENCFF851SSE.bed.gz |
| ENCFF759HTS | pseudoreplicated peaks | rep1 | 4.0 MB | `rep1-pr1_vs_rep1-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF759HTS/@@download/ENCFF759HTS.bed.gz |
| ENCFF917ZKL | pseudoreplicated peaks | rep2 | 4.0 MB | `rep2-pr1_vs_rep2-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF917ZKL/@@download/ENCFF917ZKL.bed.gz |

#### GM12878 — ENCSR095QNB

- Portal: https://www.encodeproject.org/experiments/ENCSR095QNB/ · lab Michael Snyder, Stanford · ENCODE4 · released 2020-08-05 · biosample summary: Homo sapiens GM12878 · biosamples ENCBS530ZPB, ENCBS923ZKG · bio reps [1, 2]
- Treatments: none · genetic modifications: none
- Default analysis ENCAN039YHB (/analyses/ENCAN039YHB/): GRCh38, ENCODE ATAC-seq pipeline v1.9.1, status released; other analyses: ENCAN054YRF v1.8.0 revoked
- Audits: WARNING: low read depth, mild to moderate bottlenecking, moderate library complexity

| File | Output type | Replicates | Size | Pipeline file name | URL |
|---|---|---|---|---|---|
| ENCFF180ZAY | signal p-value | rep1+2 | 502.1 MB | `rep.pooled.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF180ZAY/@@download/ENCFF180ZAY.bigWig |
| ENCFF487LOB | fold change over control | rep1+2 | 869.8 MB | `rep.pooled.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF487LOB/@@download/ENCFF487LOB.bigWig |
| ENCFF225HCI | signal p-value | rep1 | 323.4 MB | `ENCFF411MVZ.trim.srt.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF225HCI/@@download/ENCFF225HCI.bigWig |
| ENCFF672ZZF | signal p-value | rep2 | 398.9 MB | `ENCFF196VAZ.trim.srt.nodup.no_chrM_MT.tn5.pval.signal.bigwig` | https://www.encodeproject.org/files/ENCFF672ZZF/@@download/ENCFF672ZZF.bigWig |
| ENCFF316PVI | fold change over control | rep1 | 457.2 MB | `ENCFF411MVZ.trim.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF316PVI/@@download/ENCFF316PVI.bigWig |
| ENCFF745TKD | fold change over control | rep2 | 550.4 MB | `ENCFF196VAZ.trim.srt.nodup.no_chrM_MT.tn5.fc.signal.bigwig` | https://www.encodeproject.org/files/ENCFF745TKD/@@download/ENCFF745TKD.bigWig |
| ENCFF945SYZ | IDR thresholded peaks | rep1+2 | 3.5 MB | `pooled-pr1_vs_pooled-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF945SYZ/@@download/ENCFF945SYZ.bed.gz |
| ENCFF346CZA | conservative IDR thresholded peaks | rep1+2 | 3.2 MB | `rep1_vs_rep2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF346CZA/@@download/ENCFF346CZA.bed.gz |
| ENCFF470YYO | pseudoreplicated peaks | rep1+2 | 5.5 MB | `pooled-pr1_vs_pooled-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF470YYO/@@download/ENCFF470YYO.bed.gz |
| ENCFF478ZAW | replicated peaks | rep1+2 | 5.3 MB | `rep1_vs_rep2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF478ZAW/@@download/ENCFF478ZAW.bed.gz |
| ENCFF917REN | IDR thresholded peaks | rep1 | 2.4 MB | `rep1-pr1_vs_rep1-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF917REN/@@download/ENCFF917REN.bed.gz |
| ENCFF614SMH | IDR thresholded peaks | rep2 | 2.8 MB | `rep2-pr1_vs_rep2-pr2.idr0.05.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF614SMH/@@download/ENCFF614SMH.bed.gz |
| ENCFF851FBD | pseudoreplicated peaks | rep1 | 4.2 MB | `rep1-pr1_vs_rep1-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF851FBD/@@download/ENCFF851FBD.bed.gz |
| ENCFF158ORB | pseudoreplicated peaks | rep2 | 4.8 MB | `rep2-pr1_vs_rep2-pr2.overlap.bfilt.narrowPeak` | https://www.encodeproject.org/files/ENCFF158ORB/@@download/ENCFF158ORB.bed.gz |

## hTERT-RPE1 — 4DN 4DNES48BA2C8 detail

- Portal: https://data.4dnucleome.org/experiment-set-replicates/4DNES48BA2C8/ · lab Bing Ren, UCSD · award /awards/1UM1HG011585-01/ (CENTER FOR INTEGRATED MULTI-MODAL AND MULTI-SCALE NUCLEOME RESEARCH) · status released · public release 2025-06-12 · dataset 'bulk ATAC-seq in human eye tissue' / condition 'hTERT-RPE1 - WT' · description: Bulk ATAC-seq in human hTERT-RPE-1 cell line - Wild type (used as control)
- Publications: produced_in_pub = None; publications_of_set = [] → **unpublished dataset**.
- Sibling sets in the same dataset (hTERT-RPE1 clones, same lab/pipeline): 4DNES7YK6HTM (F6 clone), 4DNESJWVEOAK (G1 clone), 4DNES1TSS4EO (D4 clone). Use 4DNES48BA2C8 (WT) only.
- Replicates: 2 biological replicate experiments, each with 2 paired-end fastq (~7.3–7.5 GB each):
  - 4DNEX7A3MJJV: raw 4DNFIDNZZ2A9 (7.36 GB), 4DNFI14DT5NE (7.32 GB); per-replicate processed: 4DNFIM6LAD9J read positions (1.68 GB)
  - 4DNEX5R8DYKC: raw 4DNFIWGX7K9Y (7.49 GB), 4DNFIY8ZSQC1 (7.43 GB); per-replicate processed: 4DNFIIFNCJN9 read positions (1.72 GB)
- Set-level processed files (GRCh38, workflow `encode-atacseq-postaln 1.1.1` = 4DN CWL port of the ENCODE ATAC-seq pipeline, run 2025-06-04; both replicates pooled):

| File | Type | Size | Assembly | URL (portal) | Open-data S3 URL | md5 |
|---|---|---|---|---|---|---|
| 4DNFIS4IDZNE | signal fold change — ATAC-seq signal fold change | 4.41 GB | GRCh38 | https://data.4dnucleome.org/files-processed/4DNFIS4IDZNE/@@download/4DNFIS4IDZNE.bw | https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/df7363e3-3a77-40a9-bc8b-e0f07613055f/4DNFIS4IDZNE.bw | 5eff374d31e11c2ad61dce83e6811ff6 |
| 4DNFI4H1MT8U | conservative peaks — Conservative peak calls from ENCODE ATAC-Seq Pipeline | 15.2 MB | GRCh38 | https://data.4dnucleome.org/files-processed/4DNFI4H1MT8U/@@download/4DNFI4H1MT8U.bb | https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/de7ab24d-0e46-4d1a-81f2-4ee1c7956c5d/4DNFI4H1MT8U.bb | 43f76dfe837b9c5d29e7958ad4d3df9d |
| 4DNFIMNHSGU1 | peaks — Peak calls from ENCODE ATAC-Seq Pipeline | 15.2 MB | GRCh38 | https://data.4dnucleome.org/files-processed/4DNFIMNHSGU1/@@download/4DNFIMNHSGU1.bb | https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/7daaa35e-81ca-4094-9cad-f21a76d0c090/4DNFIMNHSGU1.bb | 718785329c877533e93061c6bb7d348c |
| 4DNFIKJGW79V | supplementary bigWig (no file_type; likely per-replicate signal) —  | 1.35 GB | GRCh38 | https://data.4dnucleome.org/files-processed/4DNFIKJGW79V/@@download/4DNFIKJGW79V.bw | https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/572e068d-e298-46aa-92f9-4ceea70dcd6e/4DNFIKJGW79V.bw | 4ceb07954828dfb69b6814886c9af36e |
| 4DNFID2TXF3X | supplementary bigWig (no file_type; likely per-replicate signal) —  | 1.37 GB | GRCh38 | https://data.4dnucleome.org/files-processed/4DNFID2TXF3X/@@download/4DNFID2TXF3X.bw | https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/cdc0f634-21d0-4807-adc7-7cb66af18e73/4DNFID2TXF3X.bw | 1bd59ae7acccceeae06f0deea9bcbe73 |
| 4DNFIM6LAD9J | read positions — Positions of aligned reads in bed format, one line per read mate, from ENCODE ATAC-Seq Pipeline | 1.68 GB | GRCh38 | https://data.4dnucleome.org/files-processed/4DNFIM6LAD9J/@@download/4DNFIM6LAD9J.bed.gz | https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/1e64bb00-5da7-4714-bc81-ff096702d2a8/4DNFIM6LAD9J.bed.gz | 6084c1c7180200e8e603cf907f721253 |
| 4DNFIIFNCJN9 | read positions — Positions of aligned reads in bed format, one line per read mate, from ENCODE ATAC-Seq Pipeline | 1.72 GB | GRCh38 | https://data.4dnucleome.org/files-processed/4DNFIIFNCJN9/@@download/4DNFIIFNCJN9.bed.gz | https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/6562a458-b64e-45f8-8a53-80ab83ce2511/4DNFIIFNCJN9.bed.gz | 75cac834202f3a612037091d8783aa7d |

Notes: 4DN emits peaks as **bigBed** (`peaks` = pooled-pseudoreplicate overlap/naive set; `conservative peaks` = IDR conservative set per the 4DN description 'Conservative peak calls from ENCODE ATAC-Seq Pipeline'); convert with `bigBedToBed`. 4DN provides a pooled **fold-change** bigWig only (no signal p-value track). The ENCODE pipeline version run by 4DN (1.1.1 CWL port, ~ENCODE atac-seq-pipeline v1.1) predates ENCODE's v1.9.1 default analyses; alignment (bowtie2, GRCh38), filtering, MACS2 and IDR logic are the same lineage but parameter/version drift (e.g. MACS2 p-value cap, TSS enrichment, blacklist version) is possible. To be strictly uniform, re-run the ENCODE pipeline from fastq for all four target lines (ATAC fastq for the ENCODE sets are listed on the portal; RPE1 fastq = 4 × ~7.4 GB).

### 4DN data-use policy for 4DNES48BA2C8

Per the 4DN Data Use Guidelines (https://data.4dnucleome.org/resources/data-use-guidelines, see also https://www.4dnucleome.org/outreach/data-use-policy/): all released 4DN data are freely usable **unless** the dataset has not yet been described in a publication by the producing lab (`produced_in_pub` empty → unpublished, as here, public release 2025-06-12). For such datasets, users wishing to publish analyses must **contact the data-producing lab (Bing Ren, UCSD — contact listed on the lab page https://data.4dnucleome.org/labs/bing-ren-lab/) before submitting the manuscript**, give the lab the opportunity to comment/participate, cite the dataset accession (4DNES48BA2C8) and acknowledge the 4DN Network and the award (1UM1HG011585-01). (Policy text should be re-verified on the portal when the manuscript is drafted; this is a metadata-derived summary.)

## Blockers and caveats

1. **Jurkat (clone E6-1): no ATAC-seq anywhere on ENCODE or 4DN.** ENCODE search for `biosample_ontology.term_name` = 'Jurkat, Clone E6-1' / 'Jurkat clone E6-1' and free-text 'Jurkat' with `assay_title=ATAC-seq` returns 0 experiments (any status, any assembly); 4DN has no Jurkat ATAC-seq ExperimentSet. ENCODE's Jurkat data are DNase-seq / ChIP-seq / RNA-seq only. Options: (a) download public Jurkat ATAC fastq from GEO/SRA (e.g. Corces et al. 2016 *Nat Genet* GSE74912 includes Jurkat; Buenrostro 2013 Omni-ATAC data) and run the ENCODE ATAC-seq pipeline (v1.9.1 or later Caper/WDL release) locally; (b) substitute the ENCODE T-ALL line **DND-41 (ENCSR660WSB, 3 reps, v1.9.1)** as a lineage-matched surrogate; (c) use ENCODE Jurkat DNase-seq as a non-uniform proxy. Any choice must be decided before the experiment design is locked.
2. **hTERT-RPE1: ENCODE has no RPE1 ATAC-seq** (0 hits for 'hTERT-RPE1'/'RPE1'/searchTerm 'RPE'). Only 4DN 4DNES48BA2C8 (WT) exists, processed with 4DN's CWL port of the ENCODE pipeline (1.1.1) — not the same release as ENCODE's v1.9.1 default analyses — with bigBed peaks and a fold-change-only bigWig. Strict pipeline uniformity therefore requires re-processing RPE1 fastq (≈30 GB raw) with the ENCODE pipeline, or accepting a documented pipeline-version mismatch for one of four target lines. Data-use: unpublished → contact Ren lab before publishing.
3. **SK-N-SH decoy is RA-treated** (ENCSR587TRP: all-trans-retinoic acid 6 µM 48 h); it is the only SK-N-SH ATAC on ENCODE. Treat as 'SK-N-SH + RA' or drop. HCT116 has 17 ATAC experiments, 16 of which are degron-tagged/5-Ph-IAA-treated; only ENCSR872WGW is plain HCT116.
4. **HeLa-S3 and HEK293 are not available as ENCODE ATAC-seq**; the decoy panel above substitutes GM23338, NCI-H929, DND-41, RWPE2, WTC11 (all Snyder ENCODE4, v1.9.1). Core decoy set of 7 untreated lines: GM12878, HCT116, MCF-7, A549, IMR-90, PC-3, Panc1 (+ optional SK-N-SH-RA, GM23338, NCI-H929, DND-41, RWPE2, WTC11).
5. **K562 has 64 ATAC experiments**; 62 are Greenleaf-lab perturbation series (drug treatments such as JQ1; ENCSR165JXS checked) and must not be used as untreated references. Untreated: ENCSR868FGK (3 reps, used) and ENCSR483RKN (2 reps).
6. ENCODE 'IDR thresholded peaks' for ATAC comes in several flavours per experiment (per-replicate self-pseudoreplicate IDR, pairwise true-replicate IDR, pooled-pseudoreplicate IDR); only `output_type` + `analysis_step_version` distinguish them. This inventory pins the pooled-pseudoreplicate IDR file (and the `conservative IDR` true-rep file) explicitly — do not pick by `output_type` alone.
7. Audit warnings: several experiments carry 'moderate library complexity' / 'mild to moderate bottlenecking' WARNINGs (listed per experiment above); none carry ERROR-level audits. In the TF ChIP CSVs, 30 K562 and 19 HepG2 experiments carry an ERROR audit (mostly 'extremely low read depth') — filter before use.
8. **Pipeline-version drift within ENCODE**: NCI-H929 (v1.8.0) and WTC11 (v1.10.0) default analyses differ from the v1.9.1 used by the other 16 experiments; prefer the 7-line core decoy set (all v1.9.1) or re-process.
9. Totals/sizes: pooled signal p-value + fold-change bigWigs are 0.7–4.5 GB each; a full download of the pooled tracks for the 15 ENCODE experiments + RPE1 is on the order of 60–70 GB — plan storage accordingly (sizes are per-file in the tables).
