# Perturb-seq data-access inventory (recon 01)

Date checked: 2026-10-11. All URL checks done with `curl` from this Mac (HEAD, or ranged GET of bytes 0–7 where HEAD is blocked). No multi-GB file was downloaded. Raw check log: `recon/perturbseq_url_checks.json`.

Verification legend: **HEAD 200** = `curl -I` returned 200 with Content-Length; **302→206** = redirector returns 302 to a short-lived (10 s) presigned S3 URL, and a ranged GET on it returns `206 Partial Content` with `Content-Range: bytes 0-7/<total>` and HDF5 magic bytes (`\x89HDF\r\n\x1a\n`).

---

## 1. Replogle et al. 2022 (Cell) — K562 essential, K562 genome-wide, RPE1 essential

- Paper: Replogle JM, Saunders RA, …, Norman TM, Weissman JS. *Cell* 185:2559 (2022). https://doi.org/10.1016/j.cell.2022.05.013
- Portal: https://gwps.wi.mit.edu/ (HEAD 200). Lists four access tiers: SRA FASTQ/BAM (BioProject PRJNA831566), MEX counts, processed h5ad, supplemental.
- **Primary processed data (h5ad)**: Figshare+ article 20029387, DOI https://doi.org/10.25452/figshare.plus.20029387 (DOI → 302 → article page; article page returns HTTP 202 to HEAD — Cloudflare challenge — but the Figshare REST API `https://api.figshare.com/v2/articles/20029387` returns 200 with the full file list).
- **Licence: CC BY 4.0** (Figshare API `license.name`). Supplemental article 21632564 is CC0.
- Published 2022-06-09, version 1.

### 1a. Files (Figshare 20029387) — sizes from the Figshare API, cross-checked by ranged GET

| File | Bytes | GB | Download URL | Check |
|---|---:|---:|---|---|
| `K562_essential_raw_singlecell_01.h5ad` | 10,661,879,995 | 10.66 | https://ndownloader.figshare.com/files/35773219 | 302→206, HDF5 ✓ |
| `K562_essential_normalized_singlecell_01.h5ad` | ~10.66e9 | 10.66 | https://ndownloader.figshare.com/files/35773075 | 302 ✓ |
| `K562_essential_raw_bulk_01.h5ad` | ~80e6 | 0.08 | https://ndownloader.figshare.com/files/35773070 | 302 ✓ |
| `K562_essential_normalized_bulk_01.h5ad` | ~80e6 | 0.08 | https://ndownloader.figshare.com/files/35780870 | (API-listed) |
| `K562_gwps_raw_singlecell_01.h5ad` | 65,830,941,948 | 65.83 | https://ndownloader.figshare.com/files/35775507 | 302→206, HDF5 ✓ |
| `K562_gwps_normalized_singlecell_01.h5ad` | ~65.83e9 | 65.83 | https://ndownloader.figshare.com/files/35774440 | 302 ✓ |
| `K562_gwps_raw_bulk_01.h5ad` | ~375e6 | 0.375 | https://ndownloader.figshare.com/files/35774443 | 302 ✓ |
| `K562_gwps_normalized_bulk_01.h5ad` | ~375e6 | 0.375 | https://ndownloader.figshare.com/files/35773217 | (API-listed) |
| `rpe1_raw_singlecell_01.h5ad` | 8,700,873,216 | 8.70 | https://ndownloader.figshare.com/files/35775606 | 302→206, HDF5 ✓ |
| `rpe1_normalized_singlecell_01.h5ad` | ~8.70e9 | 8.70 | https://ndownloader.figshare.com/files/35775554 | 302 ✓ |
| `rpe1_raw_bulk_01.h5ad` | ~95e6 | 0.095 | https://ndownloader.figshare.com/files/35775581 | 302 ✓ |
| `rpe1_normalized_bulk_01.h5ad` | ~95e6 | 0.095 | https://ndownloader.figshare.com/files/35775512 | (API-listed) |

Format (from the Figshare description): AnnData `.h5ad`; `raw_singlecell` = raw UMI counts per cell for genes with >0.01 UMI/cell; `normalized_singlecell` = gem-group Z-normalised; `*_bulk` = pseudobulk per perturbation. `.var` = genes, `.obs` = cells. "Data were aligned, UMI-collapsed, aggregated, and subset to cells bearing sgRNA(s) targeting a single target gene" (gwps portal).

Note: the "raw" and "normalized" single-cell files have identical byte sizes — likely dense float storage; the raw file is still the one to use for counts.

### 1b. Alternative formats
- **MEX (cellranger `filtered_feature_bc_matrix`) tarballs**, Figshare 20127869 (CC BY 4.0): `K562_gwps_MTX.tar.gz` 48.98 GB (files/36000572), `K562_essential_MTX.tar.gz` 8.93 GB (files/36000843), `rpe1_MTX.tar.gz` 8.20 GB (files/36000888), `ultima_MTX.tar.gz` 37.93 GB (files/36000884). These are per-gem-group cellranger outputs and do **not** carry the guide-assignment columns — use the h5ad.
- **SRA/GEO manifest**, Figshare 20022944 (CC BY 4.0; DOI → 302): `KD6_raw_files.csv`, `KD8_raw_files.csv`, `KD8_ultima_raw_files.csv`, `RD7_raw_files.csv` (gem-group ↔ FASTQ mapping). Raw reads: BioProject PRJNA831566. There is **no GEO GSE** for the 2022 paper — data are on SRA + Figshare only.
- **scPerturb harmonised copies** (Zenodo 13350497, v1.4, CC BY 4.0, HEAD 200 with a browser User-Agent; plain-UA requests get 403):
  `ReplogleWeissman2022_K562_essential.h5ad` 1.547 GB, `ReplogleWeissman2022_K562_gwps.h5ad` 8.805 GB, `ReplogleWeissman2022_rpe1.h5ad` 1.237 GB — URL pattern `https://zenodo.org/api/records/13350497/files/<name>/content`. Same cells, gzip-compressed, obs harmonised (`perturbation`, `guide_id`, `batch`, `percent_mito`, `nperts`, …). The `_K562_essential` file loaded through pertpy has been reported to contain only 207,324 cells (pertpy issue #916), so prefer the Figshare originals.

### 1c. Perturbations / cells (as stated by sources)

| Screen | Target genes (paper) | Cells in h5ad | Genes in h5ad | Median cells/pert | Largest group (= non-targeting) |
|---|---:|---:|---:|---:|---:|
| K562 genome-wide (day 8, `K562_gwps`) | 9,866 | 1,989,578 | 8,248 | 178 | 75,328 |
| K562 essential (day 6, `K562_essential`) | 2,057 | 310,385 | 8,563 | 121 | 10,691 |
| RPE1 essential (day 7, `rpe1`) | 2,393 (Nadig Table 1; paper says "DepMap essential + hand-picked") | 247,914 | 8,749 | 72 | 11,485 |

Sources: gene counts from the Cell paper text and Nadig et al. Table 1; cell/gene counts from the scPerturb dataset table (`website/datavzrd/scperturb_dataset_info_datavzrd_annotated.csv`, which reads the Figshare h5ads) and a Purdue summary of the GWPS files; the 310,385 × 8,563 figure for K562 essential was independently confirmed by a user reading the Figshare file (scverse discourse #3933).

### 1d. Non-targeting controls
**Yes.** Libraries included non-targeting sgRNAs at 5 % of elements (Methods: "non-targeting control sgRNAs accounting for 5% of the total library" for K562 essential; same design for gwps and RPE1). In the h5ad, control cells are labelled `obs["gene"] == "non-targeting"` (confirmed by the scPerturb processing script, which sets `nperts = 0` where `gene == "non-targeting"`). Controls were spread across all gem groups ("many control cells bearing diverse non-targeting sgRNAs … internal z-normalization").

### 1e. Cell-level metadata (`.obs`) — fields documented by downstream processing code
From `sanderlab/scPerturb/dataset_processing/scripts/ReplogleWeissman2022.py` (reads the Figshare `*_raw_singlecell_01.h5ad` directly):
- `gene` — target gene symbol, or `non-targeting`
- `sgID_AB` — dual-sgRNA construct identity (guide-level label)
- `gem_group` — 10x lane / batch integer
- `mitopercent` — mitochondrial fraction
- Users additionally report `gene_id`, `transcript`, `UMI_count`, `core_scale`, `core_adjusted_UMI_count`, `z_gemgroup_UMI` (not verified here — confirm on first load).
`.var` fields: `gene_name` (index is Ensembl ID), `chr`, `start`, `end`, `class`, `strand`, `length`, `in_matrix`, `mean`, `std`, `cv`, `fano`.

---

## 2. Nadig et al. 2025 (Nature Genetics) — Jurkat essential and HepG2 essential

- Paper: Nadig A, Replogle JM, Pogson AN, …, Weissman JS, O'Connor LJ. "Transcriptome-wide analysis of differential expression in perturbation atlases." *Nat Genet* 2025. https://doi.org/10.1038/s41588-025-02169-3 (PMID 40259084). Preprint: bioRxiv 10.1101/2024.07.03.601903.
- Data availability statement (preprint): raw reads on SRA **BioProject PRJNA1100571**; aligned data and processed single-cell populations on **GEO GSE264667**.
- GEO series: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE264667 (HEAD 200; the HTML page is behind an NCBI reCAPTCHA for scripted GETs, but the FTP mirror is open). Title "Transcriptome-wide characterization of genetic perturbations"; last updated 2025-05-07; 224 samples (GSM8225020–GSM8225243) = 56 gem groups × {mRNA, sgRNA} libraries × 2 cell lines.
- **Licence**: none stated. GEO data are public-domain by NCBI policy (no restrictions on reuse; cite the paper). The preprint text is CC BY-NC-ND, which does not apply to the data.

### 2a. Files (GEO FTP, HEAD 200)

| File | Bytes | GB | URL |
|---|---:|---:|---|
| `GSE264667_jurkat_raw_singlecell_01.h5ad` | 9,366,490,264 | 9.37 | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE264nnn/GSE264667/suppl/GSE264667_jurkat_raw_singlecell_01.h5ad (206 ranged GET, HDF5 ✓) |
| `GSE264667_hepg2_raw_singlecell_01.h5ad` | 5,614,460,941 | 5.61 | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE264nnn/GSE264667/suppl/GSE264667_hepg2_raw_singlecell_01.h5ad |
| `GSE264667_RAW.tar` (per-gem-group cellranger MEX: barcodes/features/matrix ×112) | 16,321,884,160 | 16.3 | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE264nnn/GSE264667/suppl/GSE264667_RAW.tar |
| `filelist.txt` | 24 KB | — | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE264nnn/GSE264667/suppl/filelist.txt |

Format: AnnData `.h5ad`, raw UMI counts, same naming convention and pipeline as Replogle 2022 (`*_raw_singlecell_01` = genes >0.01 UMI/cell). GEO sample metadata: Cell Ranger v4.0.0, reference `refdata-gex-GRCh38-2020-A`, "AnnData format where the .var annotation details genes while the .obs annotation details single-cells"; `cell line: Jurkat|HepG2`; `genotype: day 7 post transduction`.

### 2b. Design / counts (as stated by sources)
- Library dJR092: 2,688 dual-sgRNA elements targeting 20Q1 DepMap common-essential genes + **5 % non-targeting control sgRNAs**; dual-sgRNA vector pJR101 (Addgene #187241). FACS-selected day 3, harvested day 7, 10x 3′ v3.
- Nadig Table 1: **Jurkat-Essential 2,393 perturbations, median 83 cells/perturbation; HepG2-Essential 2,393 perturbations, median 45 cells/perturbation.**
- Total cell numbers are not stated in the paper/GEO. Derived copies report: scPertEval (trimmed scPerturb copies) Jurkat 258,202 cells (12,013 control; 2,137 single perts), HepG2 133,757 cells (4,976 control; 1,818 single perts). Treat as lower bounds until the GEO h5ad is opened.

### 2c. Non-targeting controls
**Yes.** `obs["gene"] == "non-targeting"` in both GEO h5ads (confirmed by scPerturb's `NadigOConner2024.ipynb`, which reads `GSE264667_jurkat_raw_singlecell_01.h5ad` / `..._hepg2_...` and remaps `non-targeting` → `control`). TRADE itself used "control cells (i.e. cells carrying a non-targeting guide RNA)" per gem group.

### 2d. Cell-level metadata (`.obs`), from the same notebook
`gene` (target or `non-targeting`), `sgID_AB` (guide pair ID), `gem_group`, `mitopercent`; `.var`: `gene_name`, `chr`, `start`, `end`, `class`, `strand`, `length`, `in_matrix`, `mean`, `std`, `cv`, `fano` (index = Ensembl `gene_id`). One duplicated symbol (`HSPA14`, ENSG00000187522) in `.var`.

### 2e. Alternative copies
- scPerturb Zenodo 13350497 (CC BY 4.0): `NadigOConner2024_jurkat.h5ad` 1,293,665,804 B (1.29 GB; HEAD 200 w/ browser UA, 206 ranged GET, HDF5 ✓), `NadigOConner2024_hepg2.h5ad` 0.851 GB. Obs harmonised (`perturbation` with `control`, `guide_id`, `batch`, `percent_mito`, `ncounts`, `ngenes`, `percent_ribo`).
- scPertEval GCS (log-normalised, trimmed to `X` + `obs.perturbation`; not for count-based work): `https://storage.googleapis.com/scperteval/processed/nadig25jurkat_processed_complete.h5ad` 2.00 GB (HEAD 200), `nadig25hepg2_processed_complete.h5ad` 1.24 GB (HEAD 200); also `replogle22k562` / `replogle22rpe1` (essential panels only).

---

## 3. X-Atlas/Pisces (Xaira Therapeutics) — status check

- Hugging Face: https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Pisces (HEAD 200; `gated: false`, `private: false`).
- Described as 25.6 M perturbed cells across 7 CRISPRi screens (README says "16 biologically diverse contexts" — the X-Cell data page says 7; the number 16 presumably counts the multi-lineage differentiation sub-contexts). Screen table from https://xaira-therapeutics.github.io/X-Cell/data/ (HEAD 200):

| Screen | Cell type | Perturbations | Perturbed cells | Median KD |
|---|---|---:|---:|---:|
| HCT116 | colorectal | 18,924 | 3.4 M | 70 % |
| HEK293T | kidney | 18,312 | 4.5 M | 48 % |
| HepG2 | liver | 9,735 | 2.6 M | 85 % |
| iPSC | iPSC | 10,095 | 4.2 M | 82 % |
| **Jurkat Resting** | T-ALL | 10,872 | 2.8 M | 79 % |
| **Jurkat Active** | CD3/CD28-stimulated | 10,878 | 2.8 M | 71 % |
| iPSC Multi-Diff | multi-lineage diff. | 12,175 | 5.1 M | 96 % |

- Paper: "X-Cell: Scaling Causal Perturbation Prediction Across Diverse Cellular Contexts via Diffusion Language Models", bioRxiv 10.64898/2026.03.18.712807 (March 2026); PDF https://www.cdn.xaira.com/papers/X_CELL_V1_0316_final.pdf.
- **Licence: CC BY-NC-SA 4.0** (`LICENSE.md` in the repo and `cardData.license`). Non-commercial + share-alike — compatible with academic use, but derivative data must be released under the same licence, and any commercial partner cannot use it.
- Planned contents (README, "Coming Soon"): (1) all valid dual-guide cells from X-Atlas/Orion (HCT116, HEK293T); (2) HepG2 and iPSC cells restricted to a 200-gene validation set **plus non-targeting controls**; (3) **all valid dual-guide cells from Jurkat Resting and Jurkat Active**.
- Documented format (X-Cell data page): `.h5ad`; `.X` = log1p CP10k (**not raw counts**); `.obs["perturbation"]` = target gene; `.obs["is_control"]` = boolean non-targeting flag; `.var_names` = Ensembl IDs.
- **Actual repo contents as of 2026-10-11: `.gitattributes`, `LICENSE.md`, `README.md` only (24 KB total). Zero data files.** Commit history: initial commit, README update, LICENSE upload — all on 2026-03-25; nothing since. **Total download size: 0 B; the dataset does not yet exist for download.**
- Reference point: sibling dataset X-Atlas/Orion (https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Orion, HEAD 200) *is* populated: per-batch parquet files (HCT116/HEK293T), ~19 GB total listed by the HF tree API. Expect Pisces Jurkat (2 × 2.8 M cells) to be on the order of 10–20 GB per screen in similar encoding, and larger as h5ad.

---

## Summary table

| # | Screen | Source / accession | Primary file (h5ad, raw counts) | Size | Format | NT controls | Perturbations | Cells | Guide label | Licence | URL check |
|---|---|---|---|---|---:|---|---|---:|---:|---|---|---|
| 1 | K562 essential (day 6) | Replogle 2022; Figshare+ 20029387 | `K562_essential_raw_singlecell_01.h5ad` — https://ndownloader.figshare.com/files/35773219 | 10.66 GB | AnnData h5ad | Yes (`gene=="non-targeting"`, ~10.7 k cells) | 2,057 genes | 310,385 | `sgID_AB`, `gene`, `gem_group` | CC BY 4.0 | 302→206 ✓ |
| 2 | K562 genome-wide (day 8) | Replogle 2022; Figshare+ 20029387 | `K562_gwps_raw_singlecell_01.h5ad` — https://ndownloader.figshare.com/files/35775507 | 65.83 GB | AnnData h5ad | Yes (~75 k cells) | 9,866 genes | 1,989,578 | same | CC BY 4.0 | 302→206 ✓ |
| 3 | RPE1 essential (day 7) | Replogle 2022; Figshare+ 20029387 | `rpe1_raw_singlecell_01.h5ad` — https://ndownloader.figshare.com/files/35775606 | 8.70 GB | AnnData h5ad | Yes (~11.5 k cells) | 2,393 genes | 247,914 | same | CC BY 4.0 | 302→206 ✓ |
| 4 | Jurkat essential (day 7) | Nadig 2025; GEO GSE264667 | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE264nnn/GSE264667/suppl/GSE264667_jurkat_raw_singlecell_01.h5ad | 9.37 GB | AnnData h5ad | Yes (`gene=="non-targeting"`, 5 % of library) | 2,393 genes (median 83 cells/pert) | not stated (≥258 k) | same | GEO / unrestricted | HEAD 200 ✓ |
| 5 | HepG2 essential (day 7) | Nadig 2025; GEO GSE264667 | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE264nnn/GSE264667/suppl/GSE264667_hepg2_raw_singlecell_01.h5ad | 5.61 GB | AnnData h5ad | Yes | 2,393 genes (median 45 cells/pert) | not stated (≥134 k) | same | GEO / unrestricted | HEAD 200 ✓ |
| 6 | X-Atlas/Pisces Jurkat Resting + Active | Xaira; HF `Xaira-Therapeutics/X-Atlas-Pisces` | **none uploaded yet** (repo = README + LICENSE) | 0 B now; ≥10 s of GB expected | h5ad planned (log1p CP10k, not raw) | Planned (`obs.is_control`) | 10,872 / 10,878 | 2.8 M / 2.8 M (stated) | `obs.perturbation` (planned) | CC BY-NC-SA 4.0 | HEAD 200, no data ✓ |

Total for the five usable screens (raw single-cell h5ad only): **~100 GB** (≈34 GB if the K562 genome-wide screen is skipped or replaced by the 8.8 GB scPerturb copy).

## Blockers

1. **X-Atlas/Pisces is not downloadable.** The HF repo has had no data files since creation (2026-03-25); only README/LICENSE. The Jurkat resting/activated screens exist only as described in the X-Cell preprint. No release date given ("Coming Soon"). Action: watch the repo / contact Xaira; plan without it.
2. **Pisces licence is CC BY-NC-SA 4.0** — non-commercial and share-alike. Fine for an academic project, but derived released data must carry the same licence and it cannot be mixed into commercially licensed outputs. Also the planned `.X` is log1p CP10k, not raw counts; raw counts may need to be requested separately.
3. **Figshare download URLs are presigned S3 links valid for 10 s**, and the S3 object rejects HEAD (403). Use `curl -L` / `wget` directly on `https://ndownloader.figshare.com/files/<id>` in one shot; resumable download (`-C -`) will need a fresh redirect each time. Figshare article HTML pages return 202 (Cloudflare) to scripted requests — use the API `https://api.figshare.com/v2/articles/<id>` for metadata.
4. **Zenodo (scPerturb mirror) returns 403 to default-UA requests**; works with a browser-like `User-Agent`.
5. **K562 genome-wide raw h5ad is 65.8 GB** (dense storage; raw and normalised files are the same size). Needs ~70 GB disk + enough RAM to load in backed mode, or use the 8.8 GB gzip scPerturb copy (same cells, CC BY 4.0).
6. **Cell totals for Jurkat/HepG2 are not stated** by GEO or the paper; only perturbation counts and medians. Verify by opening the h5ads (scPerturb-derived counts: ~258 k and ~134 k).
7. **RPE1 target-gene count discrepancy**: Cell paper says DepMap essential + hand-picked genes; Nadig Table 1 and scPerturb give 2,393; Purdue summary gives 2,392 "knocked out genes" in the h5ad (likely excluding `non-targeting`). Not a blocker, just note when counting.
8. GEO HTML pages are behind NCBI reCAPTCHA for scripted access; the FTP mirror (`ftp.ncbi.nlm.nih.gov`) is the reliable download path.
