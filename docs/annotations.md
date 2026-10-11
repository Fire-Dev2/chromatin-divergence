# Annotation sources and transformations

Built by `scripts/build_annotations.py` (library: `src/cdtest/annot.py`;
tests: `tests/test_annot.py`). All downloads live under `data/annot/`
(git-ignored); md5, size and date of every file are in
`data/annot/MANIFEST_md5.tsv`. All downloads retrieved 2026-10-11.
Spec: ANALYSIS_PLAN.md §1 (gene universe), §2.6 (knockdown groups), §3
(blacklist, CN correction, TSS), §4 (CN and CpG covariates).

## Coordinate conventions

* Interval columns are **0-based, half-open** (BED) unless the column name
  ends in `1` (`tss1`, `gene_start1`, `gene_end1`: 1-based, GTF style).
* `tss1` = GTF `start` of the Ensembl-canonical transcript on `+`, GTF `end`
  on `-`; `tss0 = tss1 - 1` (the same base, 0-based). A BED record for the
  TSS base is `chrom, tss0, tss0+1`.
* Promoter window for CpG/GC: `[tss0-500, tss0+500)` (1,000 bp, strand-agnostic).
* CN is **log2 relative copy number** (0 = the line's genome-wide baseline).

## Sources

| Item | Source / version | URL | Licence |
|---|---|---|---|
| Genes, TSS | GENCODE Human Release 50 (Ensembl 116, GRCh38.p14), **comprehensive** GTF `gencode.v50.annotation.gtf.gz` (md5 2d273848c6068682fabea72fae4896a9) | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/ | GENCODE/Ensembl: no restrictions |
| Cytobands / arms | UCSC hg38 `cytoBand.txt.gz` (md5 6677dd50caf19b0c3626b99ce67ca073) | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz | UCSC: free use |
| CpG islands | UCSC hg38 `cpgIslandExt.txt.gz` | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz | UCSC: free use |
| Genome sequence (CpG/GC only) | UCSC `hg38.fa.gz` (md5 1c9dcaddfa41027f17cd8f7a82c7293b, verified against UCSC md5sum.txt); **primary chromosomes chr1-22,X only were read**. File deleted after use (storage rule); fetch script `data/annot/genome/fetch_parallel.sh` kept | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz | UCSC/GRC: free use |
| Blacklist | ENCODE blacklist v2, hg38 (Amemiya, Kundaje & Boyle, *Sci Rep* 2019;9:9354), `hg38-blacklist.v2.bed.gz`, 636 regions (md5 83fe6bf8187a64dee8079b80f75ba289) | https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz | Boyle-Lab repo; free reuse with citation |
| Copy number | DepMap Public 24Q4 `OmicsCNSegmentsProfile.csv`, `Model.csv`, `OmicsProfiles.csv`, `OmicsDefaultModelProfiles.csv` (Figshare+ 10.25452/figshare.plus.27993248.v1) | see recon/04 | CC BY 4.0 |
| RPE1 karyotype | Volpe M et al. *Nat Commun* 2025, doi:10.1038/s41467-025-62428-z: 46,X,der(X)t(X;10)(Xq28;10q21.2) → 10q21.2-qter at 3 copies | — | CC BY |
| p53 targets | Fischer M. *Oncogene* 2017;36:3943, Table 1 (116 genes in ≥6/16 datasets) | https://www.nature.com/articles/onc2016502/tables/1 | gene list = facts; cite paper |
| Gene groups | HGNC gene groups `family.csv`, `gene_has_family.csv`, `hierarchy.csv`; `hgnc_complete_set.txt` (symbol history) | https://storage.googleapis.com/public-download-files/hgnc/ | CC0 |
| Tie-break | CORUM 5.3 human complexes (2026-04-14) | https://mips.helmholtz-muenchen.de/corum/ | CC BY-NC 4.0 |

Future FASTA use: the project's canonical genome is
`data/genome/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz` (ENCODE);
`annot.FASTA_PATH` points there. Its chr1-22,X sequence equals UCSC hg38.
pyfaidx needs it uncompressed or bgzipped.

## Outputs

### `data/annot/genes_gencode_v50.tsv` — 54,710 genes (20,033 protein_coding, 34,677 lncRNA)

* Genes with `gene_type` ∈ {protein_coding, lncRNA} on chr1-22, chrX; one row
  per gene, taken from the transcript tagged `Ensembl_canonical` (exactly one
  per gene; asserted). The **comprehensive** GTF is used because the basic GTF
  omits the canonical transcript of 5,973 such genes (mostly lncRNA).
  `gene_id` = Ensembl id without version; `gene_id_version` kept.
* `arm` = first letter of the UCSC cytoband containing `tss0`; `cytoband` = band.
* `cpg_oe_tss500` = n_CpG · L / (n_C · n_G) and `gc_frac_tss500` = (n_C+n_G)/L
  over the 1,000-bp window, L = non-N bases, sequence upper-cased
  (Gardiner-Garden & Frommer 1987; Saxonov 2006 HCP threshold o/e > 0.75 &
  GC > 0.55). `n_cpg_tss500`, `n_N_tss500` also given (10 windows contain Ns).
  Medians: protein_coding o/e 0.72, GC 0.61; lncRNA o/e 0.29, GC 0.48.
* `cgi_overlap_tss500`: window overlaps a UCSC CpG island (39.5 % of genes).
* `blacklist_overlap_tss1kb`: TSS ±1 kb overlaps ENCODE blacklist v2 (1,234 genes).
* `p53_fischer2017_core` / `p53_fischer2017_ndatasets`: 116/116 Table 1 symbols
  mapped (113 by current gene_name; 3 via HGNC previous symbol: FAM198B→GASK1B,
  FAM212B→INKA2, WDR63→DNAI3). Mapping table: `fischer2017_p53_core_mapped.tsv`.
* `knockdown_group`: as in `knockdown_groups.tsv`.

### `data/annot/knockdown_groups.tsv` (gene_name, ensembl_id, group, hgnc_groups_all, tiebreak)

All genes of the gene table; group from HGNC gene groups (family ids
including descendants via `hierarchy.csv`), joined by HGNC id from the GTF:

| group | HGNC families | n |
|---|---|---|
| cytosolic_ribosome | 728 S ribosomal proteins, 729 L ribosomal proteins (incl. paralogues RPL3L, RPL10L, RPL39L, RPS27L, …) | 87 |
| mitochondrial_ribosome | 646 (+1994, 1996) | 82 |
| proteasome | 690 | 45 |
| spliceosome | 1518 Major spliceosome, 1520 Minor spliceosome (+ all sub-complexes: Sm, LSm, SF3a/b, U1/U2/U4-U6/U5, NTC, EJC, A/B/Bact/C/P/E complexes) | 152 |
| rna_pol_ii | 726 RNA polymerase subunits, restricted to `POLR2*` (incl. shared Pol I/II/III subunits POLR2E/F/H/K/L) | 15 |
| mediator | 1061 (incl. CDK8 module: CDK8, CDK19, CCNC, MED12/12L/13/13L) | 33 |
| nuclear_pore | 1051 Nucleoporins | 30 |
| other | everything else | 54,266 |

Tie-break (gene in >1 group): count CORUM 5.3 complexes containing the gene
whose name matches the group pattern (`annot.CORUM_GROUP_REGEX`); highest
count wins, remaining ties by the table order above. With these definitions
**no gene falls in two groups**, so the tie-break was never invoked.

### `data/annot/gene_cn_log2.tsv` (gene_id × 12 lines) and helpers

* **DepMap convention.** `SegmentMean` in 24Q4 `OmicsCNSegmentsProfile.csv`
  is a *linear relative copy ratio* (README: "relative copy ratio for that
  segment"; all values > 0, median 1.01, range 4e-11–1064). We store
  `log2(SegmentMean)`, floored at −5 (homozygous deletions give ≈1e-11).
  Segment `Start`/`End` are treated as 1-based closed → `start-1, end`.
  Values are relative to each line's ploidy: they say nothing about absolute
  copy number (K562 near-triploid, HepG2 hyperdiploid, Jurkat hypotetraploid).
* Profiles (the DepMap default DNA profile of the model, matched by RRID;
  `cn_profile_map.tsv`): K562 ACH-000551 PR-aheaZL (WGS); HepG2 ACH-000739
  PR-ozp1kx (WGS); Jurkat ACH-000995 PR-Ds9QdK (WES; parental JURKAT, not
  E6-1); HCT116 PR-YmJ7mI, MCF7 PR-VLTXhP, A549 PR-oKCH5w, PC3 PR-cO6ZzX,
  PANC1 PR-Tm9gPE (all WGS).
* **GM12878 (CVCL_7526) and IMR-90 (CVCL_0347) are not in DepMap 24Q4: CN set
  to 0** (flag: `status = absent_from_DepMap_24Q4_CN0` in `cn_profile_map.tsv`).
* **RPE1** (karyotype track, primary per plan §3): log2(3/2) = 0.585 on chr10
  from the start of UCSC band 10q21.2 (59,400,000, 0-based) to qter; 0
  elsewhere. Xq28 loss from the der(X) is *not* modelled (plan specifies 0).
  The empirical DepMap clone breakpoint is at ≈60.8 Mb (inside 10q21.2).
* **RPE1_depmap**: median log2 CN of the three WGS single-cell clones
  RPE1-ss6 (PR-g5QF67), -ss77 (PR-Yf6hXD), -ss119 (PR-l4NGGv) on the union of
  their breakpoints (positions covered by ≥2 clones); the median removes
  clone-private aneuploidies.
* Gene value = CN at the TSS base. TSSs falling in gaps between segments
  (centromeres, segdup-masked regions, WES gaps) take the nearest segment on
  the same chromosome; which genes were filled is in
  `gene_cn_log2_gapfilled.tsv` (1 = filled; e.g. K562 273, Jurkat 517 genes;
  BCR's TSS is one of them).
* chrX: male lines (HepG2, Jurkat, HCT116) show log2 ≈ −1 on X relative to
  autosomes, which is correct for relative CN.
* `cn_segments_log2.tsv`: long table (line, chrom, start, end, log2cn) used by
  `annot.CNLookup` / `annot.region_cn(chrom, start, end, lines)`, which return
  the bp-weighted mean log2 CN over arbitrary 0-based intervals per line, with
  the same nearest-segment gap filling (`fill=None` to disable). For ATAC
  region CN correction (plan §3).

#### CN sanity checks (asserted in tests)

* K562: CDKN2A TSS log2 = −5 (floor; homozygous 9p21 deletion, whole
  chr9:20.96–23.24 Mb at floor); 9p median −0.48; ABL1→NUP214 (9q34) +2.2 to
  +2.8 and 22q11.21 +1.6 to +2.3 (BCR-ABL1 amplicon).
* RPE1_depmap: 10q21.2-qter median 0.61 vs proximal chr10 0.02 (karyotype
  value 0.585); chr1 median 0.01. Note RPE1_depmap also shows CDKN2A at −0.92
  (loss shared by the clones), which the karyotype track does not include.
* CDKN2A floor (−5) also in MCF7, A549, PANC1 and Jurkat (−3.98) — all known
  CDKN2A-deleted lines.
* Genome-wide median |log2 CN| < 0.3 for every DepMap line.

## Reproduce

```
unset PYTHONPATH; source env/venv/bin/activate
bash data/annot/genome/fetch_parallel.sh        # only if CpG/GC must be recomputed
python scripts/build_annotations.py --fasta data/annot/genome/hg38.fa
python -m pytest tests/test_annot.py
```
