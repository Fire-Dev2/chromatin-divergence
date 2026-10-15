# cdtest deposited intermediates

Processed inputs for every confirmatory analysis in the cdtest paper (cdtest 0.2.0), restricted to the genes used by at least one pair. With these files and the per-pair outcome tables in `results/aim1/<pair>/`, `python scripts/reproduce_from_intermediates.py` recomputes Table 2 and Figures 3–4 without the raw screens or the ATAC-seq pipeline.

| path | content |
|---|---|
| data/annot/genes_gencode_v50.tsv | GENCODE v50 gene annotation (chromosome, TSS, arm, promoter CpG o/e) |
| data/annot/gene_cn_log2.tsv | relative log2 copy number at the TSS per line (DepMap 24Q4; RPE1 karyotype) |
| data/features/gene_features.parquet | gene-level accessibility features, insertion CPM (primary), 4 screened + 7 decoy lines |
| data/features/gene_features_fc.parquet | the same from fold-change tracks (sensitivity) |
| data/features/build_info*.json | provenance of the feature build (inputs, md5 sums, parameters) |
| results/aim1/control_moments/<line>.tsv.gz | mean and variance (ddof 1) of normalized expression over kept control cells (covariates) |

The outcome vectors M, knockdown lists and reliabilities are in `results/aim1/<pair>/` (gene_outcome.tsv, knockdowns.tsv, summary.json) in the repository and in the tarball. Install into a clone with `python -c "from cdtest import intermediates; intermediates.install()"` (copies only missing files).
