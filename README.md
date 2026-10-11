# chromatin-divergence

Everything for this project lives in this one folder. Nothing large is stored elsewhere.

| Folder | Contents | In git? |
|---|---|---|
| `data/perturbseq/` | Perturb-seq screens (h5ad, ~34 GB; +66 GB K562 genome-wide later) | no |
| `data/atac/` | ENCODE/4DN ATAC signal + peaks, 11 lines (~45 GB) | no |
| `data/atac_fastq/` | Raw ATAC reads for Jurkat and RPE1 (~37 GB), deletable after processing | no |
| `data/genome/` | GRCh38 reference + bowtie2 index (~5 GB) | no |
| `data/annot/` | Gene annotations, copy number, p53 targets, gene groups | no |
| `env/` | Python environment | no |
| `src/`, `tests/` | Analysis code and tests | yes |
| `results/`, `figures/` | Outputs | yes (small files) |
| `recon/` | Data inventories, literature review, novelty check | yes |
| `docs/` | Method notes | yes |
| `log/RESEARCH_LOG.md` | Dated research log, including AI use | yes |
| `ANALYSIS_PLAN.md` | Frozen analysis plan (tag plan-v1.0) | yes |
| `outreach/` | Email drafts | yes |

`./status.sh` shows download and job progress.

To free space later: delete `data/atac_fastq/` once ATAC processing is done.
Everything under `data/` can be re-downloaded from the manifests in `recon/`.
