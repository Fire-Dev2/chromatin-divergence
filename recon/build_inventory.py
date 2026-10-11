#!/usr/bin/env python3
"""Build recon/02_atac_chip_accessions.md from fetched JSON/CSV."""
import json, csv, collections
D = json.load(open("encode_atac_candidates.json"))
F = json.load(open("fourdn_4DNES48BA2C8.json"))
BASE = "https://www.encodeproject.org"

PRIMARY = [("K562", "ENCSR868FGK"), ("HepG2", "ENCSR291GJU")]
PRIMARY_ALT = [("K562", "ENCSR483RKN"), ("HepG2", "ENCSR042AWH")]
DECOY = [("GM12878", "ENCSR637XSC"), ("HCT116", "ENCSR872WGW"), ("MCF-7", "ENCSR422SUG"),
         ("A549", "ENCSR032RGS"), ("IMR-90", "ENCSR200OML"), ("PC-3", "ENCSR499ASS"),
         ("Panc1", "ENCSR591PIX"), ("SK-N-SH", "ENCSR587TRP"), ("GM23338", "ENCSR485TLP"),
         ("NCI-H929", "ENCSR382LBS"), ("DND-41", "ENCSR660WSB"), ("RWPE2", "ENCSR080SNF"),
         ("WTC11", "ENCSR541KFY")]
DECOY_ALT = [("GM12878", "ENCSR095QNB")]

def gb(n): return f"{n/1e9:.2f} GB" if n >= 1e9 else f"{n/1e6:.1f} MB"
def url(acc, fmt): return f"{BASE}/files/{acc}/@@download/{acc}.{fmt}"

def pick(e):
    da = e["default_analysis"]; fs = [f for f in e["files"] if f["analysis"] == da]
    an = next(a for a in e["analyses"] if a["id"] == da)
    nrep = len(e["bio_reps"]); allr = e["bio_reps"]
    def one(cond):
        c = [f for f in fs if cond(f)]; return c
    sel = {}
    sel["pval_pooled"] = one(lambda f: f["output_type"] == "signal p-value" and f["biological_replicates"] == allr)
    sel["fc_pooled"] = one(lambda f: f["output_type"] == "fold change over control" and f["biological_replicates"] == allr)
    sel["pval_rep"] = one(lambda f: f["output_type"] == "signal p-value" and len(f["biological_replicates"]) == 1)
    sel["fc_rep"] = one(lambda f: f["output_type"] == "fold change over control" and len(f["biological_replicates"]) == 1)
    sel["idr_ppr"] = one(lambda f: f["output_type"] == "IDR thresholded peaks" and f["file_format"] == "bed" and "pooled-pseudoreplicated-idr" in (f["step"] or ""))
    sel["idr_cons"] = one(lambda f: f["output_type"] == "conservative IDR thresholded peaks" and f["file_format"] == "bed")
    sel["idr_truerep"] = one(lambda f: f["output_type"] == "IDR thresholded peaks" and f["file_format"] == "bed" and "replicated-idr" in (f["step"] or "") and "pseudo" not in (f["step"] or "") and f["biological_replicates"] == allr)
    sel["idr_rep"] = one(lambda f: f["output_type"] == "IDR thresholded peaks" and f["file_format"] == "bed" and len(f["biological_replicates"]) == 1)
    sel["pr_pooled"] = one(lambda f: f["output_type"] == "pseudoreplicated peaks" and f["file_format"] == "bed" and f["biological_replicates"] == allr)
    sel["pr_rep"] = one(lambda f: f["output_type"] == "pseudoreplicated peaks" and f["file_format"] == "bed" and len(f["biological_replicates"]) == 1)
    sel["overlap_pooled"] = one(lambda f: f["output_type"] == "replicated peaks" and f["file_format"] == "bed" and f["biological_replicates"] == allr)
    return an, sel

def fline(f, fmt):
    reps = "+".join(map(str, f["biological_replicates"]))
    nm = f["submitted_file_name"].replace(".gz", "")
    return f"| {f['accession']} | {f['output_type']} | rep{reps} | {gb(f['file_size'])} | `{nm}` | {url(f['accession'], fmt)} |"

out = []
W = out.append
W("# 02 — ATAC-seq and TF ChIP-seq accession inventory (ENCODE / 4DN)\n")
W("Generated 2026-10-11 from live ENCODE REST API (`www.encodeproject.org`) and 4DN data portal API (`data.4dnucleome.org`) metadata; no data files downloaded. Scripts: `fetch_encode_atac.py`, `fetch_encode_files.py`, `fetch_encode_tf_chip.py`, `build_inventory.py` (this directory). Raw JSON: `encode_atac_candidates.json`, `atac_probe.json`, `fourdn_atac_sets.json`, `fourdn_4DNES48BA2C8.json`, `fourdn_file_*.json`.\n")
W("All ENCODE experiments below are ENCODE4 Snyder-lab bulk ATAC-seq, processed by the ENCODE uniform ATAC-seq pipeline (`ENCPL344QWT` + `ENCPL787FUN`, lab `encode-processing-pipeline`) on GRCh38; the *released* default analysis is pipeline **v1.9.1** for all four target/decoy core experiments and most extras; exceptions: NCI-H929 ENCSR382LBS = **v1.8.0**, WTC11 ENCSR541KFY = **v1.10.0** (older v1.8.0 analyses on the others are `revoked`/`deleted`). File selection rule used throughout: files from the experiment's `default_analysis` only; bigWigs = MACS2 `signal p-value` and `fold change over control` (pooled and per replicate); peaks = `IDR thresholded peaks` from the pooled-pseudoreplicate IDR step (`pooled-pr1_vs_pooled-pr2.idr0.05.bfilt.narrowPeak`; this is the ENCODE ATAC 'IDR' set) plus the `pseudoreplicated peaks` pooled overlap set (`preferred_default=True` on the portal) and, where present, the `conservative IDR thresholded peaks` (true-replicate IDR).\n")

# ---------------- summary table
W("## Summary table\n")
W("| Role | Cell line | Source | Experiment | Lab | Bio reps | Assembly | Pipeline | Pooled signal p-value bigWig | Pooled fold-change bigWig | IDR peaks (pooled-pr) narrowPeak | Pseudoreplicated peaks narrowPeak | Notes |")
W("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
def row(role, cell, acc, note=""):
    e = D[acc]; an, s = pick(e)
    g = lambda k: (s[k][0]["accession"] + f" ({gb(s[k][0]['file_size'])})") if s[k] else "—"
    W(f"| {role} | {cell} | ENCODE | {acc} | {e['lab']} | {len(e['bio_reps'])} | GRCh38 | ENCODE ATAC v{an['pipeline_version']} ({an['accession']}) | {g('pval_pooled')} | {g('fc_pooled')} | {g('idr_ppr')} | {g('pr_pooled')} | {note} |")
row("target", "K562", "ENCSR868FGK")
W("| target | hTERT-RPE1 | 4DN | 4DNES48BA2C8 (exps 4DNEX7A3MJJV, 4DNEX5R8DYKC) | Bing Ren, UCSD | 2 | GRCh38 | 4DN run of ENCODE ATAC pipeline (`encode-atacseq-aln` 1.1.1 + `encode-atacseq-postaln` 1.1.1; 4DNWF2BBZ1EE/4DNWF2VTG8II) | — (not produced) | 4DNFIS4IDZNE (4.41 GB) | 4DNFI4H1MT8U 'conservative peaks' bigBed (15.2 MB) | 4DNFIMNHSGU1 'peaks' bigBed (15.2 MB) | Not ENCODE-uniform; pipeline version differs (4DN CWL port of ENCODE atac-seq-pipeline v1.1.x vs ENCODE v1.9.1); peaks are bigBed not bed; no p-value track; **contact Ren lab before publication** (see policy below) |")
W("| target | Jurkat (clone E6-1) | — | **none** | — | — | — | — | — | — | — | — | No ATAC-seq on ENCODE (any assembly/status) or 4DN; see blockers |")
row("target", "HepG2", "ENCSR291GJU")
for c, a in DECOY:
    note = ""
    if a == "ENCSR587TRP": note = "**treated: all-trans-retinoic acid 6 µM 48 h** (only SK-N-SH ATAC on ENCODE); use with caution or drop"
    if a == "ENCSR872WGW": note = "only untreated/unmodified HCT116 ATAC (16 others are degron-tagged ± auxin/5-Ph-IAA)"
    if a == "ENCSR485TLP": note = "iPSC line (male, from GM23248 fibroblasts)"
    if a == "ENCSR660WSB": note = "T-ALL line — same lineage as Jurkat; candidate Jurkat stand-in/control"
    if a == "ENCSR541KFY": note = "iPSC; a second WTC11 ATAC (ENCSR089KIJ, Yin Shen) exists"
    row("decoy", c, a, note)
W("| alt | K562 | ENCODE | ENCSR483RKN | Michael Snyder, Stanford | 2 | GRCh38 | ENCODE ATAC v1.9.1 | see detail | | | | second untreated Snyder K562 ATAC; 62 other K562 ATAC are Greenleaf-lab drug/perturbation series (e.g. ENCSR165JXS = JQ1 5 µM 4 h) |")
W("| alt | HepG2 | ENCODE | ENCSR042AWH | Michael Snyder, Stanford | 2 | GRCh38 | ENCODE ATAC v1.9.1 | see detail | | | | |")
W("| alt | GM12878 | ENCODE | ENCSR095QNB | Michael Snyder, Stanford | 2 | GRCh38 | ENCODE ATAC v1.9.1 | see detail | | | | |")
W("")
W("Cell lines requested as decoys but **absent** from ENCODE ATAC-seq (any status, any assembly, checked by `biosample_ontology.term_name` and free-text `searchTerm`): HeLa-S3, HEK293/HEK293T. 4DN has HEK293 ATAC only as CRISPR-perturbation sets from Taekjip Ha (JHU) (untreated control = 4DNES694VKB2, 2 exps); not ENCODE-uniform.\n")

# ---------------- TF ChIP
K = list(csv.DictReader(open("encode_tf_chip_K562.csv"))); H = list(csv.DictReader(open("encode_tf_chip_HepG2.csv")))
def stats(rows):
    tf = [r for r in rows if r["is_sequence_specific_TF"] == "True"]
    return len(rows), len({r["target"] for r in rows}), len(tf), len({r["target"] for r in tf})
kN, kT, kTFn, kTF = stats(K); hN, hT, hTFn, hTF = stats(H)
W("## TF ChIP-seq availability (ENCODE, released, GRCh38, IDR thresholded narrowPeak present)\n")
W("Query: `type=Experiment&assay_title=TF ChIP-seq&status=released&biosample_ontology.term_name=<cell>&assembly=GRCh38&files.output_type=IDR thresholded peaks&files.file_type=bed narrowPeak`, then each experiment re-checked to have ≥1 *released* GRCh38 `IDR thresholded peaks` bed narrowPeak. 'Sequence-specific TF' = ENCODE target `investigated_as` contains `transcription factor` (excludes cofactors, chromatin remodelers, RBPs, Pol II, cohesin, DNA repair/replication).\n")
W("| Cell line | Experiments (all TF-ChIP targets) | Distinct targets (all) | Experiments with TF targets | **Distinct sequence-specific TF targets** | CSV |")
W("|---|---|---|---|---|---|")
W(f"| K562 | {kN} | {kT} | {kTFn} | **{kTF}** | `recon/encode_tf_chip_K562.csv` |")
W(f"| HepG2 | {hN} | {hT} | {hTFn} | **{hTF}** | `recon/encode_tf_chip_HepG2.csv` |")
kset = {r["target"] for r in K if r["is_sequence_specific_TF"] == "True"}; hset = {r["target"] for r in H if r["is_sequence_specific_TF"] == "True"}
W(f"\nTF targets shared by K562 and HepG2: **{len(kset & hset)}**; K562-only {len(kset - hset)}; HepG2-only {len(hset - kset)}. CSV columns: target, investigated_as, is_sequence_specific_TF, experiment, lab, award_rfa, date_released, n_bio_reps, chosen IDR peak file accession/replicates/size/URL, number of GRCh38 IDR files, audit ERROR / NOT_COMPLIANT categories (e.g. 'extremely low read depth' — filter on these before use). Where a target has several experiments, every experiment is a row; dedupe on `target` for counts.\n")

# ---------------- detail
W("## Per-experiment file detail (ENCODE)\n")
for role, lst in [("Targets", PRIMARY), ("Decoys", DECOY), ("Alternates", PRIMARY_ALT + DECOY_ALT)]:
    W(f"### {role}\n")
    for cell, acc in lst:
        e = D[acc]; an, s = pick(e)
        W(f"#### {cell} — {acc}\n")
        W(f"- Portal: {BASE}/experiments/{acc}/ · lab {e['lab']} · {e['award_rfa']} · released {e['date_released']} · biosample summary: {e['biosample_summary']} · biosamples {', '.join(e['biosample_accessions'])} · bio reps {e['bio_reps']}")
        W(f"- Treatments: {e['treatments'] or 'none'} · genetic modifications: {e['genetic_modifications'] or 'none'}")
        W(f"- Default analysis {an['accession']} ({an['id']}): GRCh38, ENCODE ATAC-seq pipeline v{an['pipeline_version']}, status {an['status']}; other analyses: " + "; ".join(f"{a['accession']} v{a['pipeline_version']} {a['status']}" for a in e['analyses'] if a['id'] != an['id']))
        aud = e["audit_keys"]; W(f"- Audits: " + "; ".join(f"{k}: {', '.join(sorted(set(v)))}" for k, v in aud.items() if k != 'INTERNAL_ACTION') )
        W("")
        W("| File | Output type | Replicates | Size | Pipeline file name | URL |")
        W("|---|---|---|---|---|---|")
        for k in ["pval_pooled", "fc_pooled", "pval_rep", "fc_rep"]:
            for f in sorted(s[k], key=lambda f: f["biological_replicates"]): W(fline(f, "bigWig"))
        for k in ["idr_ppr", "idr_cons", "pr_pooled", "overlap_pooled", "idr_rep", "pr_rep"]:
            for f in sorted(s[k], key=lambda f: f["biological_replicates"]): W(fline(f, "bed.gz"))
        if s["idr_truerep"]:
            W(f"\nAdditional true-replicate pairwise IDR files (bed narrowPeak, all labelled `IDR thresholded peaks`, reps {'+'.join(map(str,e['bio_reps']))}): " + ", ".join(f"{f['accession']} (`{f['submitted_file_name'].replace('.gz','')}`, {gb(f['file_size'])})" for f in s["idr_truerep"]))
        W("")

# ---------------- 4DN
W("## hTERT-RPE1 — 4DN 4DNES48BA2C8 detail\n")
lab = F["lab"]["display_title"]
W(f"- Portal: https://data.4dnucleome.org/experiment-set-replicates/4DNES48BA2C8/ · lab {lab} · award {F['award']['@id']} ({F['award']['display_title']}) · status {F['status']} · public release {F['public_release']} · dataset '{F['dataset_label']}' / condition '{F['condition']}' · description: {F['description']}")
W(f"- Publications: produced_in_pub = {F.get('produced_in_pub')}; publications_of_set = {F.get('publications_of_set')} → **unpublished dataset**.")
W("- Sibling sets in the same dataset (hTERT-RPE1 clones, same lab/pipeline): 4DNES7YK6HTM (F6 clone), 4DNESJWVEOAK (G1 clone), 4DNES1TSS4EO (D4 clone). Use 4DNES48BA2C8 (WT) only.")
W("- Replicates: 2 biological replicate experiments, each with 2 paired-end fastq (~7.3–7.5 GB each):")
for ex in F["experiments_in_set"]:
    raws = ", ".join(f"{f['accession']} ({gb(f['file_size'])})" for f in ex.get("files", []))
    procs = ", ".join(f"{f['accession']} {f.get('file_type')} ({gb(f['file_size'])})" for f in ex.get("processed_files", []))
    W(f"  - {ex['accession']}: raw {raws}; per-replicate processed: {procs}")
W("- Set-level processed files (GRCh38, workflow `encode-atacseq-postaln 1.1.1` = 4DN CWL port of the ENCODE ATAC-seq pipeline, run 2025-06-04; both replicates pooled):\n")
W("| File | Type | Size | Assembly | URL (portal) | Open-data S3 URL | md5 |")
W("|---|---|---|---|---|---|---|")
for acc in ["4DNFIS4IDZNE", "4DNFI4H1MT8U", "4DNFIMNHSGU1", "4DNFIKJGW79V", "4DNFID2TXF3X", "4DNFIM6LAD9J", "4DNFIIFNCJN9"]:
    try: d = json.load(open(f"fourdn_file_{acc}.json"))
    except FileNotFoundError: continue
    W(f"| {acc} | {d.get('file_type') or 'supplementary bigWig (no file_type; likely per-replicate signal)'} — {d.get('description') or ''} | {gb(d['file_size'])} | {d.get('genome_assembly')} | https://data.4dnucleome.org{d['href']} | {d.get('open_data_url')} | {d.get('md5sum')} |")
W("")
W("Notes: 4DN emits peaks as **bigBed** (`peaks` = pooled-pseudoreplicate overlap/naive set; `conservative peaks` = IDR conservative set per the 4DN description 'Conservative peak calls from ENCODE ATAC-Seq Pipeline'); convert with `bigBedToBed`. 4DN provides a pooled **fold-change** bigWig only (no signal p-value track). The ENCODE pipeline version run by 4DN (1.1.1 CWL port, ~ENCODE atac-seq-pipeline v1.1) predates ENCODE's v1.9.1 default analyses; alignment (bowtie2, GRCh38), filtering, MACS2 and IDR logic are the same lineage but parameter/version drift (e.g. MACS2 p-value cap, TSS enrichment, blacklist version) is possible. To be strictly uniform, re-run the ENCODE pipeline from fastq for all four target lines (ATAC fastq for the ENCODE sets are listed on the portal; RPE1 fastq = 4 × ~7.4 GB).\n")
W("### 4DN data-use policy for 4DNES48BA2C8\n")
W("Per the 4DN Data Use Guidelines (https://data.4dnucleome.org/resources/data-use-guidelines, see also https://www.4dnucleome.org/outreach/data-use-policy/): all released 4DN data are freely usable **unless** the dataset has not yet been described in a publication by the producing lab (`produced_in_pub` empty → unpublished, as here, public release 2025-06-12). For such datasets, users wishing to publish analyses must **contact the data-producing lab (Bing Ren, UCSD — contact listed on the lab page https://data.4dnucleome.org/labs/bing-ren-lab/) before submitting the manuscript**, give the lab the opportunity to comment/participate, cite the dataset accession (4DNES48BA2C8) and acknowledge the 4DN Network and the award (1UM1HG011585-01). (Policy text should be re-verified on the portal when the manuscript is drafted; this is a metadata-derived summary.)\n")

# ---------------- blockers
W("## Blockers and caveats\n")
W(f"""1. **Jurkat (clone E6-1): no ATAC-seq anywhere on ENCODE or 4DN.** ENCODE search for `biosample_ontology.term_name` = 'Jurkat, Clone E6-1' / 'Jurkat clone E6-1' and free-text 'Jurkat' with `assay_title=ATAC-seq` returns 0 experiments (any status, any assembly); 4DN has no Jurkat ATAC-seq ExperimentSet. ENCODE's Jurkat data are DNase-seq / ChIP-seq / RNA-seq only. Options: (a) download public Jurkat ATAC fastq from GEO/SRA (e.g. Corces et al. 2016 *Nat Genet* GSE74912 includes Jurkat; Buenrostro 2013 Omni-ATAC data) and run the ENCODE ATAC-seq pipeline (v1.9.1 or later Caper/WDL release) locally; (b) substitute the ENCODE T-ALL line **DND-41 (ENCSR660WSB, 3 reps, v1.9.1)** as a lineage-matched surrogate; (c) use ENCODE Jurkat DNase-seq as a non-uniform proxy. Any choice must be decided before the experiment design is locked.
2. **hTERT-RPE1: ENCODE has no RPE1 ATAC-seq** (0 hits for 'hTERT-RPE1'/'RPE1'/searchTerm 'RPE'). Only 4DN 4DNES48BA2C8 (WT) exists, processed with 4DN's CWL port of the ENCODE pipeline (1.1.1) — not the same release as ENCODE's v1.9.1 default analyses — with bigBed peaks and a fold-change-only bigWig. Strict pipeline uniformity therefore requires re-processing RPE1 fastq (≈30 GB raw) with the ENCODE pipeline, or accepting a documented pipeline-version mismatch for one of four target lines. Data-use: unpublished → contact Ren lab before publishing.
3. **SK-N-SH decoy is RA-treated** (ENCSR587TRP: all-trans-retinoic acid 6 µM 48 h); it is the only SK-N-SH ATAC on ENCODE. Treat as 'SK-N-SH + RA' or drop. HCT116 has 17 ATAC experiments, 16 of which are degron-tagged/5-Ph-IAA-treated; only ENCSR872WGW is plain HCT116.
4. **HeLa-S3 and HEK293 are not available as ENCODE ATAC-seq**; the decoy panel above substitutes GM23338, NCI-H929, DND-41, RWPE2, WTC11 (all Snyder ENCODE4, v1.9.1). Core decoy set of 7 untreated lines: GM12878, HCT116, MCF-7, A549, IMR-90, PC-3, Panc1 (+ optional SK-N-SH-RA, GM23338, NCI-H929, DND-41, RWPE2, WTC11).
5. **K562 has 64 ATAC experiments**; 62 are Greenleaf-lab perturbation series (drug treatments such as JQ1; ENCSR165JXS checked) and must not be used as untreated references. Untreated: ENCSR868FGK (3 reps, used) and ENCSR483RKN (2 reps).
6. ENCODE 'IDR thresholded peaks' for ATAC comes in several flavours per experiment (per-replicate self-pseudoreplicate IDR, pairwise true-replicate IDR, pooled-pseudoreplicate IDR); only `output_type` + `analysis_step_version` distinguish them. This inventory pins the pooled-pseudoreplicate IDR file (and the `conservative IDR` true-rep file) explicitly — do not pick by `output_type` alone.
7. Audit warnings: several experiments carry 'moderate library complexity' / 'mild to moderate bottlenecking' WARNINGs (listed per experiment above); none carry ERROR-level audits. In the TF ChIP CSVs, {sum(1 for r in K if r['audit_error'])} K562 and {sum(1 for r in H if r['audit_error'])} HepG2 experiments carry an ERROR audit (mostly 'extremely low read depth') — filter before use.
8. **Pipeline-version drift within ENCODE**: NCI-H929 (v1.8.0) and WTC11 (v1.10.0) default analyses differ from the v1.9.1 used by the other 16 experiments; prefer the 7-line core decoy set (all v1.9.1) or re-process.
9. Totals/sizes: pooled signal p-value + fold-change bigWigs are 0.7–4.5 GB each; a full download of the pooled tracks for the 15 ENCODE experiments + RPE1 is on the order of 60–70 GB — plan storage accordingly (sizes are per-file in the tables).
""")
open("02_atac_chip_accessions.md", "w").write("\n".join(out))
print("wrote 02_atac_chip_accessions.md", len(out), "lines")
