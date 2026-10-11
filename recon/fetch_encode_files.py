#!/usr/bin/env python3
"""Fetch per-experiment processed-file inventory from ENCODE for ATAC-seq experiments.
Usage: python3 fetch_encode_files.py out.json ENCSR... ENCSR...
"""
import json, sys, subprocess, time

BASE = "https://www.encodeproject.org"
WANT_OUT = {"signal p-value", "fold change over control", "IDR thresholded peaks",
            "pseudoreplicated peaks", "conservative IDR thresholded peaks",
            "pseudo-replicated peaks", "replicated peaks"}

def get(url):
    for i in range(3):
        try:
            p = subprocess.run(["curl", "-sL", "-H", "Accept: application/json", url],
                               capture_output=True, timeout=300)
            return json.loads(p.stdout)
        except Exception as e:
            print("retry", url, e, file=sys.stderr); time.sleep(3)
    raise RuntimeError(url)

def exp_summary(acc):
    e = get(f"{BASE}/experiments/{acc}/?format=json")
    bios = set(); treats = set(); gms = set(); ages = set()
    for r in e.get("replicates", []):
        lib = r.get("library") or {}
        bs = lib.get("biosample") or {}
        bios.add(bs.get("accession"))
        for t in bs.get("treatments", []) or []:
            treats.add(f"{t.get('treatment_term_name')} {t.get('amount')}{t.get('amount_units','')} {t.get('duration')}{t.get('duration_units','')}")
        for g in bs.get("applied_modifications", []) or bs.get("genetic_modifications", []) or []:
            gms.add(g.get("@id") if isinstance(g, dict) else g)
    default_an = e.get("default_analysis")
    analyses = []
    for a in e.get("analyses", []):
        analyses.append({"id": a.get("@id"), "accession": a.get("accession"), "assembly": a.get("assembly"),
                         "pipeline_version": a.get("pipeline_version"), "status": a.get("status"),
                         "pipeline_labs": [l.get("title") if isinstance(l, dict) else l for l in a.get("pipeline_labs", [])],
                         "pipelines": [p.get("title") if isinstance(p, dict) else p for p in a.get("pipelines", [])],
                         "files": [f.split("/")[-2] if isinstance(f, str) else f.get("accession") for f in a.get("files", [])]})
    files = []
    for f in e.get("files", []):
        if f.get("status") != "released": continue
        if f.get("assembly") != "GRCh38": continue
        if f.get("output_type") not in WANT_OUT: continue
        files.append({"accession": f["accession"], "output_type": f["output_type"],
                      "file_format": f.get("file_format"), "file_format_type": f.get("file_format_type"),
                      "biological_replicates": f.get("biological_replicates"),
                      "file_size": f.get("file_size"), "href": BASE + f.get("href", ""),
                      "preferred_default": f.get("preferred_default", False),
                      "analysis": (f.get("analyses") or [{}])[0].get("@id") if f.get("analyses") else None,
                      "md5sum": f.get("md5sum"),
                      "step": (f.get("analysis_step_version") or {}).get("@id") if isinstance(f.get("analysis_step_version"), dict) else f.get("analysis_step_version"),
                      "submitted_file_name": (f.get("submitted_file_name") or "").split("/")[-1]})
    return {"accession": acc, "status": e.get("status"), "description": e.get("description"),
            "lab": e.get("lab", {}).get("title"), "award_rfa": e.get("award", {}).get("rfa"),
            "biosample": e.get("biosample_ontology", {}).get("term_name"),
            "biosample_summary": e.get("biosample_summary"),
            "date_released": e.get("date_released"),
            "bio_reps": sorted({r["biological_replicate_number"] for r in e.get("replicates", [])}),
            "biosample_accessions": sorted(x for x in bios if x),
            "treatments": sorted(treats), "genetic_modifications": sorted(x for x in gms if x),
            "default_analysis": default_an, "analyses": analyses, "files": files,
            "audit_keys": {k: [a.get("category") for a in v] for k, v in (e.get("audit") or {}).items()}}

if __name__ == "__main__":
    out = {}
    for acc in sys.argv[2:]:
        out[acc] = exp_summary(acc)
        s = out[acc]
        print(acc, s["biosample"], s["lab"], "reps", s["bio_reps"], "treat", s["treatments"], "gm", s["genetic_modifications"], "default", s["default_analysis"], "nfiles", len(s["files"]))
    json.dump(out, open(sys.argv[1], "w"), indent=1)
