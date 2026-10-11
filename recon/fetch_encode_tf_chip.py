#!/usr/bin/env python3
"""Count distinct sequence-specific TF ChIP-seq targets in a cell line on ENCODE
with released GRCh38 IDR thresholded peaks. Writes CSV of targets.
Usage: python3 fetch_encode_tf_chip.py K562 out.csv
"""
import json, sys, subprocess, csv, time, urllib.parse, collections

BASE = "https://www.encodeproject.org"

def get(url):
    for i in range(3):
        try:
            p = subprocess.run(["curl", "-sL", "-H", "Accept: application/json", url],
                               capture_output=True, timeout=600)
            return json.loads(p.stdout)
        except Exception as e:
            print("retry", url, e, file=sys.stderr); time.sleep(3)
    raise RuntimeError(url)

def main(cell, out_csv):
    # Experiments: TF ChIP-seq, released, this cell line, GRCh38 processed
    params = [("type", "Experiment"), ("assay_title", "TF ChIP-seq"), ("status", "released"),
              ("biosample_ontology.term_name", cell), ("assembly", "GRCh38"),
              ("files.output_type", "IDR thresholded peaks"), ("files.file_type", "bed narrowPeak"),
              ("limit", "all"), ("format", "json"),
              ("field", "accession"), ("field", "target.label"), ("field", "target.investigated_as"),
              ("field", "lab.title"), ("field", "date_released"), ("field", "award.rfa"),
              ("field", "replicates.biological_replicate_number"),
              ("field", "files.accession"), ("field", "files.output_type"), ("field", "files.file_type"),
              ("field", "files.assembly"), ("field", "files.status"), ("field", "files.preferred_default"),
              ("field", "files.biological_replicates"), ("field", "files.file_size"),
              ("field", "default_analysis"), ("field", "files.analyses"),
              ("field", "audit.ERROR.category"), ("field", "audit.NOT_COMPLIANT.category")]
    url = f"{BASE}/search/?" + urllib.parse.urlencode(params)
    exps = get(url).get("@graph", [])
    rows = []
    for e in exps:
        tgt = e.get("target") or {}
        inv = tgt.get("investigated_as") or []
        # keep only released GRCh38 IDR narrowPeak files
        idr = [f for f in e.get("files", []) if f.get("output_type") == "IDR thresholded peaks"
               and f.get("file_type") == "bed narrowPeak" and f.get("assembly") == "GRCh38"
               and f.get("status") == "released"]
        if not idr: continue
        da = e.get("default_analysis")
        pref = [f for f in idr if f.get("preferred_default")] or \
               [f for f in idr if da and any((a.get("@id") if isinstance(a, dict) else a) == da for a in f.get("analyses", []))] or idr
        f0 = max(pref, key=lambda f: len(f.get("biological_replicates") or []))
        rows.append({"cell_line": cell, "target": tgt.get("label"),
                     "investigated_as": ";".join(inv),
                     "is_sequence_specific_TF": "transcription factor" in inv,
                     "experiment": e["accession"], "lab": e.get("lab", {}).get("title"),
                     "award_rfa": (e.get("award") or {}).get("rfa"), "date_released": e.get("date_released"),
                     "n_bio_reps": len({r["biological_replicate_number"] for r in e.get("replicates", [])}),
                     "idr_peak_file": f0["accession"], "idr_peak_reps": "+".join(map(str, f0.get("biological_replicates") or [])),
                     "idr_peak_size": f0.get("file_size"),
                     "idr_peak_url": f"{BASE}/files/{f0['accession']}/@@download/{f0['accession']}.bed.gz",
                     "n_idr_files_grch38": len(idr),
                     "audit_error": ";".join((e.get("audit") or {}).get("ERROR", []) and [a.get("category") for a in e["audit"]["ERROR"]] or []),
                     "audit_not_compliant": ";".join((e.get("audit") or {}).get("NOT_COMPLIANT", []) and [a.get("category") for a in e["audit"]["NOT_COMPLIANT"]] or [])})
    rows.sort(key=lambda r: (not r["is_sequence_specific_TF"], r["target"] or "", r["experiment"]))
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    cat = collections.Counter()
    tf_targets = {r["target"] for r in rows if r["is_sequence_specific_TF"]}
    all_targets = {r["target"] for r in rows}
    for r in rows: cat[r["investigated_as"]] += 1
    print(f"{cell}: experiments_from_search={len(exps)} experiments_with_released_GRCh38_IDR={len(rows)} "
          f"distinct_targets_all={len(all_targets)} distinct_targets_investigated_as_TF={len(tf_targets)}")
    print("  investigated_as breakdown (experiments):", dict(cat.most_common(12)))
    return len(rows), len(all_targets), len(tf_targets)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
