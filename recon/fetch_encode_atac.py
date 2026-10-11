#!/usr/bin/env python3
"""Query ENCODE for ATAC-seq experiments of interest and dump metadata JSON."""
import json, sys, urllib.request, urllib.parse, time

BASE = "https://www.encodeproject.org"

import subprocess
def get(url):
    for i in range(3):
        try:
            p = subprocess.run(["curl", "-sL", "-H", "Accept: application/json", url],
                               capture_output=True, timeout=300)
            return json.loads(p.stdout)
        except Exception as e:
            print("retry", url, e, file=sys.stderr); time.sleep(3)
    raise RuntimeError(url)

def search(**params):
    params.setdefault("type", "Experiment")
    params.setdefault("limit", "all")
    params.setdefault("format", "json")
    q = urllib.parse.urlencode(params, doseq=True)
    return get(f"{BASE}/search/?{q}").get("@graph", [])

if __name__ == "__main__":
    out = {}
    # 1. any ATAC experiments for lines not found as GRCh38 'cell line'
    for term in ["Jurkat, Clone E6-1", "Jurkat clone E6-1", "hTERT-RPE1", "RPE1", "HeLa-S3", "HEK293", "HEK293T", "GM12878", "HepG2", "K562"]:
        r = search(assay_title="ATAC-seq", **{"biosample_ontology.term_name": term},
                   field=["accession","status","lab.title","assembly","award.rfa",
                          "replicates.biological_replicate_number","biosample_ontology.term_name",
                          "biosample_ontology.classification","internal_status"])
        out[term] = [{k: e.get(k) for k in ("accession","status","assembly","award","internal_status")}
                     | {"lab": e["lab"]["title"], "term": e["biosample_ontology"]["term_name"],
                        "class": e["biosample_ontology"]["classification"],
                        "nrep": len({x["biological_replicate_number"] for x in e.get("replicates", [])})}
                     for e in r]
        print(term, len(r), [x["accession"]+":"+x["status"] for x in out[term]])
    # free-text searchTerm for Jurkat / RPE
    for term in ["Jurkat", "RPE"]:
        r = search(assay_title="ATAC-seq", searchTerm=term, field=["accession","status","biosample_ontology.term_name"])
        print("searchTerm", term, [(x["accession"], x["biosample_ontology"]["term_name"], x["status"]) for x in r])
    json.dump(out, open(sys.argv[1], "w"), indent=1)
