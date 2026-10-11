#!/usr/bin/env python3
"""Single resumable download manager for every large project file.

Builds data/manifest_all.tsv (path, url, size, md5, priority) from the source
APIs, then downloads in priority order with limited concurrency. Each file is
resumed (curl -C -) until its size matches, then md5-checked; a bad md5 means
delete and refetch. Safe to rerun at any time: finished files are skipped.

Usage:  python data/fetch_all.py            # build manifest if missing, download
        python data/fetch_all.py --status   # show progress only
"""
import concurrent.futures as cf
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
MANIFEST = DATA / "manifest_all.tsv"
DONE = DATA / "manifest_done.txt"
CONCURRENCY = 4


def curl_json(url):
    out = subprocess.run(["curl", "-s", "-L", "-m", "60", "-H", "Accept: application/json", url],
                         capture_output=True, text=True).stdout
    return json.loads(out)


def build_manifest():
    rows = []  # (priority, relpath, url, size, md5)
    # 1. Perturb-seq screens
    fig = {"K562_essential_raw_singlecell_01.h5ad": 35773219, "rpe1_raw_singlecell_01.h5ad": 35775606}
    for name, fid in fig.items():
        d = curl_json(f"https://api.figshare.com/v2/file/{fid}") if False else None
        rows.append((1, f"perturbseq/{name}", f"https://ndownloader.figshare.com/files/{fid}", None, None))
    for line in ("jurkat", "hepg2"):
        name = f"GSE264667_{line}_raw_singlecell_01.h5ad"
        rows.append((2, f"perturbseq/{name}",
                     f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE264nnn/GSE264667/suppl/{name}", None, None))
    # figshare article metadata gives size + md5 for the Replogle files
    try:
        art = curl_json("https://api.figshare.com/v2/articles/20029387/files?page_size=100")
        meta = {f["id"]: f for f in art}
        rows = [(p, r, u, meta.get(int(u.rsplit("/", 1)[1]), {}).get("size") if "figshare" in u else s,
                 meta.get(int(u.rsplit("/", 1)[1]), {}).get("computed_md5") if "figshare" in u else m)
                for p, r, u, s, m in rows]
    except Exception as e:
        print("figshare metadata unavailable:", e)
    # 2. Genome (ENCODE GRCh38 no-alt + bowtie2 index)
    for acc, fname in (("ENCFF110MCL", "ENCFF110MCL.tar.gz"),
                       ("GRCh38_no_alt_analysis_set_GCA_000001405.15",
                        "GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz")):
        d = curl_json(f"https://www.encodeproject.org/files/{acc}/?format=json")
        rows.append((3, f"genome/{fname}",
                     f"https://www.encodeproject.org/files/{acc}/@@download/{fname}",
                     d.get("file_size"), d.get("md5sum")))
    # 3. Jurkat raw reads (ENA gives size + md5)
    for line in (ROOT / "recon/jurkat_atac_selected_runs.tsv").read_text().splitlines()[1:]:
        f = line.split("\t")
        for url, size, md5 in zip(f[6].split(";"), f[7].split(";"), f[8].split(";")):
            rows.append((4, f"atac_fastq/jurkat_GSE273644/{url.rsplit('/', 1)[1]}",
                         "https://" + url, int(size), md5))
    # 4. ATAC signal + peaks
    for line in (DATA / "atac/manifest.tsv").read_text().splitlines()[1:]:
        ln, role, acc, kind, url = line.split("\t")
        ext = ".bed.gz" if url.endswith(".bed.gz") else "." + url.rsplit(".", 1)[1]
        size = md5 = None
        if acc.startswith("ENCFF"):
            d = curl_json(f"https://www.encodeproject.org/files/{acc}/?format=json")
            size, md5 = d.get("file_size"), d.get("md5sum")
        else:
            d = curl_json(f"https://data.4dnucleome.org/files-processed/{acc}/?format=json&frame=object")
            size, md5 = d.get("file_size"), d.get("md5sum")
        rows.append((5, f"atac/{ln}__{acc}{ext}", url, size, md5))
    # 5. RPE1 raw reads (4DN)
    for line in (DATA / "atac_fastq/rpe1_4DNES48BA2C8/urls.tsv").read_text().splitlines():
        rep, acc, md5, url = line.split("\t")
        d = curl_json(f"https://data.4dnucleome.org/files-fastq/{acc}/?format=json&frame=object")
        rows.append((6, f"atac_fastq/rpe1_4DNES48BA2C8/{acc}.fastq.gz", url, d.get("file_size"), md5))
    rows.sort(key=lambda r: r[0])
    with open(MANIFEST, "w") as fh:
        fh.write("priority\tpath\turl\tsize\tmd5\n")
        for p, r, u, s, m in rows:
            fh.write(f"{p}\t{r}\t{u}\t{s or ''}\t{m or ''}\n")
    print(f"manifest: {len(rows)} files -> {MANIFEST}")


def read_manifest():
    out = []
    for line in MANIFEST.read_text().splitlines()[1:]:
        p, r, u, s, m = (line.split("\t") + [""] * 5)[:5]
        out.append(dict(priority=int(p), path=DATA / r, url=u, size=int(s) if s else None, md5=m or None))
    return out


def md5_of(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 24), b""):
            h.update(chunk)
    return h.hexdigest()


def done_set():
    return set(DONE.read_text().split()) if DONE.exists() else set()


def fetch(item, max_tries=50):
    path, url, size, md5 = item["path"], item["url"], item["size"], item["md5"]
    rel = str(path.relative_to(DATA))
    if rel in done_set():
        return rel, "already done"
    path.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(1, max_tries + 1):
        have = path.stat().st_size if path.exists() else 0
        if size and have > size:
            path.unlink(); have = 0
        if not size or have < size:
            subprocess.run(["curl", "-L", "-sS", "-C", "-", "--retry", "3", "--retry-delay", "5",
                            "--speed-limit", "20000", "--speed-time", "120",
                            "-o", str(path), url])
            have = path.stat().st_size if path.exists() else 0
            if size and have < size:
                time.sleep(min(60, 5 * attempt)); continue
        if md5:
            got = md5_of(path)
            if got != md5:
                print(f"[md5 mismatch] {rel}: {got} != {md5}; refetching", flush=True)
                path.unlink(); continue
        with open(DONE, "a") as fh:
            fh.write(rel + "\n")
        return rel, f"ok ({have/1e9:.2f} GB{', md5 ok' if md5 else ''})"
    return rel, "FAILED after retries"


def status():
    done = done_set()
    tot = got = 0
    for it in read_manifest():
        rel = str(it["path"].relative_to(DATA))
        have = it["path"].stat().st_size if it["path"].exists() else 0
        size = it["size"] or have
        tot += size; got += min(have, size)
        flag = "done" if rel in done else f"{have/1e9:6.2f}/{size/1e9:.2f} GB"
        print(f"  p{it['priority']} {rel:70s} {flag}")
    print(f"  total {got/1e9:.1f} / {tot/1e9:.1f} GB, {len(done)} / {len(read_manifest())} files verified")


if __name__ == "__main__":
    if not MANIFEST.exists():
        build_manifest()
    if "--status" in sys.argv:
        status(); sys.exit()
    items = read_manifest()
    with cf.ThreadPoolExecutor(CONCURRENCY) as ex:
        for rel, msg in ex.map(fetch, items):
            print(f"{time.strftime('%H:%M:%S')} {rel}: {msg}", flush=True)
    status()
    print("ALLDONE" if len(done_set()) == len(items) else "INCOMPLETE")
