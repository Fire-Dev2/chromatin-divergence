#!/usr/bin/env python3
"""Convert ENCODE ATAC filtered BAMs into compact Tn5 insertion files, then delete the BAM.

For every BAM listed in data/manifest_all.tsv under atac_bam/ that the download
manager has verified (data/manifest_done.txt), this:
  1. re-applies the ENCODE filter (properly paired, MAPQ >= 30, not
     unmapped/secondary/QC-fail/duplicate/supplementary: samtools -f 2 -F 3852 -q 30),
  2. drops chrM and non-primary contigs,
  3. takes each read's 5' end with the Tn5 shift (+4 on + strand, -5 on - strand),
  4. saves sorted insertion positions per chromosome as int32 arrays in
     data/atac_ins/<line>__rep<k>.npz, plus a small JSON with counts,
  5. deletes the BAM (and its index) so the 140 GB never sits on disk at once.

The local pipeline's tagAligns (Jurkat, RPE1, K562 validation) are converted the same
way by `--tagalign`. Safe to rerun: finished outputs are skipped.
"""
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
OUT = DATA / "atac_ins"
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]


def insertions_from_lines(lines):
    """Tn5 insertion sites from tagAlign lines (already shifted): chrom, start, end, name, score, strand."""
    pos = {c: [] for c in CHROMS}
    n = 0
    for ln in lines:
        f = ln.rstrip("\n").split("\t")
        c = f[0]
        if c not in pos:
            continue
        pos[c].append(int(f[1]) if f[5] == "+" else int(f[2]) - 1)
        n += 1
    return pos, n


def sam_stream(bam: Path):
    # samtools view → chrom, flag, start, end via a tiny awk on CIGAR-free end: use `samtools view` with
    # the reference end computed by samtools' `-e` expression is not available, so compute end from CIGAR.
    import re
    cig = re.compile(r"(\d+)([MIDNSHP=X])")
    p = subprocess.Popen(["samtools", "view", "-f", "2", "-F", "3852", "-q", "30", str(bam)],
                         stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    for ln in p.stdout:
        f = ln.split("\t", 6)
        c = f[2]
        if c not in CHROM_SET:
            continue
        flag, start = int(f[1]), int(f[3])
        ref_len = sum(int(n) for n, op in cig.findall(f[5]) if op in "MDN=X")
        yield f"{c}\t{flag}\t{start}\t{start + ref_len - 1}"
    p.wait()


CHROM_SET = set(CHROMS)


def save(pos, n, out: Path, meta: dict):
    arrs = {c: np.sort(np.asarray(v, dtype=np.int32)) for c, v in pos.items()}
    np.savez_compressed(out, **arrs)
    meta.update(n_insertions=int(n), saved=time.strftime("%F %T"))
    out.with_suffix(".json").write_text(json.dumps(meta, indent=1))


def convert_bam(rel: str):
    bam = DATA / rel
    name = bam.name.split("__ENCFF")[0]          # e.g. K562__rep1
    out = OUT / f"{name}.npz"
    if out.exists():
        bam.unlink(missing_ok=True)
        return f"{name}: already converted"
    t0 = time.time()
    pos = {c: [] for c in CHROMS}
    n = 0
    for ln in sam_stream(bam):
        c, flag, start, end = ln.split("\t")
        flag, start, end = int(flag), int(start), int(end)
        pos[c].append(start - 1 + 4 if not (flag & 16) else end - 1 - 5)
        n += 1
    save(pos, n, out, {"source": rel, "kind": "ENCODE filtered BAM", "filter": "-f 2 -F 3852 -q 30"})
    bam.unlink()
    Path(str(bam) + ".bai").unlink(missing_ok=True)
    return f"{name}: {n:,} insertions in {time.time()-t0:.0f}s; BAM deleted"


def convert_tagalign(path: str, name: str):
    out = OUT / f"{name}.npz"
    if out.exists():
        return f"{name}: already converted"
    import gzip
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        pos, n = insertions_from_lines(fh)
    save(pos, n, out, {"source": path, "kind": "local pipeline tagAlign (Tn5-shifted)"})
    return f"{name}: {n:,} insertions"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    if len(sys.argv) > 1 and sys.argv[1] == "--tagalign":
        print(convert_tagalign(sys.argv[2], sys.argv[3]))
        return
    targets = [l.split("\t")[1] for l in (DATA / "manifest_all.tsv").read_text().splitlines()[1:]
               if l.split("\t")[1].startswith("atac_bam/")]
    left = set(targets)
    while left:
        done = set((DATA / "manifest_done.txt").read_text().split()) if (DATA / "manifest_done.txt").exists() else set()
        ready = [t for t in sorted(left) if t in done]
        for t in ready:
            print(time.strftime("%H:%M:%S"), convert_bam(t), flush=True)
            left.discard(t)
        if left:
            time.sleep(60)
    print("ALLDONE")


if __name__ == "__main__":
    main()
