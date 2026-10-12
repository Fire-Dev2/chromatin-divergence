"""Helpers for the local ENCODE-equivalent ATAC-seq pipeline (scripts/atac_pipeline.py).

Pure-Python pieces that ENCODE ATAC v1.9.1 does with UCSC tools / ATAqC:
  * TSS enrichment from Tn5 insertion sites,
  * bedGraph -> bigWig (instead of bedClip + bedGraphToBigWig),
  * mean bigWig signal over intervals (validation).
Coordinates are 0-based half-open throughout (BED / tagAlign).
"""
from __future__ import annotations

import gzip
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

TSS_FLANK = 2000      # +/- bp around TSS
TSS_EDGE = 100        # bp at each end used as background
TSS_SMOOTH = 20       # bp boxcar used to smooth the aggregate profile
TSS_CENTRE = 250      # max is taken within +/- this many bp of the TSS


def read_chrom_sizes(path) -> dict:
    s = pd.read_csv(path, sep="\t", header=None, names=["chrom", "size"])
    return dict(zip(s["chrom"], s["size"].astype(np.int64)))


def iter_tagalign(path, chunksize=20_000_000):
    """Yield (chrom, start, end, strand) DataFrame chunks of a (gzipped) tagAlign."""
    cols = [0, 1, 2, 5]
    for ch in pd.read_csv(path, sep="\t", header=None, usecols=cols, chunksize=chunksize,
                          names=["chrom", "start", "end", "n", "s", "strand"],
                          dtype={"chrom": "category", "start": np.int64, "end": np.int64,
                                 "strand": "category"}):
        yield ch


def insertion_sites(df: pd.DataFrame) -> pd.DataFrame:
    """Tn5 insertion base of each read of a Tn5-SHIFTED tagAlign.

    ENCODE's shift moves '+' starts by +4 and '-' ends by -5, so the insertion base is
    `start` on '+' and `end - 1` on '-'."""
    pos = np.where(df["strand"].to_numpy() == "+", df["start"].to_numpy(), df["end"].to_numpy() - 1)
    return pd.DataFrame({"chrom": df["chrom"].astype(str).to_numpy(), "pos": pos})


def tss_profile(tagaligns, tss: pd.DataFrame, chrom_sizes: dict, flank: int = TSS_FLANK) -> np.ndarray:
    """Aggregate (summed over TSSs) per-bp Tn5 insertion counts, length 2*flank+1,
    oriented 5'->3' along the gene (minus-strand windows reversed).

    tss: DataFrame with chrom, tss0 (0-based TSS base), strand.
    tagaligns: list of tagAlign paths (summed)."""
    prof = np.zeros(2 * flank + 1, dtype=np.float64)
    # collect insertion positions per chromosome (int32 positions, compact)
    by_chrom: dict[str, list] = {}
    for path in tagaligns:
        for ch in iter_tagalign(path):
            ins = insertion_sites(ch)
            for c, g in ins.groupby("chrom", sort=False):
                by_chrom.setdefault(c, []).append(g["pos"].to_numpy(np.int64))
    offs = np.arange(-flank, flank + 1)
    for c, t in tss.groupby("chrom", sort=False):
        if c not in by_chrom or c not in chrom_sizes:
            continue
        n = int(chrom_sizes[c])
        pos = np.concatenate(by_chrom[c])
        pos = pos[(pos >= 0) & (pos < n)]
        cnt = np.bincount(pos, minlength=n)
        centres = t["tss0"].to_numpy(np.int64)
        ok = (centres - flank >= 0) & (centres + flank < n)
        centres, strand = centres[ok], t["strand"].to_numpy()[ok]
        mat = cnt[centres[:, None] + offs[None, :]]
        minus = strand == "-"
        prof += mat[~minus].sum(axis=0)
        prof += mat[minus][:, ::-1].sum(axis=0)
        del cnt
    return prof


def tss_enrichment(profile: np.ndarray, edge: int = TSS_EDGE, smooth: int = TSS_SMOOTH,
                   centre: int = TSS_CENTRE) -> dict:
    """ENCODE/Greenleaf-style TSS enrichment from an aggregate per-bp profile.

    Normalise by the mean of the `edge` bp at each end of the window, smooth with a
    centred `smooth`-bp boxcar, report the max within +/- `centre` bp of the TSS.
    Also reports the ENCODE ATAqC binned variant (10-bp bins, edge = 10 bins, max of
    the unsmoothed binned mean), computed on the same insertions."""
    profile = np.asarray(profile, dtype=float)
    flank = (len(profile) - 1) // 2
    bg = (profile[:edge].sum() + profile[-edge:].sum()) / (2 * edge)
    if bg <= 0:
        return {"tss_enrichment": float("nan"), "tss_enrichment_bin10": float("nan")}
    norm = profile / bg
    sm = pd.Series(norm).rolling(smooth, center=True, min_periods=1).mean().to_numpy()
    c = slice(flank - centre, flank + centre + 1)
    # binned variant: 400 bins of 10 bp over [-2000, +2000)
    p = profile[: 2 * flank]
    nb = len(p) // 10
    b = p[: nb * 10].reshape(nb, 10).mean(axis=1)
    ne = edge // 10
    bgb = (b[:ne].sum() + b[-ne:].sum()) / (2 * ne)
    return {"tss_enrichment": float(sm[c].max()),
            "tss_enrichment_bin10": float(b.max() / bgb) if bgb > 0 else float("nan"),
            "tss_profile_bg_per_bp": float(bg)}


def bedgraph_to_bigwig(bdg: str | Path, chrom_sizes: dict, out: str | Path,
                       chunksize: int = 20_000_000) -> int:
    """Clip a bedGraph to chromosome bounds (ENCODE: slop -b 0 | bedClip), sort and write a
    bigWig with pyBigWig. Returns number of intervals written."""
    import pyBigWig
    op = gzip.open if str(bdg).endswith(".gz") else open
    with op(bdg, "rt") as fh:
        skip = 1 if fh.readline().startswith(("track", "#")) else 0
    data: dict[str, list] = {}
    for ch in pd.read_csv(bdg, sep="\t", header=None, names=["chrom", "start", "end", "val"],
                          dtype={"chrom": "category", "start": np.int64, "end": np.int64,
                                 "val": np.float32},
                          chunksize=chunksize, skiprows=skip):
        for c, g in ch.groupby("chrom", sort=False, observed=True):
            c = str(c)
            if c not in chrom_sizes:
                continue
            n = int(chrom_sizes[c])
            s = g["start"].to_numpy(np.int64)
            e = np.minimum(g["end"].to_numpy(np.int64), n)
            v = g["val"].to_numpy(np.float32)
            k = (s >= 0) & (e > s)
            data.setdefault(c, []).append((s[k].astype(np.int32), e[k].astype(np.int32), v[k]))
    chroms = sorted(data)
    bw = pyBigWig.open(str(out), "w")
    bw.addHeader([(c, int(chrom_sizes[c])) for c in chroms], maxZooms=10)
    n = 0
    for c in chroms:
        s = np.concatenate([x[0] for x in data[c]]).astype(np.int64)
        e = np.concatenate([x[1] for x in data[c]]).astype(np.int64)
        v = np.concatenate([x[2] for x in data[c]]).astype(np.float64)
        o = np.argsort(s, kind="stable")
        s, e, v = s[o], e[o], v[o]
        # entries must be sorted and non-overlapping (MACS bedGraphs are)
        assert (s[1:] >= e[:-1]).all(), f"overlapping bedGraph intervals on {c}"
        step = 5_000_000
        for i in range(0, len(s), step):
            j = slice(i, i + step)
            bw.addEntries([c] * len(s[j]), s[j].tolist(), ends=e[j].tolist(), values=v[j].tolist())
        n += len(s)
        del data[c]
    bw.close()
    return n


def bigwig_mean(bw_path, chrom, start, end) -> np.ndarray:
    """Mean signal over [start, end) for each interval (bases without data count as 0,
    like bigWigAverageOverBed 'mean0')."""
    import pyBigWig
    bw = pyBigWig.open(str(bw_path))
    sizes = bw.chroms()
    out = np.full(len(chrom), np.nan)
    for i, (c, s, e) in enumerate(zip(chrom, start, end)):
        if c not in sizes:
            continue
        s = max(0, int(s)); e = min(int(e), sizes[c])
        if e <= s:
            continue
        iv = bw.intervals(c, s, e)
        if not iv:
            out[i] = 0.0
            continue
        a = np.array(iv, dtype=float)
        ov = np.minimum(a[:, 1], e) - np.maximum(a[:, 0], s)
        out[i] = float((ov * a[:, 2]).sum() / (e - s))
    bw.close()
    return out


def count_lines(path) -> int:
    """Line count of a plain or gzipped text file (uses system tools)."""
    path = str(path)
    cmd = f"gzip -dc '{path}' | wc -l" if path.endswith(".gz") else f"wc -l < '{path}'"
    return int(subprocess.run(["bash", "-c", cmd], capture_output=True, text=True, check=True).stdout.split()[0])


def read_narrowpeak(path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", header=None, compression="infer",
                       names=["chrom", "start", "end", "name", "score", "strand",
                              "fc", "p", "q", "summit"],
                       usecols=range(10), dtype={"chrom": str})
