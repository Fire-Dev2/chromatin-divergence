"""Per-gene chromatin-accessibility features (ANALYSIS_PLAN.md section 3).

Pipeline (all coordinates 0-based half-open, chr1-22 + chrX only):

1. **Union regions.** For every configured line, take its IDR (or conservative) peaks, build the
   window ``[summit - 250, summit + 250)``, pool all lines, merge overlapping / book-ended windows,
   and drop merged regions that overlap the ENCODE blacklist v2.
2. **Region signal.** Two versions (``signal=``):

   * ``"cpm"`` (primary, plan amendment v1.2): Tn5 insertions in each region per million total
     insertions of the line (chr1-22,X); replicates of a line are pooled (counts and totals summed).
     Insertion positions are 0-based; an insertion at p counts for ``[start, end)`` iff
     ``start <= p < end`` (``np.searchsorted`` on the sorted per-chrom arrays).
   * ``"fc"`` (sensitivity, the v1.0 track): mean fold-change bigWig signal over the region
     (uncovered bases count as 0).

   Then ``log2(x + 1)``, quantile normalisation across lines, then subtract the line's log2
   relative copy number over the region -> ``s_r``.
3. **Promoter feature.** Same scheme on the TSS +/-1 kb window ``[tss0 - 1000, tss0 + 1000)``
   computed directly from the insertions / bigWig; quantile normalisation uses the promoter matrix
   itself.
4. **Distal feature.** ``log2(1 + sum_r 2**s_r * exp(-d_r / 25 kb))`` over union regions whose
   centre lies 2.5-100 kb (inclusive) from the TSS; ``d_r`` = |centre - tss0|.
5. **Combined.** ``0.5 * z(promoter) + 0.5 * z(distal)``; z uses one mean and SD pooled over all
   genes x all configured lines.
6. ``pair_feature(A, B, kind)`` = feature_A - feature_B per gene.

Lines and input files come from ``config/feature_lines.tsv`` so the build is re-runnable when the
local Jurkat / RPE1 outputs appear.
CLI: ``python -m cdtest.features [--signal cpm|fc|both] [--config ...] [--procs 6]``.
Outputs: cpm -> ``gene_features.parquet``, ``region_signal(_raw).parquet``, ``build_info.json``,
``promoter_cpm_replicates.parquet``; fc -> the same names with suffix ``_fc``.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import gzip
import hashlib
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from cdtest import annot

ROOT = annot.ROOT
CONFIG_TSV = ROOT / "config" / "feature_lines.tsv"
OUT_DIR = ROOT / "data" / "features"
REGION_SIGNAL = OUT_DIR / "region_signal.parquet"
REGION_RAW = OUT_DIR / "region_signal_raw.parquet"
GENE_FEATURES = OUT_DIR / "gene_features.parquet"
BUILD_INFO = OUT_DIR / "build_info.json"
PROMOTER_REPS = OUT_DIR / "promoter_cpm_replicates.parquet"
SIGNALS = ("cpm", "fc")
DEFAULT_SIGNAL = "cpm"


def output_paths(signal: str = DEFAULT_SIGNAL, out_dir: Path | str = OUT_DIR) -> dict:
    """File paths of one signal version (cpm = primary, no suffix; fc = suffix '_fc')."""
    if signal not in SIGNALS:
        raise ValueError(f"signal must be one of {SIGNALS}")
    out_dir = Path(out_dir)
    sfx = "" if signal == "cpm" else f"_{signal}"
    return dict(gene_features=out_dir / f"gene_features{sfx}.parquet",
                region_signal=out_dir / f"region_signal{sfx}.parquet",
                region_raw=out_dir / f"region_signal_raw{sfx}.parquet",
                build_info=out_dir / f"build_info{sfx}.json",
                promoter_reps=out_dir / "promoter_cpm_replicates.parquet")

CHROMS = annot.CHROMS
SUMMIT_HALF = 250          # bp each side of the summit
PROMOTER_HALF = 1000       # TSS +/- 1 kb
DISTAL_MIN = 2_500         # bp, inclusive
DISTAL_MAX = 100_000       # bp, inclusive
DISTAL_SCALE = 25_000      # bp, exp(-d / scale)
PSEUDO = 1.0               # log2(x + PSEUDO)
KINDS = ("promoter", "distal", "combined")


# --------------------------------------------------------------------------- config
def insertion_files(pattern: str, root: Path = ROOT) -> list[Path]:
    """Sorted replicate insertion files matching a glob relative to root ('' -> none)."""
    if not pattern:
        return []
    return sorted(Path(root).glob(pattern))


def load_config(path: Path | str = CONFIG_TSV, root: Path = ROOT,
                signal: str = DEFAULT_SIGNAL) -> tuple[pd.DataFrame, list[dict]]:
    """Resolve the line config: per line, the lowest-priority row whose files for `signal` exist
    (cpm: insertions glob non-empty + peaks; fc: fc_bigwig + peaks).

    Returns (selected rows in config order, list of skipped lines with reason)."""
    if signal not in SIGNALS:
        raise ValueError(f"signal must be one of {SIGNALS}")
    cfg = pd.read_csv(path, sep="\t", comment="#", dtype=str, keep_default_na=False)
    cfg["priority"] = cfg.priority.astype(int)
    if "insertions" not in cfg.columns:
        cfg["insertions"] = ""
    order = list(dict.fromkeys(cfg.line))
    chosen, skipped = [], []

    def usable(r) -> bool:
        if not (root / r.peaks).exists():
            return False
        if signal == "cpm":
            return len(insertion_files(r.insertions, root)) > 0
        return (root / r.fc_bigwig).exists()

    for line in order:
        rows = cfg[cfg.line == line].sort_values("priority")
        hit = None
        for _, r in rows.iterrows():
            if usable(r):
                hit = r
                break
        if hit is None:
            col = "insertions" if signal == "cpm" else "fc_bigwig"
            skipped.append(dict(line=line, reason=f"no candidate row with existing files (signal={signal})",
                                candidates=[f"{r[col]} | {r.peaks}" for _, r in rows.iterrows()]))
        else:
            chosen.append(hit)
    sel = pd.DataFrame(chosen, columns=cfg.columns).reset_index(drop=True)
    if sel.line.duplicated().any():
        raise ValueError("duplicate line after selection")
    return sel, skipped


# --------------------------------------------------------------------------- peaks / regions
def read_summits(path: Path | str) -> pd.DataFrame:
    """Summits (chrom, summit0) from a narrowPeak(.gz) or narrowPeak bigBed.

    summit0 = chromStart + peak offset (col 10); if the offset is -1, the peak midpoint."""
    path = Path(path)
    if path.suffix in (".bb", ".bigBed", ".bigbed"):
        import pyBigWig
        bb = pyBigWig.open(str(path))
        rows = []
        for c, L in bb.chroms().items():
            if c not in CHROMS:
                continue
            for s, e, rest in bb.entries(c, 0, L) or []:
                off = int(rest.split("\t")[-1])
                rows.append((c, s, e, off))
        bb.close()
        d = pd.DataFrame(rows, columns=["chrom", "start", "end", "peak"])
    else:
        d = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2, 9],
                        names=["chrom", "start", "end", "peak"],
                        dtype={"chrom": str, "start": np.int64, "end": np.int64, "peak": np.int64})
        d = d[d.chrom.isin(CHROMS)]
    s = np.where(d.peak.values >= 0, d.start.values + d.peak.values, (d.start.values + d.end.values) // 2)
    return pd.DataFrame({"chrom": d.chrom.values, "summit0": s.astype(np.int64)})


def merge_intervals(df: pd.DataFrame) -> pd.DataFrame:
    """Merge overlapping or book-ended intervals (bedtools merge -d 0). Sorted by CHROMS order."""
    out = []
    for c in CHROMS:
        sub = df[df.chrom == c]
        if sub.empty:
            continue
        o = np.lexsort((sub.end.values, sub.start.values))
        st, en = sub.start.values[o], sub.end.values[o]
        run_end = np.maximum.accumulate(en)
        new = np.r_[True, st[1:] > run_end[:-1]]
        gid = np.cumsum(new) - 1
        ms = st[new]
        me = np.zeros(len(ms), dtype=np.int64)
        np.maximum.at(me, gid, en)
        out.append(pd.DataFrame({"chrom": c, "start": ms, "end": me}))
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame(columns=["chrom", "start", "end"])


def union_regions(summits: list[pd.DataFrame], blacklist: pd.DataFrame, chrom_sizes: dict | None = None,
                  half: int = SUMMIT_HALF) -> tuple[pd.DataFrame, dict]:
    """Union of summit +/- half windows over lines, merged, blacklist-overlapping regions removed."""
    s = pd.concat(summits, ignore_index=True)
    w = pd.DataFrame({"chrom": s.chrom.values, "start": np.maximum(s.summit0.values - half, 0),
                      "end": s.summit0.values + half})
    if chrom_sizes:
        w["end"] = np.minimum(w.end.values, w.chrom.map(chrom_sizes).values)
    m = merge_intervals(w)
    bl = annot.overlaps_intervals(m.chrom, m.start, m.end, blacklist)
    info = dict(n_windows=int(len(w)), n_merged=int(len(m)), n_blacklisted=int(bl.sum()))
    m = m[~bl].reset_index(drop=True)
    m.insert(0, "region_id", [f"r{i}" for i in range(len(m))])
    return m, info


# --------------------------------------------------------------------------- signal
def mean_bigwig(path: str, chrom, start, end) -> np.ndarray:
    """Exact mean bigWig value over each [start, end) (bases without data count as 0)."""
    import pyBigWig
    bw = pyBigWig.open(str(path))
    sizes = bw.chroms()
    chrom = np.asarray(chrom).astype(str)
    start = np.asarray(start, dtype=np.int64)
    end = np.asarray(end, dtype=np.int64)
    out = np.zeros(len(chrom))
    for i, (c, s, e) in enumerate(zip(chrom, start, end)):
        if c not in sizes:
            continue
        s2, e2 = max(int(s), 0), min(int(e), sizes[c])
        if e2 <= s2:
            continue
        v = bw.stats(c, s2, e2, type="sum", exact=True)[0]
        out[i] = (v or 0.0) / (e - s)
    bw.close()
    return out


def _mean_bigwig_job(args):
    return mean_bigwig(*args)


def signal_matrix(paths: dict, chrom, start, end, procs: int = 6) -> pd.DataFrame:
    """Columns = lines, mean FC over each interval; one process per line (max `procs`)."""
    lines = list(paths)
    jobs = [(str(paths[ln]), chrom, start, end) for ln in lines]
    if procs <= 1:
        res = [_mean_bigwig_job(j) for j in jobs]
    else:
        with ProcessPoolExecutor(max_workers=min(procs, len(jobs))) as ex:
            res = list(ex.map(_mean_bigwig_job, jobs))
    return pd.DataFrame(dict(zip(lines, res)))


def count_in_intervals(pos: np.ndarray, start, end) -> np.ndarray:
    """Number of sorted positions p with start <= p < end, per interval (np.searchsorted)."""
    pos = np.asarray(pos)
    return (np.searchsorted(pos, np.asarray(end), side="left")
            - np.searchsorted(pos, np.asarray(start), side="left")).astype(np.int64)


def count_insertions(path: Path | str, interval_sets: list[tuple]) -> tuple[list[np.ndarray], int]:
    """Insertion counts of one replicate file over several interval sets.

    path: npz with one sorted 0-based int array per chrom (only CHROMS are used; read lazily).
    interval_sets: list of (chrom, start, end) arrays. Returns (counts per set, total insertions on
    CHROMS)."""
    sets = [(np.asarray(c).astype(str), np.asarray(s, np.int64), np.asarray(e, np.int64))
            for c, s, e in interval_sets]
    out = [np.zeros(len(c), dtype=np.int64) for c, _, _ in sets]
    total = 0
    with np.load(path) as z:
        keys = set(z.files)
        for chrom in CHROMS:
            if chrom not in keys:
                continue
            pos = z[chrom]
            total += int(len(pos))
            for k, (c, s, e) in enumerate(sets):
                m = np.flatnonzero(c == chrom)
                if len(m):
                    out[k][m] = count_in_intervals(pos, s[m], e[m])
    return out, total


def _count_job(args):
    return count_insertions(*args)


def cpm_matrices(ins: dict, interval_sets: list[tuple], procs: int = 6) -> dict:
    """Pooled-replicate insertion CPM per line for each interval set.

    ins maps line -> list of replicate files. Counts and totals are summed over replicates, then
    CPM = count / total * 1e6. Returns dict(cpm=[DataFrame per set], rep_cpm=[DataFrame per set with
    columns = replicate file stem, e.g. 'K562__rep1'], totals={line: n}, rep_totals={col: n})."""
    jobs, keys = [], []
    for ln, files in ins.items():
        for f in files:
            jobs.append((str(f), interval_sets))
            keys.append((ln, Path(f).stem))
    if len({k for _, k in keys}) != len(keys):
        raise ValueError("replicate file stems must be unique across lines")
    if procs <= 1:
        res = [_count_job(j) for j in jobs]
    else:
        with ProcessPoolExecutor(max_workers=min(procs, len(jobs))) as ex:
            res = list(ex.map(_count_job, jobs))
    nset = len(interval_sets)
    pooled = {ln: [0] * nset for ln in ins}
    totals = {ln: 0 for ln in ins}
    rep_cpm = [dict() for _ in range(nset)]
    rep_totals = {}
    for (ln, col), (cnts, tot) in zip(keys, res):
        totals[ln] += tot
        rep_totals[col] = tot
        for k in range(nset):
            pooled[ln][k] = pooled[ln][k] + cnts[k]
            rep_cpm[k][col] = cnts[k] / tot * 1e6
    cpm = [pd.DataFrame({ln: pooled[ln][k] / totals[ln] * 1e6 for ln in ins}) for k in range(nset)]
    return dict(cpm=cpm, rep_cpm=[pd.DataFrame(d) for d in rep_cpm], totals=totals, rep_totals=rep_totals)


def quantile_normalize(x: pd.DataFrame | np.ndarray) -> np.ndarray:
    """Quantile normalisation across columns. Reference distribution = mean of the sorted columns;
    ties get the mean reference value of their tied ranks (via average ranks + interpolation)."""
    a = np.asarray(x, dtype=float)
    if np.isnan(a).any():
        raise ValueError("quantile_normalize: NaN in input")
    n = a.shape[0]
    ref = np.sort(a, axis=0).mean(axis=1)
    out = np.empty_like(a)
    pos = np.arange(n, dtype=float)
    for j in range(a.shape[1]):
        r = rankdata(a[:, j], method="average") - 1.0
        out[:, j] = np.interp(r, pos, ref)
    return out


def normalize_signal(raw: pd.DataFrame, cn: pd.DataFrame) -> pd.DataFrame:
    """log2(x+1) -> quantile normalise across lines -> subtract log2 CN (same column order)."""
    lg = np.log2(raw.values + PSEUDO)
    qn = quantile_normalize(lg)
    return pd.DataFrame(qn - cn[raw.columns].values, columns=raw.columns, index=raw.index)


def cn_matrix(cn_lines: dict, chrom, start, end, lookup: annot.CNLookup | None = None) -> pd.DataFrame:
    """log2 relative CN per interval; cn_lines maps feature line -> CN line name ('none' -> 0)."""
    lookup = annot.CNLookup() if lookup is None else lookup
    out = {}
    for ln, cl in cn_lines.items():
        out[ln] = np.zeros(len(chrom)) if cl in ("none", "", None) else lookup.region(cl, chrom, start, end)
    return pd.DataFrame(out)


# --------------------------------------------------------------------------- gene features
def distal_pairs(genes: pd.DataFrame, regions: pd.DataFrame, dmin: int = DISTAL_MIN,
                 dmax: int = DISTAL_MAX) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(gene row index, region row index, distance) for regions whose centre is within [dmin, dmax]
    of the gene TSS (tss0). Region centre = (start + end) / 2."""
    gi, ri, dd = [], [], []
    centre_all = (regions.start.values + regions.end.values) / 2.0
    for c in CHROMS:
        gm = np.where(genes.chrom.values == c)[0]
        rm = np.where(regions.chrom.values == c)[0]
        if len(gm) == 0 or len(rm) == 0:
            continue
        cen = centre_all[rm]
        o = np.argsort(cen, kind="stable")
        cen, rm = cen[o], rm[o]
        t = genes.tss0.values[gm].astype(float)
        lo = np.searchsorted(cen, t - dmax, side="left")
        hi = np.searchsorted(cen, t + dmax, side="right")
        cnt = hi - lo
        g = np.repeat(gm, cnt)
        start_rep = np.repeat(lo, cnt)
        offs = np.arange(cnt.sum()) - np.repeat(np.cumsum(cnt) - cnt, cnt)
        k = start_rep + offs
        d = np.abs(cen[k] - genes.tss0.values[g])
        keep = (d >= dmin) & (d <= dmax)
        gi.append(g[keep]); ri.append(rm[k][keep]); dd.append(d[keep])
    if not gi:
        return np.array([], int), np.array([], int), np.array([], float)
    return np.concatenate(gi), np.concatenate(ri), np.concatenate(dd)


def distal_feature(n_genes: int, gi: np.ndarray, ri: np.ndarray, d: np.ndarray, s: np.ndarray,
                   scale: float = DISTAL_SCALE) -> np.ndarray:
    """log2(1 + sum_r 2**s_r * exp(-d_r/scale)); s is (n_regions, n_lines). Returns (n_genes, n_lines)."""
    s = np.asarray(s, dtype=float)
    if s.ndim == 1:
        s = s[:, None]
    w = np.exp(-np.asarray(d, float) / scale)
    out = np.zeros((n_genes, s.shape[1]))
    for j in range(s.shape[1]):
        out[:, j] = np.bincount(gi, weights=w * np.exp2(s[ri, j]), minlength=n_genes)
    return np.log2(1.0 + out)


def pooled_z(mat: np.ndarray) -> tuple[np.ndarray, float, float]:
    """z-score with ONE mean and SD (ddof=0) pooled over all entries (genes x lines)."""
    a = np.asarray(mat, dtype=float)
    mu, sd = float(np.nanmean(a)), float(np.nanstd(a))
    return (a - mu) / sd, mu, sd


def combined_feature(prom: np.ndarray, dist: np.ndarray) -> tuple[np.ndarray, dict]:
    zp, mp, sp = pooled_z(prom)
    zd, md, sdd = pooled_z(dist)
    return 0.5 * zp + 0.5 * zd, dict(promoter_mean=mp, promoter_sd=sp, distal_mean=md, distal_sd=sdd)


# --------------------------------------------------------------------------- build
def _md5(p: Path, chunk=1 << 22) -> str:
    h = hashlib.md5()
    with open(p, "rb") as fh:
        while b := fh.read(chunk):
            h.update(b)
    return h.hexdigest()


def build(config: Path | str = CONFIG_TSV, out_dir: Path = OUT_DIR, procs: int = 6, md5: bool = True,
          signal: str = DEFAULT_SIGNAL) -> dict:
    """Build one signal version ('cpm' primary, 'fc' sensitivity) of all features."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = output_paths(signal, out_dir)
    sel, skipped = load_config(config, signal=signal)
    lines = list(sel.line)
    cn_lines = {r.line: r.cn_line for _, r in sel.iterrows()}
    print(f"[features:{signal}] lines ({len(lines)}): {lines}; skipped: {[s['line'] for s in skipped]}",
          flush=True)

    # 1. regions
    summits = []
    n_peaks = {}
    for _, r in sel.iterrows():
        sm = read_summits(ROOT / r.peaks)
        n_peaks[r.line] = int(len(sm))
        summits.append(sm)
    sizes = annot.chrom_sizes()
    regions, rinfo = union_regions(summits, annot.load_blacklist(), sizes)
    print(f"[features:{signal}] union regions: {len(regions):,} ({rinfo})", flush=True)

    genes = annot.load_genes()
    genes = genes[genes.chrom.isin(CHROMS)].reset_index(drop=True)
    ps = np.maximum(genes.tss0.values - PROMOTER_HALF, 0)
    pe = np.minimum(genes.tss0.values + PROMOTER_HALF, genes.chrom.map(sizes).values)
    rset = (regions.chrom.values, regions.start.values, regions.end.values)
    pset = (genes.chrom.values, ps, pe)

    # 2./3. raw signal over regions and promoters
    totals = rep_totals = None
    if signal == "cpm":
        ins = {r.line: insertion_files(r.insertions) for _, r in sel.iterrows()}
        res = cpm_matrices(ins, [rset, pset], procs)
        raw, praw = res["cpm"]
        totals, rep_totals = res["totals"], res["rep_totals"]
        rep = res["rep_cpm"][1].astype(np.float32)
        rep.insert(0, "gene_id", genes.gene_id.values)
        rep.to_parquet(paths["promoter_reps"], index=False)
    else:
        bw = {r.line: ROOT / r.fc_bigwig for _, r in sel.iterrows()}
        raw = signal_matrix(bw, *rset, procs)
        praw = signal_matrix(bw, *pset, procs)
    raw, praw = raw[lines], praw[lines]
    print(f"[features:{signal}] raw signal done", flush=True)

    lookup = annot.CNLookup()
    cn_r = cn_matrix(cn_lines, *rset, lookup)
    s = normalize_signal(raw, cn_r)
    pd.concat([regions, raw.astype(np.float32)], axis=1).to_parquet(paths["region_raw"], index=False)
    pd.concat([regions, s.astype(np.float32)], axis=1).to_parquet(paths["region_signal"], index=False)

    cn_p = cn_matrix(cn_lines, *pset, lookup)
    prom = normalize_signal(praw, cn_p).values

    # 4. distal
    gi, ri, dd = distal_pairs(genes, regions)
    dist = distal_feature(len(genes), gi, ri, dd, s[lines].values)
    n_distal = np.bincount(gi, minlength=len(genes))

    # 5. combined
    comb, zinfo = combined_feature(prom, dist)

    G, L = len(genes), len(lines)
    gf = pd.DataFrame({
        "gene_id": np.repeat(genes.gene_id.values, L),
        "gene_name": np.repeat(genes.gene_name.values, L),
        "line": np.tile(lines, G),
        "promoter": prom.ravel(), "distal": dist.ravel(), "combined": comb.ravel(),
        "promoter_raw": praw.values.ravel(), "promoter_cn_log2": cn_p[lines].values.ravel(),
        "n_distal_regions": np.repeat(n_distal, L),
    })
    gf["line"] = pd.Categorical(gf.line, categories=lines)
    gf.to_parquet(paths["gene_features"], index=False)

    def line_info(r):
        d = dict(line=r.line, role=r.role, priority=int(r.priority), peaks=r.peaks, cn_line=r.cn_line,
                 cn_note=("CN=0 (line absent from DepMap 24Q4)" if r.cn_line in ("GM12878", "IMR90")
                          else "RPE1 karyotype track (+log2(3/2) 10q21.2-qter)" if r.cn_line == "RPE1"
                          else "none" if r.cn_line == "none" else "DepMap 24Q4 segments"),
                 n_summits=n_peaks[r.line], md5_peaks=_md5(ROOT / r.peaks) if md5 else None)
        if signal == "cpm":
            files = insertion_files(r.insertions)
            d.update(insertions=[str(f.relative_to(ROOT)) for f in files],
                     n_insertions_total=int(totals[r.line]),
                     n_insertions_rep={f.stem: int(rep_totals[f.stem]) for f in files},
                     md5_insertions={f.stem: _md5(f) for f in files} if md5 else None)
        else:
            d.update(fc_bigwig=r.fc_bigwig, md5_fc_bigwig=_md5(ROOT / r.fc_bigwig) if md5 else None)
        return d

    info = dict(
        built=_dt.datetime.now().isoformat(timespec="seconds"),
        signal=signal,
        signal_desc=("Tn5 insertion counts per million insertions (chr1-22,X), replicates pooled"
                     if signal == "cpm" else "mean MACS2 fold-change bigWig signal"),
        config=str(Path(config).resolve().relative_to(ROOT)) if Path(config).resolve().is_relative_to(ROOT) else str(config),
        lines=[line_info(r) for _, r in sel.iterrows()],
        skipped_lines=skipped,
        n_lines=L, n_regions=int(len(regions)), region_info=rinfo,
        n_genes=int(G), n_gene_region_pairs=int(len(gi)),
        n_genes_no_distal_region=int((n_distal == 0).sum()),
        params=dict(summit_half=SUMMIT_HALF, promoter_half=PROMOTER_HALF, distal_min=DISTAL_MIN,
                    distal_max=DISTAL_MAX, distal_scale=DISTAL_SCALE, pseudocount=PSEUDO,
                    chroms="chr1-22,chrX", blacklist=str(annot.BLACKLIST_PATH.relative_to(ROOT)),
                    cn_fill="nearest segment for DepMap gaps (annot.CNLookup)"),
        z_pooled=zinfo,
        outputs={k: str(v.relative_to(ROOT)) if v.resolve().is_relative_to(ROOT) else str(v)
                 for k, v in paths.items() if signal == "cpm" or k != "promoter_reps"},
    )
    with open(paths["build_info"], "w") as fh:
        json.dump(info, fh, indent=2)
    print(f"[features:{signal}] genes {G:,} x lines {L}; written to {out_dir}", flush=True)
    return info


# --------------------------------------------------------------------------- access
_CACHE: dict = {}


def load_gene_features(path: Path | str | None = None, signal: str = DEFAULT_SIGNAL) -> pd.DataFrame:
    """Long gene-feature table of one signal version (cpm primary; fc sensitivity)."""
    path = Path(path) if path is not None else output_paths(signal)["gene_features"]
    key = (str(path), path.stat().st_mtime)
    if key not in _CACHE:
        _CACHE[key] = pd.read_parquet(path)
    return _CACHE[key]


def feature_matrix(kind: str = "promoter", features: pd.DataFrame | None = None,
                   signal: str = DEFAULT_SIGNAL) -> pd.DataFrame:
    """genes x lines matrix of one feature kind."""
    if kind not in KINDS:
        raise ValueError(f"kind must be one of {KINDS}")
    f = load_gene_features(signal=signal) if features is None else features
    m = f.pivot(index="gene_id", columns="line", values=kind)
    m.columns = [str(c) for c in m.columns]
    return m


def pair_feature(lineA: str, lineB: str, kind: str = "combined", features: pd.DataFrame | None = None,
                 signal: str = DEFAULT_SIGNAL) -> pd.Series:
    """x_g = feature_A - feature_B for every gene (index gene_id)."""
    m = feature_matrix(kind, features, signal)
    for ln in (lineA, lineB):
        if ln not in m.columns:
            raise KeyError(f"line {ln!r} not in features (have {list(m.columns)})")
    return (m[lineA] - m[lineB]).rename(f"{kind}:{lineA}-{lineB}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config", default=str(CONFIG_TSV))
    ap.add_argument("--out", default=str(OUT_DIR))
    ap.add_argument("--procs", type=int, default=6)
    ap.add_argument("--signal", choices=[*SIGNALS, "both"], default=DEFAULT_SIGNAL)
    ap.add_argument("--no-md5", action="store_true")
    a = ap.parse_args(argv)
    for sig in (SIGNALS if a.signal == "both" else (a.signal,)):
        build(a.config, Path(a.out), a.procs, md5=not a.no_md5, signal=sig)


if __name__ == "__main__":
    main()
