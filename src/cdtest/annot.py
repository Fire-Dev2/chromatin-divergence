"""Gene-level and region-level annotations (ANALYSIS_PLAN.md sections 1, 3, 4).

Coordinate conventions
----------------------
* All interval coordinates handled by this module are **0-based, half-open**
  (BED style) unless a column name ends in ``1`` (``tss1``, ``gene_start1``,
  ``gene_end1``), which are 1-based, closed (GTF style).
* ``tss0`` is the 0-based position of the first transcribed base of the
  Ensembl-canonical transcript; ``tss1 = tss0 + 1``. For ``+`` strand genes
  ``tss1`` is the GTF transcript ``start``; for ``-`` strand it is the GTF
  transcript ``end``.
* The promoter window used for CpG o/e and GC is ``[tss0 - 500, tss0 + 500)``
  (1,000 bp, strand-independent).

Copy number
-----------
All CN values are **log2 relative copy number** (0 = the line's genome-wide
baseline). DepMap 24Q4 ``SegmentMean`` is a *linear* relative copy ratio
(1.0 = baseline; README: "relative copy ratio for that segment"; every value
in the file is > 0 and the median is 1.01), so we take log2(SegmentMean),
floored at ``CN_LOG2_FLOOR``. See docs/annotations.md.
"""

from __future__ import annotations

import gzip
import math
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
ANNOT = ROOT / "data" / "annot"

GTF_PATH = ANNOT / "gencode" / "gencode.v50.annotation.gtf.gz"  # comprehensive: basic lacks the canonical tx of 5,973 genes
CYTOBAND_PATH = ANNOT / "ucsc" / "hg38_cytoBand.txt.gz"
CPG_ISLAND_PATH = ANNOT / "ucsc" / "hg38_cpgIslandExt.txt.gz"
BLACKLIST_PATH = ANNOT / "blacklist" / "hg38-blacklist.v2.bed.gz"
# Project canonical genome (ENCODE GRCh38 no-alt analysis set). The shipped CpG/GC columns were computed
# from UCSC hg38.fa (md5 1c9dcaddfa41027f17cd8f7a82c7293b), since deleted; identical sequence on chr1-22,X.
# pyfaidx needs an uncompressed or bgzipped FASTA.
FASTA_PATH = ROOT / "data" / "genome" / "GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz"
DEPMAP_DIR = ANNOT / "depmap_24q4"

GENES_TSV = ANNOT / "genes_gencode_v50.tsv"
GENE_CN_TSV = ANNOT / "gene_cn_log2.tsv"
CN_SEGMENTS_TSV = ANNOT / "cn_segments_log2.tsv"
CN_PROFILE_MAP_TSV = ANNOT / "cn_profile_map.tsv"
KNOCKDOWN_GROUPS_TSV = ANNOT / "knockdown_groups.tsv"

CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
GENE_TYPES = ("protein_coding", "lncRNA")
PROMOTER_HALF = 500
CN_LOG2_FLOOR = -5.0  # log2 ratio 1/32; DepMap homozygous-deletion segments go to ~1e-11
RPE1_GAIN_LOG2 = math.log2(3 / 2)

# Lines: main four + seven decoys (ANALYSIS_PLAN §3). RRIDs used to find DepMap models.
LINE_RRID = {
    "K562": "CVCL_0004",
    "HepG2": "CVCL_0027",
    "Jurkat": "CVCL_0065",  # DepMap has parental JURKAT, not E6-1 (CVCL_0367)
    "GM12878": "CVCL_7526",
    "HCT116": "CVCL_0291",
    "MCF7": "CVCL_0031",
    "A549": "CVCL_0023",
    "IMR90": "CVCL_0347",
    "PC3": "CVCL_0035",
    "PANC1": "CVCL_0480",
}
# hTERT-RPE1 parental line is not in DepMap; use the three WGS single-cell clones
# (Cohen-Sharir 2021) as a consensus (median) profile.
RPE1_CLONE_PROFILES = {"RPE1-ss6": "PR-g5QF67", "RPE1-ss77": "PR-Yf6hXD", "RPE1-ss119": "PR-l4NGGv"}
# Output column order of gene_cn_log2.tsv / CN lookup line names.
CN_LINES = ["K562", "HepG2", "Jurkat", "RPE1", "RPE1_depmap",
            "GM12878", "HCT116", "MCF7", "A549", "IMR90", "PC3", "PANC1"]


# --------------------------------------------------------------------------- GTF
_ATTR_RE = re.compile(r'(\S+) "([^"]*)"')


def _attrs(s: str) -> dict:
    d: dict = {}
    for k, v in _ATTR_RE.findall(s):
        if k == "tag":
            d.setdefault("tags", []).append(v)
        else:
            d.setdefault(k, v)
    return d


def parse_gencode_canonical(gtf: Path = GTF_PATH, chroms=CHROMS, gene_types=GENE_TYPES) -> pd.DataFrame:
    """One row per gene: the Ensembl_canonical transcript of each gene.

    Restricted to ``chroms`` and ``gene_types``. Raises if a gene has != 1
    canonical transcript.
    """
    chroms = set(chroms)
    gene_types = set(gene_types)
    genes: dict = {}
    tx: list = []
    with gzip.open(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] not in chroms or f[2] not in ("gene", "transcript"):
                continue
            a = _attrs(f[8])
            if a.get("gene_type") not in gene_types:
                continue
            if f[2] == "gene":
                genes[a["gene_id"]] = (int(f[3]), int(f[4]))
            elif "Ensembl_canonical" in a.get("tags", []):
                tx.append((a["gene_id"], a["gene_name"], a["gene_type"], a.get("hgnc_id", ""), f[0], f[6],
                           int(f[3]), int(f[4]), a["transcript_id"], a.get("transcript_type", "")))
    df = pd.DataFrame(tx, columns=["gene_id_version", "gene_name", "gene_type", "hgnc_id", "chrom", "strand",
                                   "tx_start1", "tx_end1", "canonical_transcript_id", "canonical_transcript_type"])
    dup = df.gene_id_version.duplicated(keep=False)
    if dup.any():
        raise ValueError(f"{dup.sum()} rows for genes with >1 Ensembl_canonical transcript")
    missing = set(genes) - set(df.gene_id_version)
    df.attrs["n_genes_without_canonical"] = len(missing)
    df["gene_id"] = df.gene_id_version.str.replace(r"\.\d+(_PAR_Y)?$", "", regex=True)
    if df.gene_id.duplicated().any():
        raise ValueError("duplicated version-stripped gene ids")
    df["tss1"] = np.where(df.strand == "+", df.tx_start1, df.tx_end1)
    df["tss0"] = df.tss1 - 1
    df["gene_start1"] = df.gene_id_version.map(lambda g: genes[g][0])
    df["gene_end1"] = df.gene_id_version.map(lambda g: genes[g][1])
    df["chrom"] = pd.Categorical(df.chrom, categories=CHROMS, ordered=True)
    df = df.sort_values(["chrom", "tss0", "gene_id"]).reset_index(drop=True)
    df["chrom"] = df.chrom.astype(str)
    return df


# --------------------------------------------------------------------------- cytoband
def load_cytoband(path: Path = CYTOBAND_PATH) -> pd.DataFrame:
    cb = pd.read_csv(path, sep="\t", header=None, names=["chrom", "start", "end", "band", "stain"])
    return cb[cb.chrom.isin(CHROMS)].reset_index(drop=True)


def chrom_sizes(path: Path = CYTOBAND_PATH) -> dict:
    cb = load_cytoband(path)
    return cb.groupby("chrom").end.max().to_dict()


def annotate_cytoband(chrom: pd.Series, pos0: pd.Series, cb: pd.DataFrame | None = None) -> pd.DataFrame:
    """Band name and arm (p/q) of 0-based positions. Arm = first letter of the band."""
    cb = load_cytoband() if cb is None else cb
    band = np.empty(len(chrom), dtype=object)
    chrom = np.asarray(chrom)
    pos0 = np.asarray(pos0)
    for c, sub in cb.groupby("chrom"):
        m = chrom == c
        if not m.any():
            continue
        sub = sub.sort_values("start")
        i = np.searchsorted(sub.end.values, pos0[m], side="right")
        i = np.clip(i, 0, len(sub) - 1)
        band[m] = sub.band.values[i]
    band = pd.Series(band)
    return pd.DataFrame({"cytoband": band.values, "arm": band.str[0].values})


def band_start(chrom: str, band: str, cb: pd.DataFrame | None = None) -> int:
    """0-based start of the first sub-band whose name starts with ``band`` (e.g. 'q21.2')."""
    cb = load_cytoband() if cb is None else cb
    s = cb[(cb.chrom == chrom) & (cb.band.str.startswith(band))]
    if s.empty:
        raise KeyError(f"{chrom}{band} not in cytoBand")
    return int(s.start.min())


# --------------------------------------------------------------------------- sequence
def promoter_cpg(chrom, tss0, fasta: Path = FASTA_PATH, half: int = PROMOTER_HALF) -> pd.DataFrame:
    """CpG observed/expected and GC fraction over [tss0-half, tss0+half).

    CpG o/e = n_CpG * L / (n_C * n_G) (Gardiner-Garden & Frommer 1987), with
    L the number of non-N bases; GC fraction = (n_C + n_G) / L. Sequence is
    upper-cased (soft-masking ignored). NaN if n_C * n_G == 0.
    """
    from pyfaidx import Fasta

    fa = Fasta(str(fasta), sequence_always_upper=True, as_raw=True)
    out = np.full((len(tss0), 4), np.nan)
    for k, (c, t) in enumerate(zip(chrom, tss0)):
        s = max(0, int(t) - half)
        e = min(len(fa[c]), int(t) + half)
        seq = fa[c][s:e]
        n_c, n_g, n_n = seq.count("C"), seq.count("G"), seq.count("N")
        L = len(seq) - n_n
        n_cg = seq.count("CG")
        oe = n_cg * L / (n_c * n_g) if n_c * n_g > 0 else np.nan
        gc = (n_c + n_g) / L if L > 0 else np.nan
        out[k] = (oe, gc, n_cg, n_n)
    return pd.DataFrame(out, columns=["cpg_oe_tss500", "gc_frac_tss500", "n_cpg_tss500", "n_N_tss500"])


def overlaps_intervals(chrom, pos_start, pos_end, iv: pd.DataFrame) -> np.ndarray:
    """Boolean: does [pos_start, pos_end) overlap any interval in iv (chrom,start,end; 0-based)."""
    chrom = np.asarray(chrom)
    pos_start = np.asarray(pos_start)
    pos_end = np.asarray(pos_end)
    res = np.zeros(len(chrom), dtype=bool)
    for c, sub in iv.groupby("chrom"):
        m = chrom == c
        if not m.any():
            continue
        sub = sub.sort_values("start")
        st, en = sub.start.values, np.maximum.accumulate(sub.end.values)
        i = np.searchsorted(st, pos_end[m], side="left") - 1  # last interval starting before query end
        ok = i >= 0
        hit = np.zeros(m.sum(), dtype=bool)
        hit[ok] = en[i[ok]] > pos_start[m][ok]
        res[m] = hit
    return res


def load_cpg_islands(path: Path = CPG_ISLAND_PATH) -> pd.DataFrame:
    d = pd.read_csv(path, sep="\t", header=None, usecols=[1, 2, 3], names=["chrom", "start", "end"])
    return d[d.chrom.isin(CHROMS)].reset_index(drop=True)


def load_blacklist(path: Path = BLACKLIST_PATH) -> pd.DataFrame:
    """ENCODE blacklist v2 (hg38), 0-based half-open."""
    return pd.read_csv(path, sep="\t", header=None, names=["chrom", "start", "end", "reason"])


# --------------------------------------------------------------------------- copy number
def depmap_profile_map(depmap_dir: Path = DEPMAP_DIR) -> pd.DataFrame:
    """Line -> DepMap ModelID / default DNA ProfileID (24Q4) for LINE_RRID."""
    model = pd.read_csv(depmap_dir / "Model.csv", low_memory=False)
    dflt = pd.read_csv(depmap_dir / "OmicsDefaultModelProfiles.csv")
    prof = pd.read_csv(depmap_dir / "OmicsProfiles.csv")
    dflt = dflt[dflt.ProfileType == "dna"]
    rows = []
    for line, rrid in LINE_RRID.items():
        mm = model[model.RRID == rrid]
        if mm.empty:
            rows.append(dict(line=line, rrid=rrid, model_id="", cell_line_name="", profile_id="", datatype="",
                             source="", status="absent_from_DepMap_24Q4_CN0"))
            continue
        mid = mm.ModelID.iloc[0]
        pid = dflt[dflt.ModelID == mid].ProfileID
        pid = pid.iloc[0] if len(pid) else ""
        pr = prof[prof.ProfileID == pid]
        rows.append(dict(line=line, rrid=rrid, model_id=mid, cell_line_name=mm.CellLineName.iloc[0], profile_id=pid,
                         datatype=pr.Datatype.iloc[0] if len(pr) else "", source=pr.Source.iloc[0] if len(pr) else "",
                         status="DepMap_24Q4_default_dna_profile" if pid else "no_dna_profile_CN0"))
    for clone, pid in RPE1_CLONE_PROFILES.items():
        pr = prof[prof.ProfileID == pid]
        rows.append(dict(line="RPE1_depmap", rrid="", model_id=pr.ModelID.iloc[0], cell_line_name=clone,
                         profile_id=pid, datatype=pr.Datatype.iloc[0], source=pr.Source.iloc[0],
                         status="RPE1 clone; consensus=median of 3 WGS clones"))
    rows.append(dict(line="RPE1", rrid="CVCL_4388", model_id="", cell_line_name="hTERT RPE-1", profile_id="",
                     datatype="karyotype", source="Volpe 2025 Nat Commun",
                     status="karyotype track: +log2(3/2) chr10 q21.2->qter, else 0"))
    return pd.DataFrame(rows)


def _segments_from_depmap(seg: pd.DataFrame) -> pd.DataFrame:
    """DepMap segments (1-based closed Start/End, linear ratio) -> 0-based half-open log2."""
    s = pd.DataFrame({
        "chrom": "chr" + seg.Chromosome.astype(str),
        "start": seg.Start.astype(int) - 1,
        "end": seg.End.astype(int),
        "log2cn": np.maximum(np.log2(seg.SegmentMean.astype(float)), CN_LOG2_FLOOR),
    })
    s = s[s.chrom.isin(CHROMS)]
    return _sort_segments(s)


def _sort_segments(s: pd.DataFrame) -> pd.DataFrame:
    s = s.copy()
    s["_o"] = s.chrom.map({c: i for i, c in enumerate(CHROMS)})
    return s.sort_values(["_o", "start"]).drop(columns="_o").reset_index(drop=True)


def consensus_segments(segs: list[pd.DataFrame]) -> pd.DataFrame:
    """Median log2 CN over profiles on the union of breakpoints (positions covered by >=2 profiles)."""
    out = []
    for c in CHROMS:
        parts = [s[s.chrom == c] for s in segs]
        bps = np.unique(np.concatenate([np.r_[p.start.values, p.end.values] for p in parts if len(p)] or [[]]))
        if len(bps) < 2:
            continue
        st, en = bps[:-1], bps[1:]
        mid = (st + en) / 2
        vals = np.full((len(parts), len(st)), np.nan)
        for k, p in enumerate(parts):
            if p.empty:
                continue
            i = np.searchsorted(p.start.values, mid, side="right") - 1
            ok = (i >= 0) & (mid < p.end.values[np.clip(i, 0, None)])
            vals[k, ok] = p.log2cn.values[i[ok]]
        n = np.sum(~np.isnan(vals), axis=0)
        with np.errstate(all="ignore"), warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            med = np.nanmedian(np.where(n >= 2, vals, np.nan), axis=0)
        d = pd.DataFrame({"chrom": c, "start": st, "end": en, "log2cn": med})
        d = d[~d.log2cn.isna()]
        # merge adjacent pieces with identical value
        if len(d):
            key = ((d.log2cn.diff() != 0) | (d.start != d.end.shift())).cumsum()
            d = d.groupby(key).agg(chrom=("chrom", "first"), start=("start", "min"), end=("end", "max"),
                                   log2cn=("log2cn", "first"))
        out.append(d)
    return _sort_segments(pd.concat(out, ignore_index=True))


def rpe1_karyotype_segments(cb: pd.DataFrame | None = None) -> pd.DataFrame:
    """RPE1 karyotype CN: +log2(3/2) on chr10 from the 10q21.2 band start to qter, 0 elsewhere."""
    cb = load_cytoband() if cb is None else cb
    sizes = cb.groupby("chrom").end.max().to_dict()
    b = band_start("chr10", "q21.2", cb)
    rows = []
    for c in CHROMS:
        if c == "chr10":
            rows += [(c, 0, b, 0.0), (c, b, sizes[c], RPE1_GAIN_LOG2)]
        else:
            rows.append((c, 0, sizes[c], 0.0))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "log2cn"])


def zero_segments(cb: pd.DataFrame | None = None) -> pd.DataFrame:
    sizes = chrom_sizes() if cb is None else cb.groupby("chrom").end.max().to_dict()
    return pd.DataFrame([(c, 0, sizes[c], 0.0) for c in CHROMS], columns=["chrom", "start", "end", "log2cn"])


def build_cn_segments(depmap_dir: Path = DEPMAP_DIR) -> tuple[pd.DataFrame, pd.DataFrame]:
    """All lines' log2 relative CN segments (long table with `line` column) and the profile map."""
    pm = depmap_profile_map(depmap_dir)
    need = set(pm.profile_id[pm.profile_id != ""])
    chunks = []
    for ch in pd.read_csv(depmap_dir / "OmicsCNSegmentsProfile.csv", chunksize=200_000):
        chunks.append(ch[ch.ProfileID.isin(need)])
    seg = pd.concat(chunks)
    cb = load_cytoband()
    out = []
    for line in CN_LINES:
        if line == "RPE1":
            s = rpe1_karyotype_segments(cb)
        elif line == "RPE1_depmap":
            s = consensus_segments([_segments_from_depmap(seg[seg.ProfileID == p])
                                    for p in RPE1_CLONE_PROFILES.values()])
        else:
            pid = pm.loc[pm.line == line, "profile_id"].iloc[0]
            s = _segments_from_depmap(seg[seg.ProfileID == pid]) if pid else zero_segments(cb)
        out.append(s.assign(line=line))
    return pd.concat(out, ignore_index=True)[["line", "chrom", "start", "end", "log2cn"]], pm


class CNLookup:
    """Region-level log2 relative CN lookup per line.

    >>> cn = CNLookup()
    >>> cn.region("K562", ["chr9"], [21_967_000], [21_996_000])

    ``region`` returns the bp-weighted mean of segment log2 CN over each
    0-based half-open interval. Bases not covered by any segment (gaps
    between WES/WGS segments, telomeres) take the value of the nearest
    segment on the same chromosome (``fill='nearest'``) or are ignored
    (``fill=None``; fully-uncovered regions become NaN). Chromosomes with no
    segment return 0.
    """

    def __init__(self, segments: pd.DataFrame | None = None, path: Path = CN_SEGMENTS_TSV):
        if segments is None:
            segments = pd.read_csv(path, sep="\t")
        self.segments = segments
        self._idx = {}
        for (line, c), s in segments.groupby(["line", "chrom"]):
            s = s.sort_values("start")
            self._idx[(line, c)] = (s.start.values.astype(np.int64), s.end.values.astype(np.int64),
                                    s.log2cn.values.astype(float))

    @property
    def lines(self) -> list[str]:
        return list(dict.fromkeys(self.segments.line))

    def region(self, line: str, chrom, start, end, fill: str | None = "nearest") -> np.ndarray:
        chrom = np.asarray(chrom).astype(str)
        start = np.asarray(start, dtype=np.int64)
        end = np.asarray(end, dtype=np.int64)
        if np.any(end <= start):
            raise ValueError("end must be > start (0-based half-open)")
        out = np.full(len(chrom), np.nan)
        for c in np.unique(chrom):
            m = np.where(chrom == c)[0]
            key = (line, c)
            if key not in self._idx:
                if line not in self.lines:
                    raise KeyError(f"unknown line {line}")
                out[m] = 0.0
                continue
            st, en, v = self._idx[key]
            qs, qe = start[m], end[m]
            i0 = np.searchsorted(en, qs, side="right")  # first segment ending after qs
            i1 = np.searchsorted(st, qe, side="left")  # segments with index < i1 start before qe
            res = np.full(len(m), np.nan)
            single = (i1 - i0 == 1)
            if single.any():
                j = i0[single]
                cov = np.minimum(en[j], qe[single]) - np.maximum(st[j], qs[single])
                full = cov == (qe[single] - qs[single])
                r = np.where(full | (fill is None), v[j], np.nan)
                res[np.where(single)[0]] = r
            todo = np.where(np.isnan(res))[0]
            for k in todo:
                a, b = i0[k], i1[k]
                L = qe[k] - qs[k]
                if b > a:
                    ov = np.minimum(en[a:b], qe[k]) - np.maximum(st[a:b], qs[k])
                    w = np.clip(ov, 0, None)
                    tot = (w * v[a:b]).sum()
                    covered = w.sum()
                else:
                    tot, covered = 0.0, 0
                if fill == "nearest" and covered < L:
                    # nearest segment to the region midpoint for uncovered bases
                    mid = (qs[k] + qe[k]) / 2
                    d = np.maximum(st - mid, 0) + np.maximum(mid - en, 0)
                    tot += (L - covered) * v[int(np.argmin(d))]
                    covered = L
                res[k] = tot / covered if covered > 0 else np.nan
            out[m] = res
        return out

    def point(self, line: str, chrom, pos0, fill: str | None = "nearest") -> np.ndarray:
        pos0 = np.asarray(pos0, dtype=np.int64)
        return self.region(line, chrom, pos0, pos0 + 1, fill=fill)

    def covered(self, line: str, chrom, pos0) -> np.ndarray:
        """True if pos0 lies inside a segment (no gap filling needed)."""
        return ~np.isnan(self.point(line, chrom, pos0, fill=None))

    def all_lines(self, chrom, start, end, lines=None) -> pd.DataFrame:
        lines = self.lines if lines is None else lines
        return pd.DataFrame({ln: self.region(ln, chrom, start, end) for ln in lines})


def region_cn(chrom, start, end, lines=None) -> pd.DataFrame:
    """Convenience: log2 relative CN for arbitrary 0-based intervals, columns = lines."""
    return CNLookup().all_lines(chrom, start, end, lines)


# --------------------------------------------------------------------------- p53 / groups
FISCHER_PATH = ANNOT / "fischer2017_oncogene_table1_116_p53_targets.tsv"


def load_hgnc(path: Path = ANNOT / "hgnc" / "hgnc_complete_set.txt") -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", low_memory=False, dtype=str)


def map_symbols(symbols, genes: pd.DataFrame, hgnc: pd.DataFrame | None = None) -> pd.DataFrame:
    """Map symbols to gene-table Ensembl ids: exact gene_name, then HGNC current
    symbol (via prev/alias) -> HGNC ensembl_gene_id / gene_name."""
    hgnc = load_hgnc() if hgnc is None else hgnc
    by_name = dict(zip(genes.gene_name, genes.gene_id))
    ids = set(genes.gene_id)
    cur = dict(zip(hgnc.symbol, hgnc.ensembl_gene_id))
    alias: dict = {}
    for col in ("prev_symbol", "alias_symbol"):
        for sym, al in zip(hgnc.symbol, hgnc[col].fillna("")):
            for a in al.split("|"):
                if a:
                    alias.setdefault((col, a), []).append(sym)
    rows = []
    for s in symbols:
        if s in by_name:
            rows.append((s, by_name[s], s, "gene_name"))
            continue
        hit = None
        for col in ("prev_symbol", "alias_symbol"):
            cands = alias.get((col, s), [])
            if len(cands) == 1:
                new = cands[0]
                eid = by_name.get(new) or (cur.get(new) if cur.get(new) in ids else None)
                if eid:
                    hit = (s, eid, new, f"hgnc_{col}")
                    break
        if hit is None and s in cur and cur[s] in ids:
            hit = (s, cur[s], s, "hgnc_ensembl")
        rows.append(hit or (s, None, None, "unmapped"))
    return pd.DataFrame(rows, columns=["symbol", "gene_id", "mapped_symbol", "method"])


def load_fischer2017(path: Path = FISCHER_PATH) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", comment="#")


# HGNC gene-group definitions (family ids incl. descendants) — ANALYSIS_PLAN §2 item 6.
KNOCKDOWN_GROUP_DEF = {
    "cytosolic_ribosome": dict(hgnc_families=[728, 729]),
    "mitochondrial_ribosome": dict(hgnc_families=[646]),
    "proteasome": dict(hgnc_families=[690]),
    "spliceosome": dict(hgnc_families=[1518, 1520]),
    "rna_pol_ii": dict(hgnc_families=[726], symbol_regex=r"^POLR2"),
    "mediator": dict(hgnc_families=[1061]),
    "nuclear_pore": dict(hgnc_families=[1051]),
}
# CORUM complex-name patterns used only to break ties between groups.
CORUM_GROUP_REGEX = {
    "cytosolic_ribosome": r"(?i)^(?!.*mitochond).*(ribosom|40S|60S|80S)",
    "mitochondrial_ribosome": r"(?i)(mitochondrial ribosom|28S ribosom|39S ribosom|55S ribosom)",
    "proteasome": r"(?i)proteasom",
    "spliceosome": r"(?i)(spliceosom|snRNP|PRP19|Prp19|SF3[ab]|U1 |U2 |U4|U5|U6|U11|U12|NTC|Bact|C complex)",
    "rna_pol_ii": r"(?i)(RNA polymerase II|Pol II|POLR2)",
    "mediator": r"(?i)mediator",
    "nuclear_pore": r"(?i)(nuclear pore|NPC|Nup)",
}
GROUP_PRIORITY = list(KNOCKDOWN_GROUP_DEF)


def _hgnc_descendants(fid: int, hier: pd.DataFrame) -> set:
    out, fr = {fid}, [fid]
    while fr:
        ch = [c for c in hier[hier.parent_fam_id.isin(fr)].child_fam_id if c not in out]
        out |= set(ch)
        fr = ch
    return out


def build_knockdown_groups(genes: pd.DataFrame, hgnc_dir: Path = ANNOT / "hgnc",
                           corum_path: Path = ANNOT / "corum" / "humanComplexes_corum5.3.txt") -> pd.DataFrame:
    """Assign every gene in ``genes`` to one knockdown group (default 'other').

    Membership: HGNC gene groups (with descendants). Genes in >1 group are
    assigned to the group with the most CORUM 5.3 complexes containing the gene
    whose name matches the group pattern; remaining ties by GROUP_PRIORITY.
    """
    hier = pd.read_csv(hgnc_dir / "hierarchy.csv")
    ghf = pd.read_csv(hgnc_dir / "gene_has_family.csv")
    hg = genes[genes.hgnc_id != ""].assign(hid=lambda d: d.hgnc_id.str.replace("HGNC:", "").astype(int))
    hid2gene = dict(zip(hg.hid, hg.gene_id))
    member: dict = {}
    for grp, d in KNOCKDOWN_GROUP_DEF.items():
        fams = set().union(*[_hgnc_descendants(f, hier) for f in d["hgnc_families"]])
        gids = {hid2gene[h] for h in ghf[ghf.family_id.isin(fams)].hgnc_id if h in hid2gene}
        if "symbol_regex" in d:
            names = genes.set_index("gene_id").gene_name
            gids = {g for g in gids if re.search(d["symbol_regex"], names[g])}
        for g in gids:
            member.setdefault(g, []).append(grp)
    corum = pd.read_csv(corum_path, sep="\t", usecols=["complex_name", "subunits_gene_name"], dtype=str)
    name_of = genes.set_index("gene_id").gene_name
    rows = []
    for g in genes.gene_id:
        grps = member.get(g, [])
        tie = ""
        if len(grps) == 1:
            grp = grps[0]
        elif len(grps) > 1:
            sym = name_of[g]
            has = corum[corum.subunits_gene_name.fillna("").str.split(";").apply(lambda s: sym in s)]
            score = {k: int(has.complex_name.str.contains(CORUM_GROUP_REGEX[k], regex=True).sum()) for k in grps}
            best = max(score.values())
            cands = [k for k in GROUP_PRIORITY if k in grps and score[k] == best]
            grp = cands[0]
            tie = ";".join(f"{k}:{score[k]}" for k in grps) + (" (priority)" if len(cands) > 1 else " (CORUM)")
        else:
            grp = "other"
        rows.append((name_of[g], g, grp, ";".join(grps), tie))
    return pd.DataFrame(rows, columns=["gene_name", "ensembl_id", "group", "hgnc_groups_all", "tiebreak"])


# --------------------------------------------------------------------------- loaders
def load_genes(path: Path = GENES_TSV) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype={"hgnc_id": str}, keep_default_na=False,
                       na_values={"cpg_oe_tss500": [""], "gc_frac_tss500": [""]})


def load_gene_cn(path: Path = GENE_CN_TSV) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", index_col="gene_id")


def load_knockdown_groups(path: Path = KNOCKDOWN_GROUPS_TSV) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", keep_default_na=False)


def load_cn_segments(path: Path = CN_SEGMENTS_TSV) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")
