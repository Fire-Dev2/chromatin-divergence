"""Control-cell expression moments for each screen, without downloading the full screens.

The covariate table C (src/cdtest/covariates.py) needs, per gene and line, the mean and variance
(ddof = 1) of normalized expression over the line's kept control cells. The streaming pass in
cdtest.screens computes these from the full single-cell h5ad files (~95 GB in total). This script
computes exactly the same numbers by reading only the metadata and the control-cell rows of each
h5ad over HTTP range requests (about 4 GB in total), with the same cell filter and normalization
as cdtest.screens.run:

  keep    = mitopercent < 99th percentile among control cells, and UMI_count > 0
  norm    = raw * (median UMI of kept control cells) / UMI_count          (float32, as in screens)
  moments = mean and variance (ddof = 1) over kept control cells          (float64)

Output: results/aim1/control_moments/<line>.tsv.gz (gene_id, mean, var) and <line>.json. These
files are part of the deposited intermediates; cdtest.covariates.control_moments reads them when
the pseudobulk control matrix (results/aim1/pseudobulk/<line>/ctrl_norm.npy) is absent.

If a screen's h5ad is present locally (data/perturbseq/...), it is read from disk instead.

Usage:  python scripts/fetch_control_moments.py [K562 RPE1 Jurkat HepG2 K562gw]
"""
from __future__ import annotations

import concurrent.futures as cf
import io
import json
import os
import ssl
import sys
import time
import urllib.error
import urllib.request
from collections import OrderedDict
from http.client import IncompleteRead
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from cdtest import screens  # noqa: E402

OUT = ROOT / "results/aim1/control_moments"
GEO = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE264nnn/GSE264667/suppl/"
URLS = {
    "K562": "https://ndownloader.figshare.com/files/35773219",    # K562_essential_raw_singlecell_01.h5ad
    "RPE1": "https://ndownloader.figshare.com/files/35775606",    # rpe1_raw_singlecell_01.h5ad
    "Jurkat": GEO + "GSE264667_jurkat_raw_singlecell_01.h5ad",
    "HepG2": GEO + "GSE264667_hepg2_raw_singlecell_01.h5ad",
    "K562gw": "https://ndownloader.figshare.com/files/35775507",  # K562_gwps_raw_singlecell_01.h5ad
}
UA = {"User-Agent": "cdtest-fetch-control-moments/0.1"}


def _ssl_context() -> ssl.SSLContext:
    """python.org builds of Python on macOS ship without root certificates; use certifi or the
    system bundle so HTTPS works without running 'Install Certificates.command'."""
    try:
        import certifi
        return ssl.create_default_context(cafile=certifi.where())
    except ImportError:
        pass
    for f in ("/etc/ssl/cert.pem", "/opt/homebrew/etc/openssl@3/cert.pem", "/usr/local/etc/openssl@3/cert.pem",
              "/etc/ssl/certs/ca-certificates.crt"):
        if os.path.exists(f):
            return ssl.create_default_context(cafile=f)
    return ssl.create_default_context()


CTX = _ssl_context()
THREADS = 16
#: NCBI (GEO) allows only a few concurrent connections per client; more are refused
THREADS_BY_HOST = {"ftp.ncbi.nlm.nih.gov": 3}
MAX_RUN_BYTES = 4 << 20   # one range request covers at most ~4 MB of consecutive rows
MAX_GAP_ROWS = 2          # rows between two control cells that are fetched anyway to merge requests


class RangeFile(io.RawIOBase):
    """Read-only, seekable file over HTTP range requests with a small LRU block cache
    (enough for h5py to read the HDF5 metadata and obs/var columns)."""

    def __init__(self, url: str, block: int = 1 << 20, cache_blocks: int = 512):
        self.src, self.block, self.cache_blocks = url, block, cache_blocks
        self.cache: OrderedDict[int, bytes] = OrderedDict()
        self.pos = 0
        self.url = self._resolve()
        self.n_bytes = 0

    def _resolve(self) -> str:
        req = urllib.request.Request(self.src, headers={**UA, "Range": "bytes=0-0"})
        with urllib.request.urlopen(req, timeout=120, context=CTX) as r:
            cr = r.headers.get("Content-Range")
            if not cr:
                raise IOError(f"server does not support range requests: {self.src}")
            self.size = int(cr.rsplit("/", 1)[1])
            return r.geturl()   # figshare redirects to a signed S3 URL

    def fetch(self, start: int, end: int) -> bytes:
        """Bytes start..end inclusive, with retries (re-resolving an expired signed URL)."""
        last = None
        for attempt in range(12):
            try:
                req = urllib.request.Request(self.url, headers={**UA, "Range": f"bytes={start}-{end}"})
                with urllib.request.urlopen(req, timeout=300, context=CTX) as r:
                    data = r.read()
                if len(data) == end - start + 1:
                    self.n_bytes += len(data)
                    return data
            except urllib.error.HTTPError as e:
                last = f"HTTP {e.code}"
                if e.code in (400, 401, 403, 410):
                    try:
                        self.url = self._resolve()
                    except Exception:  # noqa: BLE001
                        pass
            except (urllib.error.URLError, TimeoutError, ConnectionError, IncompleteRead, OSError) as e:
                last = repr(e)
            time.sleep(min(60, 2 ** attempt))
        raise IOError(f"range {start}-{end} failed after 12 attempts ({last}): {self.src}")

    def _block(self, i: int) -> bytes:
        if i in self.cache:
            self.cache.move_to_end(i)
            return self.cache[i]
        s = i * self.block
        b = self.fetch(s, min(s + self.block, self.size) - 1)
        self.cache[i] = b
        if len(self.cache) > self.cache_blocks:
            self.cache.popitem(last=False)
        return b

    # io interface used by h5py's fileobj driver
    def readable(self):
        return True

    def seekable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, off, whence=0):
        self.pos = off if whence == 0 else self.pos + off if whence == 1 else self.size + off
        return self.pos

    def readinto(self, buf):
        n = min(len(buf), self.size - self.pos)
        if n <= 0:
            return 0
        out, p, end = memoryview(buf), self.pos, self.pos + n
        k = 0
        while p < end:
            i = p // self.block
            b = self._block(i)
            o = p - i * self.block
            take = min(len(b) - o, end - p)
            out[k:k + take] = b[o:o + take]
            k += take
            p += take
        self.pos = end
        return n


def _runs(rows: np.ndarray, row_bytes: int):
    """Group sorted row indices into (first, last) runs for range requests."""
    max_rows = max(1, MAX_RUN_BYTES // row_bytes)
    runs, a, b = [], rows[0], rows[0]
    for r in rows[1:]:
        if r - b <= MAX_GAP_ROWS + 1 and r - a < max_rows:
            b = r
        else:
            runs.append((a, b))
            a = b = r
    runs.append((a, b))
    return runs


def moments(line: str) -> dict:
    local = screens.SCREEN_FILES[line]
    src = local if local.exists() else URLS[line]
    t0 = time.time()
    rf = None if isinstance(src, Path) else RangeFile(src)
    with h5py.File(src if rf is None else rf, "r") as h:
        obs = pd.DataFrame({c: screens._read_col(h["obs"], c) for c in ("gene", "UMI_count", "mitopercent")})
        gene_ids = screens._read_col(h["var"], "gene_id").astype(str)
        X = h["X"]
        if not isinstance(X, h5py.Dataset):
            raise TypeError(f"{line}: X is not a dense dataset")
        n_cells, n_genes = X.shape
        is_ctrl = obs["gene"].values == screens.CONTROL
        mito = obs["mitopercent"].astype(float).values
        umi = obs["UMI_count"].astype(float).values
        mito_cut = np.percentile(mito[is_ctrl], 99)
        keep = (mito < mito_cut) & (umi > 0)
        med_umi = float(np.median(umi[keep & is_ctrl]))
        scale = (med_umi / umi).astype(np.float32)
        rows = np.flatnonzero(keep & is_ctrl)
        print(f"[{line}] {n_cells:,} cells x {n_genes:,} genes; {rows.size:,} kept control cells; "
              f"source: {'local file' if rf is None else 'HTTP ranges'}", flush=True)
        acc = {"n": 0, "mean": np.zeros(n_genes), "m2": np.zeros(n_genes)}

        def add(block, idx):
            # per-block count, mean and sum of squared deviations (combined exactly below)
            v = (block.astype(np.float32) * scale[idx, None]).astype(np.float64)
            mb = v.mean(0)
            return len(v), mb, ((v - mb) ** 2).sum(0)

        def combine(part):  # Chan et al. parallel update
            nb, mb, m2b = part
            n0 = acc["n"]
            d = mb - acc["mean"]
            acc["n"] = n0 + nb
            acc["mean"] += d * nb / acc["n"]
            acc["m2"] += m2b + d * d * n0 * nb / acc["n"]

        row_bytes = n_genes * X.dtype.itemsize
        offset = X.id.get_offset() if rf is not None else None
        compressed = X.compression is not None or (X.chunks is not None and X.chunks[1] != n_genes)
        if rf is not None and offset is not None and not compressed:
            # contiguous dense matrix: fetch control rows directly by byte offset, in parallel
            runs = _runs(rows, row_bytes)
            want = set(rows.tolist())

            def get(run):
                a, b = run
                raw = rf.fetch(offset + a * row_bytes, offset + (b + 1) * row_bytes - 1)
                blk = np.frombuffer(raw, dtype=X.dtype).reshape(b - a + 1, n_genes)
                idx = np.array([r for r in range(a, b + 1) if r in want])
                return add(blk[idx - a], idx)

            done = 0
            host = urllib.request.urlparse(rf.url).hostname
            with cf.ThreadPoolExecutor(THREADS_BY_HOST.get(host, THREADS)) as ex:
                for part in ex.map(get, runs):
                    combine(part)
                    done += 1
                    if done % 200 == 0:
                        print(f"[{line}]   {done:,}/{len(runs):,} requests, {rf.n_bytes / 1e9:.2f} GB, "
                              f"{time.time() - t0:.0f}s", flush=True)
        else:
            # local file, chunked or compressed layout: read row slices through h5py
            for a, b in _runs(rows, row_bytes):
                idx = rows[(rows >= a) & (rows <= b)]
                combine(add(X[a:b + 1][idx - a], idx))
    n = rows.size
    assert acc["n"] == n
    mean, var = acc["mean"], acc["m2"] / (n - 1)
    OUT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"gene_id": gene_ids, "mean": mean, "var": var}).to_csv(
        OUT / f"{line}.tsv.gz", sep="\t", index=False, float_format="%.17g")
    info = {"line": line, "source": str(src) if rf is None else URLS[line], "n_cells": int(n_cells),
            "n_genes": int(n_genes), "n_ctrl_kept": int(n), "median_ctrl_umi": med_umi,
            "mito_cut_pct": float(mito_cut), "bytes_downloaded": 0 if rf is None else rf.n_bytes,
            "seconds": round(time.time() - t0, 1)}
    (OUT / f"{line}.json").write_text(json.dumps(info, indent=1))
    print(f"[{line}] done: {n:,} control cells, {info['bytes_downloaded'] / 1e9:.2f} GB, "
          f"{info['seconds']:.0f}s", flush=True)
    return info


def check():
    """Means must agree with mean_ctrlA/B stored in results/aim1/<pair>/gene_outcome.tsv."""
    worst = 0.0
    for d in sorted((ROOT / "results/aim1").glob("*_vs_*")):
        if d.name.endswith("_groupmean") or not (d / "gene_outcome.tsv").exists():
            continue
        a, b = d.name.split("_vs_")
        gs = pd.read_csv(d / "gene_outcome.tsv", sep="\t", index_col="gene_id")
        for ln, col in ((a, "mean_ctrlA"), (b, "mean_ctrlB")):
            f = OUT / f"{ln}.tsv.gz"
            if not f.exists():
                continue
            m = pd.read_csv(f, sep="\t", index_col="gene_id")["mean"].reindex(gs.index)
            rel = float(np.nanmax(np.abs(m.to_numpy() / gs[col].to_numpy() - 1)))
            worst = max(worst, rel)
            print(f"check {d.name} {ln}: max relative difference of control means = {rel:.2e}")
    ok = worst < 1e-5
    print("CHECK OK" if ok else "CHECK FAILED", flush=True)
    return ok


if __name__ == "__main__":
    for ln in sys.argv[1:] or ["K562", "RPE1", "Jurkat", "HepG2", "K562gw"]:
        if (OUT / f"{ln}.tsv.gz").exists():
            print(f"[{ln}] already present, skipped")
            continue
        moments(ln)
    sys.exit(0 if check() else 1)
