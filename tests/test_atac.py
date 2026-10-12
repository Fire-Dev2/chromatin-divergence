import numpy as np
import pandas as pd

from cdtest import atac


def test_tss_enrichment_flat_and_peak():
    flat = np.full(4001, 3.0)
    r = atac.tss_enrichment(flat)
    assert abs(r["tss_enrichment"] - 1.0) < 1e-9
    peak = flat.copy()
    peak[1990:2011] = 30.0  # 21-bp block at the TSS, 10x background
    r = atac.tss_enrichment(peak)
    assert 9.0 < r["tss_enrichment"] <= 10.0
    assert np.isfinite(r["tss_enrichment_bin10"])


def test_tss_enrichment_zero_flanks_is_nan_not_pass():
    prof = np.zeros(4001)
    prof[2000] = 5
    r = atac.tss_enrichment(prof)
    assert np.isnan(r["tss_enrichment"])


def test_tss_profile_strand_orientation(tmp_path):
    # one + TSS at 5000 and one - TSS at 20000 on chrT; insertions 50 bp downstream of each
    ta = tmp_path / "t.tagAlign"
    rows = [("chrT", 5050, 5080, "N", 1000, "+"),     # + read: insertion = start = 5050
            ("chrT", 19900, 19951, "N", 1000, "-")]   # - read: insertion = end-1 = 19950 (50 bp downstream of - TSS)
    pd.DataFrame(rows).to_csv(ta, sep="\t", header=False, index=False)
    tss = pd.DataFrame({"chrom": ["chrT", "chrT"], "tss0": [5000, 20000], "strand": ["+", "-"]})
    prof = atac.tss_profile([ta], tss, {"chrT": 100000})
    assert prof.sum() == 2 and prof[2000 + 50] == 2
