"""Unit tests for cdtest.features on small synthetic fixtures."""
import math

import numpy as np
import pandas as pd
import pytest

from cdtest import annot
from cdtest import features as F


def test_merge_and_union_blacklist():
    s1 = pd.DataFrame({"chrom": ["chr1", "chr1", "chr2"], "summit0": [1000, 1400, 5000]})
    s2 = pd.DataFrame({"chrom": ["chr1", "chr1"], "summit0": [1900, 10000]})
    bl = pd.DataFrame({"chrom": ["chr2"], "start": [5100], "end": [5200], "reason": ["x"]})
    reg, info = F.union_regions([s1, s2], bl)
    # chr1: [750,1250) [1150,1650) [1650,2150) -> one region (overlap + book-ended); [9750,10250)
    assert reg[["chrom", "start", "end"]].values.tolist() == [["chr1", 750, 2150], ["chr1", 9750, 10250]]
    assert info["n_blacklisted"] == 1 and info["n_merged"] == 3


def test_quantile_normalize_equalises_distributions():
    rng = np.random.default_rng(0)
    x = np.c_[rng.normal(0, 1, 500), rng.normal(3, 2, 500), rng.exponential(1, 500)]
    q = F.quantile_normalize(x)
    s = np.sort(q, axis=0)
    assert np.allclose(s[:, 0], s[:, 1]) and np.allclose(s[:, 0], s[:, 2])
    # ranks preserved
    for j in range(3):
        assert (np.argsort(q[:, j]) == np.argsort(x[:, j])).all()


def test_quantile_normalize_ties():
    x = np.array([[1.0, 5.0], [1.0, 6.0], [2.0, 7.0]])
    q = F.quantile_normalize(x)
    assert q[0, 0] == q[1, 0]
    ref = np.sort(x, axis=0).mean(axis=1)  # [3, 3.5, 4.5]
    assert q[0, 0] == pytest.approx((ref[0] + ref[1]) / 2)


def test_cn_subtraction():
    raw = pd.DataFrame({"A": [0.0, 1.0, 3.0], "B": [0.0, 1.0, 3.0]})
    cn = pd.DataFrame({"B": [0.5, 0.5, 0.0], "A": [0.0, 1.0, -1.0]})  # column order differs on purpose
    s = F.normalize_signal(raw, cn)
    base = np.log2(raw.values + 1)  # identical columns -> QN is the identity
    assert np.allclose(s["A"], base[:, 0] - cn["A"].values)
    assert np.allclose(s["B"], base[:, 1] - cn["B"].values)


def test_cn_matrix_uses_mapping_and_none():
    seg = pd.DataFrame({"line": ["MCF7", "MCF7"], "chrom": ["chr1", "chr1"], "start": [0, 1000],
                        "end": [1000, 2000], "log2cn": [1.0, -1.0]})
    lk = annot.CNLookup(segments=seg)
    m = F.cn_matrix({"MCF-7": "MCF7", "X": "none"}, ["chr1", "chr1"], [0, 500], [500, 1500], lk)
    assert m["MCF-7"].tolist() == [1.0, 0.0] and m["X"].tolist() == [0.0, 0.0]


def test_distal_pairs_window_and_weighting():
    genes = pd.DataFrame({"chrom": ["chr1", "chr2"], "tss0": [200_000, 50_000]})
    # centres: 199_000 (1 kb: too close), 202_500 (2.5 kb: in), 225_000 (25 kb), 300_000 (100 kb: in),
    # 300_250 (>100 kb: out); chr2 none
    regions = pd.DataFrame({"chrom": ["chr1"] * 5,
                            "start": [198_750, 202_250, 224_750, 299_750, 300_000],
                            "end": [199_250, 202_750, 225_250, 300_250, 300_500]})
    gi, ri, d = F.distal_pairs(genes, regions)
    assert sorted(ri.tolist()) == [1, 2, 3] and set(gi.tolist()) == {0}
    s = np.zeros((5, 2))
    s[:, 1] = 1.0  # 2**1 = 2x weight in line 2
    out = F.distal_feature(2, gi, ri, d, s)
    tot = math.exp(-2500 / 25000) + math.exp(-1) + math.exp(-4)
    assert out[0, 0] == pytest.approx(math.log2(1 + tot))
    assert out[0, 1] == pytest.approx(math.log2(1 + 2 * tot))
    assert out[1].tolist() == [0.0, 0.0]


def test_pooled_z_and_combined():
    prom = np.array([[1.0, 3.0], [2.0, 6.0]])
    dist = np.array([[0.0, 0.0], [0.0, 4.0]])
    z, mu, sd = F.pooled_z(prom)
    assert mu == pytest.approx(3.0) and sd == pytest.approx(np.std(prom))
    assert z.mean() == pytest.approx(0) and z.std() == pytest.approx(1)
    # pooled, not per-line: column means differ
    assert z[:, 0].mean() < 0 < z[:, 1].mean()
    comb, info = F.combined_feature(prom, dist)
    zd = (dist - dist.mean()) / dist.std()
    assert np.allclose(comb, 0.5 * z + 0.5 * zd)


def test_pair_feature():
    f = pd.DataFrame({"gene_id": ["g1", "g1", "g2", "g2"], "line": ["A", "B", "A", "B"],
                      "promoter": [3.0, 1.0, 0.0, 2.0], "distal": 0.0, "combined": 0.0})
    x = F.pair_feature("A", "B", "promoter", features=f)
    assert x.to_dict() == {"g1": 2.0, "g2": -2.0}
    with pytest.raises(KeyError):
        F.pair_feature("A", "Z", "promoter", features=f)


def test_read_summits_narrowpeak(tmp_path):
    p = tmp_path / "p.narrowPeak"
    p.write_text("chr1\t100\t300\t.\t0\t.\t1\t1\t1\t50\nchrUn_x\t0\t10\t.\t0\t.\t1\t1\t1\t5\n"
                 "chr2\t100\t300\t.\t0\t.\t1\t1\t1\t-1\n")
    s = F.read_summits(p)
    assert s.values.tolist() == [["chr1", 150], ["chr2", 200]]


def test_load_config_priority(tmp_path):
    (tmp_path / "a.bw").write_text("")
    (tmp_path / "a.bed").write_text("")
    cfg = tmp_path / "c.tsv"
    cfg.write_text("line\tpriority\trole\tfc_bigwig\tpeaks\tcn_line\tnote\n"
                   "L\t1\tt\tmissing.bw\tmissing.bed\tL\t\nL\t2\tt\ta.bw\ta.bed\tL\t\n"
                   "J\t1\tt\tnope.bw\tnope.bed\tJ\t\n")
    sel, sk = F.load_config(cfg, root=tmp_path, signal="fc")
    assert sel.line.tolist() == ["L"] and sel.fc_bigwig.iloc[0] == "a.bw"
    assert [s["line"] for s in sk] == ["J"]
    # no insertions column -> nothing usable for cpm
    sel, sk = F.load_config(cfg, root=tmp_path, signal="cpm")
    assert sel.empty and {s["line"] for s in sk} == {"L", "J"}


def test_load_config_signal_dependent(tmp_path):
    (tmp_path / "p.bed").write_text("")
    (tmp_path / "f.bw").write_text("")
    (tmp_path / "L__rep1.npz").write_text("")
    (tmp_path / "L__rep2.npz").write_text("")
    cfg = tmp_path / "c.tsv"
    cfg.write_text("line\tpriority\trole\tinsertions\tfc_bigwig\tpeaks\tcn_line\tnote\n"
                   "L\t1\tt\tL__rep*.npz\tnofc.bw\tp.bed\tL\t\n"
                   "L\t2\tt\tnone__rep*.npz\tf.bw\tp.bed\tL\t\n")
    sel, _ = F.load_config(cfg, root=tmp_path, signal="cpm")
    assert sel.priority.tolist() == [1]
    assert [f.name for f in F.insertion_files(sel.insertions.iloc[0], tmp_path)] == ["L__rep1.npz", "L__rep2.npz"]
    sel, _ = F.load_config(cfg, root=tmp_path, signal="fc")
    assert sel.priority.tolist() == [2]
    with pytest.raises(ValueError):
        F.load_config(cfg, root=tmp_path, signal="xx")


def test_output_paths():
    c, f = F.output_paths("cpm"), F.output_paths("fc")
    assert c["gene_features"].name == "gene_features.parquet" and f["gene_features"].name == "gene_features_fc.parquet"
    assert c["build_info"].name == "build_info.json" and f["build_info"].name == "build_info_fc.json"
    assert f["region_raw"].name == "region_signal_raw_fc.parquet"


def test_count_in_intervals_half_open():
    pos = np.array([5, 10, 10, 19, 20, 30], dtype=np.int32)
    # [10,20) holds 10,10,19 ; [0,5) none ; [5,6) one ; [20,31) holds 20,30 ; [100,200) none
    got = F.count_in_intervals(pos, [10, 0, 5, 20, 100], [20, 5, 6, 31, 200])
    assert got.tolist() == [3, 0, 1, 2, 0]


def _write_ins(path, d):
    np.savez_compressed(path, **{k: np.asarray(v, dtype=np.int32) for k, v in d.items()})


def test_count_insertions_and_cpm_pooling(tmp_path):
    r1, r2 = tmp_path / "L__rep1.npz", tmp_path / "L__rep2.npz"
    _write_ins(r1, {"chr1": [100, 150, 900], "chr2": [10], "chrY": [1, 2, 3]})  # chrY ignored
    _write_ins(r2, {"chr1": [120, 2000], "chrX": [5, 6, 7, 8, 9, 10]})
    regs = (["chr1", "chr2", "chrX", "chr3"], [100, 0, 0, 0], [200, 100, 8, 10])
    prom = (["chr1"], [0], [1000])
    (c1, c2), tot = F.count_insertions(r1, [regs, prom])
    assert tot == 4 and c1.tolist() == [2, 1, 0, 0] and c2.tolist() == [3]
    m1 = tmp_path / "M__rep1.npz"
    _write_ins(m1, {"chr1": [120, 2000], "chrX": [5, 6, 7, 8, 9, 10]})
    res = F.cpm_matrices({"L": [r1, r2], "M": [m1]}, [regs, prom], procs=1)
    assert res["totals"] == {"L": 12, "M": 8}
    # pooled: chr1 [100,200) = 2 + 1 = 3 of 12 insertions
    assert res["cpm"][0]["L"].tolist() == pytest.approx([3 / 12 * 1e6, 1 / 12 * 1e6, 3 / 12 * 1e6, 0])
    assert res["cpm"][1]["M"].tolist() == pytest.approx([1 / 8 * 1e6])
    rep = res["rep_cpm"][1]
    assert list(rep.columns) == ["L__rep1", "L__rep2", "M__rep1"]
    assert res["rep_totals"] == {"L__rep1": 4, "L__rep2": 8, "M__rep1": 8}
    with pytest.raises(ValueError):
        F.cpm_matrices({"L": [r1], "M": [r1]}, [prom], procs=1)
    assert rep["L__rep1"].iloc[0] == pytest.approx(3 / 4 * 1e6)
