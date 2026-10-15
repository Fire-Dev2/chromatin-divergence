"""Sanity tests for the annotation tables (require data/annot built by scripts/build_annotations.py)."""

import math

import numpy as np
import pandas as pd
import pytest

from cdtest import annot as A

from cdtest import intermediates  # noqa: E402

pytestmark = pytest.mark.skipif(not A.GENES_TSV.exists() or intermediates.is_deposited("data/annot/genes_gencode_v50.tsv"),
                                reason="full annotation tables not built (the deposited intermediates hold a gene subset)")


@pytest.fixture(scope="module")
def genes():
    return A.load_genes()


@pytest.fixture(scope="module")
def cn():
    return A.load_gene_cn()


@pytest.fixture(scope="module")
def lookup():
    return A.CNLookup()


def row(genes, name):
    r = genes[genes.gene_name == name]
    assert len(r) == 1, name
    return r.iloc[0]


def test_one_row_per_gene_and_scope(genes):
    assert genes.gene_id.is_unique
    assert not genes.gene_id.str.contains(r"\.").any()
    assert set(genes.chrom) == set(A.CHROMS)
    assert set(genes.gene_type) == {"protein_coding", "lncRNA"}
    n = genes.gene_type.value_counts()
    assert 19_000 < n["protein_coding"] < 20_500
    assert 15_000 < n["lncRNA"] < 40_000
    assert (genes.tss1 == genes.tss0 + 1).all()
    assert set(genes.arm) == {"p", "q"}


def test_tp53_tss(genes):
    r = row(genes, "TP53")
    assert r.chrom == "chr17" and r.strand == "-"
    assert 7_687_400 <= r.tss1 <= 7_687_600
    assert r.tss1 == r.gene_end1 or r.tss1 <= r.gene_end1
    assert r.arm == "p"
    assert r.p53_fischer2017_core == False  # noqa: E712  (TP53 is not its own target in Table 1)


def test_tss_strand_convention(genes):
    r = row(genes, "GAPDH")  # + strand, chr12p13.31
    assert r.chrom == "chr12" and r.strand == "+" and r.arm == "p"
    assert abs(r.tss1 - 6_534_512) < 2_000


def test_cpg_promoter(genes):
    if "cpg_oe_tss500" not in genes:
        pytest.skip("FASTA not processed")
    assert genes.cpg_oe_tss500.notna().mean() > 0.99
    pc = genes[genes.gene_type == "protein_coding"]
    # housekeeping CpG-island promoters
    for g in ["GAPDH", "ACTB", "RPL13A"]:
        r = row(genes, g)
        assert r.cpg_oe_tss500 > 0.6 and r.gc_frac_tss500 > 0.55
    # bimodality: protein-coding promoters are CGI-rich vs lncRNA
    assert pc.cpg_oe_tss500.median() > genes[genes.gene_type == "lncRNA"].cpg_oe_tss500.median()
    # CGI-overlapping promoters have high o/e
    assert genes[genes.cgi_overlap_tss500].cpg_oe_tss500.median() > 0.6
    assert genes[~genes.cgi_overlap_tss500].cpg_oe_tss500.median() < 0.4
    assert genes.gc_frac_tss500.between(0, 1).all()


def test_p53_targets(genes):
    p = genes[genes.p53_fischer2017_core]
    assert len(p) >= 110
    for g in ["CDKN1A", "MDM2", "BAX", "GADD45A", "SESN1", "RRM2B"]:
        assert row(genes, g).p53_fischer2017_core


def test_k562_cdkn2a_loss(cn, genes, lookup):
    gid = row(genes, "CDKN2A").gene_id
    assert cn.loc[gid, "K562"] < -1.0  # homozygous 9p21 deletion
    # 9p arm overall reduced relative to genome
    p9 = genes[(genes.chrom == "chr9") & (genes.arm == "p")].gene_id
    assert cn.loc[p9, "K562"].median() < 0
    # region lookup gives the same picture over the CDKN2A/B locus
    v = lookup.region("K562", ["chr9"], [21_967_000], [22_010_000])[0]
    assert v < -1.0


def test_k562_bcr_abl_gain(cn, genes, lookup):
    # BCR-ABL1 amplicon: 9q34.12 ABL1->NUP214 and 22q11.21-q11.23 (upstream of BCR) amplified in K562.
    # (BCR's own TSS sits in a DepMap segment gap, so test the flanks.)
    assert cn.loc[row(genes, "ABL1").gene_id, "K562"] > 1.0
    assert cn.loc[row(genes, "NUP214").gene_id, "K562"] > 1.0
    assert lookup.region("K562", ["chr22"], [20_400_000], [21_900_000])[0] > 1.0


def test_rpe1_10q_gain(cn, genes):
    gain = math.log2(1.5)
    b = A.band_start("chr10", "q21.2")
    g10 = genes[genes.chrom == "chr10"].set_index("gene_id")
    dist = g10.index[g10.tss0 >= b]
    prox = g10.index[g10.tss0 < b]
    assert np.allclose(cn.loc[dist, "RPE1"], gain)
    assert np.allclose(cn.loc[prox, "RPE1"], 0)
    assert np.allclose(cn.loc[genes[genes.chrom != "chr10"].gene_id, "RPE1"], 0)
    # empirical DepMap clone consensus agrees
    pten = row(genes, "PTEN").gene_id  # 10q23.31
    assert abs(cn.loc[pten, "RPE1_depmap"] - gain) < 0.2
    assert cn.loc[dist, "RPE1_depmap"].median() > 0.4
    assert abs(cn.loc[genes[genes.chrom == "chr1"].gene_id, "RPE1_depmap"].median()) < 0.15


def test_absent_lines_zero(cn):
    for line in ["GM12878", "IMR90"]:
        assert (cn[line] == 0).all()
    pm = pd.read_csv(A.CN_PROFILE_MAP_TSV, sep="\t")
    assert set(pm[pm.status.str.startswith("absent")].line) == {"GM12878", "IMR90"}


def test_cn_columns_complete(cn, genes):
    assert list(cn.columns) == ["gene_name"] + A.CN_LINES
    assert len(cn) == len(genes)
    assert cn[A.CN_LINES].notna().all().all()
    # relative CN: genome-wide median near 0 for every DepMap line
    for line in ["K562", "HepG2", "Jurkat", "HCT116", "MCF7", "A549", "PC3", "PANC1"]:
        assert abs(cn[line].median()) < 0.3, line


def test_region_lookup_semantics():
    seg = pd.DataFrame({"line": "X", "chrom": "chr1", "start": [0, 100, 300], "end": [100, 200, 400],
                        "log2cn": [0.0, 1.0, -1.0]})
    lk = A.CNLookup(seg)
    np.testing.assert_allclose(lk.region("X", ["chr1"], [50], [150]), [0.5])
    np.testing.assert_allclose(lk.region("X", ["chr1"], [100], [200]), [1.0])
    # gap 200-300 filled from nearest segment (midpoint 220 -> segment ending at 200)
    np.testing.assert_allclose(lk.region("X", ["chr1"], [210], [230]), [1.0])
    assert np.isnan(lk.region("X", ["chr1"], [210], [230], fill=None)[0])
    # partial: 150-250 = 50bp of 1.0 + 50bp gap filled with nearest (1.0)
    np.testing.assert_allclose(lk.region("X", ["chr1"], [150], [250]), [1.0])
    # unknown chromosome -> 0
    np.testing.assert_allclose(lk.region("X", ["chr2"], [0], [10]), [0.0])
    with pytest.raises(ValueError):
        lk.region("X", ["chr1"], [10], [10])


def test_knockdown_groups(genes):
    kg = A.load_knockdown_groups()
    assert len(kg) == len(genes) and kg.ensembl_id.is_unique
    expect = {"RPL13A": "cytosolic_ribosome", "RPS6": "cytosolic_ribosome", "MRPL13": "mitochondrial_ribosome",
              "PSMA1": "proteasome", "SF3B1": "spliceosome", "POLR2A": "rna_pol_ii", "MED12": "mediator",
              "NUP93": "nuclear_pore", "GAPDH": "other", "POLR1A": "other"}
    g = kg.set_index("gene_name").group
    for k, v in expect.items():
        assert g[k] == v, k
    assert set(kg.group) == set(A.KNOCKDOWN_GROUP_DEF) | {"other"}


def test_blacklist_loaded():
    bl = A.load_blacklist()
    assert len(bl) > 500 and set(bl.reason) >= {"Low Mappability", "High Signal Region"}
