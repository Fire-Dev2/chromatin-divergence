"""Covariate table C (plan §4), the group-mean variant (plan §2 item 6), and Aim 2/3 helpers."""
import numpy as np
import pandas as pd
import pytest

from cdtest import aim1, aim2, covariates as cov, crt, planted


# --------------------------------------------------------------------------- covariates
@pytest.fixture
def fake_pb(tmp_path):
    rng = np.random.default_rng(1)
    genes = [f"G{i}" for i in range(60)]
    for line, shift in (("RPE1", 0.0), ("K562", 0.5)):
        d = tmp_path / line
        d.mkdir()
        # extra gene not in the pair + shuffled order to test the lookup
        ids = genes[::-1] + ["EXTRA"]
        pd.DataFrame({"gene_id": ids, "gene_name": ids}).to_csv(d / "genes.tsv", sep="\t")
        X = rng.gamma(2.0, 1.0 + shift, size=(300, len(ids))).astype(np.float32)
        np.save(d / "ctrl_norm.npy", X)
    annot = pd.DataFrame({"gene_id": genes, "cpg_oe_tss500": rng.uniform(0.1, 1.0, len(genes)),
                          "chrom": ["chr1"] * 30 + ["chr2"] * 30, "tss1": np.arange(60) * 1000 + 1,
                          "arm": ["p"] * 15 + ["q"] * 15 + ["p"] * 15 + ["q"] * 15})
    cn = pd.DataFrame({"gene_id": genes, "K562": rng.normal(0, .3, 60), "RPE1": np.r_[np.zeros(40), np.full(20, .585)],
                       "RPE1_depmap": rng.normal(0, .3, 60)})
    return tmp_path, genes, annot, cn


def test_covariates_values(fake_pb):
    pb, genes, annot, cn = fake_pb
    C = cov.build_covariates("RPE1", "K562", genes, pb_dir=pb, annot=annot, cn=cn)
    assert list(C.columns) == cov.COLUMNS and list(C.index) == genes
    gt = pd.read_csv(pb / "RPE1" / "genes.tsv", sep="\t", index_col=0)
    X = np.load(pb / "RPE1" / "ctrl_norm.npy").astype(np.float64)
    j = list(gt.gene_id).index("G7")
    assert C.loc["G7", "log2_meanA"] == pytest.approx(np.log2(X[:, j].mean()))
    assert C.loc["G7", "log2_varA"] == pytest.approx(np.log2(X[:, j].var(ddof=1)))
    np.testing.assert_allclose(C.expr_diff, C.log2_meanA - C.log2_meanB)
    # RPE1 uses the karyotype column, not RPE1_depmap
    np.testing.assert_allclose(C.cn_diff, cn.RPE1.values - cn.K562.values)
    np.testing.assert_allclose(C.cpg_oe, annot.cpg_oe_tss500.values)


def test_covariate_design(fake_pb):
    pb, genes, annot, cn = fake_pb
    C = cov.build_covariates("RPE1", "K562", genes, pb_dir=pb, annot=annot, cn=cn)
    Z = cov.design(C, arm_effects=False)   # plan v1.2 design
    # intercept + 2 x ns(df 4) + 5 linear
    assert Z.shape == (60, 1 + 4 + 4 + 5)
    # expr_diff lies in the span of the two splines (linear part) -> rank one less
    assert crt.orthonormal_basis(Z).shape[1] == Z.shape[1] - 1
    # every covariate is in col(Z)
    Q = crt.orthonormal_basis(Z)
    for c in cov.COLUMNS:
        v = C[c].to_numpy()
        assert np.linalg.norm(v - Q @ (Q.T @ v)) < 1e-8 * np.linalg.norm(v)


def test_covariates_reject_nonfinite(fake_pb):
    pb, genes, annot, cn = fake_pb
    annot = annot.copy()
    annot.loc[3, "cpg_oe_tss500"] = np.nan
    with pytest.raises(ValueError):
        cov.build_covariates("RPE1", "K562", genes, pb_dir=pb, annot=annot, cn=cn)


def test_genome_info(fake_pb):
    _, genes, annot, _ = fake_pb
    g = cov.genome_info(genes, annot)
    assert g.arm.iloc[0] == "chr1p" and g.arm.iloc[-1] == "chr2q"


# --------------------------------------------------------------------------- group-mean variant
def test_generic_subtract_group():
    rng = np.random.default_rng(0)
    beta = rng.normal(size=(12, 7))
    mask = np.zeros_like(beta, bool)
    mask[0, 0] = mask[5, 3] = True
    groups = np.array(["a"] * 4 + ["b"] * 5 + ["other"] * 3)
    bt = aim1.generic_subtract(beta, mask, groups)
    assert np.isnan(bt[0, 0]) and np.isnan(bt[5, 3])
    for gname in np.unique(groups):
        np.testing.assert_allclose(np.nanmean(bt[groups == gname], axis=0), 0, atol=1e-12)
    # a single group is the global variant
    np.testing.assert_allclose(aim1.generic_subtract(beta, mask, np.array(["x"] * 12)),
                               aim1.generic_subtract(beta, mask), equal_nan=True)
    # group means differ from the global mean
    assert not np.allclose(np.nan_to_num(bt), np.nan_to_num(aim1.generic_subtract(beta, mask)))


def test_fake_screen_group_single_group_equals_global():
    rng = np.random.default_rng(3)
    ctrl = rng.poisson(2.0, size=(200, 6)).astype(np.float32)
    n_kd = np.array([20.0, 25.0, 30.0, 22.0])
    mask = np.zeros((4, 6), bool)
    s_glob = aim1.fake_screen_s2(ctrl, n_kd, mask, np.random.default_rng(9), k=5)
    s_grp1 = aim1.fake_screen_s2(ctrl, n_kd, mask, np.random.default_rng(9), k=5, groups=np.array(["g"] * 4))
    np.testing.assert_allclose(s_glob, s_grp1)
    s_grp = aim1.fake_screen_s2(ctrl, n_kd, mask, np.random.default_rng(9), k=5,
                                groups=np.array(["a", "a", "b", "b"]))
    assert s_grp.shape == (4, 6) and np.all(s_grp >= 0)


def test_knockdown_groups_lookup(tmp_path):
    p = tmp_path / "kg.tsv"
    pd.DataFrame({"gene_name": ["RPL3", "PSMA1", "NEWNAME"], "ensembl_id": ["E1", "E2", "E3"],
                  "group": ["cytosolic_ribosome", "proteasome", "spliceosome"]}).to_csv(p, sep="\t", index=False)
    kept = pd.DataFrame({"target": ["RPL3", "OLDNAME", "PSMA1", "XYZ"], "target_id": ["E1", "E3", None, None]})
    g = aim1.knockdown_groups(kept, p)
    assert g.tolist() == ["cytosolic_ribosome", "spliceosome", "proteasome", "other"]


# --------------------------------------------------------------------------- aim2 / planted helpers
def test_rho_bar_formula(tmp_path, monkeypatch):
    rng = np.random.default_rng(5)
    f = rng.normal(size=(500, 1))
    X = (np.exp(0.3 * f + rng.normal(size=(500, 8)))).astype(np.float32)
    d = tmp_path / "L"
    d.mkdir()
    ids = [f"G{i}" for i in range(8)]
    pd.DataFrame({"gene_id": ids}).to_csv(d / "genes.tsv", sep="\t")
    np.save(d / "ctrl_norm.npy", X)
    monkeypatch.setattr(aim2.screens, "OUT", tmp_path)
    R = np.corrcoef(X.astype(np.float64), rowvar=False)
    expect = (R.sum() - 8) / 56
    assert aim2.rho_bar_controls("L", ids) == pytest.approx(expect, rel=1e-6)


def test_partial_r2_many_matches_crt():
    rng = np.random.default_rng(2)
    G = 300
    Z = np.column_stack([np.ones(G), rng.normal(size=(G, 3))])
    Q = crt.orthonormal_basis(Z)
    x = rng.normal(size=G)
    Mb = rng.normal(size=(4, G)) + 0.2 * x
    got = aim2.partial_r2_many(Mb, x, Q)
    for i in range(4):
        assert got[i] == pytest.approx(crt.partial_r2(Mb[i], x, Z))


def test_planted_theta_gives_target_r2():
    rng = np.random.default_rng(4)
    G = 2000
    Z = np.column_stack([np.ones(G), rng.normal(size=G)])
    Q = crt.orthonormal_basis(Z)
    v = rng.normal(size=G)
    M = rng.normal(size=G)
    rv = v - Q @ (Q.T @ v)
    rM = M - Q @ (Q.T @ M)
    M = M - (rv @ rM) / (rv @ rv) * rv  # make theta0 = 0
    chrom = np.repeat(["chr1", "chr2"], G // 2)
    pos = np.tile(np.arange(G // 2), 2).astype(float)
    arm = chrom.copy()
    order = crt.genome_order(chrom, pos)
    row = planted._one((M, Q, chrom, pos, arm, Z[:, 1], np.zeros(G), 0.0, v, 10, 5, 0.01, 99, 1, order))
    assert row["theta0"] == pytest.approx(0, abs=1e-10)
    assert row["partial_r2"] == pytest.approx(0.01, rel=1e-8)
    assert row["theta_hat"] == pytest.approx(row["theta"], rel=1e-8)
