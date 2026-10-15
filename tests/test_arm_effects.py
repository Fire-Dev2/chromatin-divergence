"""Plan v1.3: chromosome-arm fixed effects in f(C); arm_effects=False reproduces plan v1.2."""
import json

import numpy as np
import pandas as pd
import pytest

from cdtest import aim2, covariates as cov, crt


def _C(G=80, seed=0):
    rng = np.random.default_rng(seed)
    C = pd.DataFrame({c: rng.normal(size=G) for c in cov.COLUMNS},
                     index=pd.Index([f"G{i}" for i in range(G)], name="gene_id"))
    C["log2_meanA"] = rng.uniform(0, 5, G)
    C["log2_meanB"] = rng.uniform(0, 5, G)
    C["expr_diff"] = C.log2_meanA - C.log2_meanB
    arm = np.repeat(["chr1p", "chr1q", "chr2p", "chr10q", "chrXq"], G // 5)
    return C, arm


def test_arm_indicators_reference_dropped():
    A, keep = cov.arm_indicators(["chr2q", "chr1p", "chr10p", "chr1p", "chrXp", "chr1q"])
    # genome order: chr1p (reference, dropped), chr1q, chr2q, chr10p, chrXp
    assert keep == ["chr1q", "chr2q", "chr10p", "chrXp"]
    assert A.shape == (6, 4)
    assert A[1].sum() == 0 and A[3].sum() == 0          # reference-arm genes: all zeros
    np.testing.assert_array_equal(A.sum(1), [1, 0, 1, 0, 1, 1])


def test_design_v12_unchanged_and_v13_adds_arms():
    C, arm = _C()
    Z12 = cov.design(C, arm_effects=False)
    np.testing.assert_array_equal(Z12, crt.build_design(C[cov.COLUMNS], cov.SPLINE_COLS, df=cov.SPLINE_DF))
    Z13 = cov.design(C, arm_effects=True, arm=arm)
    assert Z13.shape == (len(C), Z12.shape[1] + 4)
    np.testing.assert_array_equal(Z13[:, :Z12.shape[1]], Z12)
    assert crt.orthonormal_basis(Z13).shape[1] == crt.orthonormal_basis(Z12).shape[1] + 4
    assert cov.ARM_EFFECTS_DEFAULT is True
    np.testing.assert_array_equal(cov.design(C, arm=arm), Z13)


def test_design_arm_lookup_from_annot():
    C, arm = _C()
    annot = pd.DataFrame({"gene_id": C.index, "chrom": [a[:-1] for a in arm], "arm": [a[-1] for a in arm],
                          "tss1": np.arange(len(C)) * 100 + 1})
    np.testing.assert_array_equal(cov.design(C, annot=annot), cov.design(C, arm=arm))
    with pytest.raises(ValueError):
        cov.design(C, arm=arm[:-1])


def test_v13_residuals_have_zero_arm_means():
    C, arm = _C()
    rng = np.random.default_rng(3)
    shift = pd.Series(arm).map({"chr1p": 0, "chr1q": 2, "chr2p": -1, "chr10q": 3, "chrXq": 1}).to_numpy()
    x = rng.normal(size=len(C)) + shift
    M = 0.3 * x + rng.normal(size=len(C)) + 2 * shift
    p13 = crt.prepare(M, x, Q=crt.orthonormal_basis(cov.design(C, arm=arm)))
    for a in np.unique(arm):
        assert abs(p13.r_x[arm == a].mean()) < 1e-10 and abs(p13.r_M[arm == a].mean()) < 1e-10
    # with arm terms an arm-level shift in both x and M no longer drives theta-hat
    p12 = crt.prepare(M, x, Q=crt.orthonormal_basis(cov.design(C, arm_effects=False)))
    assert abs(p13.theta_hat - 0.3) < abs(p12.theta_hat - 0.3)


def test_decoy_similarity_uses_v13_design():
    C, arm = _C()
    rng = np.random.default_rng(4)
    shift = pd.Series(arm).map({"chr1p": 0, "chr1q": 2, "chr2p": -1, "chr10q": 3, "chrXq": 1}).to_numpy()
    x = rng.normal(size=len(C)) + shift
    D = np.column_stack([rng.normal(size=len(C)) + shift, rng.normal(size=len(C))])
    s12 = crt.decoy_similarity(x, D, cov.design(C, arm_effects=False))
    s13 = crt.decoy_similarity(x, D, cov.design(C, arm=arm))
    assert abs(s12[0]) > 0.4 and abs(s13[0]) < abs(s12[0])  # shared arm shift no longer counts


V12 = aim2.ROOT / "results/aim2/RPE1_vs_K562_test.json"


@pytest.mark.skipif(not V12.exists() or not (aim2.AIM1 / "RPE1_vs_K562/gene_outcome.tsv").exists()
                    or not (aim2.ROOT / "results/aim1/pseudobulk/RPE1/genes.tsv").exists(),
                    reason="development-pair results not present (regenerate with cdtest.aim1)")
def test_v12_reproduced_on_dev_pair():
    """arm_effects=False reproduces the frozen v1.2 numbers (theta-hat, k, partial R2, main p)."""
    old = json.loads(V12.read_text())
    P = aim2.load_pair("RPE1", "K562", arm_effects=False)
    assert P["Z"].shape[1] == old["design_columns"]
    lab = "combined_cpm"
    o = old["features"][lab]
    x = P["X"][lab]
    prep = crt.prepare(P["M"], x, Q=P["Q"])
    g = P["gen"]
    order = crt.genome_order(g.chrom.to_numpy(), g.pos.to_numpy())
    k = crt.choose_k(prep.r_x[order], g.arm.to_numpy()[order], rule="x4")
    assert prep.theta_hat == pytest.approx(o["theta_hat"], abs=1e-14)
    assert k == o["k_x4"]
    assert crt.partial_r2(P["M"], x, prep=prep) == pytest.approx(o["partial_r2_raw"], rel=1e-12)
    r = crt.block_permutation_test(P["M"], x, chrom=g.chrom.to_numpy(), pos=g.pos.to_numpy(),
                                   arm=g.arm.to_numpy(), k=k, B=9999, rng=aim2.job_seed("aim2", P["tag"], lab),
                                   prep=prep, order=order)
    assert r.p == o["main_block_x4"]["p_one_sided"] and r.p_two_sided == o["main_block_x4"]["p_two_sided"]
    # and v1.3 differs (arm terms change k)
    P13 = aim2.load_pair("RPE1", "K562")
    assert P13["arm_effects"] and P13["Z"].shape[1] == old["design_columns"] + g.arm.nunique() - 1
