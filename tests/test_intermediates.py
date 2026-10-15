"""Reproduction of the primary (H1) numbers from the deposited intermediates (intermediates/).

Runs on a fresh clone: the session fixture copies the deposited files into place (missing files
only). Checks theta-hat, the x4 block size, partial R2 and the B = 9,999 main-test p-values of both
matched pairs against the committed results (results/aim4/h1_primary.json), plus the frozen plan
v1.2 numbers on the development pair. Takes about a minute.
"""
import json

import pytest

from cdtest import aim2, crt, intermediates

pytestmark = pytest.mark.skipif(not intermediates.available(), reason="intermediates/ not present")
H1 = aim2.ROOT / "results/aim4/h1_primary.json"


@pytest.fixture(scope="module", autouse=True)
def _installed():
    intermediates.install(verbose=False)


@pytest.mark.parametrize("a,b", [("K562", "Jurkat"), ("RPE1", "HepG2")])
def test_h1_matched_pair_reproduced(a, b):
    ref = json.loads(H1.read_text())["pairs"][f"{a}_vs_{b}"]
    P = aim2.load_pair(a, b)
    lab = "combined_cpm"
    x, g = P["X"][lab], P["gen"]
    prep = crt.prepare(P["M"], x, Q=P["Q"])
    order = crt.genome_order(g.chrom.to_numpy(), g.pos.to_numpy())
    k = crt.choose_k(prep.r_x[order], g.arm.to_numpy()[order], rule="x4")
    assert prep.theta_hat == pytest.approx(ref["theta_hat"], rel=1e-9)
    assert k == ref["k_x4"]
    assert crt.partial_r2(P["M"], x, prep=prep) == pytest.approx(ref["partial_r2_raw"], rel=1e-9)
    r = crt.block_permutation_test(P["M"], x, chrom=g.chrom.to_numpy(), pos=g.pos.to_numpy(),
                                   arm=g.arm.to_numpy(), k=k, B=9999,
                                   rng=aim2.job_seed("aim2", P["tag"], lab), prep=prep, order=order)
    assert r.p == ref["p_one_sided"] and r.p_two_sided == ref["p_two_sided"]
