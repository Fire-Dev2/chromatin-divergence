"""cdtest quick-start: why gene-level tests need dependence-aware nulls.

Simulates genes whose response M and feature x share correlated "modules"
(contiguous runs of co-regulated genes), with NO true effect of x on M.
A naive gene-level OLS test treats genes as independent and rejects too
often; the cdtest block-permutation test permutes x in genome-ordered
blocks and holds the nominal 5% level. Then plants a real effect and
reports the block test's power, estimate and confidence interval.

About 30 s on a laptop; needs no downloaded data.

    python examples/quickstart.py              # 500 replicates per scenario
    python examples/quickstart.py --reps 100   # faster, noisier
    python examples/quickstart.py --arm-terms  # add one indicator per chromosome arm to f(C)

The simulated genome has no arm-scale artefact (arm-level copy number is a covariate), so the
default design has no arm terms. With --arm-terms the residual feature loses its long-range
autocorrelation, the feature-only x4 rule picks shorter blocks, and the false-positive rate in
this strongly shared setting rises (about 9%); arm terms are meant for real data whose residual
feature carries arm-scale structure (see the paper and results/calibration/).
"""
import argparse

import numpy as np

from cdtest import crt
from cdtest.covariates import arm_indicators
from cdtest.simulate import SPLINE_COLS, Setting, simulate_replicate


def run_one(setting, rng, B, arm_terms=False):
    sim = simulate_replicate(setting, rng, n_decoys=0)
    M, x, C, gen = sim["M"], sim["x"], sim["C"], sim["genome"]
    Z = crt.build_design(C, SPLINE_COLS, df=4)          # covariates f(C)
    if arm_terms:                                       # one indicator per chromosome arm
        Z = np.column_stack([Z, arm_indicators(gen["arm"].to_numpy())[0]])
    prep = crt.prepare(M, x, Z)
    naive = crt.naive_ols_test(M, x, prep=prep)
    block = crt.block_permutation_test(
        M, x, chrom=gen["chrom"].to_numpy(), pos=gen["pos"].to_numpy(),
        arm=gen["arm"].to_numpy(), B=B, rng=rng, prep=prep, k_rule="x4")
    return sim, naive, block


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--reps", type=int, default=500, help="replicates per scenario")
    ap.add_argument("--B", type=int, default=999, help="permutations per test")
    ap.add_argument("--seed", type=int, default=20261011)
    ap.add_argument("--arm-terms", action="store_true", help="add chromosome-arm indicators to f(C)")
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)

    # 1. Null: response and feature share module structure, true effect = 0.
    null = Setting(rho_M=0.5, rho_x=0.6, confound=0.0, noise="homo", effect="zero")
    rej_naive = rej_block = 0
    ks = []
    for _ in range(a.reps):
        _, naive, block = run_one(null, rng, a.B, a.arm_terms)
        rej_naive += naive["p_two_sided"] < 0.05
        rej_block += block.p_two_sided < 0.05
        ks.append(block.info["k"])
    print(f"Null simulation ({a.reps} replicates, 2,000 genes, true effect = 0)")
    print(f"  false-positive rate at 5%   naive OLS: {rej_naive / a.reps:.1%}"
          f"   cdtest block permutation: {rej_block / a.reps:.1%}")
    print(f"  Monte Carlo SE of a 5% rate: {np.sqrt(0.05 * 0.95 / a.reps):.1%};"
          f" block size k: median {np.median(ks):.0f} (range {min(ks)}-{max(ks)})")

    # 2. Power: same structure with a planted effect (5% partial R^2).
    eff = Setting(rho_M=0.5, rho_x=0.6, confound=0.0, noise="homo", effect="large")
    hits, ests = 0, []
    for _ in range(a.reps):
        sim, _, block = run_one(eff, rng, a.B, a.arm_terms)
        hits += block.p_two_sided < 0.05
        ests.append(block.T)
    lo, hi, _ = block.ci(sd_x=float(np.std(sim["x"])))
    print(f"\nPlanted effect ({a.reps} replicates, true theta = {sim['theta']:.3f})")
    print(f"  cdtest power at 5%: {hits / a.reps:.1%}   mean estimate: {np.mean(ests):.3f}")
    print(f"  last replicate: estimate {block.T:.3f}, 95% CI [{lo:.3f}, {hi:.3f}],"
          f" block size k = {block.info['k']}")

if __name__ == "__main__":
    main()
