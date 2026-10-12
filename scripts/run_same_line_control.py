"""Same-line control (proposal 'Three controls'): no accessibility difference should
explain the response difference between two screens of K562.

M is K562 essential minus K562 genome-wide (results/aim1/K562_vs_K562gw). There is no
true accessibility difference within one line, so each of the six target-pair
accessibility differences (combined CPM) is tested against this same-line M with the
frozen main test (block permutation, x4 rule, arm effects, B = 9,999). Any rejection
beyond chance would mean the test or the features pick up something technical.
Planned control; run after Aim 4 because it was omitted from the Aim 4 job.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from cdtest import aim2, covariates as cov, crt, features

ROOT = Path(__file__).resolve().parents[1]
TAG = "K562_vs_K562gw"
PAIRS = [("K562", "Jurkat"), ("RPE1", "HepG2"), ("RPE1", "K562"), ("K562", "HepG2"),
         ("Jurkat", "RPE1"), ("Jurkat", "HepG2")]


def main():
    gs = pd.read_csv(ROOT / f"results/aim1/{TAG}/gene_outcome.tsv", sep="\t")
    genes = gs.gene_id.tolist()
    try:
        C = cov.build_covariates("K562", "K562gw", genes, TAG)
    except Exception as e:  # CN lookup has no 'K562gw' column: same line, so CN difference is 0
        print("build_covariates fallback:", repr(e))
        C = cov.build_covariates("K562", "K562", genes, TAG) if False else None
        raise
    Z = cov.design(C)
    P = dict(a="K562", b="K562gw", tag=TAG, genes=genes, M=gs.M.to_numpy(float), gs=gs, C=C, Z=Z,
             gen=cov.genome_info(genes), Q=crt.orthonormal_basis(Z),
             summary=json.loads((ROOT / f"results/aim1/{TAG}/summary.json").read_text()), X={})
    for a, b in PAIRS:
        x = features.pair_feature(a, b, "combined", signal="cpm").reindex(genes).to_numpy(float)
        ok = np.isfinite(x)
        if not ok.all():
            print(f"{a}-{b}: {np.sum(~ok)} genes without feature; dropping them is not supported here")
            x = np.where(ok, x, np.nanmean(x))
        P["X"][f"{a}_vs_{b}"] = x
    out = {}
    for lab in P["X"]:
        r = aim2.test_feature(P, lab, B=9999)
        m = r["main_block_x4"]
        out[lab] = dict(theta_hat=r["theta_hat"], k=r["k_x4"], p_two_sided=m.get("p_two_sided"),
                        p_one_sided=m.get("p_one_sided"), partial_r2_share=r["partial_r2_share_reproducible"])
        print(lab, {k: (round(v, 5) if isinstance(v, float) else v) for k, v in out[lab].items()}, flush=True)
    d = ROOT / "results/aim4"
    (d / "same_line_control.json").write_text(json.dumps(out, indent=1))
    ps = [v["p_two_sided"] for v in out.values()]
    lines = ["# Same-line control (planned; run 2026-10-12 after Aim 4)", "",
             "M = K562 essential − K562 genome-wide (same line, two screens). Features: the six target-pair",
             "combined CPM accessibility differences. Frozen main test (x4 blocks, arm effects, B = 9,999).", "",
             "| feature | θ̂ | k | p two-sided | pR²/rel |", "|---|---|---|---|---|"]
    for lab, v in out.items():
        lines.append(f"| {lab} | {v['theta_hat']:.5f} | {v['k']} | {v['p_two_sided']} | {100*v['partial_r2_share']:.3f}% |")
    lines += ["", f"Rejections at 5%: {sum(p is not None and p <= 0.05 for p in ps)} of {len(ps)}."]
    (d / "same_line_control.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
