"""EXPLORATORY (post hoc) RPE1 diagnostic for Aim 4 -- NOT part of the confirmatory analysis.

H1 was run once on the frozen pipeline (pipeline-v1.1, plan v1.3). This script asks whether the
negative accessibility-response association seen only in RPE1-containing pairs (RPE1_vs_K562,
Jurkat_vs_RPE1, RPE1_vs_HepG2) is an RPE1 processing artefact or biology. It changes NO primary
result file and no frozen code; it re-uses the frozen pieces (aim2.load_pair, covariates.design,
crt.block_permutation_test with the x4 rule, B = 9,999, aim2 seeds, test-inversion CI) and only swaps
inputs.

Stages:
  python scripts/aim4x_rpe1_diagnostic.py features   variant RPE1 feature matrices -> data/features_x/
  python scripts/aim4x_rpe1_diagnostic.py tests      all variant tests -> results/aim4_exploratory/rpe1_diagnostic_raw.json
  python scripts/aim4x_rpe1_diagnostic.py report     md + json + figure

Feature variants are rebuilt from the saved raw signal of the primary build
(data/features/region_signal_raw.parquet, promoter_raw in gene_features.parquet) with the frozen
feature steps (features.normalize_signal / distal_feature / combined_feature) and only the RPE1
column (and/or the RPE1 CN track) replaced; the rebuild is checked to reproduce the primary feature.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from cdtest import aim2, annot, crt, features
from cdtest import covariates as cov

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/aim4_exploratory"
FX = ROOT / "data/features_x"
RAW_JSON = OUT / "rpe1_diagnostic_raw.json"
PAIRS = (("RPE1", "K562"), ("Jurkat", "RPE1"), ("RPE1", "HepG2"))
SCREEN_LINES = ("K562", "RPE1", "Jurkat", "HepG2")
B = 9999
LABEL = "EXPLORATORY (post hoc)"
JURKAT_INS = 115_840_964          # Jurkat pooled insertions on chr1-22,X (build_info.json)
DS_TARGETS = {"jurkat_depth_116M": JURKAT_INS, "58M_insertions": 58_000_000}
DS_SEEDS = (1, 2, 3)


# =========================================================================== features
def _base():
    info = json.loads((ROOT / "data/features/build_info.json").read_text())
    lines = [l["line"] for l in info["lines"]]
    cn_lines = {l["line"]: l["cn_line"] for l in info["lines"]}
    reg = pd.read_parquet(ROOT / "data/features/region_signal_raw.parquet")
    gf = features.load_gene_features()
    genes = annot.load_genes()
    genes = genes[genes.chrom.isin(features.CHROMS)].reset_index(drop=True)
    praw = gf.pivot(index="gene_id", columns="line", values="promoter_raw")
    praw.columns = [str(c) for c in praw.columns]
    praw = praw.loc[genes.gene_id.values, lines]
    sizes = annot.chrom_sizes()
    ps = np.maximum(genes.tss0.values - features.PROMOTER_HALF, 0)
    pe = np.minimum(genes.tss0.values + features.PROMOTER_HALF, genes.chrom.map(sizes).values)
    rset = (reg.chrom.values, reg.start.values, reg.end.values)
    pset = (genes.chrom.values, ps, pe)
    return dict(lines=lines, cn_lines=cn_lines, reg=reg, raw=reg[lines].astype(np.float64).reset_index(drop=True),
                praw=praw.reset_index(drop=True), genes=genes, rset=rset, pset=pset)


def rebuild(base, raw, praw, cn_lines, lookup, cn_cache):
    """Frozen feature steps 2-5 from raw region / promoter signal. Returns {kind: genes x lines}."""
    lines = base["lines"]

    def cnm(which, iv):
        out = {}
        for ln in lines:
            key = (which, cn_lines[ln])
            if key not in cn_cache:
                cn_cache[key] = features.cn_matrix({"x": cn_lines[ln]}, *iv, lookup)["x"].to_numpy()
            out[ln] = cn_cache[key]
        return pd.DataFrame(out)

    s = features.normalize_signal(raw[lines], cnm("r", base["rset"]))
    prom = features.normalize_signal(praw[lines], cnm("p", base["pset"])).values
    if "_pairs" not in base:
        base["_pairs"] = features.distal_pairs(base["genes"], base["reg"])
    gi, ri, dd = base["_pairs"]
    dist = features.distal_feature(len(base["genes"]), gi, ri, dd, s[lines].values)
    comb, _ = features.combined_feature(prom, dist)
    idx = pd.Index(base["genes"].gene_id.values, name="gene_id")
    return {k: pd.DataFrame(v, index=idx, columns=lines) for k, v in
            (("promoter", prom), ("distal", dist), ("combined", comb))}


def _count(arrs_by_chrom: dict, iv) -> np.ndarray:
    c, s, e = (np.asarray(iv[0]).astype(str), np.asarray(iv[1], np.int64), np.asarray(iv[2], np.int64))
    out = np.zeros(len(c), np.int64)
    for ch, pos in arrs_by_chrom.items():
        m = np.flatnonzero(c == ch)
        if len(m):
            out[m] = features.count_in_intervals(pos, s[m], e[m])
    return out


def _load_ins(path):
    with np.load(path) as z:
        return {c: z[c] for c in features.CHROMS if c in z.files}


def stage_features():
    t0 = time.time()
    FX.mkdir(parents=True, exist_ok=True)
    base = _base()
    lookup = annot.CNLookup()
    cache = {}
    variants, meta = {}, {}
    cn0 = dict(base["cn_lines"])

    def add(name, raw_r=None, raw_p=None, cn_rpe1=None, desc=""):
        raw, praw = base["raw"].copy(), base["praw"].copy()
        if raw_r is not None:
            raw["RPE1"] = raw_r
            praw["RPE1"] = raw_p
        cnl = dict(cn0)
        if cn_rpe1 is not None:
            cnl["RPE1"] = cn_rpe1
        variants[name] = rebuild(base, raw, praw, cnl, lookup, cache)
        meta[name] = desc
        print(f"[features] {name} done ({time.time() - t0:.0f}s)", flush=True)

    add("primary_rebuild", desc="primary CPM feature rebuilt from saved raw (reproduction check)")
    gf = features.load_gene_features()
    chk = {}
    for kind in ("promoter", "distal", "combined"):
        ref = features.feature_matrix(kind).loc[variants["primary_rebuild"][kind].index]
        d = (variants["primary_rebuild"][kind][ref.columns] - ref).abs().to_numpy()
        chk[kind] = float(np.nanmax(d))
    print("[features] max |rebuild - primary|:", chk, flush=True)

    add("cn_rpe1_depmap", cn_rpe1="RPE1_depmap", desc="RPE1 CN track = DepMap RPE1-clone consensus (median of 3 WGS clones)")
    add("cn_rpe1_none", cn_rpe1="none", desc="no CN correction for RPE1")

    # replicates alone
    tot = {}
    for rep in ("rep1", "rep2"):
        ins = _load_ins(ROOT / f"data/atac_ins/RPE1__{rep}.npz")
        n = sum(len(v) for v in ins.values())
        tot[rep] = n
        add(f"rpe1_{rep}_only", _count(ins, base["rset"]) / n * 1e6, _count(ins, base["pset"]) / n * 1e6,
            desc=f"RPE1 CPM from {rep} only ({n:,} insertions)")
    # downsampled pooled insertions (binomial thinning of each insertion)
    pooled = {}
    for rep in ("rep1", "rep2"):
        for c, v in _load_ins(ROOT / f"data/atac_ins/RPE1__{rep}.npz").items():
            pooled.setdefault(c, []).append(v)
    pooled = {c: np.sort(np.concatenate(v)) for c, v in pooled.items()}
    ntot = sum(len(v) for v in pooled.values())
    for dname, target in DS_TARGETS.items():
        for seed in DS_SEEDS:
            rng = np.random.default_rng(np.random.SeedSequence([20261012, seed, target]))
            p = target / ntot
            sub = {c: v[rng.random(len(v)) < p] for c, v in pooled.items()}
            n = sum(len(v) for v in sub.values())
            add(f"rpe1_ds_{dname}_s{seed}", _count(sub, base["rset"]) / n * 1e6, _count(sub, base["pset"]) / n * 1e6,
                desc=f"RPE1 pooled insertions thinned to {n:,} (target {target:,}; seed {seed})")
    del pooled
    # FC bigWigs swapped into the CPM matrix (QN makes the marginal distributions identical)
    fc4 = str(ROOT / "data/atac/RPE1__4DNFIS4IDZNE.bw")
    r4 = features.mean_bigwig(fc4, *base["rset"])
    p4 = features.mean_bigwig(fc4, *base["pset"])
    add("rpe1_4dn_fc_hybrid", r4, p4, desc="RPE1 column = mean 4DN-processed FC bigWig (4DNFIS4IDZNE) over the same "
                                         "union regions / promoters; other lines CPM; same QN/CN/distal/combined steps")
    rfc = pd.read_parquet(ROOT / "data/features/region_signal_raw_fc.parquet")
    assert (rfc[["chrom", "start", "end"]].values == base["reg"][["chrom", "start", "end"]].values).all()
    gfc = features.load_gene_features(signal="fc")
    pfc = gfc[gfc.line == "RPE1"].set_index("gene_id").promoter_raw.loc[base["genes"].gene_id.values].to_numpy()
    add("rpe1_local_fc_hybrid", rfc["RPE1"].to_numpy(np.float64), pfc,
        desc="control for the 4DN hybrid: RPE1 column = local-pipeline FC bigWig, other lines CPM")

    # store the four screen lines of each variant
    rows = []
    for name, d in variants.items():
        for kind, m in d.items():
            sub = m[list(SCREEN_LINES)].copy()
            sub["variant"], sub["kind"] = name, kind
            rows.append(sub.reset_index())
    pd.concat(rows, ignore_index=True).to_parquet(FX / "rpe1_variants.parquet", index=False)
    # simple concordance of the RPE1 promoter/combined feature with the primary
    conc = {}
    prim = variants["primary_rebuild"]
    for name, d in variants.items():
        conc[name] = {k: float(np.corrcoef(prim[k]["RPE1"], d[k]["RPE1"])[0, 1]) for k in ("promoter", "combined")}
    (FX / "rpe1_variants_meta.json").write_text(json.dumps(
        {"label": LABEL, "desc": meta, "rebuild_max_abs_diff": chk, "rpe1_feature_corr_with_primary": conc,
         "rep_insertions": tot, "pooled_insertions": ntot}, indent=1))
    print(f"[features] written ({time.time() - t0:.0f}s)", flush=True)


def stage_fc4dn_build():
    """All-FC build with RPE1 = 4DN FC bigWig + 4DN conservative peaks (config row priority 2)."""
    cfg = pd.read_csv(features.CONFIG_TSV, sep="\t", comment="#", dtype=str, keep_default_na=False)
    cfg = cfg[~((cfg.line == "RPE1") & (cfg.priority == "1"))]
    FX.mkdir(parents=True, exist_ok=True)
    p = FX / "feature_lines_rpe1_4dn.tsv"
    cfg.to_csv(p, sep="\t", index=False)
    features.build(p, FX / "fc_4dn", procs=10, md5=False, signal="fc")


# =========================================================================== tests
_VAR = None


def variant_matrix(variant, kind):
    global _VAR
    if _VAR is None:
        _VAR = pd.read_parquet(FX / "rpe1_variants.parquet")
    v = _VAR[(_VAR.variant == variant) & (_VAR.kind == kind)].set_index("gene_id")
    return v[list(SCREEN_LINES)]


def cn_table(rpe1_col):
    cn = pd.read_csv(cov.CN_TSV, sep="\t")
    if rpe1_col == "none":
        cn["RPE1"] = 0.0
    elif rpe1_col != "RPE1":
        cn["RPE1"] = cn[rpe1_col]
    return cn


def run_test(P, x, *, mask=None, C=None, seed_label="combined_cpm", with_ci=True):
    """Main test (frozen settings) on pair P with feature x, optional gene mask / covariate table."""
    M, gen = P["M"], P["gen"]
    C = P["C"] if C is None else C
    mask = np.ones(len(M), bool) if mask is None else np.asarray(mask, bool)
    Mx, xx, g, Cx = M[mask], np.asarray(x, float)[mask], gen.iloc[np.flatnonzero(mask)], C.iloc[np.flatnonzero(mask)]
    Z = cov.design(Cx, arm_effects=True, arm=g.arm.to_numpy())
    Q = crt.orthonormal_basis(Z)
    prep = crt.prepare(Mx, xx, Q=Q)
    chrom, pos, arm = g.chrom.to_numpy(), g.pos.to_numpy(), g.arm.to_numpy()
    order = crt.genome_order(chrom, pos)
    k = crt.choose_k(prep.r_x[order], arm[order], rule="x4")
    rel = P["summary"]["reliability"]
    out = {"G": int(mask.sum()), "theta_hat": prep.theta_hat, "k_x4": int(k),
           "partial_r2": crt.partial_r2(Mx, xx, prep=prep), "reliability": rel,
           "pr2_rel": crt.partial_r2(Mx, xx, prep=prep) / rel,
           "sd_x_resid": float(prep.r_x.std()),
           "theta_std": prep.theta_hat * float(prep.r_x.std()) / float(prep.r_M.std()),
           "partial_corr": float(np.corrcoef(prep.r_x, prep.r_M)[0, 1])}
    if with_ci:
        rng = aim2.job_seed("aim2", P["tag"], seed_label)
        res = crt.block_permutation_test(Mx, xx, chrom=chrom, pos=pos, arm=arm, k=k, B=B, rng=rng,
                                         prep=prep, order=order)
        grid = crt.theta_grid(prep.theta_hat, float(np.std(xx)))
        flo, fhi, _ = aim2.ci_refined(res, grid)
        pg, pl = res.p, float(crt.pvalue(res.T, res.Tstar, "less"))
        out.update({"ci95": [float(flo), float(fhi)], "p_two_sided": res.p_two_sided,
                    "p_two_sided_equal_tailed": float(min(1.0, 2 * min(pg, pl))), "p_one_sided": res.p,
                    "null_sd": float(np.nanstd(res.Tstar))})
    return out


def _pair_x(m: pd.DataFrame, a, b, genes):
    return (m[a] - m[b]).reindex(genes).to_numpy(float)


def p53_sets(genes_names: pd.Series):
    fis = pd.read_csv(ROOT / "data/annot/fischer2017_p53_core_mapped.tsv", sep="\t")
    hall = json.loads((ROOT / "data/annot/msigdb/HALLMARK_P53_PATHWAY.v2025.1.Hs.json").read_text())
    hsym = set(hall["HALLMARK_P53_PATHWAY"]["geneSymbols"])
    fis_ids = set(fis.gene_id.dropna())
    hall_ids = set(genes_names.index[genes_names.isin(hsym)])
    return fis_ids, hall_ids, len(fis), len(hsym)


def stage_tests():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    meta = json.loads((FX / "rpe1_variants_meta.json").read_text())
    gname = features.load_gene_features().drop_duplicates("gene_id").set_index("gene_id").gene_name
    res = {"label": LABEL, "B": B, "pairs": {}, "feature_meta": meta}
    cn_full = pd.read_csv(cov.CN_TSV, sep="\t").set_index("gene_id")
    for a, b in PAIRS:
        P = aim2.load_pair(a, b, arm_effects=True)
        tag, genes = P["tag"], P["genes"]
        R = {}

        def log(name, r):
            R[name] = r
            print(f"[{tag}] {name}: theta={r['theta_hat']:.5g} ci={r.get('ci95')} p2={r.get('p_two_sided')} "
                  f"pr2/rel={r['pr2_rel']:.4%} G={r['G']} ({time.time() - t0:.0f}s)", flush=True)

        x0 = P["X"]["combined_cpm"]
        log("primary", run_test(P, x0))
        # ---- (1) copy number: features AND covariate cn_diff swapped together; also each part alone
        for v, col in (("cn_rpe1_depmap", "RPE1_depmap"), ("cn_rpe1_none", "none")):
            C2 = cov.build_covariates(a, b, genes, tag, cn=cn_table(col))
            xv = _pair_x(variant_matrix(v, "combined"), a, b, genes)
            log(f"{v}", run_test(P, xv, C=C2))
            log(f"{v}__features_only", run_test(P, xv))
            log(f"{v}__covariate_only", run_test(P, x0, C=C2))
        # ---- (2) accessibility source
        rb = _pair_x(variant_matrix("primary_rebuild", "combined"), a, b, genes)
        R["_rebuild_max_abs_diff_x"] = float(np.max(np.abs(rb - x0)))
        for v in ("rpe1_rep1_only", "rpe1_rep2_only", "rpe1_4dn_fc_hybrid", "rpe1_local_fc_hybrid"):
            log(v, run_test(P, _pair_x(variant_matrix(v, "combined"), a, b, genes)))
        for dname in DS_TARGETS:
            for s in DS_SEEDS:
                v = f"rpe1_ds_{dname}_s{s}"
                log(v, run_test(P, _pair_x(variant_matrix(v, "combined"), a, b, genes)))
        # all-FC builds: primary FC (local RPE1 FC) vs RPE1 = 4DN FC + 4DN peaks
        log("fc_all_lines_primary", run_test(P, P["X"]["combined_fc"], seed_label="combined_fc"))
        f4 = features.load_gene_features(FX / "fc_4dn/gene_features_fc.parquet")
        log("fc_all_lines_rpe1_4dn", run_test(P, features.pair_feature(a, b, "combined", f4).reindex(genes).to_numpy(float),
                                              seed_label="combined_fc"))
        # ---- (3) promoter / distal
        log("promoter", run_test(P, P["X"]["promoter_cpm"], seed_label="promoter_cpm"))
        log("distal", run_test(P, P["X"]["distal_cpm"], seed_label="distal_cpm"))
        # ---- (4) chromosome exclusions
        chrom = P["gen"].chrom.to_numpy()
        log("excl_chr10", run_test(P, x0, mask=chrom != "chr10"))
        log("excl_chrX", run_test(P, x0, mask=chrom != "chrX"))
        log("excl_chr10_and_chrX", run_test(P, x0, mask=(chrom != "chr10") & (chrom != "chrX")))
        # ---- (5) leave-one-chromosome-out jackknife (point estimates)
        chs = sorted(set(chrom), key=lambda c: (len(c), c))
        jk = {c: run_test(P, x0, mask=chrom != c, with_ci=False)["theta_hat"] for c in chs}
        th = np.array(list(jk.values()))
        n = len(th)
        R["loco_jackknife"] = {"theta_drop": jk, "n": n, "mean": float(th.mean()),
                               "jackknife_se": float(np.sqrt((n - 1) / n * ((th - th.mean()) ** 2).sum())),
                               "min": [min(jk, key=jk.get), float(th.min())], "max": [max(jk, key=jk.get), float(th.max())],
                               "n_chrom_positive": int((th > 0).sum())}
        # ---- (6) p53 genes and expression terciles
        fis, hall, nf, nh = p53_sets(gname)
        gset = pd.Index(genes)
        inF, inH = gset.isin(list(fis)), gset.isin(list(hall))
        R["p53_sets"] = {"fischer_116_listed": nf, "hallmark_symbols": nh, "fischer_in_pair": int(inF.sum()),
                         "hallmark_in_pair": int(inH.sum()), "union_in_pair": int((inF | inH).sum())}
        log("excl_fischer_p53", run_test(P, x0, mask=~inF))
        log("excl_hallmark_p53", run_test(P, x0, mask=~inH))
        log("excl_p53_union", run_test(P, x0, mask=~(inF | inH)))
        # share of the theta-hat numerator carried by the p53 genes (full model)
        prep = crt.prepare(P["M"], x0, Q=P["Q"])
        num = prep.r_x * prep.r_M
        R["p53_numerator_share"] = {"union": float(num[inF | inH].sum() / num.sum()),
                                    "expected_if_average": float((inF | inH).mean())}
        rp = "log2_meanA" if a == "RPE1" else "log2_meanB"
        e = P["C"][rp].to_numpy()
        q1, q2 = np.quantile(e, [1 / 3, 2 / 3])
        log("rpe1_expr_bottom_tercile", run_test(P, x0, mask=e <= q1))
        log("rpe1_expr_middle_tercile", run_test(P, x0, mask=(e > q1) & (e <= q2)))
        log("rpe1_expr_top_tercile", run_test(P, x0, mask=e > q2))
        # ---- (7) which side: each line's own feature (oriented as in x), residualised on C
        m = features.feature_matrix("combined")
        fa, fb = m[a].reindex(genes).to_numpy(float), m[b].reindex(genes).to_numpy(float)
        rpe_x = fa if a == "RPE1" else -fb
        par_x = -fb if a == "RPE1" else fa
        log("side_rpe1_alone", run_test(P, rpe_x))
        log("side_partner_alone", run_test(P, par_x))
        rM = P["M"] - P["Q"] @ (P["Q"].T @ P["M"])
        r1 = rpe_x - P["Q"] @ (P["Q"].T @ rpe_x)
        r2 = par_x - P["Q"] @ (P["Q"].T @ par_x)
        Xj = np.column_stack([r1, r2])
        coef, *_ = np.linalg.lstsq(Xj, rM, rcond=None)
        R["side_joint"] = {"coef_rpe1": float(coef[0]), "coef_partner": float(coef[1]),
                           "corr_resid_rpe1_partner": float(np.corrcoef(r1, r2)[0, 1]),
                           "note": "joint OLS of residualised M on both residualised (oriented) line features; "
                                   "descriptive (no permutation)"}
        res["pairs"][tag] = R
        RAW_JSON.write_text(json.dumps(res, indent=1, default=aim2._json_default))
    print(f"tests done {time.time() - t0:.0f}s", flush=True)


def stage_extra():
    """Follow-up controls (EXPLORATORY): (a) partner line's CPM swapped for its own FC (RPE1 kept CPM);
    (b) RPE1 CPM feature split into the part shared with RPE1 FC and the CPM-specific remainder;
    (c) the 'one line alone' decomposition on the three pairs WITHOUT RPE1 (is a negative own-line
    association RPE1-specific?); (d) insertion-based TSS enrichment of every line (same formula)."""
    t0 = time.time()
    base = _base()
    lookup, cache = annot.CNLookup(), {}
    rfc = pd.read_parquet(ROOT / "data/features/region_signal_raw_fc.parquet")
    gfc = features.load_gene_features(signal="fc")

    def fc_hybrid(lines_to_fc):
        raw, praw = base["raw"].copy(), base["praw"].copy()
        for ln in lines_to_fc:
            raw[ln] = rfc[ln].to_numpy(np.float64)
            praw[ln] = gfc[gfc.line == ln].set_index("gene_id").promoter_raw.loc[base["genes"].gene_id.values].to_numpy()
        return rebuild(base, raw, praw, base["cn_lines"], lookup, cache)["combined"]

    raw_res = json.loads(RAW_JSON.read_text())
    ex = {}
    m_cpm = features.feature_matrix("combined")
    m_rfc = variant_matrix("rpe1_local_fc_hybrid", "combined")
    for a, b in PAIRS:
        P = aim2.load_pair(a, b, arm_effects=True)
        tag, genes = P["tag"], P["genes"]
        partner = b if a == "RPE1" else a
        sgn = 1.0 if a == "RPE1" else -1.0
        R = {}
        R["partner_fc_hybrid"] = run_test(P, _pair_x(fc_hybrid([partner]), a, b, genes))
        R["both_fc_hybrid"] = run_test(P, _pair_x(fc_hybrid([partner, "RPE1"]), a, b, genes))
        # (b) RPE1 CPM = shared-with-FC + CPM-specific (OLS of the RPE1 CPM feature on the RPE1 FC
        # feature over the pair's genes); each part enters as the oriented RPE1-side feature
        fc_ = m_rfc["RPE1"].reindex(genes).to_numpy(float)
        cp_ = m_cpm["RPE1"].reindex(genes).to_numpy(float)
        X1 = np.column_stack([np.ones_like(fc_), fc_])
        beta, *_ = np.linalg.lstsq(X1, cp_, rcond=None)
        shared, spec = X1 @ beta, cp_ - X1 @ beta
        R["rpe1_cpm_shared_with_fc"] = run_test(P, sgn * shared)
        R["rpe1_cpm_specific_part"] = run_test(P, sgn * spec)
        R["rpe1_cpm_vs_fc_r"] = float(np.corrcoef(fc_, cp_)[0, 1])
        for k, v in R.items():
            if isinstance(v, dict):
                print(f"[{tag}] {k}: theta={v['theta_hat']:.5g} p2={v['p_two_sided']} ci={v['ci95']}", flush=True)
        ex[tag] = R
    # (c) one-line-alone decomposition in the pairs without RPE1 (context)
    for a, b in (("K562", "Jurkat"), ("K562", "HepG2"), ("Jurkat", "HepG2")):
        P = aim2.load_pair(a, b, arm_effects=True)
        genes = P["genes"]
        fa, fb = m_cpm[a].reindex(genes).to_numpy(float), m_cpm[b].reindex(genes).to_numpy(float)
        ex[P["tag"]] = {"primary": run_test(P, P["X"]["combined_cpm"]), f"side_{a}_alone": run_test(P, fa),
                        f"side_{b}_alone": run_test(P, -fb)}
        for k, v in ex[P["tag"]].items():
            print(f"[{P['tag']}] {k}: theta={v['theta_hat']:.5g} p2={v['p_two_sided']}", flush=True)
    # (d) TSS enrichment from insertions: mean insertions per bp in TSS +/-50 bp over mean per bp in
    # the flanks 1,900-2,000 bp each side; protein-coding genes on chr1-22,X; replicates pooled
    g = base["genes"]
    if "gene_type" in g.columns:
        g = g[g.gene_type == "protein_coding"]
    tss = {c: np.sort(g.tss0.values[g.chrom.values == c]) for c in features.CHROMS}
    tse = {}
    for p in sorted((ROOT / "data/atac_ins").glob("*__rep*.npz")):
        line = p.stem.split("__")[0]
        c_in = c_fl = 0
        with np.load(p) as z:
            for c in features.CHROMS:
                if c not in z.files or not len(tss[c]):
                    continue
                pos, t = z[c], tss[c]
                n = lambda s, e: int((np.searchsorted(pos, e) - np.searchsorted(pos, s)).sum())
                c_in += n(t - 50, t + 50)
                c_fl += n(t - 2000, t - 1900) + n(t + 1900, t + 2000)
        a_ = tse.setdefault(line, [0, 0])
        a_[0] += c_in; a_[1] += c_fl
    ex["tss_enrichment_insertions"] = {ln: (v[0] / 100) / (v[1] / 200) for ln, v in tse.items()}
    print("TSS enrichment:", {k: round(v, 2) for k, v in ex["tss_enrichment_insertions"].items()}, flush=True)
    raw_res["extra"] = ex
    RAW_JSON.write_text(json.dumps(raw_res, indent=1, default=aim2._json_default))
    print(f"extra done {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    st = sys.argv[1]
    if st not in ("joint", "sameline"):
        {"features": stage_features, "fc4dn": stage_fc4dn_build, "tests": stage_tests, "extra": stage_extra}[st]()


def stage_joint():
    """Descriptive joint OLS (residualised M on both residualised, oriented line features) for all six
    pairs, so the RPE1 pairs' conditional coefficients can be compared with the pairs without RPE1."""
    raw_res = json.loads(RAW_JSON.read_text())
    m = features.feature_matrix("combined")
    out = {}
    for a, b in (("RPE1", "K562"), ("Jurkat", "RPE1"), ("RPE1", "HepG2"), ("K562", "Jurkat"), ("K562", "HepG2"),
                 ("Jurkat", "HepG2")):
        P = aim2.load_pair(a, b, arm_effects=True)
        Q, genes = P["Q"], P["genes"]
        rM = P["M"] - Q @ (Q.T @ P["M"])
        xa, xb = m[a].reindex(genes).to_numpy(float), -m[b].reindex(genes).to_numpy(float)
        ra, rb = xa - Q @ (Q.T @ xa), xb - Q @ (Q.T @ xb)
        coef, *_ = np.linalg.lstsq(np.column_stack([ra, rb]), rM, rcond=None)
        out[P["tag"]] = {f"coef_{a}": float(coef[0]), f"coef_{b}": float(coef[1]),
                         "corr_resid_oriented": float(np.corrcoef(ra, rb)[0, 1])}
        print(P["tag"], out[P["tag"]], flush=True)
    raw_res.setdefault("extra", {})["joint_two_line_ols_all_pairs"] = out
    RAW_JSON.write_text(json.dumps(raw_res, indent=1, default=aim2._json_default))


if __name__ == "__main__" and sys.argv[1] == "joint":
    stage_joint()


def stage_sameline():
    """EXPLORATORY: link to the planned same-line control (results/aim4/same_line_control.md).
    M = K562 essential - K562 genome-wide (technical between-screen divergence, no biology difference).
    Test single-line features and the RPE1 CPM / FC decomposition against it (frozen main test)."""
    tag = "K562_vs_K562gw"
    gs = pd.read_csv(aim2.AIM1 / tag / "gene_outcome.tsv", sep="\t")
    genes = gs.gene_id.tolist()
    C = cov.build_covariates("K562", "K562gw", genes, tag)
    gen = cov.genome_info(genes)
    Z = cov.design(C, arm_effects=True, arm=gen.arm.to_numpy())
    P = dict(tag=tag, genes=genes, M=gs.M.to_numpy(float), C=C, gen=gen, Z=Z, Q=crt.orthonormal_basis(Z),
             summary=json.loads((aim2.AIM1 / tag / "summary.json").read_text()))
    m = features.feature_matrix("combined")
    m_rfc = variant_matrix("rpe1_local_fc_hybrid", "combined")
    m_4dn = variant_matrix("rpe1_4dn_fc_hybrid", "combined")
    f = lambda mm, ln: mm[ln].reindex(genes).to_numpy(float)
    cp_, fc_ = f(m, "RPE1"), f(m_rfc, "RPE1")
    X1 = np.column_stack([np.ones_like(fc_), fc_])
    beta, *_ = np.linalg.lstsq(X1, cp_, rcond=None)
    X = {"x_RPE1_vs_K562 (CPM; reproduces same_line_control)": f(m, "RPE1") - f(m, "K562"),
         "x_RPE1_vs_HepG2 (CPM; reproduces same_line_control)": f(m, "RPE1") - f(m, "HepG2"),
         "x_RPE1_vs_K562, RPE1 = local FC hybrid": f(m_rfc, "RPE1") - f(m_rfc, "K562"),
         "x_RPE1_vs_K562, RPE1 = 4DN FC hybrid": f(m_4dn, "RPE1") - f(m_4dn, "K562"),
         "x_RPE1_vs_HepG2, RPE1 = 4DN FC hybrid": f(m_4dn, "RPE1") - f(m_4dn, "HepG2"),
         "RPE1 CPM alone": cp_, "RPE1 local FC (hybrid) alone": fc_,
         "RPE1 CPM part shared with FC": X1 @ beta, "RPE1 CPM part NOT shared with FC": cp_ - X1 @ beta,
         "K562 CPM alone": f(m, "K562"), "HepG2 CPM alone": f(m, "HepG2"), "Jurkat CPM alone": f(m, "Jurkat")}
    sd_M = float((P["M"] - P["Q"] @ (P["Q"].T @ P["M"])).std())
    out = {"label": LABEL, "M": "K562 essential - K562 genome-wide (same line)", "G": len(genes), "tests": {}}
    for k, x in X.items():
        r = run_test(P, x, seed_label=k)
        out["tests"][k] = r
        print(f"{k}: theta={r['theta_hat']:.3g} p2={r['p_two_sided']} std={r['theta_std']:.3f} ci={r['ci95']}", flush=True)
    raw_res = json.loads(RAW_JSON.read_text())
    raw_res.setdefault("extra", {})["same_line_link"] = out
    RAW_JSON.write_text(json.dumps(raw_res, indent=1, default=aim2._json_default))


if __name__ == "__main__" and sys.argv[1] == "sameline":
    stage_sameline()
