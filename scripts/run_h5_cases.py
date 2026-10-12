#!/usr/bin/env python
"""H5 mechanism case studies (ANALYSIS_PLAN.md §7 H5), K562_vs_HepG2 only.

Per sequence-specific TF that is a kept knockdown and has ENCODE IDR ChIP in K562 and/or HepG2:
  feature_g = (# ChIP peaks overlapping a K562-accessible union region, summit within 50 kb of TSS)
            - (same for HepG2), using each line's own ChIP where available;
  outcome   d_g = (btA^2 - s2A) - (btB^2 - s2B) for that knockdown row (A = K562);
  rho_own   = Spearman(d_g, feature) over genes (own gene excluded);
  swap      = 20 random other TFs' features (same list), rank of rho_own, empirical p.
No pooled test. Mixed pair: illustrative only.
"""
import json, subprocess, hashlib, os, sys, datetime
import numpy as np, pandas as pd
from scipy.stats import spearmanr

PAIR = "K562_vs_HepG2"; LINES = ["K562", "HepG2"]; WIN = 50_000; NSWAP = 20; SEED = 20261012
os.makedirs("data/chip", exist_ok=True); os.makedirs("results/aim4", exist_ok=True)

z = np.load(f"results/aim1/{PAIR}/terms.npz", allow_pickle=True)
targets = list(z["targets"]); genes = z["genes"]
kd = pd.read_csv(f"results/aim1/{PAIR}/knockdowns.tsv", sep="\t")
kept = set(kd.loc[kd.kept, "target"]) & set(targets)

# ---- (1) select ChIP experiments ----
RFA_RANK = {"ENCODE4": 0, "ENCODE3": 1, "ENCODE2": 2}
chip = {}
for line in LINES:
    d = pd.read_csv(f"recon/encode_tf_chip_{line}.csv")
    d = d[(d.is_sequence_specific_TF == True) & d.audit_error.isna() & d.idr_peak_file.notna()
          & d.target.isin(kept)].copy()
    d["rfa"] = d.award_rfa.map(RFA_RANK).fillna(9)
    d["nc"] = d.audit_not_compliant.notna().astype(int)
    d = d.sort_values(["target", "rfa", "nc", "n_bio_reps", "date_released"],
                      ascending=[True, True, True, False, False])
    chip[line] = {r.target: r for r in d.drop_duplicates("target").itertuples()}

tfs = sorted(set(chip["K562"]) | set(chip["HepG2"]))
print(len(tfs), "TFs", file=sys.stderr)

def curl_json(url):
    out = subprocess.run(["curl", "-sfL", "-H", "Accept: application/json", url],
                         capture_output=True, check=True).stdout
    return json.loads(out)

def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()

files = {}
for line in LINES:
    for tf, r in chip[line].items():
        acc = r.idr_peak_file; path = f"data/chip/{line}__{tf}__{acc}.bed.gz"
        meta = curl_json(f"https://www.encodeproject.org/files/{acc}/?format=json&frame=object")
        if not os.path.exists(path) or md5(path) != meta["md5sum"]:
            subprocess.run(["curl", "-sfL", "-o", path, r.idr_peak_url], check=True)
        ok = md5(path) == meta["md5sum"]
        if not ok: raise SystemExit(f"md5 mismatch {path}")
        files[(line, tf)] = dict(path=path, accession=acc, experiment=r.experiment, award=r.award_rfa,
                                 lab=r.lab, date=r.date_released, output_type=meta.get("output_type"),
                                 assembly=meta.get("assembly"), md5=meta["md5sum"], md5_ok=ok,
                                 not_compliant=None if pd.isna(r.audit_not_compliant) else r.audit_not_compliant)

# ---- regions, accessibility ----
reg = pd.read_parquet("data/features/region_signal.parquet")
acc = {l: (reg[l] > reg[l].median()).to_numpy() for l in LINES}
reg_by_chr = {}
for c, g in reg.groupby("chrom"):
    g = g.sort_values("start")
    reg_by_chr[c] = (g.start.to_numpy(), g.end.to_numpy(), g.index.to_numpy(),
                     np.maximum.accumulate(g.end.to_numpy()))

def accessible_summits(path, line):
    p = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1, 2, 9])
    p.columns = ["chrom", "start", "end", "summit"]
    out = {}
    for c, g in p.groupby("chrom"):
        if c not in reg_by_chr: continue
        rs, re_, ri, cummax = reg_by_chr[c]
        ps, pe = g.start.to_numpy(), g.end.to_numpy()
        hit = np.zeros(len(g), bool)
        # candidate regions: start < pe and end > ps
        hi = np.searchsorted(rs, pe, "left")
        for k in range(len(g)):
            j = hi[k] - 1
            while j >= 0 and cummax[j] > ps[k]:
                if re_[j] > ps[k] and acc[line][ri[j]]: hit[k] = True; break
                j -= 1
        summ = np.where(g.summit.to_numpy() >= 0, ps + g.summit.to_numpy(), (ps + pe) // 2)
        out[c] = np.sort(summ[hit])
    return out, len(p)

gt = pd.read_csv("data/annot/genes_gencode_v50.tsv", sep="\t").set_index("gene_id")
gt = gt.reindex(genes)
gchr, gtss = gt.chrom.to_numpy(), gt.tss0.to_numpy()

def count_near(summits):
    cnt = np.zeros(len(genes), int)
    for c in np.unique(gchr[pd.notna(gchr)]):
        idx = np.where(gchr == c)[0]; s = summits.get(c, np.array([]))
        t = gtss[idx].astype(np.int64)
        cnt[idx] = np.searchsorted(s, t + WIN, "right") - np.searchsorted(s, t - WIN, "left")
    return cnt

feat, info = {}, {}
for tf in tfs:
    srcs = {l: (l if (l, tf) in files else [o for o in LINES if (o, tf) in files][0]) for l in LINES}
    cnts, npk = {}, {}
    for l in LINES:
        s, n = accessible_summits(files[(srcs[l], tf)]["path"], l)
        cnts[l] = count_near(s); npk[l] = n
    feat[tf] = cnts["K562"] - cnts["HepG2"]
    src = "both" if srcs["K562"] == "K562" and srcs["HepG2"] == "HepG2" else f"{srcs['K562']} only"
    info[tf] = dict(chip_source=src, peaks_used={l: srcs[l] for l in LINES}, n_peaks=npk,
                    n_genes_sites_K562=int((cnts["K562"] > 0).sum()),
                    n_genes_sites_HepG2=int((cnts["HepG2"] > 0).sum()),
                    n_genes_with_sites=int(((cnts["K562"] + cnts["HepG2"]) > 0).sum()),
                    n_genes_feature_nonzero=int((feat[tf] != 0).sum()))
    print(tf, src, info[tf]["n_genes_with_sites"], file=sys.stderr)

# ---- (3)/(4) per-knockdown correlation and swap control ----
lam = pd.read_csv("data/annot/lambert2018_humantfs_v1.01.csv")
LAMBERT = set(lam.loc[lam["Is TF?"] == "Yes", "HGNC symbol"])
evaluable = [t for t in tfs if np.std(feat[t]) > 0 and info[t]["n_genes_with_sites"] > 0]
not_eval = [t for t in tfs if t not in evaluable]
rng = np.random.default_rng(SEED)
D = np.where(z["mask"], np.nan, (z["btA"] ** 2 - z["s2A"]) - (z["btB"] ** 2 - z["s2B"]))
rows = []
for tf in evaluable:
    i = targets.index(tf)
    own_id = kd.loc[kd.target == tf, "target_id"].iloc[0]
    ok = np.isfinite(D[i]) & (genes != own_id) & pd.notna(gchr)
    d = D[i, ok]
    rho_own = spearmanr(d, feat[tf][ok]).statistic
    others = [t for t in evaluable if t != tf]
    swap = rng.choice(others, size=min(NSWAP, len(others)), replace=False)
    rs = np.array([spearmanr(d, feat[s][ok]).statistic for s in swap])
    rank = 1 + int((rs > rho_own).sum())
    p_emp = (1 + int((rs >= rho_own).sum())) / (len(rs) + 1)
    # analytic Spearman p (single test, informative only)
    p_rho = spearmanr(d, feat[tf][ok]).pvalue
    rows.append(dict(tf=tf, **info[tf], n_genes_tested=int(ok.sum()), rho_own=float(rho_own),
                     p_rho_own=float(p_rho), rho_swap_median=float(np.median(rs)),
                     rho_swap_min=float(rs.min()), rho_swap_max=float(rs.max()),
                     rank_own=rank, n_compared=len(rs) + 1, p_rank_emp=float(p_emp),
                     swap_tfs=list(map(str, swap)), rho_swap=list(map(float, rs)),
                     expected_pattern=bool(rho_own > 0 and rank == 1),
                     lambert_tf=tf in LAMBERT))
allres = pd.DataFrame(rows).sort_values("rho_own", ascending=False)
res = allres[allres.lambert_tf]; supp = allres[~allres.lambert_tf]

# ---- figure ----
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(max(8, 0.38 * len(res)), 4.6))
res = allres  # plot all evaluable; non-Lambert marked with '+'
for x, r in enumerate(res.itertuples()):
    ax.scatter(np.full(len(r.rho_swap), x) + rng.uniform(-.15, .15, len(r.rho_swap)), r.rho_swap,
               s=9, c="0.65", zorder=1)
    ax.scatter(x, r.rho_own, s=55, zorder=3, marker="D",
               c="tab:red" if r.expected_pattern else ("tab:orange" if r.rho_own > 0 else "tab:blue"),
               edgecolor="k" if r.chip_source == "both" else "none")
ax.axhline(0, c="k", lw=.6)
ax.set_xticks(range(len(res)))
ax.set_xticklabels([f"{t}{'' if s == 'both' else '*'}{'' if lt else '+'}" for t, s, lt in zip(res.tf, res.chip_source, res.lambert_tf)],
                   rotation=70, fontsize=8)
ax.set_ylabel("Spearman ρ (d_g vs Δ accessible bound sites)")
ax.set_title("H5 case studies, K562 − HepG2 (mixed pair; illustrative). ◆ own TF, grey = 20 swap TFs\n"
             "red: own ρ>0 and rank 1; orange: own ρ>0; blue: ρ≤0; * ChIP in one line only; + not a Lambert TF", fontsize=9)
fig.tight_layout(); fig.savefig("figures/aim4_h5_tf_cases.png", dpi=150)
res = allres[allres.lambert_tf]

# ---- outputs ----
meta = dict(analysis="H5 mechanism case studies", pair=PAIR, date=str(datetime.date.today()),
            window_bp=WIN, accessible="region signal (cpm, normalized, CN-corrected) > line median over union regions",
            site="IDR peak overlapping >=1 accessible union region; position = summit",
            feature="count_K562 - count_HepG2", outcome="(btA^2-s2A)-(btB^2-s2B), A=K562, own gene excluded",
            n_swap=NSWAP, seed=SEED,
            tf_definition="ENCODE is_sequence_specific_TF flag AND Lambert 2018 'Is TF? = Yes' (main cases); "
                          "ENCODE-flagged non-Lambert factors reported as supplementary; swap pool = all evaluable "
                          "ENCODE-flagged factors with ChIP (self excluded)",
            n_encode_flagged_with_chip=len(tfs), not_evaluable={t: info[t] for t in not_eval},
            swap_pool=evaluable, pooled_test=None, n_tfs=len(res),
            n_tfs_both=int((res.chip_source == "both").sum()),
            n_expected_pattern=int(res.expected_pattern.sum()),
            n_rho_pos=int((res.rho_own > 0).sum()),
            chip_files={f"{l}:{t}": v for (l, t), v in files.items()}, cases=res.to_dict("records"), supplementary_non_lambert=supp.to_dict("records"))
json.dump(meta, open("results/aim4/h5_tf_cases.json", "w"), indent=1, default=str)

L = ["# H5 mechanism case studies — K562 − HepG2", "",
     f"Run {meta['date']} (`scripts/run_h5_cases.py`). Factor by factor, **no pooled test** (plan §7 H5). "
     "K562 − HepG2 is a mixed pair, so these cases illustrate a mechanism and do not establish it.", "",
     "- TFs: sequence-specific, kept knockdowns of K562_vs_HepG2, ENCODE IDR thresholded peaks (GRCh38, "
     "no audit ERROR; ENCODE4 > ENCODE3 > ENCODE2, then fewer NOT_COMPLIANT, more reps, newest).",
     f"- Site: IDR peak overlapping a union region accessible in that line (signal > line median); counted if "
     f"its summit is within {WIN // 1000} kb of the TSS. One-line ChIP: that peak set used for both lines, gated "
     "by each line's accessibility.",
     "- Feature: count_K562 − count_HepG2. Outcome: d_g = (β̃A² − s²A) − (β̃B² − s²B) for that knockdown; "
     "Spearman ρ over genes, own gene excluded.",
     f"- Swap: {NSWAP} random other TFs from the same list (seed {SEED}); rank 1 = own ρ highest of 21; "
     "p_emp = (1 + #swap ρ ≥ own)/21. Expected pattern = own ρ > 0 and rank 1.", "",
     "## Main cases (Lambert 2018 sequence-specific TFs)", "",
     "| TF | ChIP source | genes with sites | genes Δ≠0 | n genes | ρ own | p(ρ) | ρ swap median [min, max] | rank /21 | p_emp |",
     "|---|---|---|---|---|---|---|---|---|---|"]
def trow(r):
    return (f"| {r.tf}{' ✓' if r.expected_pattern else ''} | {r.chip_source} | {r.n_genes_with_sites} | "
             f"{r.n_genes_feature_nonzero} | {r.n_genes_tested} | {r.rho_own:+.3f} | {r.p_rho_own:.2g} | "
             f"{r.rho_swap_median:+.3f} [{r.rho_swap_min:+.3f}, {r.rho_swap_max:+.3f}] | {r.rank_own} | {r.p_rank_emp:.3f} |")
hdr = L[-2:]
L += [trow(r) for r in res.itertuples()]
L += ["", "## Supplementary: ENCODE-flagged factors that are not Lambert TFs (cofactors/RBPs; swap-pool members)", ""]
L += hdr + [trow(r) for r in supp.itertuples()]
L += ["", f"Not evaluable (no accessible sites near any tested gene): {', '.join(not_eval) or 'none'}."]
exp = res[res.expected_pattern]
L += ["", f"**Summary.** {len(res)} TFs ({meta['n_tfs_both']} with ChIP in both lines). Own ρ > 0 for "
      f"{meta['n_rho_pos']}/{len(res)}; expected pattern (own ρ > 0 and own sites best of 21) for "
      f"{len(exp)}: {', '.join(exp.tf) if len(exp) else 'none'}. No pooled claim is made.", "",
      "Figure: `figures/aim4_h5_tf_cases.png`. ChIP accessions and md5s: `h5_tf_cases.json`.", "",
      "**Reading.** Only PREB shows the expected pattern (own ρ = +0.028, best of 21, p_emp 0.048 — one of 20 "
      "factor-level looks, so not beyond chance). MYBL2, ATF4, CTCF (rank 4) and THAP1 (rank 2) have positive "
      "own ρ but are not uniquely best. All |ρ| ≤ 0.04: site differences explain < 0.2% of rank variance. "
      "Swap ρ medians are mostly slightly positive, i.e. any TF's K562-minus-HepG2 accessible-site count tracks "
      "a generic accessibility difference that weakly co-varies with d_g; own-vs-swap rank is the informative "
      "comparison. THAP11 and E2F6 (one-line ChIP) and the cofactors SFPQ, BCLAF1 sit at the bottom (negative ρ).",
      "", "**Caveats.** (1) Mixed pair: line identity, not chromatin alone, differs. (2) Two-line ChIP mixes "
      "binding differences with ChIP depth/quality differences (e.g. MYBL2 758 vs 24,426 peaks; BCLAF1 1,987 vs 50). "
      "(3) One-line ChIP features vary only through accessibility gating, so few genes have Δ ≠ 0. "
      "(4) The recon table lists only 'IDR thresholded peaks' (all GRCh38); conservative vs optimal is not "
      "distinguished there. (5) The ENCODE sequence-specific flag includes cofactors/RBPs; main cases are "
      "restricted to Lambert 2018 TFs (`data/annot/lambert2018_humantfs_v1.01.csv`)."]
open("results/aim4/h5_tf_cases.md", "w").write("\n".join(L) + "\n")
print("\n".join(L))
