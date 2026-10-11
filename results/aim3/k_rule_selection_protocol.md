# Block-size (k) rule selection for plan v1.1: protocol

Written and committed **before** the selection simulation was run (2026-10-11).
Simulation only; no real-data test is run for this selection.

## Candidate rules (`crt.choose_k(rule=...)`)

All rules use only the x residuals in genome order, autocorrelations within arm.
All are bounded to [5, 100] except v1.0 ([5, 50], as frozen).

| rule | definition |
|---|---|
| `v1.0` | smallest lag with ACF < 0.05 (plan v1.0) |
| `acf01` | smallest lag with ACF < 0.01 |
| `x2`, `x3`, `x4` | 2×, 3×, 4× the v1.0 k |
| `tau1`, `tau2` | ceil(c·τ), c = 1, 2; τ = 1 + 2·Σ ACF(h) for h = 1..L, L = first lag with ACF < 0.05 |
| `pw` | Politis–White (2004) automatic block length, circular-bootstrap version with the Patton–Politis–White (2009) correction |

Reference (not a candidate; needs the unknown truth): `oracle` = true modules.

## Grid

ρ_M ∈ {0.3, 0.5} × ρ_x ∈ {0.3, 0.6} × confounding ∈ {0, 0.3} × noise ∈ {homo,
hetero} × module design ∈ {geometric mean 5, geometric mean 10, geometric mean
25, mixed lognormal (median 8, σ = 1, clipped to [1, 150])} = 64 cells.
Modules are contiguous and truncated at arm ends. G = 2,000, B = 999.

- Null θ = 0: 1,000 replicates per cell.
- Power: partial R² 0.5% and 2%: 300 replicates per cell.
- Seeds: `SeedSequence([20261011, 11, cell index, replicate])`.

## Pre-set criterion (applied mechanically by `cdtest.k_rule_sim.select_rule`)

1. Eligible: worst-case FPR over the 64 null cells ≤ 0.065 at α = 0.05 **and**
   ≤ 0.02 at α = 0.01 (point estimates).
2. Among eligible rules, the one with the highest mean power at matched FPR at
   partial R² 0.5% (mean over the 64 cells; matched threshold = the 5% lower
   quantile of the rule's null p-values in the same cell).
3. Tie: any eligible rule within 0.02 of the best power counts as tied; the
   simplest tied rule wins. Simplicity order, fixed now:
   v1.0 < acf01 < x2 < x3 < x4 < tau1 < tau2 < pw.
4. If no rule is eligible: report so, report the rule with the lowest worst-case
   FPR at 5%, and make a recommendation (no automatic choice).

Monte Carlo SEs are reported for every rate.
