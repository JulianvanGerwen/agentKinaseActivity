# B2.1 + B3 — refined and calibrated Bayesian core (DEV folds 0–4, pooled OOF)

GT_FREEZE sha256 `1e6f791c6b85f15dfe51d7eb045bf392a6b25be576b10974f7bba6c5d887baac`. DEV cache only (`cache/dev`, `results/ground_truth` DEV files).
Code: `src/kinact/explore/bayes2.py`; scripts `scripts/exp_b21_{common,grid,eval,extra}.py`; tests `tests/test_bayes2.py` (7, all
114 repo tests pass). Raw: `results/b21/grid_{out,src,kern}.tsv`, `b21_eval.json`, `b21_extra.json`, `b21_tables.md`.
Ledger: 10,001 rows, branch `B2.1`. 264 of them come from a first `out` run that I aborted after an edge-case fix in `fit_rho_null`.
Grid sizes: out 2,880, src 5,880, kern 945 evaluations. Wall time ≈ 1.5 h.

## What changed vs B2 (model unchanged: y = Wa + ε, a ~ N(0, τ²K), closed form, median-centred, borrow on, no noise model)
- **Outputs** (same posterior; A = Λ⁻¹W'M):
  - `mean` = μ;
  - `z_post` = μ/√diag Λ⁻¹ (the B2 output);
  - `z_freq` = μ/√q with q = diag(AA'). This is the frequentist SD under exchangeable sites;
  - `hybrid_λ` = μ/q^{λ/2};
  - `z_vif` = z_freq/√(1+ρ_e(m_eff−1)), where m_eff = (ΣA_k)²/ΣA_k² is computed analytically. ρ_e is one global
    compound-symmetry correlation, fitted per config on the 7 DEV role=null conditions only (fitted values 0.005–0.02);
  - `empnull` = per-condition empirical null. z_freq is scaled by s0_c√(1+ρ_c(m_eff−1)), with the fit made on that
    condition's own scored kinases (robust MAD in m_eff quartiles; label-free).
- **Scale/robustness**: `wsd` = MAD-standardise, winsorise |y| ≤ c, then divide by the RMS of the clipped y.
  - Per-condition scale never changes M1 (checked: MAD and wsd give identical M1). It only affects calibration.
  - Winsorising does change rankings.
- **Library ensemble**: + r3 / cantley global top-k edges (k ∈ {100, 200}), weight α_src (≤2 extra HPs). An edge that
  several sources predict gets the sum of their weights.
- **Kernel**: + ρ_s·S, where S is the cosine² of top-500 membership of r2np or cantley (label-free; PSD by Schur). ρ_f ≤ 3 everywhere.
- **Pre-specified criterion** (written in `exp_b21_common.py` before any B2.1 result): the output is a hyperparameter
  chosen by nested inner M1-sel. The reported z is a calibrated statistic regardless of which output wins.

## Results (nested OOF; Δ = paired group bootstrap vs the B2 reference nested OOF 0.655, B=2000, C8 basic CI, one-sided p)
| variant (nested over) | M1-sel | M1-all | M2 vs PSP-z (0.633) | harmony 0 / 1-4 / 5-19 / ≥20 | T-unseen | ΔM1-sel vs B2 [CI] p | |z|>1.96 DEV nulls: all (0/1-4/5-19/≥20) | site-perm: all (0/1-4/5-19/≥20) |
|---|---|---|---|---|---|---|---|---|
| B2 reference (B2 grid, ρ_f up to 10) | 0.655 | 0.604 | 0.667 | .582/.664/.741/.646 | 0.589 | ref | 0.025 (–/0/.023/.045)ᵃ | 0.046 (.012/.028/.047/.076)ᵃ |
| **B2 model, ρ_f ≤ 3, z_post (no winsor)** | **0.670** | 0.612 | 0.695 | .583/.691/.749/.656 | 0.643 | **+0.015 [+0.006, +0.035] 0.001** | 0.006 (0/0/.007/.009)ᵇ | 0.002 (0/.001/.002/.006)ᵇ |
| out: mean | 0.668 | 0.617 | 0.689 | .598/.691/.740/.643 | 0.651 | +0.013 [+0.003, +0.036] 0.030 | (empnull z) 0.071 | 0.066 |
| out: z_freq | 0.659 | 0.614 | 0.682 | .583/.685/.732/.637 | 0.644 | +0.004 [−0.010, +0.027] 0.49 | 0.072 (.019/.055/.081/.070) | 0.063 (.067/.062/.062/.069) |
| out: hybrid (λ free) | 0.665 | 0.609 | 0.690 | .585/.693/.738/.643 | 0.644 | +0.011 [−0.000, +0.033] 0.10 | – | – |
| out: output free (pre-specified) | 0.663 | 0.610 | 0.689 | .583/.690/.739/.638 | 0.641 | +0.008 [−0.002, +0.031] 0.22 | 0.059 | 0.050 |
| winsor 3: mean | 0.669 | 0.616 | – | .567/.691/.748/.651 | 0.648 | +0.014 [−0.003, +0.043] 0.091 | (empnull z) 0.063 | 0.052 |
| winsor 3: z_freq | 0.669 | 0.628 | – | .577/.687/.745/.651 | 0.655 | +0.014 [−0.005, +0.045] 0.12 | 0.078 (.037/.053/.089/.078) | 0.061 (.064/.061/.059/.065) |
| winsor 3: empnull | 0.670 | 0.627 | – | .582/.690/.748/.644 | 0.663 | +0.015 [−0.003, +0.045] 0.091 | 0.068 (.130/.068/.079/.045) | 0.053 (.045/.057/.053/.047) |
| **winsor 3: z_vif (ρ_e on DEV nulls)** | **0.672** | **0.629** | – | .584/.690/.751/.652 | **0.664** | +0.018 [−0.001, +0.046] 0.051 | **0.064 (.037/.048/.076/.054)** | **0.050 (.060/.056/.048/.043)** |
| winsor free, output free | 0.668 | 0.623 | 0.678 | .576/.692/.740/.642 | 0.644 | +0.013 [−0.006, +0.042] 0.14 | 0.056 | 0.039 |
| src (r3/cantley ensemble): mean | 0.664 | 0.613 | 0.680 | .577/.680/.744/.645 | 0.643 | +0.009 [−0.009, +0.031] 0.18 | 0.067 | 0.053 |
| src: output free | 0.648 | 0.609 | 0.667 | .574/.660/.735/.611 | 0.635 | −0.007 [−0.028, +0.017] 0.74 | 0.055 | 0.038 |
| src: r2np-only subset, output free | 0.651 | 0.600 | 0.677 | .568/.661/.743/.621 | 0.639 | −0.004 | 0.063 | 0.042 |
| kern (+S): mean | 0.673 | 0.617 | 0.688 | .574/.694/.749/.653 | 0.656 | +0.018 [+0.001, +0.046] 0.038 | 0.060 | 0.052 |
| kern: z_freq | 0.670 | 0.624 | 0.687 | .576/.690/.741/.662 | 0.668 | +0.016 [−0.003, +0.047] 0.11 | 0.081 | 0.061 |
| ALL grids, output free | 0.642 | 0.609 | 0.662 | .567/.657/.727/.604 | 0.624 | −0.013 [−0.035, +0.014] 0.83 | 0.056 | 0.039 |
| prior-only control (B2 report) | 0.687 | 0.543 | – | .625/.705/.735/.840 | 0.324 | – | – | – |

- ᵃ From the B2 report: MAD scale, pooled-best config.
- ᵇ wsd scale. z_post is very conservative because its prior-dominated SD does not scale with the signal.
- "(empnull z)": calibration of the companion empirical-null z of the same config.
- Calibration is fold-wise: the config chosen for fold f scores the 7 DEV nulls (≈1.4k kinase scores; bin 0 has only
  54 cells) and the site-permuted GT conditions of fold f.
- Pairwise: winsor-3 empnull − mean = +0.001 [−0.007, +0.010]. Each winsor-3 variant vs "B2 model, ρ_f ≤ 3": −0.001 to +0.002, all n.s.
- Runtime: CPU single-sample path (`SingleSample2`, incl. standardisation, posterior, q, m_eff and empirical null;
  hybrid config) 18.4 ms/cond at 1 thread and 7.9 ms at 16 threads. B2 was 8.9 ms at 1 thread; the new path forms Λ⁻¹
  explicitly for q. GPU batch: posterior + q ≈ 1 ms/cond.

## Findings
1. **Most of the gain over the B2 headline is the ρ_f cap, not a new idea.** Re-running the B2 model on a ρ_f ≤ 3 grid
   gives nested 0.670 (Δ +0.015, p = 0.001). The B2 nested 0.655 was depressed by folds that picked ρ_f = 10 and lost
   decryptM (D 0.77 → 0.91 here). This passes the hard adoption test, and the point non-inferiority checks hold (M1-all
   +0.008, every harmony bin ≥ B2, T-unseen +0.054, H +0.002). It is ≥ 0.02 on no endpoint except T-unseen, so the audit
   rule is borderline; recommend a light audit.
2. **Output statistic (item 1).**
   - Without winsor: mean ≈ z_post > hybrid > z_freq. The pre-specified "output free" nested pick is 0.663, slightly
     worse than either fixed output (selection noise across 8 outputs).
   - With winsor 3 all outputs converge to 0.669–0.672. The posterior-mean advantage seen in B2 (0.672) is not robust:
     here mean − empnull = −0.001.
3. **Calibration (item 2).**
   - B2's z_post grows with substrate count (1.2% → 7.6%) and its absolute level is arbitrary (it depends on the scale
     convention).
   - z_freq with the wsd scale is about 6% flat across bins under site permutation (MAD scaling alone gives 24%:
     log2FC are heavy-tailed, so MAD underestimates SD).
   - On real nulls, z_freq still inflates with m (bins 5–19 and ≥20 at 8–9%). This is correlated site noise.
   - z_vif (one ρ_e ≈ 0.01, fitted on DEV nulls only) gives 5.0% under permutation (4.3–6.0% per bin) and 6.4% on DEV
     nulls (3.7–7.6% per bin).
   - The per-condition empirical null is similar: 5.3% permutation, 6.8% on nulls (bin 0 13%, but n = 54).
   - Both are at no M1 cost; z_vif is +0.002 vs the z_post config, n.s. Calibration is now close to the 5% target and
     roughly independent of substrate count.
4. **Robustness (item 5).** MAD vs SD has no effect on M1 (scale invariance within a condition) but decides calibration
   (above). Winsorising at 3 adds about +0.004 at the grid optimum for every output, and M1-all +0.015. Within nested
   CV it is n.s. but consistent. Winsor 2 is worse than 3.
5. **Library ensemble (item 3): negative.**
   - r3 / cantley edges never beat r2np-only at the grid optimum (0.680 vs 0.679 pooled-best). Adding them at α ≥ 0.5
     lowers M1 (cantley@200·α1: 0.661).
   - r2np@200 < r2np@100 here, and α_r2np = 2 < 1.
   - The large grid (5,880) costs about 0.02–0.03 in nested selection noise (0.648).
   - Do not adopt.
6. **Kernel (item 4): flat.**
   - S at ρ_s 0.5–1 is pooled-best ±0.001 (cantley S 0.679, r2np S 0.679, none 0.679); ρ_s = 3 hurts.
   - Nested "kern: mean" 0.673 is within noise of "winsor 3: mean" 0.669, and it costs an extra HP.
   - The family/group tree stays; ρ_f = 3 is chosen in 24/25 fold picks.
7. **Selection noise is the dominant hazard.** Nested estimates fall as grids grow: ALL grids 0.642 vs pooled-best
   0.680. Keep the final search space small and fixed.

## Recommended config (≤5 HPs; DEV nested 0.672, harmony .584/.690/.751/.652, T-unseen 0.664)
- Library: curated_union (w = 1) + r2np global top-100, binary, α = 1 (fixed). No r3/cantley, no S-kernel.
- Standardisation: median-centre, MAD-standardise, winsorise |y| ≤ 3, rescale to unit RMS (`scale="wsd", winsor=3`).
- Prior and borrowing: τ = 0.125, ρ_f = 3, ρ_g = 1 (folds split between 0.25 and 1: 3/5 picked 1); borrow on.
- Output: z_vif with ρ_e = 0.01 (fitted on DEV role=null conditions). It is a calibrated z (5.0% permutation / 6.4%
  DEV nulls, flat across bins). Report the per-condition empirical-null z as a sensitivity check.
- HPs: τ, ρ_f, ρ_g, winsor c, ρ_e (the last is fitted on nulls, not on GT).
- vs B2 reference: +0.018 [−0.001, +0.046], p = 0.051.
- vs prior-only: still below on M1-sel (0.672 vs 0.687) and H, far above on T-unseen (0.664 vs 0.324). Stacking (B6)
  remains the route to "beats".

## Caveats
- Sequential design: winsor 3 was fixed for the src/kern grids after seeing the out grid (pooled), which carries mild
  optimism. The out-grid nested numbers do not have this problem.
- ρ_e rests on 7 null conditions from one study (decryptM 0-vs-0 nM). Bin-level null rates have wide uncertainty
  (≈1.4k cells; bin 0 has 54).
- D-source numbers are inflated by residue class and PSP memorisation (r2np audit); do not read them as method quality.
- Runtime (18 ms/cond at 1 thread) misses the ≤5 ms target. q only needs diag(Λ⁻¹GΛ⁻¹), so a triangular-solve
  implementation could recover about 2×. Not done.
