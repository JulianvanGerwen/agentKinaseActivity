# B7b — late fusion of perturbation-response signatures (DEV folds 0–4, nested OOF; DESIGN §12 C10)

## PRE-DECLARED DESIGN (written 2026-10-27 before any B7b number was computed; not changed afterwards)
Follow-up to B7 (results/b7/b7_report.md). Signatures leave the Bayesian kernel and are added as a separate score.
- **Signatures**: B7 PRIMARY variant `clean+contrast`: positives only from conditions with n_targets ≤ 10, generic
  response subtracted, same fit code (`kinact.explore.signatures.fit_signatures`) and same training-set rules / leakage guards.
  Shrinkage κ and n_top are NOT re-tuned. For outer fold f they are fixed to B7's nested picks for fold f
  (fold0 100/κ0, fold1 200/κ0, fold2 50/κ4, fold3 200/κ4, fold4 200/κ1). The same (n_top, κ) is used for fold f's inner
  scoring. Those picks used only folds ≠ f.
- **Signature score** (for kinases with a signature in the relevant training fit only):
  z_sig(k,c) = Σ_s w_ks·y_cs / sqrt(Σ_s w_ks²).
  - y_cs = core-v1 input (wsd-standardised log2FC), clipped at |y| ≤ 4 as in the fit.
  - The sums run over k's top-n_top signature sites measured in c. At least 5 measured sites are required, else NaN.
  - Then z_sig is divided by one scalar ρ_sig per signature fit: ρ_sig = RMS of z_sig over the DEV role=null conditions
    (the analogue of core-v1's ρ_e). A scalar divisor keeps it odd in y and adds no intercept or kinase offset (C10 ii).
- **Fusion**: z = (z_core + α·z_sig)/sqrt(1+α²) where both are defined. z = z_core when the kinase has no signature
  or z_sig is NaN. z = z_sig if z_core is NaN (rare).
  - z_core is the core-v1 z_vif with core-v1 nested core HPs per fold, i.e. exactly the B7 α_sig = 0 reference.
  - The fusion assumes independence. Calibration is re-checked; there is no further recalibration.
- **HP**: α ∈ {0.25, 0.5, 1} only, nested.
  - inner_f(α) = M1-sel over folds g ≠ f. Fold g is scored with signatures fitted on folds ∉ {f, g} and z_core from
    the core-v1 inner scores (B7 grid store).
  - α*_f = argmax inner_f. Fold f is then scored with signatures fitted on the 4 folds ≠ f.
- **ONE primary comparison**: B7b nested vs core-v1 nested. Paired group bootstrap (B = 2000, seed 0, C8 inference).
  **Hard test: one-sided p < 0.05 for ΔM1-sel > 0.**
- **C10 adoption checks**, all required:
  1. T-unseen point non-inferiority, ΔM1-sel(T-unseen) ≥ −0.01.
  2. Each harmony bin point ΔM1-sel ≥ −0.02.
  3. Shuffled-signature control: site labels permuted within each signature (B7 `Signatures.shuffled`), same fold-wise
     α* / n_top / κ, 3 seeds. Passes if the mean shuffled ΔM1-sel vs core-v1 lies in [−0.01, +0.01] (collapses to core).
     Real − shuffled is reported.
  4. Drug-held-out within fold: grouped folds are drug-disjoint (shared drugs / groups across folds = 0), and the fit
     raises on a shared drug group. Confirmed by re-running the B7 disjointness check.
  5. Calibration: |z| > 1.96 rate on the DEV nulls (fold-f nulls scored with the fold-f design) and on site-permuted GT
     conditions. Passes if the overall rate is ≤ core-v1 + 1 pp on both (core-v1: 6.4% / 5.4%). Per-substrate-bin rates
     are reported.
- **Secondary (descriptive)**: T-seen; kinases with / without a training signature (without: identical z by
  construction); M1-all; harmony; sources H/D/O; M2 vs core-v1. Not part of the decision.
- Ledger: branch `B7b`, 3 per-α OOF rows + 1 nested row. Code: `src/kinact/explore/signatures_late.py`,
  `scripts/exp_b7b_eval.py`. Time box about 1.5 h.

## Results (appended after the run; design above unchanged)
- Run: `scripts/exp_b7b_eval.py`, about 25 s. Raw output: `results/b7b/{b7b_eval.json, grid.tsv, eval.log}`.
- GT_FREEZE `1e6f791c…`. Ledger branch `B7b`: 3 per-α rows + 1 nested row.
- ρ_sig ranges 0.68–1.00 across fits, so z_sig is slightly over-dispersed on the 7 DEV nulls before scaling.

**α grid** (OOF, descriptive): 0.25 → 0.677 · 0.5 → 0.680 · 1 → 0.681.
- Nested picks: α = 1 in folds 0–3 and α = 0.5 in fold 4.
- α = 1 is the top of the pre-declared grid. The grid was not extended.

**Primary comparison**: B7b nested vs core-v1 nested, paired group bootstrap B = 2000.

| endpoint | core-v1 | B7b | Δ [95% basic CI] | one-sided p |
|---|---|---|---|---|
| **M1-sel (hard test)** | 0.663 | 0.680 | **+0.018 [−0.000, +0.023]** | **0.0085 → PASS** |
| M1-all | 0.615 | 0.623 | +0.009 [+0.001, +0.011] | 0.010 |
| harmony | 0.663 | 0.683 | +0.020 [−0.000, +0.026] | 0.006 |
| M2 (matched, 872 pairs) | 0.662 | 0.679 | +0.017 | 0.010 |

The basic CI's lower bound touches 0 while the shifted-null p is 0.0085. The two C8 forms disagree slightly here; the pre-declared rule uses p.

**C10 checks**
1. **T-unseen** (272 pairs): 0.648 → 0.648, Δ +0.0006 [−0.0001, +0.0013]. **PASS** (margin −0.01). The TYK2-style coupling from B7 is gone, because signatures stay outside the kernel.
2. **Harmony bins**, ΔM1-sel by bin. **PASS**, all ≥ −0.02:

   | bin | 0 | 1–4 | 5–19 | ≥20 |
   |---|---|---|---|---|
   | ΔM1-sel | +0.028 | +0.024 | +0.008 | +0.021 |

3. **Shuffled signatures** (3 seeds): Δ −0.002 / −0.006 / −0.008, mean −0.005. **PASS**: it collapses to within ±0.01 of core-v1.
   - Real − shuffled = +0.023.
   - Seed 2's CI excludes 0 on the negative side, so random edges cost a little.
4. **Drug-held-out within fold**: 0 shared drugs and 0 shared groups between each eval fold and its training folds. The fit's SignatureLeak guards were active for all 15 fits. **PASS**.
5. **Calibration**, |z| > 1.96. **PASS**:

   | | DEV nulls, 3,519 cells: all (0 / 1–4 / 5–19 / ≥20) | site-permuted GT: all (0 / 1–4 / 5–19 / ≥20) |
   |---|---|---|
   | core-v1 | 6.4% | 5.4% |
   | B7b | 6.2% (1.9 / 5.5 / 6.8 / 5.9) | 5.5% (6.6 / 5.8 / 5.3 / 5.3) |

   z_sig alone on permuted GT is over-dispersed (8.8%, SD 1.15). Fusion with weight 1/√2 dilutes this.

**Secondary (descriptive)**
- T-seen: 0.668 → 0.692, Δ +0.024 [+0.002, +0.032].
- Kinases with a training signature (366 pairs): 0.734 → 0.836, Δ +0.102 [+0.068, +0.154], p 0.0005.
- Kinases without one (511 pairs): Δ −0.0002. Their z is identical, which the script asserts; only the ranks shift.
- Sources: H +0.025 (p .002), D −0.001, O +0.023.
- Sign: −1 +0.018, +1 +0.012.
- M2 vs PSP-z: 0.708 vs 0.633 (core-v1 0.687).

## Verdict
All pre-declared criteria pass: hard test p = 0.0085, T-unseen, harmony bins, shuffled control, drug-held-out and calibration.
- ΔM1-sel +0.018 is below the 0.02 audit trigger.
- The gain is concentrated in signature kinases (+0.10). An independent audit is still recommended before adoption.
- Points for the audit:
  - α sits at the grid edge;
  - (n_top, κ) were inherited from B7's nested picks, which used inner folds only;
  - the CI lower bound is ≈ 0;
  - B7b is one look after B7 on the same DEV folds (garden of forking paths: B7 → B7b). DEV-confirm is the natural next check.
