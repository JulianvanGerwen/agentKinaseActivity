# B6 — stacking / ensembling of score matrices (DEV folds 0–4, pooled nested OOF)

GT_FREEZE sha256 `1e6f791c6b85f15dfe51d7eb045bf392a6b25be576b10974f7bba6c5d887baac`. DEV only (`cache/dev`, DEV files of
`results/ground_truth`, existing DEV artefacts of B4/B4b/baselines). Code: `src/kinact/explore/stacking.py`,
`scripts/exp_b6_build.py` (components → `cache/b6/components.npz`), `scripts/exp_b6_eval.py` (`--stage comp|stack`).
Raw: `results/b6/complementarity.json`, `predeclared.json`, `stack_eval.json`, `oof_*.parquet`, logs. Ledger: 4 rows, branch `B6`.
**Verdict: no stacker passes the hard test. The pre-declared primary (0-parameter average of core-v1 z_vif and B4
co-regulation wz) gives ΔM1-sel +0.012 [−0.001, +0.030], p = 0.057, with significant gains on M1-all, harmony and
T-unseen and *better* calibration. The learned pairwise stacker gives −0.001 (selection noise eats the inner gains).**

## Components (OOF validity checked for each)
Every component is cond × kinase and label-free given the fold. Fold-specific components are stored per outer fold f
(config chosen for f by the producing branch's nested CV on folds ≠ f, applied to all conditions).
| comp | source | OOF validity |
|---|---|---|
| zvif | core-v1 z_vif, `results/b4b/fit_<f>.npz` `core` | τ/ρ_g nested per fold (D7); ρ_e from 7 DEV nulls (label-free); recompute max\|Δ\| 0 |
| mean | core-v1 posterior mean (same model, `BCtx.score(out="mean")`) | same as zvif |
| prom | core with prom^−0.25 edge weights (`prom_beta0.25`) | β fixed a priori (B4b recommendation), no fitting |
| wzexp | B4 weighted z + co-regulation expansion, B4 nested pick per fold | expansion of a fold-g row fitted on allowed conditions excl. fold g + shared controls; recomputed, max\|Δ\| 4e-6 vs B4 grid |
| phosx / kinlib | PhosX-like top15 / KinLib enrich top15 thr .5 (`cache/scores`) | fixed Cantley atlas, no fitting; recompute max\|Δ\| 0 |
| curz | z-score curated_union min1 | fixed library |
| self | B5 PSP self-site signal (RMS-standardised per condition) | label-free table |
Each component is also scored on the 7 DEV null conditions and on site-permuted DEV GT conditions (seed 0, the
*same* permuted input for every component) for normalisation, recalibration and calibration.

## 1. Complementarity (OOF, 877 M1-sel pairs; M2 = matched universe vs zvif)
| comp | M1-sel | M1-all | cov | harmony bins 0/1-4/5-19/≥20 | H/D/O | T-unseen | M2 comp vs zvif |
|---|---|---|---|---|---|---|---|
| zvif | 0.663 | 0.615 | 0.994 | .582/.688/.737/.645 | .629/.917/.675 | 0.648 | – |
| mean | 0.667 | 0.612 | 0.994 | .580/.689/.745/.656 | .641/.915/.677 | 0.647 | .667/.662 |
| prom | 0.668 | 0.629 | 0.994 | .601/.684/.736/.656 | .632/.932/.677 | 0.652 | .667/.662 |
| wzexp | 0.660 | 0.631 | 1.000 | **.645**/.662/.691/.658 | .631/.913/.648 | **0.684** | .657/.662 |
| phosx | 0.557 | 0.539 | 0.880 | .545/.557/.566/.553 | .527/.644/.531 | 0.579 | .594/.687 |
| kinlib | 0.550 | 0.541 | 0.880 | .534/.537/.572/.549 | .526/.587/.527 | 0.575 | .584/.687 |
| curz | 0.563 | 0.543 | 0.860 | .521/.573/.594/.572 | .557/.693/.552 | 0.575 | .650/.688 |
| self | 0.505 | 0.502 | 0.347 | ≈.50 | ≈.51 | 0.506 | .726/.785 (304 pairs) |

Median per-condition Spearman between components (cells scored by both):
zvif–mean 0.98, zvif–prom 0.96, zvif–wzexp **0.56**, zvif–curz 0.46, zvif–phosx 0.22, zvif–kinlib 0.20, zvif–self 0.16;
phosx–kinlib 0.79; wzexp–curz 0.53.
Correlation of per-pair AUROC with zvif: mean .99, prom .96, wzexp .73, curz .70, self .39, phosx .36, kinlib .30.
- Pairs better / worse than zvif by >0.1 AUROC: wzexp 20% / 18%; phosx 20% / 67%; curz 15% / 70%.
- **Only wzexp is both strong and complementary**, strongest exactly where core-v1 is weakest (harmony bin 0, unseen
  kinases). mean and prom are near-duplicates of zvif; motif and curated baselines are weak on the matched universe.
- Null SD on permutations (raw scale): zvif 0.99, prom 0.99, curz 1.16, phosx 1.03. wzexp is 2.0–2.4, so it is not
  calibrated alone: correlated expansion sites inflate it.

## 2. Stackers (pre-declared in `predeclared.json` after §1, before any stacker was run; 4 in total)
- **Normalisation.** Each component is divided by its null SD σ_j; NaN → 0 (no evidence on the z scale).
- **Missingness.** No additive missingness indicator. It is not sign-symmetric: it would add a fixed offset to kinases
  lacking a component, which acts as a label-sign prior (DESIGN C5/C10 spirit).
- **Combination.** s = Σ w_j z_j.
- **Recalibration.** Divide by the per-cell null SD √(w_oᵀ R_o w_o), where R is the null correlation of the
  normalised components on site-permuted data and o is the set of components observed in that cell.
  - zvif–wzexp R = 0.49–0.59 (0.50 on the DEV null conditions).
  - σ_j and R are **cross-fitted**: for outer fold f they come from permuted rows not in f, so calibration on fold f
    is out-of-sample.
- **A1 (PRIMARY):** mean of {zvif, wzexp}, 0 params. **A2:** mean of {zvif, wzexp, phosx, curz}.
  **A3:** rank-normal (per condition) mean of {zvif, wzexp}.
- **L1:** linear stacker over {zvif, prom, wzexp, phosx, curz, self} (6 params, shared across kinases).
  - Loss: smooth-AUROC pairwise loss within training M1-sel conditions, with pair weights that reproduce the M1 macro
    (pairs → (kinase, group) → kinase), trained on the recalibrated output.
  - Penalty: ridge λ‖w − w_core‖², so λ = ∞ is exactly core-v1. λ ∈ {∞, .1, .01, .001} is chosen by inner CV (4 folds).

| variant | M1-sel | M1-all | cov | harmony (0/1-4/5-19/≥20) | H/D/O | sign −/+ | T-uns | nsig | M2 vs ref |
|---|---|---|---|---|---|---|---|---|---|
| core-v1 z_vif (ref) | 0.663 | 0.615 | 0.994 | .663 (.582/.688/.737/.645) | .629/.917/.675 | .637/.805 | 0.648 | 49 | – |
| **A1 avg zvif+wzexp** | **0.675** | **0.635** | 1.000 | .678 (.640/.687/.726/.658) | .649/.932/.675 | .649/.819 | **0.683** | 57 | .671/.662 |
| A2 avg +phosx+curz | 0.675 | 0.630 | 1.000 | .678 (.647/.688/.721/.655) | .654/.917/.657 | .646/.833 | 0.682 | 47 | .672/.662 |
| A3 rank-avg zvif+wzexp | 0.672 | 0.634 | 1.000 | .675 (.634/.689/.719/.657) | .648/.929/.671 | .646/.817 | 0.681 | 55 | .669/.662 |
| L1 linear, λ nested | 0.662 | 0.616 | 1.000 | .663 (.585/.685/.732/.652) | .632/.923/.672 | .635/.810 | 0.645 | 49 | .664/.662 |

Paired group bootstrap vs ref (B = 2000, seed 0, C8 basic CI, one-sided p); T-unseen on the same resamples:
| variant | ΔM1-sel | ΔM1-all | Δharmony | ΔT-unseen |
|---|---|---|---|---|
| **A1** | **+0.012 [−0.001, +0.030], p = 0.057** | +0.020 [+0.010, +0.034], p<.001 | +0.015 [+0.001, +0.030], p .016 | +0.036 [+0.017, +0.074], p<.001 |
| A2 | +0.012 [−0.007, +0.037], p .14 | +0.015 [+0.001, +0.031] | +0.015 [−0.006, +0.035] | +0.035 [−0.001, +0.089] |
| A3 | +0.010 [−0.004, +0.028], p .10 | +0.019 [+0.010, +0.032] | +0.012 [−0.003, +0.026] | +0.033 [+0.012, +0.068] |
| L1 | −0.001 [−0.010, +0.007], p .61 | +0.001 | +0.000 | −0.003 |

Calibration, |z| > 1.96 (target 5%), overall with bins 0 / 1-4 / 5-19 / ≥20:
| variant | DEV nulls (3.5k cells) | site-permuted (cross-fitted) |
|---|---|---|
| ref | 6.4% (1.9/5.5/7.2/5.9) | 5.4% (6.5/5.8/5.2/4.9) |
| **A1** | **4.6%** (1.8/3.5/5.0/4.8) | **5.2%** (3.0/4.7/5.2/6.5) |
| A2 | 4.4% | 5.0% |
| A3 | 5.8% | 5.3% |
| L1 | 6.5% | 5.2% |

## Findings
1. **A1 fails the hard test narrowly** (p = 0.057 vs < 0.05). Point non-inferiority passes everywhere: M1-all +0.020;
   H +0.020, D +0.015, O ±0.000; bins 0 +0.058, 1-4 −0.001, 5-19 **−0.011** (within −0.02), ≥20 +0.013; M2 +0.009;
   T-unseen +0.036. Calibration improves 6.4% → 4.6% on nulls (bin ≥20 permutation rate 6.5%).
2. **The gain is fold-heterogeneous and TK-heavy.** Per-fold ΔM1-sel +0.041 / −0.016 / −0.014 / +0.022 / +0.050;
   of 106 kinases, 54 improve and 33 lose (by >0.01).
   - Largest gains are on RTKs: GAK +.33, INSR +.24, IGF1R +.23, EPHA2 +.18, EPHA6 +.16, ROS1, KIT.
   - Largest losses are on AGC kinases: PRKX −.25, PRKCB −.21, PRKD2, PRKCZ; also CDK12 and CSNK2A2.
   - The T-unseen gain is not the SRC-family effect: excluding SRC-family kinases and ABL1/2, T-unseen still goes
     0.569 → 0.610 (112 pairs) and M1-sel 0.646 → 0.658. Fifteen of 28 unseen kinases improve and 6 lose.
3. **The learned stacker does not generalise (L1).** Inner-CV gains ≤ 0.01, λ picks .01/.1/.01/.01/.1, weights
   unstable across folds (wzexp 0 → 0.44, phosx 0 → 0.23, curz −0.09 → 0.14, self 0.04 → 0.20). With ~70 groups,
   6 weights already overfit (same lesson as B4b).
4. **A2 (weak baselines added)** changes nothing on M1-sel and widens the CI (helps bin 0 and s = +1, costs on O).
5. **A3 (rank-normal)** is slightly worse than z; its "calibration" is per condition by construction, not a null one.
6. **mean, prom, self add nothing beyond zvif** (correlations 0.96–0.98; self coverage 0.35).

## Recommendation / caveats
- Do not adopt (no hard-test pass). A1 is the only ensemble worth carrying: 0 params and better calibration, so a
  pre-registered single confirmation (DEV-confirm look, or TEST secondary) is possible. Quote it as significant on
  M1-all, harmony and T-unseen, not on M1-sel.
- **Selection disclosure.** The A1 subset was chosen after seeing §1 on the same folds. wzexp itself carries B4's
  nested config selection (valid OOF). The L1 inner CV uses fold-f component configs that saw inner-fold labels.
  This is the same nuance as B2/B4, and it is minor.
- **Deployment cost.** wzexp needs the co-regulation fit (~8 s once on allowed conditions) plus a sparse wz per
  condition (<1 ms); G2 must exclude DEV conditions sharing a control with TEST from that fit (B4 note).
- **Runtime.** Component build 58 s; stacking eval 14 min (mostly L1 inner fits on CPU).
