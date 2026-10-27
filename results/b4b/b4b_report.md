# B4b — edge weights as a shared function of label-free edge/site features (DEV folds 0–4, pooled nested OOF)

GT_FREEZE sha256 `1e6f791c6b85f15dfe51d7eb045bf392a6b25be576b10974f7bba6c5d887baac`. DEV only (`cache/dev`, `results/ground_truth` DEV
files, `cache/b4/allowed.json`). Code `src/kinact/explore/edgeweights.py` (imports bayes2, unmodified); scripts
`scripts/exp_b4b_{common,fit,eval}.py`; tests `tests/test_edgeweights.py` (5 pass). Raw: `results/b4b/fit_fold*.{npz,json}`,
`b4b_eval.json`, `b4b_tables.md`. Ledger branch `B4b` (14 rows). Wall time: fits 6–14 min per fold (4 GPUs, parallel), eval 2 min.
**Verdict: [-] for the learned function (nested CV picks θ = 0 in 5/5 folds). The 1-parameter promiscuity variant gives
+0.007 M1-sel (n.s., fails the hard test) and +0.019 M1-all (significant), but costs calibration. Not adopted.**

## Reference (per PI update D7)
core-v1 = curated_union + r2np@100 (α 1), wsd, no winsor, ρ_f 3, **z_vif output**. τ and ρ_g are nested per fold: fold0 τ .125 / ρ_g 0;
fold2 τ .06 / ρ_g .25; others τ .125 / ρ_g .25. These are the picks over the 45-config B2.1 z_vif grid; they match the audit.
ρ_e is fitted on the 7 DEV nulls. Nested OOF M1-sel = **0.663**, which reproduces the audit's 0.663.
At θ = 0, `DiffCore` reproduces `bayes2` z_vif to within 5e-6 (max |Δ|) and gives the same NaN pattern. The test checks this.

## Method
- **Weights.** W_ks = base_ks·exp(5·tanh(θ·φ/5)). The tanh is a soft clamp for numerical safety; for |θ·φ| < 2 it is ≈ exp(θ·φ).
  θ = 0 is exactly core-v1. There are 11 parameters, with no per-kinase or per-edge parameters and no GT-derived features.
- **φ.** Each feature is centred and scaled over edges:
  - `is_r2np`, `psp_vivo`, `psp_vitro`;
  - `n_src` = log #curated source tags;
  - r2np / r3 / cantley within-kinase percentile (from the top-500 lists);
  - `log_prom` = log #kinases with the site in the library;
  - `res_mismatch` = TK group XOR Y residue;
  - `log_sd` and `meas_freq` of the standardised site value. These two use allowed conditions only: they exclude the
    outer fold, conditions sharing a control with it, and the inner validation fold.
- **Objective.** A smooth-AUROC pairwise loss, sigmoid(−(t_k − t_j)/(0.5·sd_c)) with t = sign·z_vif. It runs within each
  training M1-sel condition over scored evaluable negatives, with pair weights that reproduce the harness macro
  (pairs → (kinase, group) → kinase). Penalty λ‖θ‖². Full-batch L-BFGS (20 iterations) on GPU, with autograd through Cholesky.
  During training ρ_e is fixed at the fold's core value. After fitting it is refitted on the DEV nulls (label-free).
- **Nested CV.** λ ∈ {∞ (= core-v1), 0.1, 0.01, 0.001} is chosen per outer fold by pooled inner-CV M1-sel. Each inner fit
  trains on 3 folds and scores the 4th: 60 inner fits + 20 outer fits.
  Promiscuity variant: w = prom^−β, β ∈ {0, .25, .5, .75, 1}, chosen per outer fold by pooled M1-sel over the other 4 folds.
- **Disclosed look.** Before the nested run, one look at held-out fold1 using z_post decided the loss.
  - The pairwise **logistic** loss lowered M1-sel even on the training folds (fold1 0.684 → 0.677; fold0, a training fold, 0.696 → 0.675).
  - The sigmoid loss raised it (0.684 → 0.701).
  - Sigmoid was therefore pre-specified, and logistic is kept as a side variant.

## Results (nested OOF, 372 GT conds / 877 M1-sel pairs; Δ = paired group bootstrap vs core-v1 z_vif, B = 2000, C8 basic CI, one-sided p)
| variant | M1-sel | M1-all | harm (0/1-4/5-19/≥20) | H/D/O | sign −/+ | T-uns | ΔM1-sel [CI] p | ΔM1-all [CI] | ΔT-unseen [CI] | |z|>1.96 nulls / perm |
|---|---|---|---|---|---|---|---|---|---|---|
| **core-v1 z_vif (ref)** | **0.663** | 0.615 | .663 (.582/.688/.737/.645) | .629/.917/.675 | .637/.805 | 0.648 | ref | ref | ref | 0.064 / 0.054 |
| learned, λ nested (PRIMARY) | 0.663 | 0.615 | identical (θ = 0 chosen 5/5) | | | 0.648 | +0.000 | +0.000 | +0.000 | 0.064 / 0.054 |
| learned, λ=0.1 fixed | 0.659 | 0.614 | .661 (.588/.682/.734/.638) | .625/.906/.675 | .634/.806 | 0.645 | −0.003 [−.007,+.001] .93 | −0.000 | −0.002 [−.010,+.005] | 0.061 / 0.053 |
| learned, λ=0.01 fixed | 0.663 | 0.615 | .667 (.605/.677/.732/.654) | .618/.903/.679 | .638/.812 | 0.655 | +0.001 [−.012,+.012] .44 | +0.001 | +0.008 [−.003,+.027] | 0.068 / 0.056 |
| learned, λ=0.001 fixed | 0.662 | 0.615 | .668 (.614/.665/.730/.666) | .603/.903/.681 | .635/.817 | 0.657 | −0.000 [−.019,+.016] .51 | −0.000 | +0.010 [−.008,+.040] | 0.082 / 0.066 |
| learned, logistic, λ=0.01 | 0.651 | 0.617 | .656 (.605/.644/.726/.648) | .607/.908/.638 | .622/.808 | 0.650 | −0.012 [−.034,+.014] .83 | +0.002 | +0.002 | 0.036 / 0.045 |
| **prom^−β, β nested** | **0.670** | **0.634** | .672 (.606/.677/.738/.669) | .634/.939/.678 | .648/.798 | 0.653 | +0.007 [−.007,+.019] .16 | **+0.019 [+.013,+.031]** p<.001 | +0.005 [−.017,+.035] | 0.076 / 0.068 |
| prom^−β, β=0.25 fixed | 0.668 | 0.629 | .669 (.601/.684/.736/.656) | .632/.932/.677 | .644/.807 | 0.652 | +0.005 [−.006,+.014] .18 | +0.015 [+.011,+.022] | +0.004 | 0.072 / 0.061 |
| prom^−β, β=0.5 fixed | 0.666 | 0.636 | .670 | .626/.915/.678 | .644/.800 | 0.649 | +0.004 [−.016,+.020] .35 | +0.021 | +0.001 | 0.082 / 0.081 |
| prom^−β, β=1 fixed | 0.662 | 0.632 | .669 | .616/.875/.675 | .641/.790 | 0.639 | −0.000 | +0.017 | −0.008 | 0.106 / 0.107 |
| (core-v1 z_post ranking) | 0.666 | 0.615 | .666 | .644/.910/.674 | .643/.798 | 0.640 | +0.003 [−.005,+.011] .25 | +0.000 | −0.007 | z_post ≈10× conservative |

- Coverage is 0.994 for every row (fold0 uses ρ_g 0). M2 vs ref, matched universe: prom nested 0.670 vs 0.662; learned λ .01 0.663 vs 0.662.
- β picks: .5/.25/.25/.5/.25. The inner gain of β > 0 is small: +0.000 to +0.006 per fold.
- ΔT-unseen = M1-sel on pairs whose kinase was never a positive in the training folds (as in B2.1), bootstrapped on the same
  resamples. Calibration is fold-wise for the chosen config: 7 DEV nulls (3.5k cells) and site-permuted GT conditions (186k cells).

## Findings
1. **The learned function does not generalise.** Inner CV prefers θ = 0 in every fold. Learned weights lose 0.003–0.017 inner
   M1-sel, even though the training loss falls by 5–10%. With 11 parameters fitted on about 700 M1-sel pairs from about 70
   groups, the smooth-AUROC objective overfits group-specific patterns.
   - Fixed-λ OOF rows are within ±0.003 of core-v1.
   - Small T-unseen and bin-0 gains appear at λ ≤ 0.01: bin 0 .582 → .605/.614, T-unseen +0.008/+0.010, both n.s.
   - H drops at small λ (.629 → .603) and calibration on nulls worsens (8.2% at λ .001).
2. **θ is sane and consistent across folds**, mean ± sd over 5 outer fits at λ = 0.01:
   - log_prom −0.53±0.14 (the dominant term; λ .001: −0.83±0.15), so hub sites are down-weighted, as the B4 lesson predicts;
   - r2np_pct +0.21±0.04 (better-ranked predictions weigh more);
   - meas_freq +0.15±0.05 and log_sd +0.19±0.09 (frequently measured, responsive sites up);
   - n_src +0.08;
   - source indicators ≈ 0 (psp_vivo −0.08±0.06, is_r2np −0.05). r3/cantley percentiles are noisy (cantley −0.11±0.16)
     and res_mismatch is ≈ 0 (−0.06±0.07).
   Promiscuity is essentially the only stable signal, so the 1-parameter variant captures what the model learns.
3. **Promiscuity^−β (one parameter) is the useful part, but it does not pass the hard test.**
   - ΔM1-sel +0.007, p = 0.16. Harmony +0.009 [−0.005, +0.021]. All point non-inferiority checks pass:
     bins .606/.677/.738/.669 vs .582/.688/.737/.645 (bin 1-4 −0.011, within −0.02); H/D/O all up; T-unseen +0.005.
   - M1-all +0.019 [+0.013, +0.031] is robust for every β > 0. Down-weighting hubs helps most in multi-target conditions.
   - The D gain (+0.022) is in the D source, which is flagged as inflated by residue/memorisation effects.
4. **Calibration cost of re-weighting.**
   - Under the nested β, z_vif goes from 6.4% to 7.6% on nulls and from 5.4% to 6.8% under permutation. At β = 1 both are
     about 10.6%, and bin 0 on nulls is about 31% (n = 54 cells).
   - The ρ_e refitted on the 7 nulls shrinks toward 0 as β grows.
   - Likely cause: unequal weights cut the effective number of sites, and without winsorising the heavy-tailed log2FC then
     give non-Gaussian tails. This is not tested.
   - A fix would need winsor-3 (re-opening D7) or an empirical null.
5. **The logistic loss is misaligned with M1.** It is dominated by grossly misranked pairs (−0.012 OOF).

## Recommendation
- Do not adopt B4b into core-v1: the hard test fails and calibration degrades.
- Keep `prom^−β` (β = 0.25 fixed, 1 parameter) as a pre-registered candidate for B6 or a winsor-3 revisit (re-check calibration). It is the only
  edge-weight idea with a consistent sign, mainly on M1-all, and it costs nothing at runtime.
- A >50-parameter variant is not warranted: the 11-parameter model already overfits.

## Caveats
- One exploratory look at fold1 (z_post) chose the loss family before the nested run (disclosed above). The nested
  λ-selection used no other look.
- ρ_e rests on 7 null conditions from one study.
- Within inner CV the ρ_e refit uses the same nulls for all folds (label-free). The site statistics exclude the outer fold
  and the inner validation fold.
- Not covered: DEV conditions that share a control with TEST conditions (G2 rule, as in B4).
