# Independent audit — B7b late-fusion signatures (DEV, nested OOF)

- Auditor: independent agent. I did not produce B7b.
- Date: 2026-10-27. GT_FREEZE `1e6f791c6b85f15d`.
- Inputs read:
  - RESEARCH_TREE §5.7 and DESIGN §7 and §12 (C8, C10);
  - `results/b7/b7_report.md` and `results/b7b/b7b_report.md`;
  - `src/kinact/explore/{signatures,signatures_late,bayes2}.py`;
  - `scripts/exp_b7b_eval.py` (and its helpers `exp_b7_common.py`, `exp_b7_eval.py`, `exp_b21_common.py`);
  - `tests/test_signatures.py` and `docs/harness.md`.
- I did not touch `cache/sealed/` or the history of `results/ground_truth`.
- My runs had `ledger.append` disabled, so `results/dev_ledger.tsv` is unchanged.
- Audit scripts, in the session scratchpad:
  - `run_rep.py`: the original `main()`, with output redirected;
  - `audit.py`, `audit2.py`, `audit3.py`: the extra checks;
  - raw output in `audit.json` and `audit.log`.

## Verdict: ADOPT, conditionally
1. The pre-declared decision rule is met, and the result reproduces exactly.
2. No leakage path was found. Direct label and data randomisation of fold f leaves the fold-f scores bit-identical.
3. α = 1 sits on a plateau (0.75–2), not on a rising edge.
4. Two caveats should be recorded with the adoption:
   - **Concentration.** The top 5 kinases carry 71% of Δ. Without them, Δ = +0.005.
   - **Over-dispersion.** Signature kinases get over-dispersed and false-positive scores:
     - FP rate as negatives: 8.1% → 9.8%;
     - |z| > 1.96 on site-permuted data, signature kinases only: 5.6% → 7.0% (bin 0: 9.3%; single kinases up to 11%).
5. Required before G2:
   - one DEV-confirm look, as the report suggests;
   - a note in DESIGN/RESEARCH_TREE that signature-kinase z values are not calibrated at 5% (see §4).
   - Any ρ_sig fix must be declared before it is evaluated. Because α sits on a plateau, a rescale of up to about 1.5× is not expected to change the ranking.

## 1. Reproduction and adherence to the declared design
- **Exact reproduction.** I re-ran `exp_b7b_eval.main()` with the ledger disabled. The whole `results` block of `b7b_eval.json` is byte-identical to the committed file.
  - M1-sel: 0.6627 → 0.6803, Δ +0.0176, basic CI [−0.0002, +0.0231], p = 0.0085.
  - Nested picks: α = 1 in folds 0–3 and α = 0.5 in fold 4.
  - T-unseen Δ: +0.0006.
  - Harmony bins: +0.028 / +0.024 / +0.008 / +0.021.
  - Shuffled controls: −0.002 / −0.006 / −0.008.
  - Calibration: 6.2% / 5.5%; z_sig alone 8.8%.
- **Cached outer fits.** I refitted the five outer signature fits from scratch (`fit_signatures` on the 4 training folds). They are identical to the cached `cache/b7/sig_clean+contrast_fold*.pkl` used by the script (sites and weights equal).
- **Declaration vs code.** Every declared element is implemented as written:
  - signature fit: clean+contrast;
  - (n_top, κ) taken per fold from B7's nested picks;
  - z_sig formula, with clip 4 and ≥ 5 measured sites;
  - ρ_sig = RMS over the DEV nulls, one value per fit;
  - fusion formula and NaN fallbacks;
  - α grid {0.25, 0.5, 1}, with inner_f over folds g ≠ f, using fits that exclude {f, g} and core-v1 inner z;
  - one primary comparison;
  - the 5 C10 checks with their declared thresholds.
  - I found no undeclared knob.
- **Ledger.** It shows exactly one B7b run: 3 per-α rows and 1 nested row, `n_configs_tried` = 4. There are no other `late-sig` rows, so there was no hidden re-running within B7b.
- **Caveat: pre-registration timing cannot be verified.** The declaration was committed in the same commit as the results (f5648f2). The report file's mtime is after `b7b_eval.json`. Content and code agree, but the timing rests on the producer's word.
- **Garden of forking paths.** B7b is a second look at the same DEV folds. It was motivated by B7's diagnosis that kernel coupling hurt T-unseen. The producer acknowledges this, and DEV-confirm is the right guard.

## 2. Leakage
| check | result |
|---|---|
| Fold f's fit uses only GT conditions of the other folds (`training_conditions` filters on split, qc and group; the guards raise) | code reviewed; guards active in all 15 fits |
| Drug and group sharing between the eval fold and its training folds | 0 shared drugs and 0 shared groups in every fold (re-run) |
| Generic-response subtraction | computed inside `fit_signatures` from the training columns only (`Y[:, train]`) |
| Standardisation | per-condition wsd (column-wise, label-free); site σ from training groups only; ρ_sig from the 7 DEV null conditions (label-free scalar) |
| Nested α | inner_f restricts pairs to folds ≠ f; its signatures exclude f |
| **Direct test A**: fold-f GT labels replaced by random kinases and signs (1.6k–6.6k rows per fold), n_targets randomised | refitted signatures identical; fold-f z_sig bit-identical, all 5 folds |
| **Direct test B**: fold-f data columns replaced by Gaussian noise before fitting | signatures identical, all 5 folds |
| **Direct test C**: Evaluator rebuilt with randomised fold-f direct labels; inner_f(α) recomputed | inner_f identical for every α and fold, so α*_f is unchanged |
| Popularity via a generic response (late-fusion generic control: every signature kinase gets the generic-response signature) | Δ +0.002 [−0.009, +0.008], p = 0.29: none |
| Push toward the "target" direction on signature-kinase negatives (mean s_c·z_sig) | −0.03 to −0.06 overall, with no positive shift in shared-control studies (largest per-study mean +0.47, n = 176; most are ≤ 0) |

- **Shared controls.** 45% of M1-sel pairs come from conditions that share a (PMID, Control) with a training condition, for example Hijazi.
  - Δ is +0.029 in those pairs, against +0.013 in the others (p 0.017).
  - The gain is larger where controls are shared. There is no mean push on negatives, so this is not a popularity leak.
  - It may still reflect study-specific readouts. The not-shared subset alone still gains +0.013.

## 3. Robustness
**Fixed-α OOF curve** (outer fits; Δ vs core-v1, C8 inference). The values 1.5–1e6 are descriptive only, not for selection.

| α | 0.1 | 0.25 | 0.5 | 0.75 | **1** | 1.5 | 2 | 3 | 5 | z_sig only |
|---|---|---|---|---|---|---|---|---|---|---|
| ΔM1-sel | .009 | .014 | .017 | .018 | **.019** | .018 | .018 | .018 | .017 | .013 |
| p | .020 | .014 | .007 | .004 | .005 | .008 | .010 | .016 | .026 | .064 |

- α = 1 is the top of a broad plateau (0.75–2), not a rising edge.
- A nested grid extended to {.., 1.5, 2} (descriptive) would have picked α = 1.5–2 and given Δ +0.016 (p 0.016). Extending the grid would not have helped.

**Basic CI vs shifted-null p.** These do not actually contradict each other.
- The replicates are biased upward: mean Δ* = 0.0221 against Δ̂ = 0.0176 (bias +0.0046). The skew is +0.52, and the percentile CI is [0.012, 0.035].
- The basic interval reflects the replicates around Δ̂. Its lower bound is 2Δ̂ − q97.5 = −0.0002, so the bias is subtracted.
- The shifted-null p uses only the spread around mean(Δ*): SE 0.0061, so Δ̂/SE ≈ 2.9 and p ≈ 0.0085.
- The 95% two-sided basic CI corresponds to a one-sided 2.5% test. The matching one-sided 5% test is the basic 90% lower bound: 2Δ̂ − q95 = **+0.0022 > 0**.
- The p-value consistent with the basic interval, P(Δ* ≥ 2Δ̂), is **0.028**.
- So both C8 forms reject at the pre-declared one-sided 0.05. The basic form is weaker (about 0.03, not 0.0085).

**Breakdowns.**
- Per fold (within-fold M1-sel, core → B7b):

  | fold | 0 | 1 | 2 | 3 | 4 |
  |---|---|---|---|---|---|
  | core-v1 → B7b | 0.632 → 0.657 | 0.705 → 0.723 | 0.691 → 0.717 | 0.756 → 0.796 | 0.665 → 0.753 |
  | Δ | +0.024 | +0.017 | +0.026 | +0.040 | +0.088 |

  The sign is positive in 5 of 5 folds.
- Sources:
  - H: +0.025 [+0.006, +0.037], p .002;
  - D: −0.001 [−0.015, +0.010], p .57;
  - O: +0.023 [−0.003, +0.035], p .017.
- Sign:
  - −1: +0.018, p .016;
  - +1: +0.012, p .15.
- Other endpoints:
  - M1-all: +0.009, p .009;
  - harmony: +0.020, p .005;
  - M2: +0.017, p .010;
  - T-unseen: +0.0006 (p .051), so non-inferior.

## 4. Calibration (|z| > 1.96)
| set | core-v1 | B7b | z_sig alone |
|---|---|---|---|
| Site-permuted GT, all kinases | 5.37% | 5.50% | – |
| Site-permuted GT, **signature kinases only** (16.2k cells), bins 0 / 1–4 / 5–19 / ≥20 | 5.6% (5.6 / 5.9 / 6.1 / 5.0) | **7.0% (9.3 / 6.1 / 7.4 / 6.9)** | 8.8% (11.5 / 8.3 / 9.5 / 8.2); SD 1.15–1.22 |
| Site-permuted GT, per signature kinase (≥ 50 cells): median / 90th percentile / max | 5.6 / 7.3 / 8.1% | 6.8 / 10.1 / 11.3% (MAP2K1, MAPK1/3, ROS1, MTOR, PRKACA) | – |
| DEV nulls, signature kinases (only 312 cells, 7 conditions) | 6.1% | 4.2% | 4.8% |

- The pre-declared check (overall rate ≤ core + 1 pp) passes only because signature kinases are a small share of the cells.
- For the kinases that change, the fused z is over-dispersed: +1.5 pp overall and up to 2× nominal for single kinases.
- Likely cause: ρ_sig is estimated on only 7 null conditions and comes out at 0.68–1.00. Dividing by ρ < 1 inflates z_sig, which already has SD ≈ 1.15 on permuted data.
- The independence assumption in the fusion adds to this.
- This does not affect M1 ranking much, because the α plateau absorbs it. It does affect the "calibrated z" deliverable (desideratum 3) for about 40 kinases.

## 5. Sanity: where the gain comes from and false positives
**Concentration.**
- Per-kinase contributions to Δ (106 GT kinases): 56 gain and 44 lose.
- The top 5 carry **71%** of Δ: MAP3K7 +0.0034, ROCK2 +0.0026, MET +0.0023, ROCK1 +0.0023, IKBKB +0.0018.
- **Without these 5 kinases, Δ = +0.0054.**
- Largest losers: JAK2 −0.0011, ATM −0.0009, EGFR −0.0008, CDK2 −0.0006.
- By group (leave-one-group-out):
  - top contributors: H1152 (ROCK) +0.0039, cobimetinib+pictilisib +0.0010, MEKi +0.0006, CHIR99021 +0.0006, dinaciclib +0.0005;
  - Δ without the top 5 groups: +0.011;
  - largest negative groups: Y27632 −0.0019, ripasudil −0.0016 (both also ROCK inhibitors), dasatinib −0.0014, U0126 −0.0009.
- The ROCK signature generalises unevenly across ROCK inhibitors.

**False positives: signature kinases as evaluable negatives** in fold-f conditions (13,316 cells). t = s_c·z, with s_c the majority target sign.

| | core-v1 | B7b |
|---|---|---|
| FP rate, t > 1.96 (all conditions) | 8.1% | **9.8%** (Δ +1.7 pp; group-bootstrap 95% [+1.1, +2.4]) |
| FP rate, M1-sel conditions | 8.5% | 10.0% |
| FP rate, H / D / O | 3.1 / 10.5 / 9.9% | **6.3** / 11.3 / 11.1% |
| \|z\| > 1.96 | 14.7% | 18.8% |
| mean t | 0.07 | 0.02 (no directional push) |
| share of (negative, target) comparisons where the signature negative outranks the target | 31.2% | **28.8%** (M1-sel: 26.7 → 22.7%) |

- Per kinase, the largest FP increases are ERBB2 (9 → 24%), MTOR (11 → 21%), ALK, MAP2K2, MAPK1/3, PRKACA and MAP2K1.
- The biggest single false positives are pathway neighbours of the true target:
  - PRKCB t 10.7 under the GSK3 inhibitor SB216763 (study 32034161; its signature was trained on PKC inhibitors of the same study);
  - MAPK8 (JNK) under the p38 inhibitor SB203580;
  - GSK3A/B under U0126;
  - MET under AZD3759 (EGFR);
  - ROCK2 under ceralasertib.
- Masked or downstream kinases are, correctly, excluded. Signature kinases that are masked in a condition have mean t ≈ +0.8 to +1.0: the signatures are pathway readouts, not kinase-specific, which matches B7's 32/2,400 direct-substrate finding.
- Net effect: signature kinases are more often falsely "active" in absolute terms. But targets gain more, so negatives overtake targets less often. That is why M1 rises.
- Any downstream use that thresholds z for single kinases (for example "kinase X is inhibited") will see more pathway-neighbour false calls.

## 6. Recommendations
1. Adopt B7b as the next DEV candidate. The hard test and every pre-declared C10 check pass, they reproduce, and leakage is ruled out by direct tests.
2. Spend one DEV-confirm look on it. If DEV-confirm Δ vs core-v1 is ≤ 0, revert: the gain sits in a few kinases and is a second look at the same folds.
3. Record the calibration and FP caveats in the tree and the final report.
4. Any fix must be pre-declared before its effect is seen, for example ρ_sig estimated on site-permuted training conditions instead of 7 nulls, or a correlation-aware fusion variance.
5. At TEST, report T-seen, T-unseen and has-signature separately (C10 iv), and report the FP rate of signature kinases.
