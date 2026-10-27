# Review r2 — statistics / benchmark design (DESIGN.md v1)

Reviewer: statistics agent, 2026-10-27.

Evidence: a quick check on candidate_gt.tsv, using direct, high/medium pairs with expected_direction ≠ 0 and ≤10 targets per condition.
- Scores were plain PSP z-scores (n≥3 or n≥1), computed from the compendium.
- The neutral ½ rule was applied with U = the PSP kinases.
- Groups were approximated as drug, or as study|treatment for other sources.
- Script: /tmp/claude-588982/-workspace/4cf38946-739b-40b2-907e-e597b4a480f7/scratchpad/pw*.py.

These numbers are indicative only.

## 1. Round-1 must-items

| r1 item | Status |
|---|---|
| 1 sign t = s·score | Addressed (A1). Make the oracle unit test (oracle −1 on inhibitor targets ⇒ M1 = 1) explicit in §9. |
| 2 unscored = ½ | Addressed (A2), including M1-scored and pair coverage. |
| 3 averaging order | Addressed (A3). |
| 4 coverage metric | Partly addressed. The permutation null (target rank uniform within each condition, group structure kept) and the "expected false discoveries" line are not written down. State them in §5. |
| 5 inference | Mostly addressed. Missing: the p formula (1 + #{Δ* ≤ 0})/(B + 1); that the whole metric, including the kinase set and bin membership, is recomputed on each resample; how groups with mixed sources are stratified. Add these to §5. |
| 6 DEV comparisons | Addressed in structure (pooled OOF, group bootstrap, ledger, nested CV). The rule itself changed and is now mis-calibrated; see N5. |

## 2. New or remaining problems

**N1 [must] §5: prior-only control must be in the confirmatory claim.**
- The quick check gives M1-sel ≈ 0.51 for PSP z-score (n≥3) and ≈ 0.55 (n≥1). The r1 prior-only score reached ≈ 0.63 on Hijazi.
- So "beats z-score/KSEA" on M1 can be won by a kinase-popularity prior, for example one that enters through τ²K, library size, or a learned offset.
- §5 says "must beat prior-only" but gives no test, no metric, and no definition of how sign is handled.
- Fix:
  - Add "final vs prior-only on M1-sel" as test 6 in the Holm family.
  - Add a second control, **shift-prior**: score_k = freq_k(DEV-train) × mean(x_c). It uses the condition's direction but no kinase-specific data.
  - Define prior-only as freq_activated − freq_inhibited in the training folds. It is blind to condition type and is computed on DEV only for TEST.

**N2 [must] §4/§5: split targets and power must refer to the primary endpoint.**
- §4 balances "≈30% of GT pairs". Hijazi accounts for 2,286 of 3,647 pairs overall but only 399 of 1,760 M1-sel pairs.
- Fix: apply every split-selection criterion to M1-sel pairs and kinases, by source and by bin. Keep M1-all as a secondary check.
- Pre-register how the MDE is computed: MDE = (z_{1−α/6} + z_{0.8}) · SE_boot(Δ) on TEST. SE_Δ comes from two DEV-available baselines with different coverage (z-score vs a full-coverage motif baseline), projected to the chosen TEST groups.
- Check (30% of M1-sel pairs; about 89 kinases and about 100 groups):
  - SE_Δ = 0.007 for two z-score variants that differ in coverage.
  - SE_Δ ≈ 0.017 for a simulated full-coverage noisy method vs z-score.
  - That gives **MDE ≈ 0.02–0.05**. δ_min = 0.01 therefore never binds.
- State in §8: if the Holm test fails, the result is reported as an estimate with its CI and no superiority claim. Do not enlarge TEST after looking.

**N3 [must] §5: the M1-sel inclusion rule is underspecified.**
- Define n_targets(c) as the number of annotated direct positives in U, counted after the decryptM Kd/dose rule and before off-target/paralog masking. It must be label-only and frozen with the GT hash.
- Exclude conditions with 0 positives, i.e. the specificity controls.
- Assert expected_direction ∈ {−1, +1}: candidate_gt has 3 direct pairs with direction 0.
- Decide whether a decryptM drug can move in and out of M1-sel across doses, and write the decision down.

**N4 [should] §5: M1-sel mostly measures coverage.**
- With the ½ rule, AUROC_ck = ½ + f_c(A_ck − ½), where f_c is the fraction of negatives the method scores.
- In the check, 50% of M1-sel kinases sit at exactly ½ for z-score and the per-kinase SD is 0.031. A full-coverage, weakly accurate method will beat an accurate PSP method.
- This matches desideratum 2, but the accuracy claim then rests on a single M2 test.
- Fix:
  - Report the decomposition (pair coverage, M1-scored) next to every M1.
  - For "beats", require the M2 point Δ ≥ 0 against every baseline in the family. This is descriptive and adds no tests.

**N5 [should] §7: the adoption rule is both too strict and too lax.**
- Too strict, per idea:
  - A 95% two-sided lower bound > 0 is one-sided α = 0.025.
  - DEV SE_Δ is about 0.004–0.012, so each idea needs Δ ≈ 0.01–0.03.
  - Real increments of 0.005 can never be adopted.
  - On top of that, 4 conjunctive conditions apply, several with low power. Source D has only 24 kinases and 13 drugs.
- Too lax overall: no multiplicity control across K ideas.
- Fix:
  - Use one-sided p < 0.05 on ΔM1-sel.
  - Turn ΔM1-all, harmony bins and T-unseen into point non-inferiority checks with margin −0.02 (agrees with r2_pi N1).
  - Make source-sign agreement descriptive only.
  - Every 5 adoptions, test the cumulative bundle vs the last checkpoint on DEV-confirm. This counts against the 5-use budget.
- The current harmony rule (bin CI upper bound < 0) almost never fires, so it gives no protection.

**N6 [should] §5: assign harmony bins per (kinase, group).**
- With bins per pair, a kinase can fall into several bins.
- Fix: bin = the median number of quantified PSP substrates over the kinase's conditions within the group. Bin-macro = mean over bins of the kinase means within each bin. Freeze this with the GT hash.

**N7 [should] §5: define the bootstrap strata on the primary endpoint.**
- Resample only groups that contain M1-sel pairs for M1-sel, and all GT groups for M1-all.
- Use separate resamples per endpoint, shared across methods, with seed 0.
- Groups with mixed sources go into the stratum of the majority source by M1-sel pairs.

**N8 [could] §5 Holm family.**
- Define "best-other-baseline" as the best baseline not already in the family.
- If native JUMPsem fails, the validated reimplementation is the comparator and is labelled as such.
- With N1 the family has 6 tests.

**N9 [could] §5 coverage-significant.**
- 60 of 151 M1-sel kinases come from a single group, which caps their power.
- Report the count per bin and the expected false discoveries. Do not use it for adoption.
