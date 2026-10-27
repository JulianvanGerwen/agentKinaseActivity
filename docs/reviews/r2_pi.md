# Review R2 — external PI (strategy, scope, priorities)

Reviewer: external PI agent, 2026-10-27. Read: PROMPT.md, DESIGN.md v1 (§12), r1_pi.md, results/ground_truth/*. Read-only check [chk] below.

## 1. Round-1 must-items
- R1-1 (metric dominated by promiscuous drugs): **addressed** (A4, M1-sel with ≤10 targets; M1-all with 1/n weighting).
- R1-2 (curated and in vitro KS sources): **addressed** (§2). Phosformer is on disk but not mentioned: either list it as predicted (it needs a leakage audit) or state why it is cut.
- The r1 should-items are addressed, with three gaps:
  - Specificity controls appear in §3 but not in the §5 calibration report.
  - The giant-group check (R1-10) is not stated.
  - The kinase–kinase analysis has no gold set or metric (see N4).

[chk] On candidate_gt (direct, high/medium):
- M1-sel keeps 151 of 320 kinases and 1,760 of 3,647 pairs.
- 64 kinases come from a single treatment, so their per-kinase AUROC rests on one drug group.
- Only 254 pairs have s=+1.
- PIK3CD is among the top-10 positives, so the U filter must remove it.

None of this blocks the design, but it is why the harmony and coverage metrics carry the understudied-kinase claim.

## 2. New or remaining problems
N1. **§7 adoption rule is too conjunctive and will block real gains** (must).
- Five simultaneous conditions are required: CI lower bound > 0, ΔM1-all ≥ 0, sign agreement in 2 of 3 sources, no harmful harmony bin, T-unseen gain.
- With ~150 M1-sel kinases (~100 in DEV), many of them from a single drug, every condition beyond the first is underpowered. Requiring all of them at once biases the project toward the baseline. That works against "best possible estimator".
- Fix: primary = 95% CI lower bound of ΔM1-sel > 0. All other conditions become non-inferiority checks with margin −0.02 on the point estimate, not CI-based. T-unseen: point estimate ≥ −0.01.
- Keep the ledger and the audit of large gains.

N2. **No prioritised build order or time boxes for method branches** (must).
- §7 lists A–E as equals. LITERATURE and benchmarKIN show that the library is the dominant lever, so effort should follow it.
- Fix: add an ordered plan with a time box for each step:
  1. Library builder, with z-score/mlm/ulm on each library variant. This step also answers "which library is best" for all methods.
  2. Bayesian core B on the best library.
  3. Calibration.
  4. C: learned edge-feature weights and co-regulation expansion.
  5. D: self-sites.
  6. E: stacking, only if 1–4 leave ≥2 complementary methods with low score correlation.
- Kill rule: a branch with no DEV gain after its time box is marked [-].

N3. **Metric zoo is over-engineered for the harness v1** (should).
- §5 lists about 15 metrics. Building and verifying all of them before any method runs delays G1.
- Fix: frozen harness v1 = M1-sel, M1-all, harmony, coverage-significant, M2, calibration (N-eval + specificity controls), M3, and the controls.
- Defer to the post-G2 report (descriptive only, never used for selection): by-sign, M1-family, downstream, M4, GT-A, M1-scored. Defining them is still allowed before TEST.

N4. **Kinase–kinase deliverable is undefined** (should).
- Fix: gold set = PSP edges where the substrate site is a kinase's "enzymatic activity, induced" site, plus a canonical cascade list.
- Metric = AUROC of |Spearman| between atlas scores over non-TEST conditions.
- Compare against z-score and against a substrate-overlap-only control. The control is essential, because shared substrates create spurious correlation.
- Keep it to one analysis and one figure.

N5. **GT-A definition is inconsistent** (should). §3 describes an activity proxy, while §5 says "GT-A AUROC".
- Fix: for each kinase, Spearman between the inferred score and the mean log2FC of its own activating sites, across non-inhibitor conditions within study, averaged over studies. Report the macro mean against a site-permutation null.
- The kinase's own sites must be removed from the library, and so must sites on the same peptide.

N6. **Redundant site evidence inflates z and n for some kinases** (should, §2 shared builder). Two sources:
- Mouse/rat orthology maps several edges to the same human site.
- Multi-site peptides report the same log2FC several times.

Fix: collapse edges per (kinase, human site), and down-weight sites from identical-vector peptides (weight 1/n_sites with the same value in the condition). This affects every method equally and is cheap.

N7. **The ≤10-target filter is asymmetric across sources** (could, §5).
- Hijazi n_targets comes from the DiscoverX profile, while decryptM comes from a curated primary list.
- State that n_targets = direct positives in U after masking, as built by the GT script, and freeze it.

N8. **Giant-group check** (could, §4): log the size of the largest union-find group, and treat co-stimuli (EGF) as non-linking if the largest group holds >25% of conditions.

No other changes needed. §4 sealing, §8 TEST protocol, §9 and §10 are adequate.
