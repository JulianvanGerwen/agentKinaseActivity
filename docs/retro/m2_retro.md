# M2 process retro (fresh reviewer, 2026-10-28)

Scope: process only. Read PROMPT.md, RESEARCH_TREE.md, LESSONS.md, m1_retro.md, git log (68 commits), results/b1..b8 reports, results/g2_prep/g2_prep.md, results/test/LEDGER.md, docs/audit/*, docs/retro/ideas_m2.md. Did not read ground_truth history or cache/sealed.

## Stop doing

1. **Running method branches whose best case sits below the MDE.** Preregistered TEST MDE is 0.06 (G1.json). Of 8 branches, B4, B4b, B5, B6 and B7 ended at -0.010 to +0.012 or n.s. (b4 -0.010, b4b +0.007 p 0.16, b5 +0.003, b6 +0.012 p 0.057, b7 +0.012 p 0.15). Each cost a 100-line report, a tree commit and a D-note. The adopted gains were small: B2.1 +0.015, B7b +0.018, E1 -0.0007. The final model beat core-v1 by only +0.021 on TEST [0.001, 0.031], and it was worse on the external Hernandez set (0.552 vs 0.613). The fresh-reviewer verdict "diminishing returns, remaining ideas <= +0.01 << 0.06" (D10, ideas_m2.md) came after B7b. Ask for that verdict after the first two null branches, and set an adoption floor of at least MDE/3 before opening a branch.
2. **Opening a box with a calibration fix that did not show up on real data.** B8-E1 moved site-permuted calibration (signature kinases 7.0% to 5.1%, dM1-sel -0.0007). On TEST real nulls, core-v3 and core-v1 are identical at 8.1% (LEDGER calibration table). Site-permuted nulls looked fine while real nulls were off by 3 points (8.1% vs the 5% target). Validate calibration on the real-null set that will be reported, or stop describing permuted-null calibration as a result.
3. **Freezing the harness without edge-case tests.** run_test_review.md found that `metric_replicates` returns 0 (not NaN) for empty subsets, which gave CIs like "+0.494 [+0.494, +0.987]" (J7). The harness was G1-frozen, so it could only be patched around. It also found that frozen inputs were never re-hashed. Both surfaced only in the pre-TEST review of the runner. Run an empty-stratum and re-hash test at G1, not at G2.
4. **Leaving external-set validation to the sealed run.** `run_external` was dry-run with random/oracle only (2e2d719). The first real signal that core-v3 loses to core-v1 on Hernandez (-0.061) came at TEST, after development was closed (D14). A DEV-time external run on the non-sealed sets would have shown whether the signature fusion transfers.
5. **Tiny "tree:" commits.** 24 of 68 commits are "tree ..." (c.f. 7 of 22 in M1, so the habit grew in M2). M1 retro asked to fold them into work commits. Do it, or batch tree updates once per branch.
6. **Duplicated per-branch outputs.** b1 has both b1_report.md and b1_report_auto.md; b4, b4b, b21 have separate *_tables.md; the tree restates every branch number in four places (summary, tree line, decision log, report). Keep one report per branch and have the tree link and give only the verdict and delta.
7. **Fighting the prior-only control for the primary claim.** The control was already known from M1 to be about 0.69 on DEV (D5 banned the only fix). TEST gave +0.002, p 0.50, and no claim. The honest outcome was correct, but the endpoint was kept as the Holm family member while three branches tried to close the gap (B4, B6, B7). The evidence that mattered was T-unseen (0.557 vs 0.253) and harmony. Pre-declare at G1 which stratum is the decision metric when a popularity control is near ceiling, and stop tuning the pooled number.
8. **Post-development descriptive analyses.** 7c3aaec (M3, M4, GT-A, kinase-kinase, residue null) was run after B8 and not used for any decision. Some are useful for REPORT, but queue them as REPORT-time tasks, not milestone tasks.

## Worked well

1. **Pre-declared boxes written before any number.** B7b and B8 each open with "written before any number was computed; not changed afterwards" plus adoption criteria and a one-comparison inference. B8 closed in one box with 1 adoption and 2 rejections (E2 p 0.80, E3 p 0.73). C10 was pre-registered before B7 results (b83a3cd). The final model has no post hoc component except the stated ρ_f cap.
2. **Refusing near-misses.** D9 declined to adopt A1 (p 0.057) by relaxing the rule, and registered it as a secondary test instead. On TEST A1 is 0.686, below core-v3 (0.714) and core-v1 (0.693; -0.007 vs core-v1). A relaxed rule would have shipped a loser.
3. **Conditional adoption with a confirm look.** D11 adopted core-v2 "conditional on DEV-confirm look #2, revert if delta <= 0". Look #2 gave +0.039 [0.006, 0.069], and TEST gave a consistent +0.021 with the gain on signature kinases (+0.050; no-signature -0.001). Only 2 of 5 looks were used.
4. **Pre-declared secondary criteria catching a trap.** The T-unseen non-inferiority criterion caught B7 (TYK2 borrowed through the JAK kernel, T-unseen -0.026) and forced the late-fusion redesign (B7b: T-unseen +0.0006, shuffled control -0.002 to -0.008).
5. **Re-split stability check.** Seeds 1-4 reproduced core-v3 0.6800 +/- 0.0037 vs 0.6796 on seed 0, and the reviewer flagged the weakness (ARI 0.68-0.75, so the partitions overlap). This ruled out fold-selection inflation cheaply (g2_prep.md).
6. **Independent audit and reproduction at every adoption.** core_v1_audit (exact reproduction, no leakage), b7b_audit (byte-identical; found that 71% of the gain comes from 5 kinases, and over-dispersed z), run_test_review (5 fixes, input re-hash, 162 tests), and TEST reproduced bit-identically from 03ab343 by a second agent. The b7b audit's caveats directly produced the E1 box.
7. **Freeze file pinned to a commit, with run-once guard.** FREEZE.json pins code_commit 03ab343, config, manifest and GT hashes; the run refuses if the partial file exists; the result and freeze were committed together (be2ef15). There was no sealed-data incident.
8. **Null controls shaped the claims.** prior-only, shift-prior, random, shuffled signatures and the residue-matched random library (0.529 vs real 0.668) kept the report honest: claims were restricted to four Holm comparators and an explicit "no claim vs prior-only".
9. **A fresh-reviewer go/no-go ranking before the last box** (ideas_m2.md) gave ranked ideas with P(pass) and cost, producing a single time-boxed B8 (<= 3 h) rather than open-ended exploration.
10. **Small nested search spaces.** The "bigger grids lose ~0.03 to selection noise" lesson (B2.1) kept later grids small; the nested M1-sel figures (0.680) were close to the TEST result.

## Proposed LESSONS.md (12 rules)

```
# LESSONS — process rules (one line each, with evidence)

- Have a different agent independently reproduce or audit every headline number, adoption and the final runner before relying on it (TEST reproduced bit-identically from 03ab343; b7b audit found 71% of gain from 5 kinases; run_test_review fixed 5 bugs; GT audit found 12 rule errors, 4 high).
- Audit label/rule tables BEFORE building splits and caches on them (split seed 79 -> 90 and three audit re-runs, a4b8814); have builders list "judgement calls from memory" and audit only those.
- Put a trivial null control (prior-only, random, shuffled, residue-matched) next to every metric and in the test family, and restrict claims to what beats it (prior-only 0.712 vs core-v3 0.714 on TEST -> no claim; residue-matched library 0.529 vs 0.668).
- Pre-declare each experiment box (hypothesis, one comparison, adoption rule, time box) before any number exists, and never relax a rule after seeing a near-miss (B7b/B8 headers; A1 p 0.057 not adopted and later 0.686 < 0.714 on TEST).
- Adopt marginal gains only conditionally, with a pre-declared confirm look and revert rule (D11: core-v2 conditional; look #2 +0.039; TEST +0.021); add non-inferiority gates on the strata that expose failure (T-unseen caught B7).
- Stop branching when the best remaining idea is below MDE/3: get a fresh-reviewer go/no-go after two null branches, not after five (B4-B7 gave -0.010..+0.012 vs MDE 0.06; ideas_m2 verdict came at D10).
- Freeze with a hash file pinned to a code commit, re-hash every input file the run reads, and make the run-once guard refuse a re-run (FREEZE.json at 03ab343; input re-hash added only in run_test_review).
- Check selection inflation by re-splitting with several seeds and keep search grids small (seeds 1-4: 0.6800+/-0.0037 vs 0.6796; bigger grids lose ~0.03 to selection noise).
- Validate calibration and transfer on the same real data you will report, before freezing (E1 fixed permuted-null calibration but real nulls stayed 8.1%; Hernandez core-v3 0.552 < core-v1 0.613 seen only at TEST).
- Unit-test the frozen harness for empty strata, NaN vs 0 and re-hash before G1 (metric_replicates returned 0 for empty subsets, J7; zero_policy="raise" prevented the silent-zero class).
- Seal by mechanism and repeat the "do not read history/sealed" ban in every agent prompt; never commit draft label tables before the split exists (candidate_gt.tsv in history; no sealed incident in M2).
- Keep one report per branch and one source of truth per fact; fold tree-only updates into the work commit (24 of 68 commits are "tree:"; b1 has duplicate reports; same finding stated in 4 places).
```

Dropped from the old list (low value or merged): "look for sibling dirs" (one-time, done), "tell concurrent agents about cache regeneration" (not triggered in M2), "statistician reviews metric formulas" and "parallel specialist reviewers, one round" (design-phase only, merged into the audit rule), "NaN not 0" (merged into the harness-test rule), "freeze with hash and gate" (merged into the pinned freeze rule).
