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
