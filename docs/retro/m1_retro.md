# M1 process retro (fresh reviewer, 2026-10-27)

Evidence: PROMPT.md, RESEARCH_TREE.md, LESSONS.md, git log (22 commits), docs/reviews (11 files, 650 lines), docs/audit (2 files, 334 lines), GT_REPORT.md, docs/harness.md. Process only, not science.

## 1. What should we stop doing

1. **Building the GT and split more than once.** The split went from seed 79 (3c639c3) to seed 90 (a4b8814). The audit re-ran on three GT versions, and the bootstrap was fixed after the harness was "frozen" (C8, ac7c7cc). The audit found 12 rule errors, 4 of them high, in the builder's from-memory tables. Auditing before building the split would have avoided the rebuild.
2. **Spending the whole milestone on design and infrastructure before any method result.** M1 produced 11 reviews, DESIGN v0 to v3, a 594-line LITERATURE.md and a 683-line data profile. Method work (branch 5) has not started and baselines (branch 4) are still [~]. Review round 1 alone cost 6 agents and produced A1-A10. The prompt caps design review at three rounds, but the result is far more machinery than the 2-3 key decisions needed. Cap the next gate at one reviewer pair and one round, with a hard budget.
3. **Writing documents nobody will reread.** data_profile.md is 683 lines. Harness verification is 182 lines, and the same facts appear again in the tree and in GT_REPORT (14 KB). The tree's status block repeats the finding "z-score about random, prior-only 0.687 is the bar" four times (summary, 2.2, 3.2c, 3.2 [!]). Write one source of truth and link to it.
4. **Tiny "tree" and "tree update" commits** (7 of 22). Fold them into the work commit that caused them.
5. **Leaving sealing by convention.** candidate_gt.tsv with future TEST labels is in git history (already in LESSONS). It is only a risk if an agent runs `git log -p`. Squash or rewrite history before any method agent starts, or state the ban in every method-agent prompt. A mechanism beats a promise (r1_aiagent says the same).
6. **Reviewing everything by default.** Expert audits (GT audit, harness verification) found real errors. Generic breadth reviews (r1_priorwork, r2_priorwork, r2_pi) mostly confirmed the top issues. Spawn specialist reviewers only where a builder lists from-memory judgement calls or a formula is novel.
7. **Over-engineered endpoint stack for a barely-above-random baseline.** Many endpoints (M1-sel, M1-all, scored variants, harmony bins, M2, bootstrap bias corrections) plus the ½ rule mean the PSP z-score sits at 0.515, which is 0.01 above random. The harness verification traced the gap from 0.78 to 0.52 to the metric design (the ½ rule alone costs 0.15). Keep the pre-registered primary, but report the coverage-vs-ranking split as a single headline table so decisions are not made on one compressed number.

## 2. What worked well (keep)

1. **Independent verification by someone who did not build it.** The from-scratch reimplementation matched kinact to 6 digits and reproduced the benchmarKIN reference of 0.775. The GT audit found 12 fixes (a 10.3% sample failure rate). These were the highest-value steps.
2. **Reviewers with different expertise in parallel.** Round 1 caught the M1 sign error, the downstream-mask problem and the prior-only control. The prior-only control reframed the whole objective.
3. **Pre-registration, a freeze hash and a frozen interface.** G1.json, the GT_FREEZE sha, `zero_policy="raise"` (NaN is not 0) and the sealed TEST cache stopped silent drift. The dev ledger and MDE (0.036 on DEV, about 0.06 on TEST) set an honest bar before any method existed.
4. **Reusing existing assets.** The kinaseSubMain audit (which sets are safe, which leak) and the benchmarKIN reference check saved weeks of rebuild.
5. **The research tree as memory,** with a status summary and numbered decisions D0-D4. Keep it, but shorter.
6. **Builders list their from-memory judgement calls.** The audit then targets those and does not re-review everything.

## 3. Proposed pruned LESSONS.md (12 lines)

```
# LESSONS — process rules (one line each, with evidence)

- Have a different agent independently re-implement or audit every headline number and label table before relying on it (reimplementation matched kinact to 6 digits and reproduced benchmarKIN 0.775; GT audit found 12 rule errors, 4 high, 10% of sampled conditions wrong).
- Audit label/rule tables BEFORE building splits and caches on them (split seed 79 -> 90 and three audit re-runs after the GT fixes, a4b8814).
- Ask builders to list "judgement calls made from memory", then audit only those (GT agent flagged 10 -> targeted audit instead of re-reviewing everything).
- Add a trivial null control (prior-only, random) to every metric before trusting a baseline (prior-only 0.687 vs z-score 0.515 reframed the objective; r2_stats).
- Have a statistician review metric formulas and resampling before coding them (r1 caught a sign error; C8 later fixed a bootstrap bias that put the CI bound at the point estimate).
- Run diverse specialist reviewers in parallel, one round (6 r1 reviews in ~7 min converged on the same top issues); use more rounds only for items that change labels or metrics.
- Never commit draft label tables before the split exists (candidate_gt.tsv with TEST labels is in history); seal by mechanism, not convention, and repeat the "do not read history/sealed" ban in every agent prompt.
- Freeze with a hash and a gate file, and write amendments before seeing their effect (GT_FREEZE sha 1e6f791c, G1.json at 7de8eef).
- Unscored must be NaN, never 0; make the harness raise on suspicious zeros (zero_policy="raise").
- Look for sibling project dirs and existing audited resources before building (kinaseSubMain at $KINASESUB_DIR, benchmarKIN local; their audit marked which sets leak).
- Tell concurrent agents when a shared cache is being regenerated, and require each to report the GT hash behind its numbers.
- Timebox infrastructure: get to a first end-to-end method number before the next design round; keep one source of truth per fact and fold tree-only updates into work commits (M1 spent all effort on design and GT, with no method result and baselines still in progress).
```
