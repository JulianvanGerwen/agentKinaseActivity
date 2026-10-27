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
- Look for sibling project dirs and existing audited resources before building (kinaseSubMain at /cluster/project/beltrao/..., benchmarKIN local; their audit marked which sets leak).
- Tell concurrent agents when a shared cache is being regenerated, and require each to report the GT hash behind its numbers.
- Timebox infrastructure: get to a first end-to-end method number before the next design round; keep one source of truth per fact and fold tree-only updates into work commits (M1 spent all effort on design and GT, with no method result and baselines still in progress).
