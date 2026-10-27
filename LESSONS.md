# LESSONS — process rules (one line each, with evidence)

- Look for sibling project dirs outside the workspace before assuming a resource is missing (kinaseSubMain was at /cluster/project/beltrao/..., not /workspace).
- Have a statistician review metric formulas before coding them (r1 caught a sign error in M1 that would have rewarded wrong-sign scores).
- Parallel independent reviewers converge quickly: 6 r1 reviews in ~7 min wall-clock flagged the same top issues (downstream masks, promiscuous-drug labels) from different angles.
- Never commit draft label tables before the split exists: candidate_gt.tsv (incl. future TEST labels) is in git history (commit before G1); sealing is therefore by convention — agents must not read git history of results/ground_truth.
- Ask builders to list "judgement calls made from memory" explicitly (GT agent flagged 10 → targeted expert audit instead of re-reviewing everything).
