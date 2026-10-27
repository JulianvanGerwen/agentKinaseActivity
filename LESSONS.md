# LESSONS — process rules (one line each, with evidence)

- Look for sibling project dirs outside the workspace before assuming a resource is missing (kinaseSubMain was at /cluster/project/beltrao/..., not /workspace).
- Have a statistician review metric formulas before coding them (r1 caught a sign error in M1 that would have rewarded wrong-sign scores).
- Parallel independent reviewers converge quickly: 6 r1 reviews in ~7 min wall-clock flagged the same top issues (downstream masks, promiscuous-drug labels) from different angles.
