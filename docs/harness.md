# Evaluation harness `kinact` (frozen at G1) — DESIGN v3 §5

Install: `/home/plm_env/bin/python -m pip install -e /workspace`. Tests: `python -m pytest -q tests/test_eval.py`.

## Inputs (docs/interfaces.md)
- `scores`: DataFrame conditions × kinases, larger = more active. **Unscored = NaN, never 0.**
  `Evaluator(..., zero_policy="raise")` (the default) raises if a condition has no NaN and more than 50% exact 0.0 (a likely NaN→0 bug).
  It also warns when a kinase's score is constant across conditions. Controls that emit legitimate zeros (oracle, prior-only) use `zero_policy="ignore"`.
- `gt_pairs` (conditionID, kinase, cls, sign, confidence), `gt_conditions` (conditionID, source, group_id, n_targets, in_M1sel, role, split), universe U.
  `load_gt()` reads `results/ground_truth/` (DEV + devconfirm) via `read_gt_tsv` (keep_default_na=False: role "null" stays a string). It accepts `gt_pairs.tsv` or `gt_pairs_dev.tsv`, plus `harmony_bins.tsv` if present (`harmony_bins_group.tsv` can be passed as `bins` too). The harness never reads `cache/sealed/draft_candidate_v0/`.
- Harmony bins:
  - `harmony_bins.tsv` keyed by (kinase, group_id) with a frozen `bin` column; or
  - per (conditionID, kinase) with `n_psp_quantified`. The bin is then computed per (kinase, group) as the median over the group's positive conditions, with bins 0 | 1–4 | 5–19 | ≥20.

## Formulas
- t = s·score. Evaluable negatives N_c = U minus every kinase listed for c in gt_pairs (direct of any confidence, downstream, masked). S_c = the kinases scored in c.
- Pair AUROC for target k:
  - if k ∈ S_c: [Σ_{j∈N_c∩S_c}(1[t_j<t_k] + ½·1[t_j=t_k]) + ½·|N_c∖S_c|] / |N_c|
  - if k ∉ S_c: ½
  - `auroc_scored` = wins / |N_c∩S_c| (scored target, scored negatives only).
- Eligible conditions: role=='gt', n_targets>0, and in split_filter. Targets: cls=='direct', confidence ∈ {high, medium}, sign ∈ {−1,+1} (asserted).
- Averaging: pairs → (weighted) mean per (kinase, group) → mean over the kinase's groups → macro over kinases.
- Endpoints (string grammar `base[:filter,...]`):
  - `M1-sel` (PRIMARY): pairs from conditions with in_M1sel, unweighted.
  - `M1-all`: all pairs, weight 1/n_targets(c) inside the (kinase, group) mean.
  - `M1-sel-scored` / `M1-all-scored`: `auroc_scored` over pairs with a scored target and ≥1 scored negative. Pair coverage = fraction of pairs with a scored target.
  - `harmony`: M1-sel per bin (bin per (kinase, group)), kinase means within each bin, then mean over non-empty bins (bin-macro).
  - Filters: `source=H|D|O`, `sign=-1|1`, `unseen` (kinases ∉ `train_pos_kinases`), `bin=1-4`.
  - `M2` (A vs B): mask every score outside S_c(A)∩S_c(B), then compute M1-sel-scored for both on the surviving pairs.
- Coverage-significant (per kinase on M1-sel):
  - Null: the target's rank is uniform on {0..|N_c∩S_c|}; unscored negatives stay at ½ and unscored targets at ½. The null goes through the same hierarchy, with n_perm=2000 and seed 0.
  - p = (1 + #{null ≥ obs}) / (n_perm + 1). Kinases with no scored pair get p = 1.
  - BH across GT kinases; significant = q ≤ 0.10 and AUROC > ½. Also reported: expected false discoveries (q·n_sig) and the count per bin.
- Calibration: `calibration(scores, null_conditions, n_substrates=None)` returns:
  - the |z| > 1.96 rate overall, per kinase, and per substrate bin (taken from the method's n_substrates frame);
  - the KS distance of the pooled two-sided p-values from U(0,1).

## API
```python
from kinact.eval import Evaluator, load_gt, calibration, load_sealed_test
from kinact import bootstrap as bs, controls as ctl, ledger
P, C, U, bins = load_gt()
ev = Evaluator(P, C, U, split_filter=["fold0", "fold1", "fold2", "fold3", "fold4"],
               bins=bins, train_pos_kinases=train_pos)  # raises TestSplitLocked if 'test' selected w/o FREEZE.json
res = ev.evaluate(scores)          # dict: M1-sel, M1-all, M1-*-scored, coverage_*, tie_share_M1-sel,
                                   # by_source, by_source_all, by_sign, harmony, T-unseen, coverage_significant,
                                   # per_kinase (DataFrame), pairs (DataFrame)
ev.m2(scores_a, scores_b)          # matched-universe M2 point estimate + n_pairs
plan = bs.BootstrapPlan.from_evaluator(ev, B=2000, seed=0)   # build ONCE, reuse for all methods/endpoints
d = bs.paired_delta(ev, scores_final, scores_z, "M1-sel", plan)  # estimate, se, ci_low/high, p_one_sided
bs.paired_delta(ev, a, b, "M2", plan); bs.bootstrap_ci(ev, a, "harmony", plan)
bs.holm([p1, ..., p6]); bs.beats(d, holm_p, dM2, dUnseen)
ledger.append("zscore", config_dict, "dev", ledger_metrics, branch="1-library", folds=range(5))
```
- `evaluate(scores, gt_pairs, gt_conditions, split_filter, universe, bins, ...)` is the one-call contract wrapper.
- `Evaluator(..., positive_cls="downstream")` gives the secondary downstream-positive AUROC (descriptive).

## Bootstrap (`kinact.bootstrap`)
- Unit = GT group. Stratum = the group's majority source by M1-sel pairs (by all pairs if it has none).
- One count matrix W (B × groups) with seed 0, shared by all methods and endpoints (B6). Groups without pairs for an endpoint contribute nothing.
- Because whole groups are resampled, (kinase, group) means are fixed. Each replicate is then exact and vectorised (checked against brute-force duplication in the tests):
  - kinase_b = Σ_g W_bg·v_kg / Σ_g W_bg, and a kinase is present iff the denominator > 0;
  - metric_b = macro over present kinases (for harmony: within each bin, then over bins).
- Paired Δ = A − B on identical resamples. Inference is bias-corrected (DESIGN §12 C8). The group bootstrap shifts the kinase composition of the macro average, because single-group kinases drop out of about 37% of resamples:
  - CI (`ci_low`/`ci_high`) = basic interval [2Δ̂ − q97.5, 2Δ̂ − q2.5];
  - one-sided p (`p_one_sided`) = shifted-null (1 + #{Δ*_b − mean(Δ*) ≥ Δ̂}) / (B + 1);
  - diagnostics: `ci_low_pct`/`ci_high_pct` (raw percentile interval), `p_one_sided_pct` (old (1 + #{Δ* ≤ 0}) / (B + 1)), `mean_replicates`, and `bias` = mean(Δ*) − Δ̂;
  - SE = sd(Δ*), which the location correction does not change.
- `holm`: step-down adjusted p. `beats` = Holm p < .05 AND ΔM1-sel ≥ 0.01 AND ΔM1-sel(T-unseen) > 0 AND ΔM2 ≥ 0.

## Controls (`kinact.controls`)
- `random_scores(conds, kinases, seed, coverage)`: expected M1 = ½.
- `oracle(P, conds, kinases, anti=False)`: M1 = 1 (anti-oracle: 0).
- `training_frequency(P, C, train_conds, U)`: freq_k = fraction of training GT conditions with k as a direct positive. Refuses TEST conditions.
- `prior_only(P, C, train_conds, eval_conds, U)`: s_c·freq_k, where s_c = majority sign of c's direct positives (a tie gives +1).
- `shift_prior(X, freq)`: freq_k·mean_sites(x_c).
- `site_permutation(X, seed, mode="within_condition"|"rows")`.
- `shuffled_library(lib, seed, mode="relabel")`: a random bijection on site labels, which preserves per-kinase counts, weights and overlaps.
- `sign_flip_check(fn, X, lib)`: requires score(−X) = −score(X) and the same NaN pattern.

## Guards / ledger / runtime
- TEST guard: any Evaluator whose selected conditions include split=='test' raises `TestSplitLocked` unless FREEZE.json exists (`KINACT_FREEZE_PATH` or `/workspace/FREEZE.json`).
- `load_sealed_test(cache/sealed)` additionally requires FREEZE.json `gt_hash` == `sealed_gt_hash()`. That hash is sha256 over the sorted lines `name\tsha256(file)` of `cache/sealed/*.tsv`.
- Ledger `results/dev_ledger.tsv`:
  - columns: timestamp, git_commit (+dirty / nogit), method, config_hash (sha256[:16] of sorted JSON), branch, split, folds, n_configs_tried, metrics_json;
  - n_configs_tried = distinct config hashes so far in the branch;
  - written with O_APPEND + flock and never rewritten; it refuses TEST.
- `runtime.bench_runtime(method, X, lib, n_rep=5, threads=(1, 16))`:
  - one warm-up, then the median of n_rep runs, under threadpoolctl limits;
  - if the method has `prepare` / `score` (/`null`), load/mapping, scoring and null are timed separately;
  - reports ms/condition and an extrapolation to 1,937 conditions. This is reported, not gated.

## Timings (64-core node; single process)
- 1000 conditions × 500 kinases (synthetic): GT indexing + plan + paired ΔM1-sel at B=2000 takes 0.3 s.
- Synthetic-scale GT (780 conditions, 3,644 pairs, 371 groups; measured on the retired draft layout): full `evaluate` with n_perm=2000 takes 0.18 s; paired Δ (M1-sel, M1-all or M2) at B=2000 takes about 0.06 s each.
