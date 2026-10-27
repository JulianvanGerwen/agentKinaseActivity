# G2 prep — fold-sensitivity (winner's-curse) check + DEV-confirm look #2 (2026-10-27)

Scripts: `scripts/g2_prep_folds.py` (re-split), `g2_prep_resplit.py --seed s` (full nested pipeline), `g2_prep_aggregate.py`, `g2_prep_devconfirm.py` (look #2). Raw outputs in `results/g2_prep/`: `resplit_seed{0-4}.json`, `resplit_per_seed.tsv`, `resplit_deltas.tsv`, `resplit_summary.json`, `devconfirm_look2{_frozen,}.json`.

GT_FREEZE `1e6f791c…`. I never opened `cache/sealed/`, the git history of `results/ground_truth`, or `external/`. I made no commit.

## 1. Re-split check (DESIGN §4)

**Folds.** I reimplemented the fold step of `make_split.py` (step 3) from the agent-visible GT files. With `rng = default_rng(s)`, the DEV-confirm greedy shuffles are replayed and their results discarded, so membership is unchanged. Groups are then shuffled, stable-sorted by (M1-sel pairs·1000 + size), and assigned to the fold with the least load in their stratum. Seed 0 reproduces the frozen folds exactly (asserted); TEST and DEV-confirm are untouched.

**Caveat: seeds 1–4 are only partly different partitions.** The size sort is deterministic, so the seed only reorders groups with tied keys. The big groups always land in the same fold: G0005 (194 pairs) always goes to fold0, which is why fold0 has 323 M1-sel pairs in every seed. Adjusted Rand index vs seed 0 is 0.68–0.75 on M1-sel pairs and 0.39–0.58 on GT conditions. This is the pre-registered logic, but it reduces how much the check can detect.

**Pipeline per seed (all nesting redone):** core-v1 z_vif grid (45 configs) → B7 joint grid (clean+contrast) gives (n_top, κ) → B7b late-fusion α → E1 permutation scales give core-v3. Seed 0 reproduces B8 bit-exactly: core-v1 0.66274, core-v2 0.680308, core-v3 0.679572.

### M1-sel (nested OOF, 877 pairs)
| method | s0 | s1 | s2 | s3 | s4 | mean±sd (0–4) | s0 − mean(1–4) |
|---|---|---|---|---|---|---|---|
| **core-v3** | 0.6796 | 0.6810 | 0.6751 | 0.6791 | 0.6854 | **0.6800±0.0037** | −0.0006 |
| core-v2 | 0.6803 | 0.6809 | 0.6749 | 0.6789 | 0.6855 | 0.6801±0.0038 | +0.0002 |
| core-v1 | 0.6627 | 0.6631 | 0.6642 | 0.6618 | 0.6668 | 0.6637±0.0019 | −0.0013 |
| B7 joint (nested; rejected form) | 0.6747 | 0.6622 | 0.6587 | 0.6646 | 0.6755 | 0.6672±0.0076 | +0.0094 |
| PSP z-score min3 | 0.5153 | = | = | = | = | 0.5153 | 0 |
| curated z-score min1 (= nested z lib/min, curated min1 in every fold) | 0.5634 | = | = | = | = | 0.5634 | 0 |
| prior-only (OOF freq) | 0.6875 | 0.6960 | 0.6861 | 0.6820 | 0.6768 | 0.6857±0.0071 | +0.0023 |

### Other endpoints (mean±sd over seeds 0–4; seed 0 in brackets)
| method | M1-all | harmony | T-unseen |
|---|---|---|---|
| core-v3 | 0.6279±0.0037 [0.6229] | 0.6819±0.0050 [0.6822] | 0.6423±0.0157 [0.6482] |
| core-v1 | 0.6211±0.0038 [0.6147] | 0.6640±0.0018 [0.6631] | 0.6422±0.0163 [0.6478] |
| PSP z min3 | 0.5122 | 0.5169 | 0.5127±0.0029 |
| curated z min1 | 0.5429 | 0.5652 | 0.5690±0.0077 |
| prior-only | 0.5438±0.0085 | 0.7245±0.0086 | 0.3031±0.0172 |

T-unseen varies by seed because "seen" depends on the folds.

### Δ core-v3 − comparator, M1-sel per seed
Paired C8 bootstrap, B = 2000. Cells show the estimate with the one-sided p in parentheses.

| vs | s0 | s1 | s2 | s3 | s4 | mean±sd |
|---|---|---|---|---|---|---|
| core-v1 | +0.0168 (.010) | +0.0179 (.018) | +0.0109 (.071) | +0.0173 (.014) | +0.0186 (.006) | +0.0163±0.0031 |
| PSP z min3 | +0.164 | +0.166 | +0.160 | +0.164 | +0.170 | +0.165±0.004 (all p .0005) |
| curated z min1 | +0.116 | +0.118 | +0.112 | +0.116 | +0.122 | +0.117±0.004 (all p .0005) |
| prior-only | −0.008 | −0.015 | −0.011 | −0.003 | +0.009 | −0.006±0.009 (all n.s.) |

Other endpoints, mean Δ over seeds 0–4: vs core-v1 M1-all +0.007, harmony +0.018, T-unseen +0.000 (range −0.0006 to +0.0015); T-unseen vs PSP-z +0.130, vs curated-z +0.073; vs prior-only M1-all +0.084, harmony −0.043, T-unseen +0.339.

**Reading.**
- **No measurable winner's curse on the core-v3 path.** Seed-0 M1-sel (0.6796) sits at the mean of seeds 1–4 (0.6802); the sd over fold seeds is 0.004.
- **The gain over core-v1 replicates in direction in every seed** (+0.011 to +0.019). One seed (s2) misses p < 0.05 at 0.071.
- **The joint-posterior B7, which was rejected, did show optimism** (+0.009 at seed 0). That is consistent with its rejection.
- **Limits of this check:** the baselines are fixed configs with no training part, so their M1-sel cannot change with folds (nesting the z-score library/min choice picks curated min1 in every fold of every seed). Design choices made on seed-0 folds were not re-done (ρ_f cap, signature variant, late vs joint fusion, E1). Fold overlap with seed 0 is high.
- **Seed-to-seed spread** of Δ vs PSP-z is about 0.004. That is much smaller than the TEST MDE of 0.06.

## 2. DEV-confirm look #2 (logged; 3 looks remain)

**Frozen before loading DEV-confirm** (`devconfirm_look2_frozen.json`, hash `46caff727af577b3`):
- **Core:** the look-1 core-v1 config (hash `04b69b65…`, τ .125, ρ_f 3, ρ_g .25); ρ_e refitted on the 7 DEV-CV nulls = 0.01, identical to look 1. Prior-only frequencies from DEV-CV.
- **Signatures:** clean+contrast, fitted on all DEV-CV GT conditions (463 training conditions, 55 signature kinases); groups disjoint from DEV-confirm (asserted).
- **HPs:** (n_top, κ) = (100, 0), argmax of pooled seed-0 OOF in the B7 grid; α = 1, argmax of pooled seed-0 late-fusion OOF (0.6820 vs 0.6791 at α .5); ρ_sig = 0.898 from the DEV nulls; E1 v_k from site-permuted DEV-CV GT conditions (median 2.16, range 1.62–2.96).
- **Dry run:** the post-load code path was first tested with DEV fold0 standing in for DEV-confirm; it logged nothing (`dryrun/`).

**DEV-confirm set:** 179 M1-sel pairs, 37 kinases, 22 groups. 143 of the 179 pairs are signature kinases.

| method | M1-sel | M1-all | M1-sel-scored / cov | harmony | H / D / O | T-unseen |
|---|---|---|---|---|---|---|
| **core-v3** | **0.842** | 0.699 | 0.844 / 1.00 | 0.902 | .605/.923/.898 | 0.994 (1 pair) |
| core-v2 | 0.844 | 0.700 | 0.846 / 1.00 | 0.904 | .603/.929/.901 | 0.994 |
| core-v1 | 0.803 | 0.668 | 0.805 / 1.00 | 0.872 | .544/.910/.816 | 0.994 |
| PSP z-score min3 | 0.545 | 0.524 | 0.792 / 0.76 | 0.540 | .532/.596/.532 | 0.5 |
| curated z-score min1 | 0.618 | 0.566 | 0.756 / 0.96 | 0.609 | .485/.785/.621 | 0.5 |
| KSEA curated min1 | 0.618 | 0.565 | 0.754 / 0.96 | 0.609 | .493/.786/.618 | 0.5 |
| prior-only (DEV-CV freq) | 0.856 | 0.769 | 0.856 / 1.00 | 0.728 | .876/.981/.858 | 0.259 |

The core-v1, PSP-z, curated-z and prior-only rows reproduce look 1 exactly.

Paired C8 bootstrap, Δ M1-sel [95% basic CI], one-sided p:
- core-v3 − core-v1: **+0.039 [+0.006, +0.069], p 0.016** (M1-all +0.031, harmony +0.030, M2 +0.039); core-v2 − core-v1: **+0.041 [+0.008, +0.070], p 0.012**; core-v3 − core-v2: −0.002 [−0.003, −0.000] (E1 cost, same size as on DEV).
- core-v3 − PSP-z +0.298 [+0.237, +0.352] (M2 +0.029, 136 pairs); − curated-z +0.224 [+0.185, +0.266] (M2 +0.078); − KSEA curated +0.225 [+0.181, +0.272]; all p .0005.
- core-v3 − prior-only: −0.014 [−0.071, +0.068], p 0.67 (M1-all −0.070; harmony +0.174).
- Gain over core-v1 comes from signature kinases (0.793 → 0.852; non-signature 0.823 → 0.823) and the Hijazi source (+0.06).
- Calibration |z| > 1.96, core-v1 → core-v3: DEV-confirm nulls (4 conditions, 2,012 cells) 6.5% → 5.9%; site-permuted 5.9% → 5.8%; signature kinases, permuted 6.6% → 5.2%.

**D11 condition: PASSED.** The point ΔM1-sel vs core-v1 is above 0 for both core-v2 (+0.041) and core-v3 (+0.039). **core-v3 is kept; no revert to core-v1.**

Logged: `results/devconfirm_looks.tsv` 7 rows (look = 2, bundle core-v3, hash `46caff727af577b3`); `results/dev_ledger.tsv` 3 rows in branch DEVCONFIRM (core-v3, v2, v1; look = 2) plus 10 rows in branch G2PREP (core-v3 and core-v1 per fold seed).

## Notes for G2
- Quote the DEV number as **core-v3 M1-sel 0.680 ± 0.004 (fold seeds 0–4)**, with Δ vs PSP-z +0.165, vs curated-z +0.117, and vs core-v1 +0.016.
- **Still below prior-only on M1-sel.** On DEV the mean Δ is −0.006 and on DEV-confirm it is −0.014, both n.s. The "beats prior-only" Holm test on TEST is not expected to pass.
- **T-unseen carries the case against prior-only** (+0.34 on DEV). On DEV-confirm it cannot be assessed: there is only 1 unseen pair.
- **The DEV-confirm gain over core-v1 is concentrated** in signature kinases (MAP2K1/2, MTOR and others). That is consistent with the B7b audit caveat.
- `cache/g2_prep/core_v3_alldev_fit.pkl` holds the all-DEV signature fit and the v_k scales. Before TEST, these must be refit excluding conditions that share a control with TEST; that exclusion needs metadata the DEV agents don't have.
