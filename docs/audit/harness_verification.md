# Harness verification (independent audit), 2026-10-27

**Final GT:** GT_FREEZE.json sha256 `1e6f791c6b85f15dfe51d7eb045bf392a6b25be576b10974f7bba6c5d887baac`, split seed 90. All numbers below were re-run on this version.

Scope: DESIGN §5, docs/harness.md, `src/kinact/{eval,bootstrap,controls,data}.py`, `scripts/smoke_dev.py`.

I never read `cache/sealed/` or the git history of results/ground_truth. The audit scripts live outside the repo, in the auditor's scratchpad: `indep.py`, `spot.py`/`spot2.py`, `bk.py`, `cmp.py`, `boot.py`, `leak.py`.

The GT was rebuilt twice during the audit. The conclusions held on all three GT versions.

## Verdict
**No code bug found, so `src/kinact` is unchanged and no regression test was added.**

The z-score M1-sel of about 0.52 comes from the metric design, mainly the ½ rule, plus differences in the GT definition. It is not caused by a sign, alignment, library, averaging or masking error.

There are two design-level issues for the PI. Changing either needs a dated §12 amendment.
1. The ½ rule compresses M1-sel for methods with partial coverage, and it decides the outcome of z-score vs prior-only.
2. The group bootstrap is biased for endpoints that macro-average over kinases. This explains why the CI upper bound sits at the point estimate.

## 1. z-score M1-sel vs benchmarKIN ≈0.78

### (a) Independent reimplementation
`indep.py` uses no kinact imports. It:
- reads the raw `quant_dev.npy`;
- maps peptides to sites with a pandas `groupby` mean;
- scores PSP z with mean·√n/sd (ddof 1, n ≥ 3);
- applies its own reading of the eligibility, negative, ½-rule and hierarchy rules.

| final GT | independent | kinact `smoke_dev` |
|---|---|---|
| M1-sel | 0.515348 | 0.515 |
| M1-sel-scored | 0.633203 | 0.633 |
| pair coverage | 0.578107 | 0.578 |

The results are identical. They were also identical on both earlier GT versions: 0.521237 / 0.666031 and 0.523181 / 0.677297. The z matrices match element-wise (max abs diff 0).

### (b) Spot checks: z recomputed by hand from the raw substrate values
Hand z, independent z and kinact z agree in every case.

| condition | target (s) | z | n sub | rank among scored negatives | AUROC scored | AUROC with ½ rule |
|---|---|---|---|---|---|---|
| Ochoa MEKi GSK-1120212 (trametinib) | MAP2K1 (−) | −6.82 | 7 | 3/87 | 0.977 | 0.585 |
| Torin1 vs insulin (21659604) | MTOR (−) | −8.51 | 29 | 1/46 | 1.000 | 0.546 |
| MCF7 MK2206 * | AKT1 (−) | −4.25 | 82 | 1/111 | 1.000 | 0.613 |
| HL60 MK2206 * | AKT1 (−) | +0.19 | 73 | 67/105 | 0.365 | 0.471 |
| Ochoa EGFRi PD-168393 * | EGFR (−) | −1.59 | 4 | 7/80 | 0.924 | 0.570 |
| NTERA2 AZD3759 | EGFR (−) | NaN | 2 (< 3) | unscored | — | 0.5 |
| EGF 5 min (17081983_1) | MAPK1 (+) | +2.34 | 37 | 5/42 | 0.902 | 0.535 |
| EGF vs NoSerum (32034161_1_1) * | MAPK1 (+) | +3.48 | 104 | 3/93 | 0.978 | 0.595 |

\* Checked on the previous GT version; the condition is no longer in DEV.

Signs, condition alignment, site alignment and kinase alignment are all correct.

**The last column is the key finding.** About 80% of the negatives are unscored by PSP, and under the ½ rule each counts ½. A perfectly ranked target therefore scores only about 0.5 + 0.5·(m/|N_c|) ≈ 0.55–0.62, where m is the number of scored negatives.

Across DEV the metric follows the approximation:

M1-sel ≈ 0.5 + coverage_target × frac_neg_scored × (M1-scored − 0.5) = 0.5 + 0.578 × 0.21 × 0.133 ≈ 0.516

The observed value is 0.515.

### (c) benchmarKIN-style metric on DEV Hijazi
Setup:
- **Conditions:** 72 DEV Hijazi conditions map to benchmarKIN experiments (HL60 and MCF7). benchmarKIN has no NTERA2 and no U73122 conditions. Mapping is by name, and each condition's best Spearman match among all benchmarKIN Hijazi columns is its own name for 68 of 68 (checked on the earlier GT version). Only DEV conditions were used.
- **Metric:** the reference `ksp/downstream/evaluate.py::perturbation_eval`, loaded read-only. It uses per-experiment RMS scaling, rm_bg, pooled AUROC, and 1000 draws of n_pos negatives with replacement.

| step (DEV Hijazi, final GT) | value | Δ |
|---|---|---|
| benchmarKIN data + benchmarKIN targets (102 pairs) + min 5 substrates, pooled MC-AUROC | **0.775** | (reference on the full Hernandez+Hijazi set: 0.782) |
| → our compendium values for the same conditions | 0.707 (min 3: 0.716) | −0.06 to −0.07 |
| → our GT targets (kinhib matrix: 824 pairs vs 102; 82 shared) | 0.617–0.632 | −0.08 to −0.10 |
| → per-experiment scored AUROC, kinase-macro, our masks | 0.662–0.681 | +0.04 |
| harness on the same 66 conditions: M1-sel-scored / M1-all-scored | 0.667 / 0.661 | ≈ |
| → ½ rule (M1-sel); negatives scored 21%, target coverage 0.44 | **0.519** | **−0.15** |

Notes on the table:
- **Subsampling:** pooled vs full AUROC differ by 0.001 or less, so the balanced subsampling changes nothing.
- **Pooling:** pooled vs per-experiment differs by 0.04 or less.
- **Reproduction:** benchmarKIN's own inputs on DEV Hijazi reproduce the published value (0.775).
- **Data difference:** the benchmarKIN Hijazi values are heavy-tailed (median column SD 6.9 vs our 0.92), with Pearson 0.35 and Spearman 0.72 against our values. This is a preprocessing difference, not a misalignment.

**Explanation of the gap:**

| source | effect |
|---|---|
| our input processing | ≈ −0.06 |
| broader multi-target Hijazi GT | ≈ −0.09 |
| per-experiment / kinase-macro / masks | ≈ +0.04 |
| ½ rule | ≈ −0.15 |

Over all DEV, M1-sel-scored = 0.633 is the number that is comparable to benchmarKIN-style metrics.

**Consequence: prior-only vs z-score is dominated by coverage.** Prior-only scores all 503 kinases, so it escapes the ½-rule compression. The scored-only and matched comparisons give different pictures:

| comparison | z-score | prior-only |
|---|---|---|
| M1-sel | 0.515 | 0.687 |
| M1-sel-scored | 0.633 | 0.687 |
| M2 (matched universe) | 0.633 | 0.562 |

On M2, z-score beats prior-only by +0.072. With M1-sel as the primary endpoint, coverage of U outweighs ranking quality. This follows from design rule A2, but the PI should confirm it is intended.

## 2. Bootstrap

**Code is correct.**
- **Brute force:** recomputing replicates 0, 1 and 7 by duplicating groups reproduces the harness to 1e-5.
- **Percentiles and p-value:** the percentile CI and the one-sided p are computed as specified.
- **Shared W:** a single W matrix is reused for all methods and endpoints.

**Why the CI upper bound equals the point estimate: the replicate distribution is biased (final GT).**
- 51 of the 106 M1-sel kinases occur in only one GT group. Each is absent from about 37% of resamples, while kinases with several groups are almost always present.
- The per-kinase Δ(z − prior) is −0.004 for single-group kinases and −0.328 for multi-group kinases.
- The replicates therefore over-weight the multi-group kinases:

| | point estimate | mean of replicates | replicates above point |
|---|---|---|---|
| prior-only M1-sel | 0.6875 | 0.7184 | 95% |
| Δ(z − prior) | −0.172 | −0.202 | 6.6% |

The percentile CI [−0.255, −0.166] barely contains the point estimate. A basic (pivot) interval would be about [−0.179, −0.089].

The same bias affects:
- p-values for any Δ that correlates with how many groups a kinase has;
- SE(Δ), and through it the MDE (now 0.0605).

**Recommendation (amendment, not a silent fix):** report the bootstrap bias mean(Δ*) − Δ̂ with every paired Δ, and pre-register a bias-corrected interval and p (BC percentile, or centring Δ* on Δ̂). Alternatively, hold the kinase set fixed at the point-estimate set and resample groups within each kinase.

**Random control:**
- Over 20 seeds M1-sel averages 0.495, with SD 0.018 and range 0.430–0.525. It is unbiased; the earlier 0.476 was about −1.3 SD.
- Seed 0 gives 0.505 on the final GT, and the bootstrap SE of a single run is 0.012.

## 3. Prior-only is out of fold
- No GT group spans more than one split.
- For each fold, the training conditions and the evaluated conditions are disjoint.
- Dropping all of the evaluated fold's label rows, or permuting their kinase labels, changes freq_k and the prior-only scores by exactly 0.
- The only use of a condition's own labels is the oracle sign s_c, as designed.
- For contrast, an in-fold (leaky) frequency would give 0.861 vs 0.687 out of fold.

## Tests
`pytest tests/`: 94 passed on the final GT, with no code changes.

## Final smoke (`scripts/smoke_dev.py`, re-run last on GT_FREEZE sha256 1e6f791c…)
372 eligible conditions, 89 M1-sel groups.

| method | M1-sel | M1-all | M1-sel-scored | pair coverage | harmony |
|---|---|---|---|---|---|
| zscore | 0.515 | 0.512 | 0.633 | 0.578 | 0.517 |
| random | 0.505 | 0.501 | 0.505 | 1.000 | 0.510 |
| prior_only | 0.687 | 0.543 | 0.687 | 1.000 | 0.726 |
| shift_prior | 0.546 | 0.514 | 0.546 | 1.000 | 0.561 |
| oracle | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |

Paired bootstrap (B = 2000):

| comparison | Δ | 95% CI | p |
|---|---|---|---|
| z vs prior-only | −0.172 | [−0.255, −0.166] (upper bound biased, see §2) | — |
| z vs random | +0.010 | [−0.015, +0.039] | 0.17 |

MDE = 0.0605.

## C8 implemented (after PI decision)
`src/kinact/bootstrap.py` now reports:
- the basic CI [2Δ̂ − q97.5, 2Δ̂ − q2.5];
- the shifted-null one-sided p;
- the raw percentile CI and p, the mean of the replicates and the bias, as diagnostics.

New tests in `tests/test_eval.py`:
- `test_c8_single_group_kinase_bias_basic_ci`: a toy case that reproduces the bias;
- `test_c8_null_pvalues_uniform`: 200 equal-method nulls, where the fraction with p < 0.05 must fall in [0.01, 0.10].

`pytest tests/`: 96 passed.

Smoke results on the same GT (`1e6f791c…`):

| comparison | Δ | basic CI | shifted-null p | bias |
|---|---|---|---|---|
| z vs prior-only | −0.172 | [−0.179, −0.089] | 1.0 | −0.030 |
| z vs random | +0.010 | [−0.018, +0.035] | 0.217 | — |

The SE is unchanged at 0.0239 because the correction is a pure location shift. The recomputed MDE is therefore also unchanged at 0.0605 (`results/smoke/mde.json`).
