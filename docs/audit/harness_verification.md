# Harness verification (independent audit), 2026-10-27

Scope: DESIGN §5, docs/harness.md, `src/kinact/{eval,bootstrap,controls,data}.py`, `scripts/smoke_dev.py`.
I never read `cache/sealed/` or the git history of results/ground_truth. The audit scripts live outside the repo, in the auditor's scratchpad: `indep.py`, `spot.py`, `bk.py`, `cmp.py`, `boot.py`, `leak.py`.
The GT and DEV cache were rebuilt and re-frozen at 13:01–13:03 during this audit. All numbers below come from the frozen GT. The final smoke run's GT_FREEZE hash is given in the section at the end.

## Verdict
**No code bug found. No change to `src/kinact`.** The z-score DEV M1-sel of about 0.52 comes from the metric design, mainly the ½ rule, plus differences in the ground truth. It is not caused by a sign, alignment, library, averaging or masking error.

There are two design-level issues for the PI. Changing either needs a dated §12 amendment.
1. The ½ rule compresses M1-sel for methods with partial coverage.
2. The group bootstrap is biased for endpoints that macro-average over kinases.

## 1. z-score M1-sel 0.52 vs benchmarKIN ≈0.78

### (a) Independent reimplementation
`indep.py` uses no kinact imports. It:
- reads the raw `quant_dev.npy`;
- maps peptides to sites with a pandas `groupby` mean;
- scores PSP z with mean·√n/sd (ddof 1, n ≥ 3);
- applies its own reading of the eligibility, negative, ½-rule and hierarchy rules.

| | independent | kinact `smoke_dev` |
|---|---|---|
| M1-sel | 0.523181 | 0.523181 |
| M1-sel-scored | 0.677297 | 0.677297 |
| pair coverage | 0.651252 | 0.651252 |

The results are identical. The z-score matrix also matches element-wise (max abs diff 0) on the spot-checked rows. On the pre-freeze GT both implementations also agreed, at 0.521237 / 0.666031.

### (b) Spot checks: z recomputed by hand from the raw substrate values

Hand z, independent z and kinact z agree in every case.

| condition | target (s) | z | n sub | rank among scored negatives | AUROC scored | AUROC with ½ rule |
|---|---|---|---|---|---|---|
| Ochoa MEKi GSK-1120212 (trametinib) | MAP2K1 (−) | −6.82 | 7 | 3/87 | 0.977 | 0.585 |
| Selumetinib 27422710_1 | MAP2K1 (−) | −1.17 | 6 | 6/128 | 0.961 | 0.621 |
| MCF7 MK2206 | AKT1 (−) | −4.25 | 82 | 1/111 | 1.000 | 0.613 |
| HL60 MK2206 | AKT1 (−) | +0.19 | 73 | 67/105 | 0.365 | 0.471 |
| Ochoa EGFRi PD-168393 | EGFR (−) | −1.59 | 4 | 7/80 | 0.924 | 0.570 |
| Ochoa PI-103 | MTOR (−) | −2.31 | 33 | 7/75 | 0.919 | 0.568 |
| EGF vs NoSerum 32034161_1_1 | MAPK1 (+) | +3.48 | 104 | 3/93 | 0.978 | 0.595 |
| EGF vs NoSerum 32034161_1_1 | EGFR (+) | +1.46 | 7 | 17/93 | 0.826 | 0.565 |

Signs, condition alignment, site alignment and kinase alignment are all correct.

**The last column is the key finding.** About 80% of the negatives (≈375 of ≈480) are unscored by PSP. Under the ½ rule they all count ½. A perfectly ranked target therefore scores only about 0.5 + 0.5·(m/|N_c|) ≈ 0.60, where m is the number of scored negatives.

Across DEV the metric follows the approximation:

M1-sel ≈ 0.5 + coverage_target × frac_neg_scored × (M1-scored − 0.5) = 0.5 + 0.651 × 0.22 × 0.177 ≈ 0.525

The observed value is 0.523.

### (c) benchmarKIN-style metric on DEV Hijazi
Setup:
- **Conditions:** 78 DEV Hijazi conditions map to benchmarKIN experiments (HL60 and MCF7). benchmarKIN has no NTERA2 and no U73122 conditions. Mapping is by name, and each condition's best Spearman match among all benchmarKIN Hijazi columns is its own name for 68 of 68. Only DEV conditions were used.
- **Metric:** the reference `ksp/downstream/evaluate.py::perturbation_eval`, loaded read-only. It uses per-experiment RMS scaling, rm_bg, pooled AUROC, and 1000 draws of n_pos negatives with replacement.

| step (DEV Hijazi) | value | Δ |
|---|---|---|
| benchmarKIN data + benchmarKIN targets (124 pairs) + min 5 substrates, pooled MC-AUROC | **0.744** | (reference on the full Hernandez+Hijazi set: 0.782) |
| → our compendium values for the same conditions | 0.703 (min 3: 0.716) | −0.03 to −0.04 |
| → our GT targets (kinhib matrix: 1,076 pairs vs 124; 94 shared) | 0.589 | −0.12 |
| → per-experiment scored AUROC, kinase-macro, our masks | 0.625 | +0.04 |
| harness on the same 72 conditions: M1-all-scored / M1-sel-scored | 0.624 / 0.654 | ≈ |
| → ½ rule (M1-sel); negatives scored 22%, target coverage 0.39–0.65 | **0.517** | **−0.14** |

Notes on the table:
- **Subsampling:** pooled vs full AUROC differ by less than 0.001, so the balanced subsampling changes nothing.
- **Pooling:** pooled vs per-experiment differs by 0.03 or less.
- **Data difference:** the benchmarKIN Hijazi values are heavy-tailed (median column SD 6.9 vs our 0.92), with Pearson 0.35 and Spearman 0.72 against our values. This is a preprocessing difference, not a misalignment.

**Explanation of the gap:**

| source | effect |
|---|---|
| condition subset (Hijazi DEV vs full benchmark) | ≈ −0.04 |
| input processing | ≈ −0.04 |
| broader multi-target Hijazi GT | ≈ −0.12 |
| hierarchy / pooling | ≈ +0.03 |
| ½ rule | ≈ −0.14 |

Over all DEV, M1-sel-scored = 0.677 is the number that is comparable to benchmarKIN-style metrics.

**Consequence: prior-only vs z-score is dominated by coverage.** Prior-only scores all 503 kinases, so it escapes the ½-rule compression. The scored-only and matched comparisons give different pictures:

| comparison | z-score | prior-only |
|---|---|---|
| M1-sel | 0.523 | 0.707 |
| M1-sel-scored | 0.677 | 0.707 |
| M2 (matched universe) | 0.677 | 0.562 |

On M2, z-score beats prior-only by +0.115. M1-sel as the primary endpoint rewards coverage of U more than ranking quality. This is by design (A2), but the PI should confirm it is intended.

## 2. Bootstrap

**Code is correct.**
- **Brute force:** recomputing replicates 0, 1 and 7 by duplicating groups reproduces the harness exactly.
- **Percentiles and p-value:** the percentile CI and the one-sided p are computed as specified.
- **Shared W:** a single W matrix is reused for all methods and endpoints.

**Why the CI upper bound equals the point estimate: the replicate distribution is biased.**
- 65 of the 108 M1-sel kinases occur in only one GT group. Each is absent from about 37% of resamples, while kinases with several groups are almost always present.
- The per-kinase Δ(z − prior) is −0.084 for single-group kinases and −0.335 for multi-group kinases.
- The replicates therefore over-weight the multi-group kinases:

| | point estimate | mean of replicates | replicates above point |
|---|---|---|---|
| prior-only M1-sel | 0.707 | 0.733 | 97.5% |
| Δ(z − prior) | −0.184 | −0.208 | 4% |

The percentile CI [−0.241, −0.181] therefore barely contains the point estimate. A basic (pivot) interval would be about [−0.187, −0.127]. p-values for any Δ that correlates with how many groups a kinase has are biased in the same way.

**Recommendation (amendment, not a silent fix):** report the bootstrap bias mean(Δ*) − Δ̂ with every paired Δ, and pre-register a bias-corrected interval and p (BC percentile, or centring Δ* on Δ̂). Alternatively, hold the kinase set fixed at the point-estimate set and draw groups within each kinase.

**Random control:**
- Over 20 seeds M1-sel averages 0.496, with SD 0.016 and range 0.466–0.525. It is unbiased; 0.476 was a −1.5 SD draw.
- The bootstrap SE of a single random run is 0.012.
- On the frozen GT, seed 0 gives 0.500.

## 3. Prior-only is out of fold
- No GT group spans more than one split.
- For each fold, the training conditions and the evaluated conditions are disjoint.
- Dropping all of the evaluated fold's label rows, or permuting their kinase labels, changes freq_k and the prior-only scores by exactly 0.
- The only use of a condition's own labels is the oracle sign s_c, as designed.
- For contrast, an in-fold (leaky) frequency would give 0.837 vs 0.707 out of fold.

## Tests
`pytest tests/`: 64 passed, with no code changes.

## Final smoke (re-run last)
See the section appended below.
