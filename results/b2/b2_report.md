# Branch 2 — Bayesian linear-Gaussian core (DEV, folds 0-4 pooled OOF)

GT: `GT_FREEZE.json` sha256 `1e6f791c6b85f15dfe51d7eb045bf392a6b25be576b10974f7bba6c5d887baac`; data `cache/dev/quant_dev.npy`
(`d3babe0d…`). Only DEV cache + `results/ground_truth/` DEV files were used. Code: `src/kinact/explore/bayes.py`,
`scripts/exp_b2_grid.py` (grid), `scripts/exp_b2_eval.py` (nested CV, ablations, calibration, runtime), tests
`tests/test_bayes.py` (6 pass: planted activities recovered, closed form = numpy = CPU path, family borrowing,
sign-flip/shift invariance, toy null calibration, noise model). Raw numbers: `results/b2/b2_eval.json`, grids
`results/b2/grid_noise-both{,_lowtau,_rhof10}.tsv`. Ledger: 6,464 grid rows + nested/reference rows, branch `B2`.

## Model (single-sample, closed form)
y_s = (x_s − med_c)/σ_c = Σ_k W_sk a_k + ε_s, ε_s ~ N(0, v_s); a ~ N(0, τ²K); posterior mean μ = Λ⁻¹W'Dy,
Λ = W'DW + K⁻¹/τ², SD = √diag Λ⁻¹, output z = μ/SD. Only measured sites enter (per-condition mask).
- σ_c = 1.4826·MAD of (x−med)/√v over all measured sites of c (label-free, within sample).
- v_s (noise model): robust per-site variance of standardised values over allowed DEV conditions (qc ok,
  **excluding the scored condition's own fold**), median(r²)/0.455 shrunk to 1 with 10 pseudo-obs; off = v ≡ 1.
- K = (I + ρ_f·[same Manning family] + ρ_g·[same group, not Other/Atypical]) normalised to unit diagonal.
- Library: W = curated_union (weight 1) + α·r2np@100 (top-100 predicted sites per kinase, binary; edges not
  already curated). α = 0 → curated only. Kinase scored iff ≥1 measured library substrate; `borrow` also reports
  kinases with no measured substrate but a K-relative with one (scored purely from relatives).
- Hyperparameters (5): τ, ρ_f, ρ_g, α, intercept (median-centring vs joint flat-prior intercept); plus the two
  design switches asked for (noise model on/off, borrow) and output (z vs posterior mean) — all chosen by nested CV
  inside each variant. No GT-derived per-kinase information; no per-kinase/per-edge parameters.

## Nested CV protocol
Each condition is scored with v from the other folds, so a config's score of c does not depend on the outer split.
For each outer fold f the config with best pooled M1-sel on the other 4 folds is chosen; OOF scores are pooled and
evaluated once (frozen harness, n_perm 2000). Grid: τ {0.03,0.06,0.125,0.25,0.5,1,2,4} × ρ_f {0,0.5,1,3,10}
× ρ_g {0,0.25,1} × α {0,0.25,0.5,1} × intercept × noise × borrow × output = 5,888 evaluations (τ<0.25 and
ρ_f=10 added after the first grid hit those edges; nested result on the first grid alone: M1-sel 0.670).

## Results (nested OOF, DEV folds 0-4; 372 GT conditions, 89 M1-sel groups)
| variant | M1-sel | M1-all | M1-sel-scored | cov | harmony | H | D | O | sign− | sign+ | T-unseen | n sig | same-fam |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **B2 (z; borrow, noise, all HPs free)** | **0.655** | 0.604 | 0.655 | 1.000 | 0.658 | 0.644 | 0.768 | 0.674 | 0.630 | 0.800 | 0.589 | 51 | 0.588 |
| B2 no-borrow | 0.648 | 0.594 | 0.655 | 0.983 | 0.650 | 0.635 | 0.767 | 0.662 | 0.625 | 0.784 | 0.584 | 51 | 0.590 |
| abl: posterior mean instead of z | 0.672 | 0.617 | 0.672 | 1.000 | 0.671 | 0.644 | 0.911 | 0.679 | 0.649 | 0.802 | 0.654 | 53 | 0.581 |
| abl: no K (ρ_f=ρ_g=0) | 0.617 | 0.581 | 0.624 | 0.983 | 0.626 | 0.595 | 0.854 | 0.617 | 0.593 | 0.763 | 0.593 | 39 | 0.590 |
| abl: noise model forced on | 0.646 | 0.608 | 0.646 | 1.000 | 0.645 | 0.618 | 0.897 | 0.653 | 0.623 | 0.758 | 0.640 | 44 | 0.584 |
| abl: curated only (α=0) | 0.607 | 0.610 | 0.615 | 1.000 | 0.612 | 0.575 | 0.903 | 0.603 | 0.577 | 0.742 | 0.597 | 41 | 0.562 |
| abl: no K, α=0 | 0.553 | 0.536 | 0.623 | 0.860 | 0.557 | 0.552 | 0.697 | 0.538 | 0.543 | 0.617 | 0.565 | 37 | 0.588 |
| abl: no K, α=0, τ≥1 (≈ ridge mlm) | 0.550 | 0.523 | 0.615 | 0.860 | 0.558 | 0.539 | 0.634 | 0.551 | 0.538 | 0.624 | 0.546 | 37 | 0.563 |
| z-score psp_all min3 (smoke ref) | 0.515 | 0.512 | 0.633 | 0.578 | 0.517 | 0.518 | 0.560 | 0.507 | 0.513 | 0.536 | 0.512 | 24 | 0.577 |
| z-score curated min1 (best baseline) | 0.563 | 0.543 | 0.650 | 0.860 | 0.565 | 0.557 | 0.693 | 0.552 | 0.554 | 0.627 | 0.575 | 40 | 0.634 |
| prior-only (OOF freq, oracle sign) | 0.687 | 0.543 | 0.687 | 1.000 | 0.726 | 0.753 | 0.536 | 0.709 | 0.695 | 0.769 | 0.324 | 53 | 0.635 |

T-unseen = pooled OOF M1-sel over pairs whose kinase is not a DEV positive in the other folds (28 kinases, 272 pairs).
same-fam = scored-only AUROC of the target vs evaluable negatives of its own Manning family (r1_algo diagnostic).
Harmony per bin (0 | 1-4 | 5-19 | ≥20): B2 0.582 | 0.664 | 0.741 | 0.646; z-score curated 0.521 | 0.573 | 0.594 | 0.572;
prior-only 0.625 | 0.705 | 0.735 | 0.840.

Paired group bootstrap (B=2000, seed 0, C8 basic CI, one-sided shifted-null p), ΔM1-sel of B2:
| comparison | Δ | 95% CI | p |
|---|---|---|---|
| vs z-score psp_all min3 | +0.140 | [+0.078, +0.180] | 0.0005 |
| vs z-score curated min1 | +0.092 | [+0.037, +0.123] | 0.0005 |
| vs prior-only | −0.033 | [−0.068, +0.057] | 0.84 |
| vs no K | +0.038 | [−0.008, +0.069] | 0.032 |
| vs curated only (α=0) | +0.048 | [+0.004, +0.079] | 0.012 |
| vs noise model forced on | +0.008 | [−0.018, +0.023] | 0.25 |
| vs no-borrow | +0.007 | [−0.000, +0.011] | 0.012 |
| vs posterior-mean output | −0.017 | [−0.038, −0.008] | 1.0 |
M2 (matched universe, M1-sel-scored): vs z-score psp_all 0.667 vs 0.633 (Δ +0.033, 507 pairs); vs z-score curated
min1 0.674 vs 0.650 (Δ +0.023, 754 pairs) — i.e. the gain is not only coverage.

## Selected hyperparameters
Every fold chose: α=1, intercept=median, noise model OFF, borrow ON, τ=0.125 (0.25 in fold0), ρ_f ∈ {3, 10},
ρ_g ∈ {0.25, 1}. Pooled-best (optimistic, non-nested) config: τ=0.125, ρ_f=3, ρ_g=0.25, α=1, median, borrow,
no noise → M1-sel 0.673. **Instability:** folds that picked ρ_f=10 lose most of decryptM (D 0.91 → 0.77) at an
almost identical inner score; the nested z estimate dropped 0.670 → 0.655 when ρ_f=10 was added to the grid. The
posterior-mean output is more stable (0.672, D 0.911); grid maxima are equal (mean 0.674 vs z 0.673).
Recommendation: fixed config above (ρ_f capped at 3), output posterior mean for ranking, z for calibrated reporting.

## Findings
1. Strong shrinkage wins: optimum τ≈0.125 (σ_c units), interior of the grid. With small τ the posterior mean is ≈
   τ²·K·W'Dy, a family-smoothed **sum** of substrate signals (little deconvolution); τ≥1 (≈ ridge mlm) is 0.55.
2. Two independent gains, each ≈ +0.05: family/group kernel K (no-K 0.617) and predicted edges (α=0: 0.607);
   together vs neither: +0.10 (0.655 vs 0.553). K helps H/O (Hijazi multi-target drugs); D is not helped (no-K D 0.85).
3. The label-free site-noise model **hurts** (−0.008 nested, −0.017 at the pooled optimum, 0.673 vs 0.656): high-
   variance sites are the responsive ones, so down-weighting them removes signal. It also worsens calibration.
4. Borrowing (scoring edge-less kinases from relatives) is a small positive (+0.007, pair coverage 0.983 → 1.0;
   bin-0 harmony 0.574 → 0.582). The price of K is within-family discrimination: same-family AUROC 0.588 vs 0.634
   for z-score curated (r1_algo warning confirmed; prior-only also 0.635).
5. vs prior-only: B2 is below on H (0.64 vs 0.75) and ≥20 bin, far above on D (0.77–0.91 vs 0.54) and T-unseen
   (0.59 vs 0.32) — complementary signals (stacking candidate, §7.6), not yet a "beat" (Δ −0.03, n.s.).

## Calibration (|z|>1.96 rate; target 0.05)
| scorer | DEV nulls all | 1-4 | 5-19 | ≥20 | site-perm all | 0 | 1-4 | 5-19 | ≥20 |
|---|---|---|---|---|---|---|---|---|---|
| B2 pooled-best (τ .125, ρ_f 3, no noise) | 0.025 | 0.000 | 0.023 | 0.045 | 0.046 | 0.012 | 0.028 | 0.047 | 0.076 |
| B2 fold0 choice (τ .25, ρ_f 10, ρ_g 1) | 0.058 | 0.027 | 0.069 | 0.063 | 0.088 | 0.039 | 0.070 | 0.095 | 0.107 |
| B2 noise-on fold0 choice | 0.048 | 0.009 | 0.055 | 0.063 | 0.100 | 0.042 | 0.068 | 0.106 | 0.144 |
| z-score psp_all min3 | 0.067 | 0.093 | 0.046 | 0.074 | 0.100 | – | 0.072 | 0.093 | 0.177 |
| z-score curated min1 | 0.054 | 0.052 | 0.051 | 0.073 | 0.082 | – | 0.064 | 0.090 | 0.171 |
DEV nulls = 7 decryptM 0-vs-0 nM conditions (one study, 36926954; small). Site-perm = within-condition permutation of
the 372 GT conditions. Bins = measured library substrates (incl. predicted; borrowed kinases = bin 0). B2's z is
better calibrated than z-score overall and much less √m-inflated at ≥20 (0.076 vs 0.17), but still increases with
substrate count (small τ makes z sum-like) — an effective-m / empirical-null correction (branch 3) is still needed.

## v leakage check (DESIGN §4)
log v estimated without fold f vs from all DEV: r = 0.92/0.75/0.89/0.92/0.93 (folds 0-4). At the best noise-on
config, M1-sel with leave-own-fold-out v 0.657 vs all-DEV v 0.659; Spearman of scores 0.986 → negligible (and the
selected configs do not use v at all).

## Runtime (200 random DEV conditions, warm, median of 3; p = 503 kinases, 27,601 library sites)
CPU single-sample path (`SingleSample.score_one`; sparse W'DW + 503×503 Cholesky, incl. median/MAD over all
measured sites): **8.9 ms/cond at 1 thread, 4.5 ms at 16 threads** (≈17 s / 9 s for 1,937 conditions).
GPU batched (incl. CPU standardisation): 4.7 ms/cond; the posterior solve alone is 0.3 ms/cond on one RTX 4090
(grid of 5,888 evaluations in 15 min). One-off noise-model fit 5.2 s. numpy `inv` of K with 64 BLAS threads took
2.4 s (oversubscription) — use torch / limit threads.

## Caveats
- DEV only; the grid is large (5,888 evals) — nested CV protects the estimate but folds disagree on ρ_f (above).
- r2np@100 percentile weights are ≈1 for all top-100 edges (weight=percentile is uninformative) → binary edges.
- α at its upper grid edge (1.0, i.e. predicted edges weighted like curated); α>1 not tried.
- Joint intercept never selected (median-centring equivalent within 0.003).
- T-unseen uses a pooled-OOF definition (not in the frozen harness); smoke/baselines did not report it.
