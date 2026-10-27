# B7 — perturbation-derived kinase response signatures (DEV folds 0–4, nested OOF; DESIGN §12 C10)
GT_FREEZE `1e6f791c6b85f15d…` (checked by `Ctx`). Data read: `cache/dev` + `results/ground_truth` DEV files only.
- Code: `src/kinact/explore/signatures.py` (design pre-declared in its docstring before any result); `scripts/exp_b7_{common,grid,eval}.py`.
- Tests: `tests/test_signatures.py` (7, pass).
- Raw: `results/b7/{b7_eval.json, grid_<variant>.tsv, fits_<variant>.json, grid.log, eval.log}`; score/fit archives in `cache/b7/`.
- Ledger: 196 rows, branch `B7` (4 variants × 48 configs per-config OOF + 4 nested rows).
- Wall time ≈ 35 min (grid 9 min, eval 25 min).
## Method (as pre-declared)
**Signature fit.** For each training set:
- Data: wsd-standardised log2FC, clipped at |y| ≤ 4.
- Positives: direct, high/medium.
- Sign-aligned s·y, averaged within each drug group, then across groups.
- Generic response subtracted ("contrast"): the mean over all training groups of the condition-sign-aligned response.
- Standardised by the site's between-group SD.
- Shrinkage: w = (m/σ)·n/(n+κ).
- Kinase eligible if it has ≥ 2 training drug groups; edges are its top `n_top` sites by |w| (signed weights).
**Scoring.** W = curated + r2np@100 + α_sig·W_sig, in the core-v1 posterior and z_vif. ρ_e is refitted on the DEV nulls for each design.
- Core HPs (τ, ρ_f, ρ_g) are the core-v1 nested picks for each outer fold.
- α_sig = 0 reproduces core-v1: M1-sel 0.6627 (audit 0.6626), T-unseen 0.648; max |Δz| 1e-6 in fold 0.
**HPs (3, nested):** α_sig {0, .25, .5, 1, 2} × n_top {20, 50, 100, 200} × κ {0, 1, 4}.
**Nesting.**
- Inner score for outer fold f / eval fold g uses signatures fitted on the folds ∉ {f, g}.
- The outer score for fold f uses signatures fitted on the 4 other folds.
**Variants.**
- `clean`: positives only from conditions with n_targets ≤ 10. This is the pre-declared PRIMARY, run with contrast on.
- `all`: every GT condition.
- Each of `clean` and `all` is run with and without contrast.
**C10 safeguards.**
- (i) Leakage guards raise on: an evaluated condition in the training set; a training condition sharing a drug group with an evaluated condition; a positives row from a non-training condition. Drug groups are disjoint across folds (0 shared groups, 0 shared drug names).
- (ii) No kinase offset: signatures are design columns only. Fit and score are sign-symmetric (tested).
## Kinases with signatures per fit
| fit | clean+contrast | all+contrast |
|---|---|---|
| outer (4 training folds), folds 0–4 | 48 / 40 / 40 / 43 / 45 | 87 / 109 / 90 / 121 / 116 |
| inner (3 training folds) | 28–40 | 50–107 |
## Results (nested OOF; Δ = paired group bootstrap vs core-v1, B = 2000, C8 basic CI, one-sided p)
| variant | M1-sel | M1-all | harmony | T-unseen | T-seen | has-sig. kin. | no-sig. kin. |
|---|---|---|---|---|---|---|---|
| core-v1 (α = 0) | 0.663 | 0.615 | 0.663 | 0.648 | 0.668 | 0.733 | 0.646 |
| **clean+contrast (PRIMARY)** | **0.675** | 0.613 | 0.680 | **0.622** | 0.694 | 0.809 | 0.646 |
| all+contrast | 0.669 | 0.606 | 0.672 | 0.644 | 0.678 | 0.726 | 0.629 |
| clean, no contrast | 0.670 | 0.603 | 0.674 | 0.612 | 0.691 | 0.819 | 0.639 |
| all, no contrast | 0.660 | 0.598 | 0.667 | 0.637 | 0.668 | 0.720 | 0.618 |
- The has-sig. / no-sig. columns refer to each variant's own kinase set.
- The core-v1 entries for T-seen and the two kinase-set columns are back-computed (B7 value − Δ), so they are approximate.
PRIMARY Δ vs core-v1:

| endpoint | Δ [95% CI] | p |
|---|---|---|
| **M1-sel** | **+0.012 [−0.017, +0.026]** | **0.15 (fails the hard test)** |
| M1-all | −0.002 [−0.020, +0.006] | 0.59 |
| harmony | +0.017 [−0.013, +0.033] | 0.07 |
| M2 (matched, 872 pairs) | +0.011 [−0.020, +0.026] | 0.18 |
| **T-unseen** (272 pairs) | **−0.026 [−0.064, +0.007]** | **fails the −0.01 margin** |
| T-seen (605 pairs) | +0.026 [−0.006, +0.044] | 0.024 |
| kinases with a training signature (366 pairs) | +0.076 [+0.033, +0.133] | 0.003 |
| kinases without one (511 pairs) | −0.000 [−0.029, +0.016] | 0.50 |
| H | +0.036 [+0.003, +0.060] | 0.015 |
| D | −0.007 | – |
| O | +0.013 | – |
| sign −1 | +0.021 | 0.05 |
| sign +1 | −0.020 [−0.080, +0.007] | – |
- Harmony bins (point non-inferiority holds in all four): 0 +0.003, 1–4 +0.014, 5–19 +0.018, ≥20 +0.033.
- M2 vs PSP-z (507 pairs): 0.728 vs 0.633 (core-v1: 0.687).
- The other variants are weaker: all+contrast +0.006 (p 0.32); clean, no contrast +0.007 (p 0.27; M1-all −0.012 and T-unseen −0.036, both significant); all, no contrast −0.003.
- **Contrast helps** (+0.005 / +0.009 in each pair). **Clean > all** on seen kinases, but clean loses more on T-unseen.
- **Selection noise is large.** Pooled-best fixed configs reach 0.676 (clean+contrast) and 0.684 / 0.686 (all variants), but nested picks scatter: α 0.25–2 and n_top 20–200 across folds.
## Controls (PRIMARY, fold-wise chosen HPs)
- **Shuffled signatures** (site labels permuted within each signature; weights and sizes kept): M1-sel 0.651 / 0.647 / 0.647 (3 seeds), vs core −0.011 [−0.026, +0.002]. T-unseen 0.637.
  - The gain does not just collapse to core-v1; it turns into a loss, because random signed edges add noise.
  - Real − shuffled ≈ +0.025, so the gain comes from site identity.
- **Generic-response signature** (every signature kinase gets the same generic-response top sites, same norm): M1-sel 0.658, vs core −0.005 (n.s.). There is no popularity leak through a generic response.
- **Calibration** (|z| > 1.96; fold-f nulls scored with the fold-f design):

  | | DEV nulls (3,519 cells), all (0 / 1–4 / 5–19 / ≥20) | site-permuted GT, all (0 / 1–4 / 5–19 / ≥20) |
  |---|---|---|
  | PRIMARY | 5.9% (3.8 / 5.8 / 6.2 / 5.5) | 5.2% (6.1 / 5.7 / 5.0 / 4.8) |
  | core-v1 | 6.4% | 5.4% |
  Calibration is unchanged or slightly better.
## Why T-unseen drops
The signature kinases do not outrank unseen targets: the share of signature-kinase negatives beating an unseen target goes 0.216 → 0.209. Instead, unseen targets' own scores move.
- Indirect effect: signature columns explain shared sites (explaining-away), and the family kernel couples the posteriors.
- Of 28 unseen kinases, 12 lose and 5 gain more than 0.01.
- **TYK2 alone: −0.53**, about −0.019 of the −0.026. Its JAK-family relatives JAK2 and JAK3 have signatures.
- Next: MAP4K2 −0.15, MAP4K5 −0.12, PRKCZ −0.09. Gains: MAPK11 +0.21, CDK19 +0.10.
- So the C10 requirement that kinases without a signature stay "unchanged" holds on average (−0.000) but is violated for single unseen kinases that sit in a signature kinase's family.
## Signature sanity (clean+contrast, fold-0 outer fit; top 8 sites, w; PSP = annotated kinases)
- **MTOR (7 groups):** mTORC1 readouts — EIF3B S81, **LARP1 S766 (PSP: MTOR)**, SNX2 S117, KMT2B S844, ESYT2 S743, **RRAGC S15**, **EIF4EBP2 T46**. In the `all` variant: ULK1 S450, EIF4EBP2 S65, CAD S1859 (S6K site).
- **MAPK1 (5):** GIGYF2 S593, SAMD4B S271, ZNF185 S465, **STAT3 S727** (ERK readout), PLIN3 S130 (−), EPN1, AMFR, TJP1. Proline-directed, plausible.
- **MAP2K1 (5):** BAZ2A S135, SPAG1 S347 (−), PDLIM2, RIF1; `all` adds **BRAF S151 (ERK feedback)** and CIC S173.
  - Jaccard(MAPK1, MAP2K1) top-50 = 0.01. In this GT, MEK inhibitors make MAP2K1/2 the direct targets; ERK is downstream and masked, so the overlap is small.
  - Paralogs that are always co-targeted get identical signatures: MAPK1 = MAPK3 (Jaccard 1.0), MAP2K1 ≈ MAP2K2 (0.96).
- **EGFR (10):** weak, generic-looking (NCAPD2, ACACA, CTNND1). In the `all` variant: **SHC1 Y427, EGFR Y1172 (autophosphorylation)**, DYRK1A Y321.
- **PRKACA (3):** VIM S73 and CAD S1406 (PSP PRKACA), ITPR3 S1832, ARPP19 S104 (PSP PRKACA, `all` variant).
- **Enrichment of known substrates:**
  - 14% of top-50 signature sites are PSP substrates of any kinase, vs 7.6% of eligible sites (about 2×).
  - Only 32 of 2,400 are PSP substrates of the kinase itself: signatures are mainly pathway readouts, not direct substrates.
- Mean pairwise top-50 Jaccard across kinases is 0.012, so signatures are distinct.
## Verdict: not adopted ([-] at this time box; C10 iii fails)
- Hard test fails: +0.012, p = 0.15.
- T-unseen non-inferiority fails: −0.026 against a −0.01 margin, driven by TYK2 / MAP4K.
- Shuffled control: real − shuffled ≈ +0.025, which supports signal from site identity, but the shuffled score goes below core instead of collapsing to it.
- What does hold:
  - a strong, specific gain for kinases with a training signature (+0.076, p 0.003; H +0.036);
  - calibration unchanged;
  - no generic-response or popularity leak;
  - all four harmony bins point non-inferior.
- Best use: as a stacking feature (B6) or with a decoupling fix, i.e. signatures as a separate additive score not entering the joint posterior (that should leave unseen kinases exactly unchanged). That fix was not tried here: it was outside the pre-declared design and time box.
## Caveats
- ρ_e is refitted on 7 nulls for each design. In the inner grid the nulls' groups may overlap a training fold; this affects only the global scalar ρ_e. The calibration table uses fold-clean scoring.
- Only 2–11 training groups per kinase. Paralogs that are always co-targeted cannot be separated.
- Core HPs were taken from the core-v1 nested picks rather than re-nested jointly. This is equivalent to the reference and adds no outer-fold information.
- T-unseen and has-signature sets are within-CV definitions: positives in other folds, and kinases in the outer fit, respectively.
