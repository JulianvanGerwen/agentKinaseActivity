# Independent audit of core-v1 (D6) and DEV-confirm look #1 — 2026-10-27

Auditor: independent (did not produce B2/B2.1). I did not open `cache/sealed/` or the git history of `results/ground_truth`.
GT_FREEZE sha256 `1e6f791c…` (checked by `Ctx`).

Scripts:
- `scripts/audit_core_v1_dev.py` → `results/audit_core_v1/dev_audit.json`, `dev_grid_rescored.tsv`, log `results/audit_core_v1_dev.log`.
- `scripts/audit_core_v1_devconfirm.py` → `results/audit_core_v1/devconfirm_look1.json`, `devconfirm_frozen_config.json`,
  `results/devconfirm_looks.tsv`, plus 2 rows in `results/dev_ledger.tsv` (branch `DEVCONFIRM`).

## Verdict
- **Reproduced, no code bug, no leakage found.** The adoption evidence holds: "B2 model, ρ_f ≤ 3, z_post, no winsor"
  gives nested M1-sel 0.670, Δ vs the B2 reference +0.015 [+0.006, +0.035], p = 0.001. Every number matches.
- **Inconsistency in D6 (main finding).** The 0.670 evidence uses **z_post** for ranking. z_post is *not* calibrated:
  - |z| > 1.96 occurs in 0.6% of DEV-null cells and 0.25% of site-permuted cells.
  - The calibrated statistic named in D6 (z_vif, ρ_e ≈ 0.01, no winsor) ranks differently and does worse:
    - τ, ρ_g nested: 0.663, Δ vs B2 ref +0.008 [−0.005, +0.030], p = 0.26;
    - τ fixed at 0.125: 0.666, Δ +0.011 [−0.002, +0.033], p = 0.09.
  - **Without winsorising, z_vif fails the hard adoption test.**
  - The tree's claims "calibrated z … 5.0% overall, 4.3–6.0% per bin, no M1 cost" come from the **winsor-3** z_vif row,
    which was *not* adopted. The no-winsor numbers are below.
  - The PI has to choose between:
    - (a) ranking by z_post (adopted evidence holds, but the output is not a calibrated z, which breaks DESIGN §1);
    - (b) z_vif without winsor (calibrated, but DEV evidence n.s.);
    - (c) re-opening winsor-3 z_vif (0.672; the threshold was chosen post hoc).
  - Shipping z_vif and quoting z_post's DEV numbers is not allowed.
- **The T-unseen jump is real but fragile and concentrated.**
  - Size: 0.589 → 0.643, Δ +0.054 [+0.001, +0.124], p = 0.040.
  - It comes almost entirely from fold 0, where B2 picked τ 0.25 / ρ_f 10 / ρ_g 1. Unseen M1-sel in fold 0 went 0.612 → 0.691.
  - Per kinase, the gains are SRC-family kinases plus ABL1: LCK +0.65, FYN +0.27, YES1 +0.26, LYN +0.23, SRC +0.22,
    ABL1 +0.05. That is about 1.67 of the 28 × 0.054 ≈ 1.5 total. This looks like one dasatinib-like group.
  - **Do not cite it as a general unseen-kinase improvement.**
- **The ρ_f-cap gain is a fold effect, not a bug.** Folds 1 and 2 are identical (both grids pick ρ_f = 3). The gain comes from:

  | fold | B2 ref | core-v1 |
  |---|---|---|
  | 0 | 0.586 | 0.615 |
  | 3 | 0.763 | 0.770 |
  | 4 | 0.666 | 0.677 |

  Caveat: dropping ρ_f = 10 was decided after seeing the outer-fold OOF drop. The +0.015 is therefore a comparison
  between two grids chosen with outer knowledge, and p = 0.001 overstates a pre-specified test. Mitigation: ρ_f ≤ 3 *was*
  the original pre-extension B2 grid (0.670 there too).

## 1. Reproduction (DEV folds 0–4, 372 GT conditions, 877 M1-sel pairs / 89 groups)
- **Grid.** I re-scored the core-v1 space with `bayes2`: curated_union + r2np@100 (α = 1), wsd scale, no winsor, borrow;
  τ {.06, .125, .25, .5, 1} × ρ_f {0, 1, 3} × ρ_g {0, .25, 1}; outputs z_post / z_vif / mean / z_freq (180 rows).
  It is identical to `results/b21/grid_out.tsv`: max |Δ| = 0 for pooled M1-sel, all 5 inner-fold values and ρ_e.
- **Nested selection.** I re-implemented it independently: inner_f = M1-sel on pairs whose condition ∉ fold f, argmax,
  OOF pooled. Same picks as the producer: τ 0.125, ρ_f 3, ρ_g 0 in fold 0 and 0.25 in folds 1–4.
- **Minimal scorer, written from scratch.** Dense numpy, built from the library parquet files: own W, own Manning kernel,
  own posterior / q / m_eff / VIF. Run on 50 random GT conditions at τ .125, ρ_f 3, ρ_g .25, ρ_e .01, against `bayes2`:
  - same scored set;
  - Spearman pooled 1.000000, minimum per condition 0.99997 (z_post) and 0.99998 (z_vif);
  - max |Δ| 2e-6.
- **Leakage checks.**
  - (a) For each fold f, I replaced fold f's scores with N(0, 1) noise in all 45 z_post configs and re-ran the selection.
    The pick for f was unchanged in 5/5 folds, so selection never sees the held-out fold.
  - (b) ρ_e is fitted only on the 7 conditions with `role == null`: decryptM 36926954, 0 vs 0 nM, in folds 0/2/3/4.
    None of them is a GT condition. Refitting ρ_e leave-outer-fold-out gave ρ_e 0.005–0.02 and OOF M1-sel 0.6619
    vs 0.6626, which is negligible.
  - (c) z_vif = μ/√(q(1 + ρ_e(m_eff − 1))) uses only posterior quantities and ρ_e. `empirical_null` uses only the
    condition's own scores. There is no label input; checked in the code.
  - (d) T-unseen within CV: "seen" = direct positives in the other folds only. Checked; the mask is identical for both methods.

### Core-v1 DEV metrics (nested OOF)
| variant | M1-sel | M1-all | M1-sel-scored / cov | M2 vs PSP-z (507 pairs) | harmony bin-macro (0 / 1-4 / 5-19 / ≥20) | H / D / O | sign − / + | T-unseen (272 pairs) | n sig |
|---|---|---|---|---|---|---|---|---|---|
| B2 reference | 0.655 | 0.604 | 0.655 / 1.000 | 0.667 vs 0.633 | 0.658 (.582/.664/.741/.646) | .644/.768/.674 | .630/.800 | 0.589 | 51 |
| **core-v1, z_post rank** (τ, ρ_g nested) | **0.670** | 0.612 | 0.670 / 0.994 | 0.695 vs 0.633 | 0.670 (.583/.691/.749/.656) | .646/.909/.678 | .647/.800 | 0.643 | 54 |
| **core-v1, z_vif** (τ, ρ_g nested) | **0.663** | 0.615 | 0.662 / 0.994 | 0.687 vs 0.633 | 0.663 (.581/.688/.737/.645) | .629/.917/.675 | .637/.805 | 0.647 | 49 |
| core-v1, z_vif, τ .125 fixed (ρ_g nested) | 0.666 | 0.612 | 0.665 / 0.994 | 0.693 vs 0.633 | 0.666 | .629/.915/.679 | .640/.808 | 0.648 | 52 |
| core-v1, z_vif, ρ_g = 1 fixed | 0.666 | 0.630 | 0.666 / 1.000 | 0.688 vs 0.633 | 0.667 | .626/.902/.682 | .641/.807 | 0.652 | 51 |
| core-v1, z_post, ρ_g = 1 fixed | 0.671 | 0.626 | 0.671 / 1.000 | 0.690 vs 0.633 | 0.672 | .643/.906/.679 | .646/.809 | 0.652 | 53 |
| PSP z-score min3 | 0.515 | 0.512 | 0.633 / 0.578 | – | 0.517 | .518/.560/.507 | .513/.536 | 0.512 | 24 |

Note: coverage 0.994 comes from fold 0 picking ρ_g = 0, which leaves family-less kinases unborrowed. Fixing ρ_g = 1 restores 1.000.

Paired bootstrap vs the B2 reference (B = 2000, C8 basic CI, one-sided p):

| Δ | z_post | z_vif |
|---|---|---|
| M1-sel | +0.015 [+0.006, +0.035], p 0.001 | +0.008 [−0.005, +0.030], p 0.26 |
| M1-all | +0.008 [+0.005, +0.018] | +0.010 [+0.006, +0.022] |
| harmony | +0.012 [+0.004, +0.027] | +0.005 [−0.007, +0.023] |
| M2 | +0.016 [+0.007, +0.035] | +0.009 [−0.005, +0.031] |
| T-unseen | +0.054 [+0.001, +0.124], p 0.040 | +0.058 [−0.006, +0.138], p 0.048 |

### Calibration (|z| > 1.96; fold-wise chosen config)

| statistic | data | all | bin 0 | bin 1–4 | bin 5–19 | bin ≥20 |
|---|---|---|---|---|---|---|
| z_vif, no winsor | DEV nulls (3,519 cells; bin 0 has 54) | 0.064 | 0.019 | 0.055 | 0.072 | 0.059 |
| z_vif, no winsor | site-permuted (186k cells) | 0.054 | 0.065 | 0.058 | 0.052 | 0.049 |
| z_post | DEV nulls | 0.006 | 0.000 | 0.000 | 0.007 | 0.009 |
| z_post | site-permuted | 0.002 | – | – | – | – |

z_vif without winsor is close to the target. The tree should quote these figures, not the winsor-3 ones. z_post is about 10× conservative.

## 2. DEV-confirm look #1 (of at most 5)
The configuration was frozen **before** any DEV-confirm value was loaded (`devconfirm_frozen_config.json`, hash `04b69b6513e3927d`):
- curated_union + r2np@100 (α = 1), wsd, no winsor, borrow;
- τ 0.125, ρ_f 3;
- ρ_g = 0.25: argmax of pooled all-DEV M1-sel for z_vif over {0, 0.25, 1}; pooled DEV 0.668;
- ρ_e = 0.01, fitted on the 7 DEV nulls;
- pre-declared PRIMARY = z_vif, SECONDARY = z_post (same posterior).

Data loaded: only the 97 `split == devconfirm` columns of `cache/quant_wide_f32.npy` (memmap column subset) plus the
7 DEV nulls from `cache/dev`. No TEST. Prior-only freq = all DEV-fold GT conditions.

Evaluation set: 179 M1-sel pairs, 37 kinases, 22 groups (O 16 / H 3 / D 3). MAP2K1/2 and MTOR make up 62 pairs.
**Note:** the tree says "devconfirm 119"; the harness counts 179 high/medium M1-sel pairs. This is unresolved;
the numbers here are harness counts.

| method | M1-sel | M1-all | M1-sel-scored / cov | harmony | H / D / O |
|---|---|---|---|---|---|
| **core-v1 z_vif (primary)** | **0.803** | 0.668 | 0.805 / 1.00 | 0.872 | .544/.910/.816 |
| core-v1 z_post (secondary) | 0.809 | 0.667 | 0.810 / 1.00 | 0.876 | .522/.930/.826 |
| PSP z-score min3 | 0.545 | 0.524 | 0.792 / 0.76 | 0.540 | .532/.596/.532 |
| curated_union z-score min1 | 0.618 | 0.566 | 0.756 / 0.96 | 0.609 | .485/.785/.621 |
| prior-only (DEV freq) | 0.856 | 0.769 | 0.856 / 1.00 | 0.728 | .876/.981/.858 |

Δ for core-v1 z_vif (paired group bootstrap, B = 2000, C8):
- **vs PSP z-score:** M1-sel +0.259 [+0.192, +0.312], p 0.0005; M1-all +0.145 [+0.019, +0.176]; M2 0.801 vs 0.792 (+0.008, 136 pairs).
- **vs curated z-score:** M1-sel +0.185 [+0.142, +0.225], p 0.0005; M1-all +0.103; M2 +0.041 (172 pairs).
- **vs prior-only:** M1-sel −0.052 [−0.123, +0.032], p 0.91; M1-all −0.101 [−0.217, −0.030]; harmony +0.144
  (bin 0 is a single kinase).
- z_post gives the same picture: +0.264 vs PSP-z, −0.047 vs prior-only.

**DEV-confirm verdict:** pass (Δ vs z-score > 0, by a wide margin). The positive M2 shows the gain is mostly coverage;
on the matched universe it is small vs PSP-z.

Other observations:
- **Below prior-only.** Core-v1 is still below the prior-only control. This replicates DEV.
- **Not a T-unseen test.** T-unseen has only 1 pair on DEV-confirm.
- **Harmony bins are small.** Bin 0 holds one kinase, so the per-bin values are uninformative.
- **Calibration on DEV-confirm:**
  - 4 DEV-confirm nulls (2,012 cells): 6.5% overall (bin 1–4 9.5%, ≥20 4.8%);
  - site-permuted DEV-confirm GT conditions: 5.9%;
  - z_post on the same nulls: 0.4%.

Logged:
- `results/dev_ledger.tsv`: 2 rows, branch DEVCONFIRM, split devconfirm, look = 1 (primary and secondary).
- `results/devconfirm_looks.tsv`: created, look 1, 5 method rows.
- Remaining looks: 4.

## Recommendations
1. Correct D6 and §5.3 of the tree:
   - state which statistic ranks;
   - replace the calibration figures with the no-winsor z_vif ones (6.4% nulls / 5.4% permutation);
   - drop "no M1 cost" (it is −0.008 without winsor).
2. Pre-register the ranking/output choice before B6. If z_vif is the deliverable, its DEV evidence vs B2 ref is n.s.
3. Fix ρ_g = 1 (or 0.25) instead of nesting it. Nesting it gives coverage 0.994 in fold 0 for no gain.
4. Report the T-unseen gain as SRC-family-driven (fold 0), not as a general unseen-kinase improvement.
