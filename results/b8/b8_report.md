# B8 — final development box: E1 calibration, E2 shared-factor removal, E3 hyperparameter bagging (DEV folds 0–4)

## PRE-DECLARED DESIGN (written 2026-10-27 before any B8 number was computed; not changed afterwards)
Source of ideas: docs/retro/ideas_m2.md (ideas 1 and 3), docs/audit/b7b_audit.md §4/§6.4 (calibration fix). DEV only
(cache/dev + results/ground_truth DEV files). Code: `src/kinact/explore/final_box.py`, `scripts/exp_b8_*.py`.
Ledger branch `B8`. Time box ≤ 3 h; after this box development stops.

**Reference (core-v2, D11)** = B7b nested exactly as in results/b7b: core-v1 z_vif with per-fold core-v1 nested core picks,
late-fused with `clean+contrast` signatures (B7 per-fold (n_top, κ)), z_sig/ρ_sig (RMS over DEV nulls),
z = (z_core + α z_sig)/√(1+α²), per-fold α*_f = B7b nested picks (fold0–3: 1, fold4: 0.5). The script first reproduces
B7b nested M1-sel 0.6803 (must match to 1e-6). In E1–E3, α*_f, (n_top, κ) and the signature fits are INHERITED, not
re-tuned: each experiment changes exactly one thing.

**Inference (all three)**: ONE comparison each, variant (nested OOF) vs core-v2 (nested OOF), paired group bootstrap
B = 2000, seed 0, C8 forms (shifted-null one-sided p, basic CI). Hard test (E2, E3): one-sided p < 0.05 for ΔM1-sel > 0.
Point non-inferiority (E2, E3): M1-all, each source H/D/O, each harmony bin and M2 (variant vs core-v2) ≥ −0.02;
T-unseen ≥ −0.01. Everything else descriptive.

### E1 — per-signature permutation-null recalibration of the fused z (calibration fix; not accuracy-motivated)
- For outer fold f and each kinase k with a signature in fit f: F_k = z_core + α*_f·z_sig (core-v2 quantities).
  On SITE-PERMUTED GT conditions of the training folds ≠ f (label-free; `controls.site_permutation` seed 0, the same
  permuted matrix B7b used, but only the columns of folds ≠ f), v_k = mean(F_k²) over cells where both are defined.
  Fused output z = F_k / √v_k (replaces the fixed √(1+α²)). Cells with only z_sig defined: z_sig/√mean(z_sig²) over the
  same permuted cells. Kinases with < 30 defined permuted cells use the pooled v over all signature kinases of fit f.
  Kinases without a signature are unchanged. No HP.
- Evaluation conditions for calibration are disjoint from the fit conditions: fold-f site-permuted GT conditions
  (scored with fold-f design, as in B7b) and the DEV role=null conditions (fold-f nulls with fold-f design).
- **Adoption**: (a) |rate(|z|>1.96) − 5%| for SIGNATURE kinases on fold-f site-permuted data decreases vs core-v2, AND
  (b) signature-kinase rate on DEV nulls stays in [2%, 8%] (only ~300 cells), AND (c) point ΔM1-sel ≥ −0.01.
  No hard test (calibration fix). Reported: overall / per bin / signature kinases (+ per bin, per-kinase median/p90/max).

### E2 — remove top-r shared response factors from the input (idea 1)
- Allowed set for held-out folds H: `coreg.allowed_mask` (DEV, qc ok, not in H, not sharing (PMID, Control, Sample)
  with H; GT, non-GT and null conditions alike — labels never used). (TEST-shared-control exclusion cannot be applied from
  the DEV cache; at G2 the PI applies it from metadata.)
- Matrix: core input y (per-condition wsd standardisation), clipped |y| ≤ 4, sites measured in ≥ 20% of the allowed
  conditions, missing → 0, uncentred; V_r = top-r right singular vectors (site space).
- Projection of any scored condition y: on its measured V-sites β = least squares of clip(y) on V_r; y' = y − V_r β on
  V-sites (other sites unchanged); y' rescaled to RMS 1 over measured sites (wsd convention). Conditions with < 50
  measured V-sites are not projected. Applied to the full input: core z (h only; G unchanged), z_sig, DEV nulls (ρ_e
  and ρ_sig refit on projected nulls with the same rules), permuted conditions (calibration).
- HP r ∈ {1, 2, 3}, nested: inner_f(r) = M1-sel over folds g ≠ f, fold g scored with V fitted on allowed({f, g}),
  B7b signature fits excluding {f, g}, core picks of f, α*_f. Fold f scored with V fitted on allowed({f}).
- Pre-declared check (descriptive, not an HP): same nested r* with V fitted on NON-GT allowed conditions only.
- Adoption: hard test + non-inferiority above.

### E3 — hyperparameter bagging of the core (idea 3)
- Core grid as core-v1: ρ_f = 3 (picked in every fold), τ ∈ {0.06, 0.125, 0.25, 0.5, 1}, ρ_g ∈ {0, 0.25, 1} (15 configs),
  each z_vif with its own ρ_e (DEV nulls). For outer fold f, rank by inner_f (core-v1 z_vif inner M1-sel over folds ≠ f,
  recomputed; top-1 must equal the core-v1 pick) and take the top 3 (n_bag = 3 fixed, no HP).
- z_bag = Σ_i z_i / √(Σ_ij c_ij), c_ij = uncentred correlation of z_i, z_j pooled over site-permuted GT conditions of folds
  ≠ f (label-free; keeps unit null scale). z_bag replaces z_core in the core-v2 fusion (same z_sig, α*_f).
- Adoption: hard test + non-inferiority above.

### Final bundle
core-v2 + every adopted experiment, combined in the order E2 (input) → E3 (core) → E1 (fusion calibration, refitted on the
bundle's permuted training conditions). If more than one is adopted, the combination is evaluated descriptively vs
core-v2 (no further test, no further choice). Full DEV table: M1-sel, M1-all, M1-sel-scored/coverage, M2 vs PSP-z,
harmony bins, source, sign, T-unseen/T-seen, calibration (overall/per bin/signature kinases), runtime.

## Results (appended after the runs; design above unchanged)
- Runs: `scripts/exp_b8_eval.py --only E1|E2|E3` (about 25–60 s each), then `scripts/exp_b8_bundle.py`.
- Raw output: `results/b8/{b8_eval.json, b8_bundle.json, e2_grid.tsv, e3_grid.tsv, eval_E*.log, bundle.log}`.
- GT_FREEZE `1e6f791c…`. core-v2 reproduced exactly (M1-sel 0.680308).
- Ledger branch `B8`: 6 rows (E1; E2 3 per-r + nested; E3). No experiment was re-run.

**Summary**: one comparison each vs core-v2 (nested OOF). Paired group bootstrap B = 2000, C8.

| | ΔM1-sel [95% basic CI] | one-sided p | non-inferiority | decision |
|---|---|---|---|---|
| E1 perm-null calibration | −0.0007 [−0.0016, +0.0006] | (none needed) | all pass | **ADOPT** (calibration criteria a–c all pass) |
| E2 factor removal (r nested) | −0.0016 [−0.0058, +0.0015] | 0.80 | all pass | reject (hard test fails) |
| E3 top-3 core bagging | −0.0012 [−0.0056, +0.0018] | 0.73 | all pass (D −0.014) | reject (hard test fails) |

**E1**
- Per-kinase scale v_k replaces 1+α² (2.0; 1.25 in fold 4). Median v_k per fold was 2.32 / 1.92 / 3.41 / 2.18 / 1.26, range 0.99–5.45.
  So the fixed √(1+α²) under-normalised most signature kinases. The cause is site correlation within signatures plus z_core–z_sig correlation (the independence assumption fails).
- Calibration, |z| > 1.96, on held-out fold-f conditions (fit used folds ≠ f):

  | set | core-v2 | E1 |
  |---|---|---|
  | site-permuted, all kinases (0 / 1–4 / 5–19 / ≥20) | 5.50% (6.6 / 5.8 / 5.3 / 5.3) | 5.33% (6.5 / 5.7 / 5.2 / 4.9) |
  | site-permuted, **signature kinases** (16.2k cells) | **7.04%** (9.3 / 6.1 / 7.4 / 6.9) | **5.10%** (5.4 / 5.0 / 5.5 / 4.8) |
  | per signature kinase (n = 55, ≥ 50 cells): median / p90 / max | 6.8 / 10.1 / 11.3% | 5.1 / 7.3 / 9.2% |
  | DEV nulls, all kinases (3.5k cells) | 6.22% | 6.11% |
  | DEV nulls, signature kinases (312 cells, 7 conditions) | 4.2% | 2.9% |

- Adoption checks:
  - (a) signature-kinase rate 7.0% → 5.1%: PASS;
  - (b) DEV-null signature-kinase rate 2.9%, inside [2, 8]%: PASS. It is low, but 312 cells give a binomial SE of ±1.0 pp;
  - (c) ΔM1-sel −0.0007 ≥ −0.01: PASS.
- Cost: every endpoint is within 0.003. Source H −0.0011, with a CI that just excludes 0.

**E2**
- Inner M1-sel per r is flat: OOF 0.681 / 0.682 / 0.683 for r = 1 / 2 / 3. Nested picks: r = 1, 3, 1, 2, 2.
- Δ by source and bin:
  - H −0.003, D +0.003, O −0.001;
  - bin 0 +0.003, bins 1–4 −0.003, ≥20 −0.006.
- The pre-declared non-GT-factor check gives the same result: Δ −0.0016, p 0.83.
- The top singular vector is dominated by folds 1–4 (σ₁ ≈ 710 against ≈ 320 when fold 0 is in the allowed set). The generic response is either already removed by the wsd/median standardisation and the signature contrast, or it carries target signal.
- The hypothesis (Hijazi +0.01 to +0.02) is not supported: H −0.003.

**E3**
- The top-3 configs are near-duplicates: uncentred null correlations 0.96–0.99. So bagging barely changes z.
- Results: ΔM1-sel −0.001; M1-all +0.004 (p .08); D −0.014 (worst endpoint, still inside −0.02); T-unseen −0.005.
- Selection noise in the core is not the binding constraint once (τ, ρ_g) is restricted to its plateau.

## Final recommended bundle: core-v3 = core-v2 + E1 (DEV nested OOF)
| endpoint | core-v1 (ref) | core-v2 | **core-v3 (bundle)** |
|---|---|---|---|
| M1-sel (primary) | 0.663 | 0.6803 | **0.6796** |
| M1-all | 0.615 | 0.6235 | 0.6229 |
| M1-sel-scored / pair coverage | – | 0.6796 / 0.998 | 0.6789 / 0.998 |
| M2 vs PSP-z (507 pairs, 67 kinases) | 0.687 vs 0.633 | 0.708 vs 0.633 | 0.709 vs 0.633 (Δ +0.075) |
| harmony (bin-macro) | – | 0.6831 | 0.6822 |
| bins 0 / 1–4 / 5–19 / ≥20 | .581/.688/.737/.645 (b21) | .610/.712/.744/.666 | .609/.712/.744/.664 |
| source H / D / O | – | .654 / .916 / .698 | .653 / .916 / .697 |
| sign −1 / +1 | – | .655 / .817 | .655 / .814 |
| T-unseen (272 pairs) / T-seen | 0.648 / 0.668 | 0.648 / 0.692 | 0.648 / 0.691 |
| kinases with / without a signature | 0.734 / – | 0.836 / 0.646 | 0.832 / 0.646 |
| coverage-significant kinases (BH q ≤ .10) | – | 52 | 52 |
| calibration, DEV nulls (all) | 6.4% | 6.2% | 6.1% |
| calibration, site-permuted (all; bins 0 / 1–4 / 5–19 / ≥20) | 5.4% | 5.5% (6.6/5.8/5.3/5.3) | 5.3% (6.5/5.7/5.2/4.9) |
| calibration, site-permuted, signature kinases | 5.6% | 7.0% | **5.1%** |

Runtime (CPU, 1 thread, 200 DEV GT conditions, median of 5):
- core-v1 scorer (`KinaseActivityModel.score`): 7.3 ms/condition. These conditions are denser than the 4.6 ms compendium average.
- Late fusion with the E1 scale (z_sig projection + calibrated fusion): 0.04 ms/condition.
- Total: about 7.35 ms/condition, or about 14 s extrapolated to 1,937 conditions. The E1 scale is precomputed per kinase, so it costs nothing at scoring time.

Notes for G2:
- At TEST, the E1 scales v_k are fitted on site-permuted DEV-CV conditions together with the signatures (all 5 folds).
- The TEST-shared-control exclusion for unsupervised fits does not affect the bundle, because E2 was rejected.
- DEV-confirm look #2 should evaluate core-v3. E1 changes M1 by < 0.001, so D11's conditional applies unchanged.
- Development stops here (D10).
