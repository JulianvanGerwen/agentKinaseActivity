# Pre-freeze review of the one-shot TEST runner (`scripts/run_test.py`)

Reviewer: independent agent (did not write the runner). Date 2026-10-27. Scope: DESIGN §5/§8/§12 C7–C10, RESEARCH_TREE
D11–D13, `scripts/run_test.py`, `src/kinact/{secondary,core,eval,bootstrap,external_eval}.py`, `config/core_v3.json`,
`tests/test_run_test.py`, `results/test_dryrun/`, `results/baselines/baselines_dev.md`. Not done: real-mode run,
reading `cache/sealed/` contents (only the file names were listed, to check the loader's keys), creating FREEZE.json.
Line numbers refer to `scripts/run_test.py` after the fixes below.

**Verdict: the runner is fit to freeze.** The statistical core (Holm family, "beats", bootstrap, comparators) matches
the pre-registration. I fixed 5 robustness/reporting bugs (none changes a Holm-family number; the dry-run Holm table
is byte-identical before and after). There are 8 judgement items for the PI, listed below; none blocks the freeze.

## 1. Holm family, "beats", inference: correct
- Family = exactly 6 tests (l. 435–441): core-v3 vs z-score / KSEA / JUMPsem / best-other / prior-only on M1-sel, and
  core-v3 vs z-score on M2 (matched universe via `paired_delta('M2')`, i.e. M1-sel-scored on S(A)∩S(B)). Matches §5 + D13.
- Holm step-down `bs.holm` is correct (monotone, capped). One-sided p is the C8 shifted-null p for Δ>0; the CI is
  the C8 basic CI. B = 2000 and seed 0 (l. 81). One `BootstrapPlan` is shared by every method, endpoint and masked
  subset (l. 422).
- "Beats" (l. 444–452, `bs.beats`): Holm p < 0.05 AND ΔM1-sel ≥ 0.01 AND point ΔM1-sel(T-unseen) > 0 AND point
  ΔM2 ≥ 0. ΔM2 is the comparator-specific `ev.m2`, the same quantity as the M2 test's point estimate. A NaN in any
  criterion fails the rule (conservative).
- MDE: the pre-registered value comes from `results/G1.json` (C7, 0.060). The realised 2.5·SE is labelled
  descriptive.

## 2. T-unseen / T-family-unseen: consistent with DESIGN, with a caveat (J1)
- Real mode: DEV positives = direct, high/medium positives in role=gt conditions of folds 0–4 **and DEV-confirm**
  (l. 252–264). This matches DESIGN §4 ("never a DEV positive"; DEV = folds + DEV-confirm) and the G1 count
  (`make_split.py:295`, 15 unseen). Dry run: folds 0–4 only, which is right because DEV-confirm is the evaluated split.
- The family map falls back to the kinase itself when Manning family is empty.

## 3. Leakage: none found
- Signatures: the training list is 463 conditions, all in folds 0–4 (none in DEV-confirm). l. 755–757 asserts no
  shared group or condition with the evaluated split. The ρ_e and ρ_sig fits use the 7 fold-0–4 nulls (N-fit) only.
- Prior-only freq uses `training_frequency` over folds 0–4 conditions (l. 269, 310). The sign is the oracle s_c.
- A1 co-regulation is fitted on DEV cache columns with split ∈ folds, qc ok, and no (PMID, Control, Sample) key
  shared with any evaluated condition (l. 338–340). Standardisation is per column (`coreg.standardise`), so excluded
  columns do not leak. σ and R come from site-permuted fold GT conditions.
- N-eval nulls (TEST-group, ~5) are used only in `calib_block` (l. 484), never for fitting. Core scaling is per
  condition (single-sample), so the composition of the TEST batch does not affect core-v3 scores.
- JUMPsem runs on the TEST qc-ok conditions only, in per-study units built from metadata (l. 300–303).

## 4. Comparators: match D13 / B14 / baselines_dev.md
z-score curated_union min1 (0.563, DEV best), KSEA curated_union min1 (0.562), JUMPsem curated_union min3 aligned
(0.507 > native 0.506), and best-other = PhosX-like cantley site-top15 min3 (0.557, best method excluding
z/KSEA/JUMPsem). PSP z min3 and native JUMPsem are descriptive only.

## 5. External evaluation: wired correctly
- `run_external(confirm_test=True)` runs per method. The `rec` store order (perturbation sets, then COHORTS) matches
  the call order in `run_external`.
- Per-set paired bootstrap: B = 2000, seed 0. Perturbation sets resample PMID×drug groups; the tumour set resamples
  kinases.
- `drug_in_dev` means the agent or stimulus appears in any DEV or DEV-confirm condition (`prep_external.py:47–60`,
  TEST excluded). See J3.

## 6. Determinism, hashes, outputs
- The dry run reproduces the Holm table exactly (also equal to DEV-confirm look #2: 0.842, −0.014 vs prior-only).
  Seeds: bootstrap 0, random 0, coverage-significant 0, PhosX GSEA 42.
- FREEZE check covers: sealed GT hash, GT_FREEZE sha, config sha, manifest sha plus all 6 artefacts, HEAD commit,
  and clean tracked src/scripts/config/models/tests.
- Outputs written: LEDGER.md, test_results.json, per_kinase/*.tsv, scores/*.parquet, freeze_check.json. Runtime is
  reported (bench 6.7 ms/cond at 1 thread, 2.0 at 16 threads; compendium 7 s at 1 thread). It is not gated.

## 7. Crash risks
- Memory: the runtime block loads 190k × 1937 float64 (about 3 GB; peak about 15 GB). Fine at 375 GB.
- JUMPsem on TEST: 42 units, max 26 conditions (from metadata), against DEV units of 72 conditions that took
  ≤ 1574 s. Expected time is under about 10 min.
- A1 needs CUDA (`CoregWZ` uses device="cuda"). It runs inside `score_all`, so with no GPU the run aborts before
  scores are saved. The GPU is present and the dry run passes (J8).

## Fixes applied (minimal, with tests in `tests/test_run_test.py`; 12 pass, full suite 162 pass)
1. **Frozen-input re-hash (new `check_frozen_inputs`, l. 125; enforced in `verify_freeze`).** Before this fix, only
   GT_FREEZE.json's own sha was checked. The files the run actually reads (visible gt_pairs/conditions used for DEV
   positives and prior-only freq, condition_qc, harmony bins, `cache/quant_wide_f32.npy`, the DEV cache, and the
   curated_union/r2np libraries in the manifest inputs) were never re-hashed. All of them are git-ignored or outside
   the dirty-tree check. All 34 now match. Comparator-only libraries (cantley, psp_all) and the external sets have no
   frozen hash; they are logged in `freeze_check.json["inputs"]["logged_only"]`.
2. **Masked-subset bootstrap (`delta`, l. 361–381).** `bootstrap.metric_replicates` returns 0, not NaN, for a
   resample that draws none of the subset's groups. For T-unseen/T-seen/has-signature subsets this piled Δ* = 0 into
   the distribution. In the dry run this produced CIs like "+0.494 [+0.494, +0.987]". Such resamples are now dropped
   and `n_resamples` is reported. This affects descriptive CIs and p-values only: "beats" uses the point estimate.
   The frozen harness was not changed (see J7).
3. **JUMPsem timeout** (`src/kinact/methods/jumpsem.py`). `subprocess.TimeoutExpired` used to propagate through the
   thread pool and abort the whole TEST run. Now the unit is recorded as "R timeout" and left NaN.
4. **Secondary blocks cannot kill the run.** The runtime and external blocks run under `secondary_block` (l. 707),
   which records the error and traceback. `test_results.partial.json` is written before them (l. 772). The ledger
   prints "FAILED" for a failed block instead of raising KeyError at the very end. The run-once guard now also
   refuses if the partial file exists (a re-run after a crash is an explicit PI decision).
5. **Ledger shows the external drug seen/unseen strata** (DESIGN §2 requires stratification; it was in the JSON only).

## Judgement issues for the PI (not changed)
- **J1 T-unseen reference set.** At TEST, "unseen" means not a positive in DEV-CV ∪ DEV-confirm (DESIGN-literal,
  G1-consistent). Signatures and prior-only, however, see only DEV-CV. Nine DEV-confirm-only positive kinases
  (AURKC, EIF2AK2, GRK4, MAP3K2, MKNK2, NEK2, PDGFRA, PRKG2, STK10) would count as T-seen without any fitted
  component having seen their labels. Recommendation: keep it as is, and state it in REPORT.
- **J2 Prior-only freq uses folds 0–4 only, not DEV-confirm.** This is defensible ("training-fold conditions";
  symmetric with the signatures). Including DEV-confirm would give a slightly stronger control. Either way, decide it
  explicitly in FREEZE/D13.
- **J3 External stratification is thin.** The seen/unseen split is computed only on an unweighted all-pairs macro
  (labelled "M1-all_drug_*", but it has no 1/n weights and includes n_targets > 10). There is no M1-sel stratum and
  no paired bootstrap per stratum. The tumour set is not stratified (N/A).
- **J4 Calibration binning.** core-v1 and A1 are binned with core-v3's `n_substrates` (identical core library; A1
  expansion edges are not counted). Signature-only kinases have n = 0 and so fall in bin 0 (dry: 0.25 on 8 cells).
- **J5 N-eval has ~5 TEST nulls,** so the null calibration will be noisy. The site-permuted TEST GT calibration is
  the more stable number. Say so in REPORT.
- **J6 Dirty check uses `--untracked-files=no`.** A stray untracked module under src/ would not be detected. Run
  `git status` (including untracked files) by hand before writing FREEZE.json. Optionally add the logged library and
  external hashes to FREEZE.json.
- **J7 Harness `metric_replicates` 0-for-empty behaviour** also affects descriptive `paired_delta` filters
  (`M1-sel:source=…`, `:bin=…`). It is negligible when a stratum has several groups, but bin "0" may be sparse. I left
  it unfixed because the harness is G1-frozen. An amendment is needed if the PI wants it changed.
- **J8 A1 failure is fatal.** A1 (secondary) runs inside `score_all`, so a CUDA or co-regulation failure aborts the
  run before any score is saved. Accept (dry run OK, GPU present), or make it non-fatal before the freeze.

The dry run was re-run after the fixes (`results/test_dryrun/`, 287 s, rc 0). These numbers are a smoke test, not a
DEV-confirm look.
