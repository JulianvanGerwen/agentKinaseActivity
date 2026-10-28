# REPORT — Fast, broad-coverage kinase activity inference from phosphoproteomics

Status: final for this project phase (2026-10-28). All TEST numbers come from the single sealed run in `results/test/LEDGER.md`:
- The freeze file is `FREEZE.json`, pinned to code commit `03ab343`; the results were committed in `be2ef15`.
- An independent agent reproduced the run bit-identically. Only runtime fields differed.

The research record is in `RESEARCH_TREE.md` (decisions D0–D14), the design in `DESIGN.md` (v3 plus amendments C7–C10), the data sources in `DATA.md`, and the process lessons in `LESSONS.md`.

## 1. Summary

**core-v3** takes one vector of site-level log2FC (or any signed statistic) and returns a calibrated activity z-score for all **503 human protein kinases**. It runs in **6.7 ms per condition on one CPU thread**, 2.0 ms on 16 threads and 1.7 ms on GPU. The whole 1,937-condition compendium is scored in 2.9 s.

On the sealed TEST set it scores M1-sel **0.714**. The set has 468 conditions, of which 245 have ground truth, giving 470 kinase–condition pairs over 67 kinases, 15 of them unseen in development.

| method | M1-sel (TEST) | Δ core-v3 [95% CI] | Holm p | "beats" (pre-registered) |
|---|---|---|---|---|
| **core-v3** | **0.714** | — | — | — |
| z-score (curated, min 1) | 0.584 | +0.130 [+0.090, +0.169] | 0.003 | **yes** |
| KSEA (curated, min 1) | 0.582 | +0.132 [+0.091, +0.171] | 0.003 | **yes** |
| PhosX-like (best other baseline) | 0.549 | +0.165 [+0.119, +0.202] | 0.003 | **yes** |
| JUMPsem (best DEV config) | 0.505 | +0.209 [+0.149, +0.262] | 0.003 | **yes** |
| z-score, PSP only (min 3; classic setup) | 0.527 | +0.187 [+0.136, +0.237] | descriptive | — |
| prior-only control (drug-target frequency, oracle sign) | 0.712 | +0.002 [−0.036, +0.083] | 0.50 | **no** |
| random | 0.487 | +0.227 | descriptive | — |

**Beats rule.** "Beats" requires all four of the following: Holm p < 0.05, ΔM1-sel ≥ 0.01, a gain on unseen kinases, and no loss on the matched-coverage metric M2. core-v3 meets the rule against all three named benchmarks (z-score, KSEA and JUMPsem) and against the best other baseline.

**What it does not beat.** core-v3 does not beat the trivial prior-only control. That control knows which kinases are frequent drug targets and the true sign of the perturbation. It is excellent on heavily studied kinases (TEST T-seen 0.844 vs 0.759) and useless on new ones (T-unseen 0.253 vs core-v3 0.557). We deliberately did not add any popularity prior to close this gap (D5), because the brief asks for harmony across kinases, not exploitation of annotation bias.

## 2. Problem set-up and evaluation

**Ground truth (`results/ground_truth/`, `GT_FREEZE.json`).**
- Positives:
  - Kinase-inhibitor conditions (Hijazi and others in the compendium) labelled with a dose-aware target set. The Hijazi targets come from the paper's annotation; Klaeger, Reinecke and Davis occupancy data add off-targets.
  - decryptM dose–response conditions.
  - First-tier stimulus kinases from a pre-registered table (D2).
- Expected direction: inhibited kinases should decrease.
- Masks: downstream and paradoxical cases are masked.
- An independent audit corrected 12 rule errors before the split was built.

**Splits.** Studies, drugs and controls are linked into union-find groups so that no drug or shared control crosses splits.
- TEST: seed 90, sealed in `cache/sealed` and opened once.
- DEV: 5 cross-validation folds plus a "devconfirm" set used for at most 5 pre-declared looks; 2 were used.

**Primary metric, M1-sel.**
- What it measures: for each condition, how well the kinase score separates its targets from all other kinases. The score is AUROC with the true sign applied, over conditions with ≤ 10 direct targets.
- Averaging: pair → (kinase, group) → kinase → macro over kinases, so heavily studied kinases cannot dominate.
- Unscored kinases count ½. Coverage is part of the objective, which is why many published methods sit near chance.

**Secondary metrics.**
- M1-all (all conditions, 1/n weighting).
- M2: matched-coverage, scoring only kinases both methods cover.
- Harmony: mean over bins of quantified PSP substrate count (0, 1–4, 5–19, ≥ 20).
- T-unseen: kinases never positive in DEV. T-family-unseen: no family member positive in DEV.
- Calibration: fraction of |z| > 1.96 on null and site-permuted data.

**Statistics.**
- Bias-corrected cluster bootstrap over GT groups (C8), B = 2000, seed 0.
- A pre-registered Holm family of 6 tests.
- Minimum detectable effect about 0.06 (C7).

**Controls.** Every method is reported next to the prior-only, shift-prior, random, site-permuted, shuffled-library and residue-matched-library controls.

**Verification.**
- The harness was re-implemented from scratch by a second agent, with identical results.
- It reproduces the published benchmarKIN number (0.775) on benchmarKIN's own data.
- Reference implementations agree with ours: decoupler and KSEAapp at Spearman ≈ 1.0 for the z-score and KSEA baselines; the R JUMPsem package was used directly.
- Toy-data tests have known answers (`tests/`).

## 3. Approach (core-v3)

1. **Kinase–substrate library.**
   - Curated: the union of PhosphoSitePlus with other curated sources. 405 PMIDs from compendium and high-throughput MS studies are excluded to prevent leakage. Coverage: 439 kinases.
   - Predicted: the top 100 sites per kinase from the kinaseSubMain "R2 no-phospho" model, at weight α = 1. That model was audited as leakage-safe (`docs/audit/r2np_leakage.md`) because it does not use phosphoproteomics features.
   - Every one of the 503 kinases therefore has evidence.
2. **Bayesian linear-Gaussian core.**
   - Model: x = W·a + b + ε, with prior a ~ N(0, τ²K).
   - K is a kinase-family kernel (ρ_f = 3; ρ_g = 0.25), so understudied kinases borrow strength from relatives.
   - Shared substrates are explained away jointly, unlike per-kinase z-score or KSEA.
   - Closed-form posterior; no per-kinase trained weights.
3. **Calibrated output (z_vif).**
   - The posterior is turned into a z-score with a variance-inflation term, ρ_e ≈ 0.01, fitted on DEV null conditions.
   - Result: about 5% |z| > 1.96 on site-permuted data.
4. **Late-fused perturbation signatures (B7b, amendment C10).**
   - For kinases with clean perturbation data in the DEV training conditions, a signature of their top 100 responsive sites is correlated with the input.
   - It is added to the core z at α = 1.
   - Signatures are fitted on DEV-CV ground-truth conditions only. TEST and devconfirm conditions were never used.
5. **E1 rescaling.** Each signature kinase gets a scale from a permutation null, which restores calibration for those kinases.

**Deliverables.**
- API `kinact.core.KinaseActivityModel.from_config("config/core_v3.json").score(df)` and CLI `kinact score`; see `docs/usage.md`.
- PSP-derived artefacts are licensed CC BY-NC-SA and not committed. They are rebuilt by `scripts/build_core_v3_artifacts.py`, and their hashes are checked against `models/core_v3/manifest`.

## 4. Results

### 4.1 Development (DEV, nested CV; selection data)

| step | DEV M1-sel | decision |
|---|---|---|
| best published baseline (z-score curated min 1) | 0.563 | — |
| prior-only control | 0.687 | bar |
| B1 extended library (z-score + r2np) | 0.638 (nested) | led to the core library |
| B2 Bayesian core + family kernel | 0.655 | adopted |
| B2.1 ρ_f cap, z_vif output (core-v1) | 0.663 | adopted (D6/D7) |
| B7b late signatures (core-v2) | 0.680, Δ +0.018, p 0.0085 | adopted after audit (D11) |
| B8 E1 calibration (core-v3) | 0.680 | adopted (calibration fix, D12) |
| re-split seeds 1–4 | 0.680 ± 0.004 | no selection inflation |
| devconfirm look 2 | core-v3 0.842 vs core-v1 0.803 (p 0.016) | confirmed |

Negative or not adopted:
- B4 co-regulation expansion (−0.010).
- B4b learned edge weights (nested pick θ = 0).
- B5 self-sites (+0.003, n.s.).
- B6 stacking: A1 +0.012, p 0.057, so it became a pre-registered secondary.
- B7 signatures inside the kernel (failed on unseen kinases).
- E2 shared-factor removal and E3 hyperparameter bagging.
- Winsorising and posterior-mean output.

Reports are in `results/b*/`.

### 4.2 Sealed TEST (`results/test/LEDGER.md`)

| method | M1-sel | M1-all | harmony | bins 0 / 1–4 / 5–19 / ≥20 | Hijazi / decryptM / other | T-seen | T-unseen | T-fam-unseen |
|---|---|---|---|---|---|---|---|---|
| **core-v3** | **0.714** | 0.631 | 0.725 | 0.72 / 0.66 / 0.70 / 0.82 | 0.650 / 0.909 / 0.737 | 0.759 | 0.557 | 0.661 |
| core-v1 | 0.693 | 0.628 | 0.697 | 0.69 / 0.63 / 0.69 / 0.78 | 0.653 / 0.898 / 0.700 | 0.731 | 0.559 | 0.666 |
| z-score curated | 0.584 | 0.544 | 0.590 | 0.52 / 0.57 / 0.60 / 0.68 | 0.584 / 0.692 / 0.575 | 0.608 | 0.503 | 0.498 |
| KSEA curated | 0.582 | 0.544 | 0.588 | 0.52 / 0.57 / 0.60 / 0.66 | 0.575 / 0.692 / 0.573 | 0.604 | 0.503 | 0.499 |
| PhosX-like | 0.549 | 0.529 | 0.559 | 0.56 / 0.51 / 0.52 / 0.64 | 0.519 / 0.610 / 0.541 | 0.560 | 0.511 | 0.528 |
| JUMPsem | 0.505 | 0.506 | 0.506 | 0.50 / 0.51 / 0.52 / 0.50 | 0.503 / 0.510 / 0.505 | 0.507 | 0.498 | 0.502 |
| prior-only | 0.712 | 0.678 | 0.750 | 0.61 / 0.71 / 0.81 / 0.87 | 0.690 / 0.854 / 0.780 | 0.844 | 0.253 | 0.254 |

**Per kinase (67 TEST kinases).**

| method | kinases with M1-sel > 0.5 | kinases ≥ 0.7 | median |
|---|---|---|---|
| core-v3 | 55 | 37 | 0.753 |
| z-score | 45 | 12 | 0.571 |

core-v3 beats z-score on 51 kinases and loses on 15.

**Harmony.** The largest gain over the published methods is in the kinases with no quantified PSP substrates (bin 0: 0.72 vs 0.52). That is the understudied regime the brief emphasises.

**Paired secondaries.**
- core-v3 vs core-v1: +0.021 [+0.001, +0.031]. The gain comes entirely from signature kinases: +0.050 on 311 pairs, versus −0.001 for kinases without a signature.
- The A1 stacking secondary (0.686) did not improve on core-v1.
- M2 (matched coverage) vs z-score is +0.031 [−0.013, +0.073], not significant (Holm p 0.15). Much of the advantage therefore comes from **scoring more kinases well**; on the kinases both methods cover, core-v3 is only modestly better.

**Calibration (|z| > 1.96, target 5%).**
- Site-permuted TEST conditions: 4.8% (bins 4.2–5.2%).
- The 5 TEST null conditions (control vs control): 8.1%. z-score and KSEA give a similar 7.8–8.0%.
- Real "null" experiments carry biological and batch variance that site permutation does not. Treat |z| ≈ 2 as suggestive, not significant.

**Runtime.**

| method | time per condition |
|---|---|
| core-v3, CPU 1 thread (200 DEV conditions, median of 5) | 6.7 ms |
| core-v3, CPU 16 threads | 2.0 ms |
| core-v3, GPU (RTX 4090) | 1.7 ms |
| JUMPsem (R) | 6–12 s |

core-v3 is in the same speed class as z-score and KSEA.

### 4.3 External sets (secondary, independent of the compendium)

| set | core-v3 | core-v1 | z-score | KSEA | PhosX-like | PSP z (min 3) |
|---|---|---|---|---|---|---|
| Tyrosine-kinase perturbations (18 exps), M1-sel | **0.737** | 0.737 | 0.526 | 0.524 | 0.526 | 0.504 |
| Hernández (new studies, 15 exps), M1-sel | 0.552 | **0.613** | 0.568 | 0.568 | 0.498 | 0.509 |
| CPTAC tumours (7 cohorts; gold = protein abundance, C9), M1-kinase | 0.672 | **0.677** | 0.571 | 0.564 | 0.540 | 0.551 |

The same sets scored with benchmarKIN's original metric are in the ledger, where core-v3 is 0.817, 0.568 and 0.698.

- **Tyrosine and CPTAC:** core-v3 is clearly better than every baseline (Δ +0.10 to +0.23, p < 0.001).
- **Hernández:** this set is small (7 kinases) and core-v3 is not better than z-score (−0.016, CI −0.11 to +0.14). The signature fusion *lowered* it relative to core-v1 (−0.061, n.s.), which is the one external warning sign for the signature component.
- **benchmarKIN metric on the Hernández and tyrosine sets:** PSP z-score also scores well on the original metric (0.718–0.753), because that metric ignores coverage. Our coverage-aware metric is stricter by design (D1/D3).

## 5. Limitations

1. **No gain over a popularity prior.** On TEST the prior-only control (drug-target frequency with oracle sign) ties core-v3 overall. It is better on heavily studied kinases and far worse on unseen ones. Users interested only in canonical targets get little extra over "what is usually inhibited"; the value is in breadth and in data-driven evidence.
2. **Gains are concentrated.**
   - The signature component helps about a quarter of kinases; 71% of its DEV gain came from 5 kinases (MAP3K7, ROCK1/2, MET, IKBKB).
   - Signatures depend on which kinases were perturbed in the compendium, which is itself a study bias, even though fitting was restricted to DEV conditions.
   - On the Hernández external set the signature fusion hurt.
3. **Kinase–kinase specificity.** Within-condition AUROC between kinases falls to 0.53–0.60 once kinases are stratified by substrate overlap. Much of the apparent specificity between related kinases (CDKs, MAPKs, AGC) is shared-substrate signal. Treat paralog-level calls with caution.
4. **Ground-truth limitations.**
   - Hijazi doses could not be recovered and default to 10 µM, so off-target sets are approximate.
   - Some Hijazi target kinases may not be expressed in the cell line, so some "positives" cannot change.
   - Hijazi remains the hardest source: 0.65 on TEST.
   - decryptM scores are inflated by the dasatinib / tyrosine composition: a residue-matched random library already reaches about 0.70 there.
5. **Calibration on real nulls** is about 8% rather than 5%.
6. **Predicted library dependency.** Breadth relies on one predicted kinase–substrate model (kinaseSubMain R2 no-phospho). It was audited as leakage-safe for our evaluation, but its errors propagate into every understudied kinase.
7. **Small external sets.** Tyrosine has 16 kinases and Hernández 7, so external CIs are wide.
8. **Licensing.** The PSP-derived model edges are CC BY-NC-SA and not redistributed. Users must rebuild them with their own PSP download.

## 6. Follow-ups and remaining ideas

Untested or only partly tested ideas, ranked by expected value (from `docs/retro/ideas_m2.md` plus TEST observations). TEST is now spent, so any new evaluation needs a fresh, independent test set: new external perturbation studies, or a future benchmarKIN release.

1. **Fix the signature weakness seen on Hernández.**
   - Shrink signature weight by how well a signature transfers across cell types, for example a cross-study consistency score.
   - Or require at least 2 independent studies per signature.
2. **Sparse (ARD/Laplace) prior via reweighted ridge.** This would sharpen calls between paralogs and hubs in multi-target inhibitor conditions (Hijazi, bin ≥20).
3. **Adaptive predicted-library depth.** Use a score-percentile threshold with a floor instead of a fixed top 100, to help bin 0 and unseen kinases.
4. **Context priors without popularity.** Weight the library by kinase expression in the sample's cell type, using proteomics or RNA from the dataset itself or from CCLE. This addresses the unexpressed-target problem without encoding study bias.
5. **Better null calibration.** Fit ρ_e on real control-vs-control experiments from many labs rather than site permutation, to bring 8% down to 5%.
6. **More independent ground truth.**
   - Recover the true Hijazi doses.
   - Add new dose-resolved inhibitor datasets (decryptM successors) and CRISPR/kinase-dead perturbations.
   - Add CPTAC activity proxies beyond protein abundance, such as activation-loop phosphorylation.
7. **Multi-sample mode.** For time courses and dose series, use a joint model across conditions (a JUMPsem-like idea, but fast), sharing activity structure.

## 7. Reproducibility

- **Environment and data.** Python env `/home/plm_env`, with the R environment in `external/envs/r`. Sources, versions, checksums and licenses are in `DATA.md`; seeds are in `config/seeds.json`.
- **Rebuild steps.**
  - Ground truth: `scripts/build_gt.py` → `scripts/make_split.py` → `scripts/make_dev_cache.py`.
  - Model: `scripts/build_core_v3_artifacts.py`.
  - TEST: `scripts/run_test.py`, which refuses to run unless all hashes in `FREEZE.json` and the code commit match.
- **Tests and audits.** Unit and toy-data tests are in `tests/`. Audits are in `docs/audit/`, design reviews in `docs/reviews/`.
