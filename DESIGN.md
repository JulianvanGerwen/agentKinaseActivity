# DESIGN — Kinase activity inference from phosphoproteomics (v1, 2026-10-27)

Status: v1 = v0 revised after review round 1 (docs/reviews/r1_*.md). Amendment log in §12. After gate G1 (GT + split + harness freeze) any change to metrics, labels, splits or tests requires a dated amendment with rationale written BEFORE its effect is seen.

## 1. Task (fixed by PROMPT.md)
Input: one vector of site-level directional statistics (log2FC, or t/z) keyed by `{UniProt}_{residue}{pos}` (15-mer optional).
Output: for every kinase in universe U the method can score, a signed activity-change score (positive = activated) on a calibrated z scale (0 = no evidence), plus n_substrates/coverage and a posterior SD where available.
Desiderata: (1) accuracy on GT perturbations, macro-averaged over kinases; (2) number of kinases scored with demonstrable (above-null) accuracy, incl. understudied ones; (3) calibration (null conditions: |z|>1.96 ≈5% per kinase, independent of substrate count); (4) speed: ≤5 ms per condition end-to-end, ≤10 s for the 1,937-condition compendium (single-threaded and 16-thread, warm, median of 5; load/mapping, scoring and null timed separately; `bench_runtime.py` on 200 fixed DEV conditions for all methods).
Deliverable: Python package + CLI (`kinact score in.tsv > out.tsv`), precomputed activity matrix for the whole compendium, and a kinase–kinase relationship analysis (do score correlations recover known kinase→kinase links?).
Primary mode: single-sample. Multi-sample methods (JUMPsem) run natively per study.

## 2. Data
- Compendium `data/03_*`. QC: exclude from evaluation only conditions with demonstrable scale errors (unlogged; ≈3–4 GT conditions) and exact-duplicate vectors (35688186 caffeine=quercetin; 22115753_1=_2 keep one). Other `qc_flag`s are reported, not excluded. No imputation.
- Kinase universe **U = the 503 human protein kinases (include==True) of kinaseSubMain's kinase list**, HGNC-normalised. Lipid kinases (PIK3C*, PI4K*, PIP4K*, PIKFYVE…) are outside U.
- KS libraries (all through ONE shared library builder used by every method incl. baselines):
  - PSP KS (human kinases; human substrates + mouse/rat via SITE_GRP_ID; in vivo/in vitro flags).
  - Additional curated / in-vitro sources from kinaseSubMain `external/`: OmniPath, iPTMnet, Sugiyama 2019, iKiP-DB (provenance + licence in DATA.md; edges whose only supporting reference is a compendium study are dropped).
  - Predicted (leakage-audited, `docs/kinasesub_audit.md`): kinaseSubMain R2-nophospho, R3, Cantley/Kinase Library atlas. Input as top-k *measured* sites per kinase, weighted by percentile^γ.
  - Banned: kinaseSubMain release/`clean_bench`, SELPHI2, CoPheeKSA, PTMsigDB perturbation signatures, Hijazi downstream marker sets, any kinaseSubMain setting tuned on benchmarKIN. Checksums of all used predicted files recorded in DATA.md; confirm kinaseSubMain did not select those models on Hijazi (audit note).
- Drug–target profiles for masking and dose-aware GT, precedence: Klaeger 2017 (Kd_app) > Reinecke 2024 > Davis 2011 > kinoplex > user-provided Hijazi matrix only. Klaeger parser must expand all genes of a protein group.
- External (never used for development): benchmarKIN tyrosine set (18 exps), the 3 Hernandez studies not in the compendium, CPTAC tumour benchmark (7 cohorts). Truly independent; reported as secondary, stratified by drug seen/unseen in DEV.

## 3. Ground truth (built by script from rule files; frozen with SHA-256 at G1 after expert-agent review)
GT unit = (condition c, kinase k ∈ U, sign s). Classes per (c,k): **direct positive**, **downstream positive** (secondary only), **masked** (unknown; excluded from both), **negative**.
- **GT-H (primary source)**: Hijazi 183 conditions. Direct positives = user kinhib matrix ∩ U, s = −1.
- **GT-D**: decryptM. Positives = manually curated primary target list per drug (documented rule file) at doses ≥ Kd_app(drug,k) where Kd is known; allosteric inhibitors and PIKKs exempt from the Kd rule (dose ≥ published IC50 rule instead). Lowest doses below every target Kd → negative-only. Only high/medium confidence pairs enter metrics.
- **GT-O**: other inhibitors, stimuli (receptor kinase + canonical first-tier kinases, s=+1), genetic, cell-cycle; high/medium confidence only; low/`inferred` = sensitivity analysis.
- **Masks (all sources, incl. Hijazi/Wilkes)**: (i) pre-registered downstream/feedback mask table by target class — MEK/ERK/RAF inhibitors → ERK1/2, RSKs, MNKs, MSKs, RAFs; PI3K/PDK1/AKT inhibitors → AKT, SGK, S6K, mTOR, GSK3; mTOR → S6K, AKT, ULK1; p38 → MK2/3; plus canonical cascades for stimuli and cell-cycle; (ii) paradoxical RAF activation: in BRAF-wt cells (all Hijazi lines) RAF-inhibitor RAF positives are dropped and MEK/ERK/RSK masked; (iii) off-targets: kinases with Kd_app ≤ 3× dose (any profile source, precedence above) that are not annotated targets; (iv) negative paralogs sharing ≥80% of PSP substrates with a positive.
- Lipid-kinase-only drugs (CAL101, TGX221, PIK294), Edelfosine and U73122 → no positives; used as **specificity controls** with PKC/AKT-axis masked where biologically indicated (PLC inhibition lowers PKC).
- **Null controls**: 27 decryptM "0 nM vs 0 nM" conditions, split into two halves (seed 0): N-fit (allowed for empirical-null fitting) and N-eval (calibration reporting only). Site-permutation nulls for development.
- **GT-A (secondary, label-independent)**: in non-inhibitor conditions, a kinase's own activating sites (PSP Regulatory_sites "enzymatic activity, induced") as activity proxy; the kinase's own sites removed from its library for this test. Gives evaluation for understudied kinases.
- Validation tests before G1: every GT kinase ∈ U; no bad tokens (BCR/ABL→ABL1, Q6ZSR9 dropped, ALL/no-kinase rows handled); s ∈ {−1,+1}; masks applied; summary counts per source/confidence.

## 4. Splits and leakage control
- **Grouping**: union-find over all 1,937 conditions linking: same canonical drug anywhere (treatment or control field), same study for non-drug studies, Hijazi↔Wilkes by drug; identical-vector duplicates. Groups stratified by majority GT source.
- **Split selection rule (pre-registered)**: generate 100 candidate splits (seeds 0–99) with TEST ≈ 30% of GT pairs; choose the one minimising a label-only distance to targets (TEST share 25–35% of pairs per source and per PSP-coverage bin; ≥20% of TEST GT kinases unseen in DEV). Uses no method outputs. State the minimum detectable ΔM1 on the chosen TEST at G1.
- **Physical sealing**: build a DEV-only cache (TEST conditions removed) that development agents use; pytest guard fails if any TEST condition ID appears in DEV artefacts; the evaluator refuses TEST until FREEZE.json exists.
- **DEV**: 5-fold grouped CV (seed 0), plus a DEV-confirm hold-out of ~15% of DEV groups, evaluated at most 5 times total (logged).
- **Unsupervised use of compendium conditions** (co-regulation, factorisation, site noise variances): allowed on non-TEST conditions excluding those sharing a control with a TEST condition; inside CV, additionally excluding the held-out fold.
- **Kinase-novelty tiers**: T-seen / T-unseen (kinase never a positive in DEV) / T-family-unseen.
- Seeds: `SEEDS=[0,1,2,3,4]` (config/seeds.json). Winner's-curse check before FREEZE: re-split DEV folds with seeds 1–4.

## 5. Evaluation (frozen harness `kinact.eval`, verified before methods)
Sign-adjusted score t = s·score (larger = more correct). Unscored kinases: an unscored target contributes AUROC ½; each unscored negative contributes ½ to a scored target's AUROC (shift/scale invariant).
- **Pair AUROC**: AUROC_ck = P(t_k > t_j) over negatives j in c (+½ ties/unscored).
- **Averaging order**: pairs → mean within (kinase, GT group) → mean over groups per kinase → macro over kinases.
- **M1-sel (PRIMARY)**: direct positives from conditions with ≤10 direct targets in U, all GT sources (high/medium).
- **M1-all (key secondary)**: all direct positives, pairs weighted by 1/n_targets(c) within group.
- Also reported: M1-scored (only kinases the method scores) + pair coverage; M1 by source (H/D/O) and by sign (s=−1, s=+1); M1-family (negatives restricted to same-family kinases); downstream-positive AUROC; M2 (matched universe of kinases scored by both methods); M3 benchmarKIN-compatible (pooled condition-centric AUROC with balanced subsampling + scaled rank, validated on DEV/benchmarKIN with benchmarKIN's own target lists); M4 kinase-centric (within study, scored kinases only, TEST conditions only on TEST); GT-A AUROC.
- **Harmony**: M1-sel per PSP-substrate bin {0, 1–4, 5–19, ≥20 quantified PSP substrates; bins from PSP only, frozen with GT hash} and the bin-macro mean.
- **Coverage-significant**: number of GT kinases whose mean AUROC beats a rank-permutation null (BH q ≤ 0.10).
- **Calibration**: on N-eval nulls and site-permuted data: fraction |z|>1.96 overall and per substrate-count bin.
- **Negative controls in harness**: random scores (M1≈0.5), **prior-only control** (kinase's target frequency in training conditions, no phospho data; ≈0.63 leave-drug-out) — every learned method must beat it and must also improve M4; site-permuted input (≈0.5); sign-flip test (negating input negates scores); oracle (=1).
- **Inference**: stratified (by source) cluster bootstrap over GT groups, B=2000, seed 0, shared resamples across methods; one-sided p. TEST family of 5 pre-registered tests, Holm: final vs z-score, KSEA, JUMPsem, best-other-baseline (chosen on DEV at G2) on M1-sel, and final vs z-score on M2. "Beats" = Holm p<0.05 AND ΔM1 ≥ 0.01 AND ΔM1 > 0 on T-unseen.

## 6. Baselines (shared library builder; same small DEV tuning grid as our method, e.g. library variant, min substrates)
Must (time-boxed): z-score (benchmarKIN/RoKAI form), KSEA (Casado mean-based z; KSEAapp), JUMPsem (try R + lavaan via micromamba for 30 min to run published package; else Python reimplementation validated on toy data), decoupler ulm / mlm / wmean / norm_wmean / viper (pip decoupler as reference; Spearman ≥0.99 vs our implementation, else labelled unverified), fgsea-style GSEA, RoKAI with network propagation, Kinase Library motif enrichment / PhosX-like motif-only.
Could: KSTAR, KARP, PhosR, IKAP. Cut: PTM-SEA, separate Wilcoxon/KS, INKA.

## 7. Method development (DEV only)
- **Core model (B)**: Bayesian linear-Gaussian x_c = W·a_c + b_c + ε, per-site noise variances estimated label-free on allowed conditions, prior a ~ N(0, τ²K) with kinase–kinase kernel K (family tree / motif similarity / substrate overlap). Closed-form posterior mean and z (one solve per condition; covers shrinkage, mlm-deconvolution, family borrowing, multi-target drugs). ≤5 hyperparameters.
- **Library (A)**: weighted union of curated + predicted edges; top-k and γ selected by nested CV.
- **Calibration**: effective-m analytic null and/or per-condition empirical null (fitted on N-fit/site permutations only) to remove √m inflation.
- **Learned parts (C)**: NO per-edge or per-kinase free weights. Edge weights only as a shared function of edge/site features (≤20 parameters); >50 parameters needs leave-kinase-out non-inferiority first. Data-driven substrate expansion via co-regulation / factorisation fitted only on allowed conditions (§4).
- **Kinase self-sites (D)**: activating sites as auxiliary signal, tested separately (paradoxical behaviour under inhibitors).
- **Ensembling (E)**: learning-to-rank stacking, pairwise within condition, ≤30 parameters, no kinase-specific params, pairs weighted to match kinase-macro M1.
- **Adoption rule**: pooled out-of-fold predictions, group-level paired bootstrap: 95% CI lower bound of ΔM1-sel > 0, ΔM1-all ≥ 0, sign agreement on ≥2 of 3 GT sources, no harmony bin with 95% CI upper bound of Δ < 0; gains must hold on T-unseen kinases (within-CV analogue: kinases not positive in training folds). Nested CV for internal hyperparameters.
- **Ledger**: every DEV evaluation appended to `results/dev_ledger.tsv` (config hash, commit, metrics, n configs tried per branch). Any ΔM1 ≥ 0.02 audited by a non-producing agent before adoption.

## 8. Model selection and TEST
G2: final config chosen on DEV CV + DEV-confirm, frozen in FREEZE.json (commit + config hash + GT hash). TEST run once; independent agent reproduces; report DEV→TEST drop. Post-TEST explorations labelled, cannot change the final.

## 9. Verification
pytest: toy data with known answers for every metric and method; sign-flip and permutation tests; leakage tests (drug/study/control disjointness, DEV-only cache); GT/evaluator hash check at each agent start; deterministic seeds. REPORT numbers generated by script from committed result files; LITERATURE.md PMIDs/DOIs checked by script.

## 10. Compute
64 CPU, 375 GB RAM, 4× RTX 4090. Sparse z-score over the whole compendium measured at 0.09 s (algo review). Bayesian core ≈7 ms/condition. GPUs only if needed.

## 11. Agent organisation & guards
PI designs, decides, integrates, keeps RESEARCH_TREE.md. Agents: GT+split builder; harness builder; baselines; method explorers (one tree branch each); independent auditors. Agents work from the DEV-only cache, never touch TEST, cite committed files, return ≤300-word summaries. Baselines split into must/could tiers and time-boxed.

## 12. Amendment log
- A1 (r1, stats): sign-adjusted score corrected to t = s·score (v0 had −s, which rewarded wrong signs).
- A2 (r1, stats): unscored kinases contribute ½ (shift/scale invariant) instead of score 0.
- A3 (r1, stats/data): averaging order pair→(kinase,group)→kinase→macro.
- A4 (r1, PI/priorwork/stats/data): primary metric restricted to conditions with ≤10 direct targets (M1-sel); all-target M1 with 1/n weighting as key secondary. Rationale: 44–70% of Hijazi kinases come only from promiscuous drugs profiled at 200 nM but dosed at 1–10 µM, least reliable labels. Decision (PI): ≤10 cutoff, not 1/n weighting, for the primary, because it is simpler and its labels are most trustworthy; coverage of understudied kinases is still measured by M1-all and harmony.
- A5 (r1): U = 503 kinaseSubMain protein kinases; lipid kinases and bad tokens out of positives.
- A6 (r1, priorwork/data/aiagent): downstream/feedback masks for all sources incl. Hijazi; paradoxical RAF rule; multi-source off-target masking; paralog mask; downstream pairs as secondary.
- A7 (r1, data): drug linking across all conditions incl. controls; pre-registered 100-candidate split selection; minimal QC exclusion; duplicate vectors removed; co-regulation exclusion rule.
- A8 (r1, stats/algo/aiagent): Coverage@0.6 replaced by permutation-null significant coverage; bootstrap/Holm/min-gain rule; DEV adoption rule, ledger, DEV-confirm set; prior-only negative control; split null controls.
- A9 (r1, algo): runtime targets tightened to ≤5 ms/condition, ≤10 s compendium; Bayesian linear-Gaussian core; ban on per-edge/per-kinase learned weights.
- A10 (r1, PI): extra curated KS sources; more drug-profile sources; GT-A self-site secondary GT; CLI/API deliverable and kinase–kinase analysis; baseline must/could tiers.
- Rejected/deferred: none of the must items rejected. Deferred: 1/n weighting as primary (kept as M1-all).
