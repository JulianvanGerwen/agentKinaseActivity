# DESIGN — Kinase activity inference from phosphoproteomics (v0 draft, 2026-10-27)

Status: DRAFT v0, before review round 1. Changes after review go in §12 (amendment log); nothing below is silently edited after freeze.

## 1. Task (fixed by PROMPT.md)
Input: one vector of site-level directional statistics (log2FC, or t/z) from a human phosphoproteomics experiment, keyed by `{UniProt}_{residue}{pos}` (15-mer optional).
Output: for every human protein kinase that the method can score, a signed activity-change score (positive = activated), a standardised z-like scale where 0 = no evidence of change, and a coverage flag (n substrates used).
Desiderata, in order: (1) accuracy on ground-truth (GT) perturbations, macro-averaged over kinases so well-studied kinases do not dominate; (2) number of kinases scored with non-trivial accuracy (coverage), incl. understudied ones; (3) calibration (null conditions give |z|>2 at ≈5%); (4) speed (target: <50 ms per condition for single-sample mode on CPU, whole compendium of 1,937 conditions in < 1 min).
Primary mode is **single-sample** (one condition in, scores out). Multi-sample methods (JUMPsem) are benchmarked in their native mode per study.

## 2. Data
- Compendium `data/03_*` (1,937 conditions, 178 studies, 190k sites, 4.6% dense; centred log2FC). QC: exclude conditions flagged `qc_flag` (unlogged/scaled; `results/data_profile/conditions_profile.tsv`) from evaluation; methods must be robust to missingness (no imputation).
- PSP kinase–substrate (KS) dump (human kinases; human substrates + mouse/rat mapped via SITE_GRP_ID; in vivo/in vitro flags), PSP Regulatory_sites, Phosphorylation_site_dataset.
- Predicted KS (leakage-audited, `docs/kinasesub_audit.md`): kinaseSubMain **R2-nophospho** (DEV-trained, no phospho quantification features), **R3** (sequence/structure neural, no phospho features), **Cantley/Kinase Library** atlas (motif only). **Banned**: kinaseSubMain deployment/release and `clean_bench` variants (they use Hijazi + compendium co-regulation/inhibitor deconvolution), SELPHI2, CoPheeKSA, PTMsigDB perturbation signatures, Hijazi's own downstream marker sets, NetworKIN versions trained on phospho data (if used, only after audit).
- Inhibitor target profiles: `kinhib_inhib_comb__1.csv` (user-provided Hijazi target matrix: primary GT); Klaeger 2017 kinobeads Kd_app (29/61 Hijazi drugs + decryptM drugs) for off-target masking and for dose-aware GT on decryptM.
- External (never used for development): benchmarKIN tyrosine set (18 exps), benchmarKIN Hernandez experiments not in the compendium, benchmarKIN CPTAC tumour benchmark (7 cohorts) — copies in kinaseSubMain `external/`, to be re-copied to our `external/` with checksums in DATA.md.

## 3. Ground truth
GT unit = (condition c, kinase k, expected sign s∈{−1,+1}). Built by `scripts/build_candidate_gt.py` from rule files; **frozen with a SHA-256 before any method is tuned (gate G1)**.
- **GT-H (primary)**: Hijazi 177 conditions, targets = kinhib matrix, s = −1.
- **GT-D**: decryptM dose series (340 conds): kinase k is a target at dose d iff Klaeger Kd_app(drug,k) ≤ d (pre-registered rule; confidence by d/Kd). Lowest doses with no target → negative-only conditions.
- **GT-O**: other inhibitors (high/medium confidence), stimuli (growth factor/cytokine/GPCR/DNA damage; direct receptor kinase and canonical first-tier kinases, s=+1), genetic perturbations, cell-cycle. Low-confidence and `inferred` rows reported only as sensitivity analysis.
- **Negatives** in condition c: all kinases in the evaluation universe U (§5) except: targets of c, kinases with Klaeger Kd_app ≤ 3× dose for that drug (masked, "unknown"), and kinases annotated as known downstream of the target (`kinases_downstream`; masked — their activity legitimately changes).
- **Null controls**: 27 decryptM "0 nM vs 0 nM" conditions → calibration only.
- A domain-expert agent reviews the rule files before G1 freeze.

## 4. Splits and leakage control
Leakage pathways considered: (L1) same drug in train and test (also across studies: Hijazi/Wilkes/decryptM/RPE1 share drugs); (L2) same lab/cell line/drug batch (Wilkes ≈ Hijazi); (L3) predicted KS libraries trained on test phospho data (handled by §2 bans); (L4) any learned component (site weights, co-regulation expansions, hyperparameters) fitted on conditions whose GT kinases appear in test; (L5) evaluator choices tuned on test.
- **Grouping**: union-find over conditions linking (a) the same canonical drug in any study, (b) same study for non-drug studies, (c) Hijazi↔Wilkes 2015 per drug. Each group goes wholly to TEST or DEV.
- **TEST (sealed, ~30% of GT conditions)**: stratified random selection of groups (seed 0), stratified by GT source (H/D/O) and by target-kinase PSP-coverage bin; plus all external sets (§2). TEST labels and conditions hashed; the evaluator refuses TEST until FREEZE.json is written (single, logged TEST run; independent reproduction).
- **DEV**: 5-fold grouped CV over the remaining groups (seed 0). All learned components are fit inside folds only on training-fold conditions. Unsupervised use of compendium conditions (e.g. site co-regulation) may use any non-TEST, non-held-out-fold conditions — never TEST conditions or TEST-group studies.
- **Kinase-novelty tiers** on TEST: T-seen (target kinase is a target in some DEV condition) vs T-unseen. Reported separately; a learned method must not win only on T-seen.
- Seeds: fixed list `SEEDS = [0, 1, 2, 3, 4]` (config/seeds.json). Split uses seed 0 only.

## 5. Evaluation (frozen harness `kinact.eval`, built and verified before methods)
Evaluation universe U = all human protein kinases (HGNC/Manning + PSP + kinaseSubMain list ≈ 520). Each method returns scores for a subset; **unscored kinase ⇒ score 0 (neutral) with tie-aware ranking** in the coverage-aware metrics. All scores are sign-adjusted: t = −s·score (larger = correct).
- **M1 (primary, coverage-aware, kinase-macro)**: for each GT pair (c,k), AUROC_ck = fraction of negatives j in c with t_j < t_k (+0.5 ties). Average over conditions per kinase, then macro-average over GT kinases with ≥1 evaluable pair. Kinase-macro neutralises over-representation of popular targets.
- **M2 (algorithm-fair, matched)**: same as M1 but U restricted to kinases scored by both methods in a pairwise comparison (Δ reported with paired CIs) — separates library coverage from algorithm quality.
- **M3 (benchmarKIN-compatible)**: condition-centric pooled AUROC with balanced subsampling and scaled rank of targets, as in Müller-Dott 2025, for comparison with literature numbers (harness validated by reproducing PSP+z-score on benchmarKIN ≈ 0.78–0.79).
- **M4 (kinase-centric)**: for kinase k, AUROC separating conditions where k is a target vs DEV/TEST conditions where it is a negative, within the same studies; macro over k with ≥2 target conditions. Penalises kinase-specific offsets.
- **Harmony**: M1 stratified by PSP-substrate bins of the kinase {0, 1–4, 5–19, ≥20 quantified substrates}; report each bin and the bin-macro mean. **Coverage@0.6**: number of GT kinases with per-kinase AUROC ≥0.6.
- **Calibration**: on null conditions and on label-shuffled/site-permuted data, fraction |z|>1.96 per kinase (target ≈5%).
- **Runtime**: wall-clock per condition and full-compendium, same hardware, single process + note threads.
- **Inference**: cluster bootstrap over GT groups (B=2000, seed 0) for CIs; paired bootstrap for Δ. Pre-registered superiority tests on TEST: final vs {z-score, KSEA, JUMPsem, best other baseline} on M1 (and M2 vs z-score), Holm-adjusted; "beats" = Holm p<0.05 and Δ>0.
- **Controls in harness tests**: toy data with planted kinase effects (known answer), random-score negative control (M1≈0.5), oracle positive control (M1=1), shuffled-library control.

## 6. Baselines (implemented in Python, each verified on toy data; R reference outputs where feasible)
z-score (RoKAI/benchmarKIN form), KSEA (Casado mean-based z; KSEAapp variant), JUMPsem (faithful reimplementation of the R package, native per-study multi-sample), decoupler-style ulm / mlm / wmean / norm_wmean / viper, fgsea-style GSEA, Wilcoxon/KS rank tests, PTM-SEA-style ssGSEA, RoKAI with network propagation, Kinase Library motif enrichment (Johnson 2023) / PhosX-like motif-only scoring, INKA-lite if time. Each with PSP library (default), and where natural with the extended libraries.

## 7. Method development plan (DEV only; ranked by expected value)
A. **Library**: PSP (in vivo ∪ in vitro, orthologue-mapped) + predicted KS (R2-nophospho, R3, Cantley) as *weighted* edges; weights via calibrated probability or rank; per-kinase cap/top-k; choose mixing α and cut-offs by DEV CV. Kinase families: hierarchical borrowing for understudied kinases.
B. **Statistic**: weighted, variance-aware z (per-condition robust scale), shrinkage for small n, empirical-null calibration; multivariate deconvolution (ridge/elastic-net mlm over the weighted KS matrix) to separate kinases sharing substrates.
C. **Learned components (with grouped CV)**: low-dimensional learned hyperparameters (α, top-k, λ); site reliability/responsiveness priors learned from non-TEST compendium variance; data-driven substrate expansion via site co-regulation across non-TEST compendium conditions (RoKAI-like propagation but data-driven); per-library weights.
D. **Kinase self-sites**: activation-loop/regulatory sites of the kinase itself (PSP Regulatory_sites) as an auxiliary signal — tested separately because inhibitors can paradoxically increase them.
E. **Ensembling** of complementary statistics/libraries; final speed optimisation (sparse matrix ops).
Each idea is an experiment node in RESEARCH_TREE.md with DEV-CV M1/M2/harmony/runtime, compared to the current best by paired fold-level bootstrap. An idea is adopted only if it improves M1 without degrading harmony bins by >0.01.

## 8. Model selection and TEST protocol
Gate G2: the final configuration (code commit + config hash) is chosen on DEV CV only and written to FREEZE.json. TEST is run once; an independent agent reproduces it. Post-TEST explorations are labelled as such and cannot change the reported final.

## 9. Verification
`tests/` with pytest: toy datasets with known answers for every metric and every method; leakage tests (split disjointness by drug/study; no TEST condition in any fit); GT freeze hash check; deterministic seeds. Every result that enters REPORT.md is re-checked by an agent that did not produce it.

## 10. Compute
64 CPU, 375 GB RAM, 4× RTX 4090. All baselines are sparse linear algebra → seconds for the compendium; JUMPsem and viper-style permutations are the slowest (timed first). GPUs only if a learned model needs them.

## 11. Agent organisation & failure-mode guards
PI (me) designs, decides, integrates. Agents: harness+GT (one), baselines (one or two), method explorers (parallel, each owns one tree branch), independent reviewers/auditors. Guards: agents never touch TEST; all numbers cite committed files; summaries ≤300 words; every metric change is an amendment with rationale before seeing its effect.

## 12. Amendment log
(empty)
