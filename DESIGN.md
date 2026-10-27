# DESIGN — Kinase activity inference from phosphoproteomics (v3, 2026-10-27)

Status: v3 = final design after three review rounds (docs/reviews/r1_*.md, r2_*.md, r3_final.md), merged into one consistent text. Amendment history in §12. After gate G1 (GT + split + harness freeze) any change to metrics, labels, splits or tests needs a dated amendment in §12, written BEFORE its effect is seen.

## 1. Task (fixed by PROMPT.md)
Input: one vector of site-level directional statistics (log2FC or t/z), keyed by `{UniProt}_{residue}{pos}`.
Output: for every kinase in U the method can score, a signed activity-change score (positive = activated) on a calibrated z scale (0 = no evidence), n_substrates, and posterior SD where available. **Unscored kinases are NaN, never 0.**
Desiderata: (1) accuracy on GT perturbations, macro-averaged over kinases; (2) number of kinases scored with above-null accuracy, incl. understudied ones; (3) calibration (|z|>1.96 ≈5% on nulls, independent of substrate count); (4) speed — reported, not a gate: target ≤5 ms/condition end-to-end, ≤10 s for the 1,937-condition compendium (single- and 16-thread, warm, median of 5; load/mapping, scoring and null timed separately; `bench_runtime` on 200 fixed DEV conditions for all methods).
Deliverables: Python package + CLI (`kinact score in.tsv > out.tsv`); precomputed activity matrix for the compendium; kinase–kinase relationship analysis (§5, descriptive).
Primary mode: single-sample. JUMPsem runs natively per study (multi-sample), on the evaluated split only.

## 2. Data
- **Compendium** `data/03_*`. Evaluation excludes only conditions with demonstrable scale errors and exact-duplicate vectors (keep one). Other QC flags are reported. No imputation. "Measured" sites are determined per input at runtime.
- **Universe U** = 503 human protein kinases (include==True in kinaseSubMain list), HGNC-normalised; Manning group/family from the same list. Lipid kinases are outside U.
- **KS libraries** — ONE shared library builder serves all methods incl. baselines:
  - PSP KS (human kinases; human + mouse/rat substrates via SITE_GRP_ID; orthologue duplicates collapsed; in vivo/in vitro flags).
  - Curated union: PSP ∪ cleaned OmniPath/SIGNOR/iPTMnet. Every edge citing a compendium PMID or a high-throughput MS paper is dropped (list built by script, committed); reference-less edges only from whitelisted manually curated sources; MIMP/ProtMapper/text-mining-only edges dropped. kinaseSubMain "iKiP-DB" files dropped (mislabelled). Sugiyama 2019 only as a separate in-vitro library if provenance is clean. Edge PMIDs kept in `evidence`.
  - Predicted (leakage-audited, `docs/kinasesub_audit.md`): kinaseSubMain R2-nophospho, R3, Cantley atlas, as top-k sites per kinase with percentile weights.
  - Banned: kinaseSubMain release / `clean_bench`, SELPHI2, CoPheeKSA, PTMsigDB perturbation signatures, Hijazi downstream marker sets, Phosformer (not audited), any kinaseSubMain setting tuned on benchmarKIN.
  - **Methods may not use per-kinase features derived from GT labels** (e.g. target frequency); such information appears only in the prior-only control.
- **Drug–target profiles** (masking and dose-aware GT), precedence Klaeger 2017 > Reinecke 2024 > Davis 2011 > kinoplex > user Hijazi matrix; Klaeger protein groups expanded to all genes. Doses in `results/ground_truth/condition_doses.tsv` (Hijazi supplement / study methods; unknown → pre-registered default 10 µM).
- **External** (never used for development): benchmarKIN tyrosine set (18 exps), the 3 Hernandez studies not in the compendium, CPTAC tumour benchmark (7 cohorts). Secondary, stratified by drug seen/unseen in DEV. All provenance in DATA.md.

## 3. Ground truth (scripts/build_gt.py from rule TSVs; frozen with SHA-256 at G1)
GT unit (c, k ∈ U, s). Classes: direct positive, downstream positive (secondary only), masked (unknown), negative (all other U kinases in a GT condition).
**Precedence per (c,k): paradox-drop > direct > downstream > masked > negative.**
- **GT-H (Hijazi, 183 conds)**: direct = user kinhib matrix ∩ U, s = −1.
- **GT-D (decryptM)**: direct = curated primary-target list per drug at doses ≥ Kd_app (Kd known) — allosteric/PIKK/covalent inhibitors use published IC50 instead; target at sub-Kd dose → masked. Each dose is its own condition within the drug group. Only high/medium confidence pairs enter metrics.
- **GT-O**: other inhibitors; stimuli ≤4 h (later time points masked): receptor kinase + canonical near-certain first-tier kinases listed in `stimulus_kinases.tsv` are direct (s=+1), the rest downstream; genetic; cell-cycle. High/medium only; low/`inferred` = sensitivity analysis. LCK/TGFBR2 rules at medium.
- **Masks (all sources)**: downstream/feedback table by target class (MEK/ERK/RAF → ERK1/2, RSKs, MNKs, MSKs, RAFs; PI3K/PDPK1/AKT → AKT, SGK, S6K, MTOR, GSK3; MTOR → S6K, AKT, ULK1; p38 → MK2/3; canonical cascades for stimuli/cell-cycle); paradoxical RAF (cell-line BRAF/RAS genotype × RAF-inhibitor type; RAF positives dropped, MEK/ERK/RSK masked); off-targets with Kd_app ≤ 3× dose (any profile source) plus literature off-target table; negative paralogs sharing ≥80% PSP substrates with a positive.
- **n_targets** = # direct positives in U after the dose/Kd rule, before masking, label-only. Conditions with 0 direct positives are not in M1. **M1-sel** set: conditions with 1 ≤ n_targets ≤ 10.
- **Specificity controls**: lipid-kinase-only drugs (CAL101, TGX221, PIK294), Edelfosine, U73122 (PKC axis masked), `negative_control` rows.
- **Null controls**: the 16 decryptM "0 nM vs 0 nM" conditions follow their drug group: DEV-group nulls = N-fit (may be used to fit empirical nulls), TEST-group nulls = N-eval (calibration report only). Site-permutation nulls for development.
- **GT-A (descriptive)**: per-kinase Spearman(score_k, own activating-site log2FC; PSP Regulatory_sites "activity, induced") over non-inhibitor conditions, permutation null, macro over kinases; self-sites removed from the input vector for every method when computing GT-A.
- Tests before G1: every GT kinase ∈ U; no bad tokens; s ≠ 0 for direct/downstream; deterministic rebuild.

## 4. Splits and leakage control
- **Grouping**: union-find over all 1,937 conditions linking same canonical drug anywhere (treatment or control field; `drug_aliases.tsv`), same study for non-drug studies, Hijazi↔Wilkes by drug, duplicates. Largest group size logged.
- **Split selection (pre-registered, label-only)**: 100 candidate splits (seeds 0–99), TEST ≈30% of M1-sel pairs; choose the minimum of a distance on M1-sel pairs/kinases per source × harmony bin (defined in `make_split.py` docstring) with a penalty if <20% of TEST GT kinases are unseen in DEV.
- **DEV**: 5 grouped folds (seed 0) + DEV-confirm ≈15% of DEV groups.
- **Sealing**: TEST labels live only in `cache/sealed/` (hash in GT_FREEZE.json); `results/ground_truth/` holds DEV + DEV-confirm rows. Development agents use the DEV-only cache (`cache/dev/`, excluding TEST and DEV-confirm). The evaluator loads sealed labels only if FREEZE.json exists and its GT hash matches; FREEZE.json and the TEST run are written only by the PI-launched TEST script. pytest fails if a TEST condition ID appears in DEV artefacts.
- **Unsupervised use of compendium conditions** (co-regulation, factorisation, site variances): non-TEST conditions excluding those sharing a control with a TEST condition; inside CV also excluding the held-out fold.
- **Tiers**: T-seen / T-unseen (never a DEV positive) / T-family-unseen (no DEV positive in its Manning family).
- **benchmarKIN M3 before FREEZE**: only on DEV conditions + the tyrosine set (benchmarKIN contains Hijazi and 22 compendium Hernandez studies).
- Seeds `SEEDS=[0,1,2,3,4]` (config/seeds.json). Winner's-curse check before FREEZE: re-split DEV folds with seeds 1–4.

## 5. Evaluation (frozen harness `kinact.eval`)
t = s·score. Pair AUROC_ck = P(t_k > t_j) over negatives j in c; ties ½; unscored target → ½; unscored negative counts ½. Averaging: pairs → (kinase, GT group) → kinase → macro.
- **Frozen core (selection may use only these)**: **M1-sel (PRIMARY)**; M1-all (1/n_targets pair weights within group); M1-scored + pair coverage (reported with every M1); M2 (universe = kinases scored by both methods); harmony = M1-sel per bin (bins per (kinase, condition): 0, 1–4, 5–19, ≥20 quantified PSP substrates; `harmony_bins.tsv`) + bin-macro; by source (H/D/O), by sign, T-unseen; coverage-significant (# kinases beating a rank-permutation null, BH q ≤ 0.10); calibration (|z|>1.96 rate overall and per bin, on N-eval / site-permuted data); runtime.
- **Controls (all reported)**: random; prior-only = s_c·freq_k (freq from training-fold conditions; oracle sign — a strong control); shift-prior = freq_k·mean(x_c); site-permuted input; shuffled library (per-kinase counts preserved); sign-flip check; oracle.
- **Descriptive (never used for selection)**: M3 benchmarKIN-compatible; M4 kinase-centric within study; M1-family; downstream-positive AUROC; GT-A; kinase–kinase analysis (gold = kinase→kinase edges onto PSP regulatory sites of the substrate kinase; AUROC of across-compendium score correlation vs random pairs, stratified by substrate Jaccard).
- **Inference**: cluster bootstrap over GT groups stratified by source, B=2000, seed 0, shared resamples for all methods and endpoints, metric recomputed per resample; one-sided p. **TEST Holm family (6)**, all on M1-sel unless stated: final vs z-score, KSEA, JUMPsem, best-other-baseline (fixed at G2 from DEV), prior-only; final vs z-score on M2. "Beats" = Holm p < 0.05 AND ΔM1-sel ≥ 0.01 AND ΔM1-sel > 0 on T-unseen AND point ΔM2 ≥ 0. MDE pre-registered at G1 as 2.5·SE(ΔM1-sel) of z-score vs prior-only, SE estimated by bootstrap on DEV and scaled by sqrt(G_DEV/G_TEST) (groups), so TEST is never touched (C7). If a test fails: estimate + CI, no claim.

## 6. Baselines (shared library builder; each run with its published default library AND on the shared-library grid incl. our final library; the better DEV config goes to TEST)
Must: z-score (benchmarKIN/RoKAI form), KSEA (Casado mean-based z; KSEAapp R as reference), JUMPsem (R package 1.1.0 natively), decoupler 2.3.1 ulm / mlm / wmean / norm_wmean / viper (used directly), GSEA (decoupler/gseapy), RoKAI network propagation, Kinase Library motif enrichment / PhosX-like motif-only. Python reimplementations, if needed for speed, must match the reference at Spearman ≥0.99. Could: KSTAR, KARP, PhosR, IKAP. Cut: PTM-SEA, separate Wilcoxon/KS, INKA.

## 7. Method development (DEV only; build order with time boxes in RESEARCH_TREE.md)
1. **Library**: variants scored with z-score and mlm (PSP, curated union, predicted top-k with percentile^γ, mixtures).
2. **Bayesian core**: x_c = W·a_c + b_c + ε; per-site noise variances estimated label-free on allowed conditions; prior a ~ N(0, τ²K), K = family-tree / motif-similarity / substrate-overlap kernel; closed-form posterior mean and z; ≤5 hyperparameters (nested CV).
3. **Calibration**: effective-m analytic null and/or per-condition empirical null (fit on N-fit / site permutations only).
4. **Learned parts**: no per-edge or per-kinase free weights, no GT-derived per-kinase features. Edge weights only as a shared function of edge/site features (≤20 params); >50 params needs a leave-kinase-out non-inferiority check first. Data-driven substrate expansion via co-regulation/factorisation on allowed conditions only.
5. **Kinase self-sites** as auxiliary signal (separately; paradoxical behaviour under inhibitors).
6. **Stacking** (learning-to-rank, pairwise within condition, ≤30 params, no kinase-specific params) only if complementary methods exist.
**Kill rule**: no DEV gain by a branch's time box → [-] abandoned.
**Adoption rule**: on pooled DEV out-of-fold predictions, the one hard test is one-sided paired group-bootstrap p < 0.05 for ΔM1-sel > 0. Point non-inferiority (margin −0.02): M1-all, each source, each harmony bin, M2 vs current best; T-unseen margin −0.01. Any ΔM1 ≥ 0.02 is audited by a non-producing agent before adoption.
**DEV-confirm**: the cumulative adopted bundle is evaluated on DEV-confirm at most 5 times (logged), at milestones and at G2. A failure (Δ vs z-score ≤ 0) triggers a review of the last adopted ideas, not tuning on DEV-confirm.
**Ledger**: every DEV evaluation appended to `results/dev_ledger.tsv` (config hash, commit, metrics, n configs per branch).

## 8. Model selection and TEST
G2: final config (commit + config hash + GT hash) chosen on DEV CV + DEV-confirm, written to FREEZE.json by the PI. TEST run once; independent agent reproduces it; report DEV→TEST drop. Post-TEST explorations are labelled and cannot change the final.

## 9. Verification
pytest: toy data with known answers for every metric and method; sign-flip, permutation and shuffled-library tests; leakage tests (drug/study/control disjointness, DEV-only cache, no library edge cites a compendium PMID); GT/evaluator hash check at each agent start; deterministic seeds. REPORT numbers generated by script from committed result files; LITERATURE.md PMIDs/DOIs checked by script.

## 10. Compute
64 CPU, 375 GB RAM, 4× RTX 4090. Sparse z-score over the compendium ≈0.09 s; Bayesian core ≈7 ms/condition. GPUs only if needed.

## 11. Agent organisation & guards
PI designs, decides, integrates, keeps RESEARCH_TREE.md, writes FREEZE.json. Agents: GT+split builder, harness builder, library builder, baselines, method explorers (one branch each), independent auditors. Agents use the DEV-only cache, never touch sealed data, cite committed files, return ≤300-word summaries.

## 12. Amendment history
- v0→v1 (round 1, 6 reviewers): A1 sign fixed (t = s·score; v0 had −s). A2 unscored → ½. A3 averaging order. A4 primary metric M1-sel (≤10 targets; PI decision over 1/n weighting — simpler, most trustworthy labels; M1-all keeps coverage view). A5 U = 503 protein kinases. A6 downstream/paradoxical/off-target/paralog masks. A7 drug-linked grouping, pre-registered split selection, minimal QC. A8 significance-based coverage, bootstrap/Holm/min-gain, ledger, DEV-confirm, prior-only control. A9 runtime targets, Bayesian core, no per-edge/per-kinase weights. A10 extra KS sources, drug profiles, GT-A, CLI, baseline tiers.
- v1→v2 (round 2, 4 reviewers): B1 adoption rule (one hard test + non-inferiority). B2 build order + kill rule. B3 prior-only in Holm family; shift-prior. B4 split on M1-sel pairs; MDE. B5 n_targets. B6 reporting/bootstrap strata. B7 frozen core vs descriptive. B8 KS-source leakage filters. B9 null controls follow groups. B10 sealing. B11 doses, literature off-targets, sub-Kd masking. B12 GT rule details. B13 GT-A. B14 baselines' default libraries. B15 kinase–kinase analysis. B16 families/aliases. B17 Phosformer cut.
- v2→v3 (round 3, final reviewer): C1 merged text (16 nulls; 6-test family; M4 descriptive; iKiP-DB dropped; sealed TEST labels). C2 precedence paradox-drop > direct > downstream > masked > negative. C3 stimulus first-tier kinases are direct via `stimulus_kinases.tsv` (decides most s=+1 pairs). C4 benchmarKIN M3 only on DEV + tyrosine before FREEZE. C5 ban GT-derived per-kinase features in methods. C6 unscored = NaN; shuffled-library control; decoupler used directly; runtime reported not gated; DEV-confirm semantics.
- C7 (pre-G1, PI, before any result seen): MDE estimated on DEV and scaled to TEST size instead of bootstrapping on TEST.
- C8 (pre-G1, PI; only controls seen, no method comparisons): the bootstrap over GT groups changes the kinase composition of the macro average (kinases in a single group drop out of ~37% of resamples), biasing replicates (docs/audit/harness_verification.md). Inference therefore uses bias-corrected forms: CI = basic (reverse-percentile) interval [2Δ̂ − q97.5, 2Δ̂ − q2.5]; one-sided p = fraction of (Δ* − mean(Δ*)) ≥ Δ̂ (shifted-null bootstrap). MDE recomputed with the bias-corrected SE. The ½ rule for unscored kinases is confirmed as intended (coverage is part of the objective); M2 ≥ 0 in "beats" protects algorithm quality.
- C9 (post-G1 factual correction, no result seen): the benchmarKIN CPTAC tumour gold standard is the top/bottom 5% of kinase PROTEIN abundance per sample (not activating-site based). External sets prepared (docs/external_eval.md): Tyr 18 exps/16 kinases, Hernandez non-compendium 15 exps/7 kinases (PMIDs 23592840, 23612710, 24204825), CPTAC 7 cohorts/367 kinases. Secondary/descriptive; paired bootstrap per set added at TEST time.
