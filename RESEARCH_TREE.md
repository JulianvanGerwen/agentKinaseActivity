# RESEARCH TREE — Kinase activity inference from phosphoproteomics

Legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding

## Status summary (2026-10-27, M0: kickoff)
Project initialised. Hardware: 64 CPU, 375 GB RAM, 4× RTX 4090; Python env `/home/plm_env/bin/python` (torch 2.6+cu124, pandas, sklearn, statsmodels; no lightgbm/decoupler yet). Internet available. Provided data in `data/` (compendium of log2FC: 1.09 GB long table; PSP KS/site/regulatory dumps; Hijazi inhibitor→kinase target matrix; Klaeger tables; papers Hijazi, KSEA, JUMPsem). The predicted-KS project lives at `/cluster/project/beltrao/JvG/Projects/hackathonKinaseSubs/kinaseSubMain` (NOT under /workspace).

Key early observation: kinaseSubMain's deployed library uses (F6b) ridge deconvolution of inhibitor conditions incl. Hijazi and (F6a) co-regulation over the compendium → leaks into our evaluation. They also built a "clean R2" library (benchmark PMIDs + CPTAC removed) and a benchmarKIN harness reproducing Müller-Dott 2025 (PSP+z-score perturbation AUROC 0.782). Phase 1 = literature + data profiling + audit of kinaseSubMain reuse/leakage, in parallel.

## Decision log
- D2: Stimulus first-tier kinases (e.g. EGF→ERK1/2) are DIRECT positives via a pre-registered table — they are the activity changes users most need detected; deeper cascades = downstream (secondary).
- D1: Primary endpoint M1-sel (≤10 direct targets) rather than 1/n weighting — simpler, most trustworthy labels; M1-all + harmony carry the coverage objective.
- D0 (M0): Work in /workspace (git). `data/`, `external/`, `cache/` git-ignored (PSP is CC BY-NC-SA; compendium large). Python env = /home/plm_env.

## Tree
- [~] 1. Phase 1 — groundwork
  - [x] 1.1 LITERATURE.md (594 lines; [V]/[M] verification tags)
    - [!] benchmarKIN (Müller-Dott 2025): z-score best on perturbation (AUROC ~0.79); library choice > algorithm; naive top-k predicted KS hurts (~0.65). JUMPsem = multi-sample one-factor SEM per kinase (not single-sample).
    - Hijazi targets profiled at 200 nM but dosed 1–10 µM → off-targets under-annotated; decryptM dose-resolved + Klaeger Kd → dose-aware GT.
  - [x] 1.2 Data profile → docs/data_profile*.md, DATA.md, scripts/prep_data.py, results/ground_truth/candidate_gt.tsv
    - 1,937 conds / 178 studies / 190k sites; 1,020 candidate GT conds (Hijazi 177 usable, decryptM 340); 52% of PSP KS sites quantified.
    - [!] Only ~1/3 of Hijazi target kinases (98/306) evaluable by PSP methods (≥3 substrates/cond); top-10 kinases = 37% of annotations.
  - [x] 1.3 kinaseSubMain audit → docs/kinasesub_audit.md
    - [!] Their data/01-03 + kinhib file are byte-identical to ours. SAFE predicted-KS: R2 'nophospho' (DEV-trained, s0, full 350k×503 matrix) + R3 (no phospho feats) + Cantley atlas. UNSAFE: deployment release & 'clean_bench' (still uses 37 compendium studies incl. decryptM via F6a/F6b).
    - [!] benchmarKIN available locally (230 perturbation exps: Hernandez 94, Hijazi 118, Tyr 18; 7 CPTAC cohorts) + validated harness. Independent of compendium: Tyr set, 3 Hernandez studies, CPTAC.
    - Reusable: site universe mapping, 692 inhibitor conditions with occupancy-based targets (Klaeger/Reinecke/Davis), kinase families/paralogs, ESM kinase embeddings.
- [x] 2. DESIGN.md + 3 review rounds (v3 final)
  - [x] 2.1 v0 → round 1 (6 reviewers, all 'major'): sign error in M1 caught; masks for downstream/paradoxical; M1-sel primary; Bayesian core; ban per-kinase weights; prior-only control → v1 (A1–A10)
  - [x] 2.2 round 2 (pi, stats, priorwork, data) → v2 (§13 B1–B17)
    - [!] prior-only control (kinase target frequency, oracle sign) ≈0.63 > PSP z-score M1 ≈0.51–0.55 → M1 is coverage-dominated; prior-only added to TEST Holm family.
    - [!] OmniPath/iPTMnet edges cite compendium PMIDs (leakage); kinaseSubMain 'iKiP-DB' mislabelled → filtered/dropped.
  - [x] 2.3 reference tools: R 4.4.3 + JUMPsem 1.1.0 + KSEAapp 2.0 (external/envs/r), decoupler 2.3.1, gseapy (docs/reference_tools.md)
  - [x] 2.4 round 3 (final reviewer) → DESIGN v3 merged (C1–C6). Design phase closed.
- [~] 3. G1 build
  - [~] 3.1 GT + split + DEV-only cache (scripts/build_gt.py, make_split.py)
  - [~] 3.2 harness src/kinact/eval.py + tests
  - [~] 3.3 KS libraries with leakage filters (scripts/build_libraries.py)
- [ ] 4. Baselines (z-score, KSEA, JUMPsem, …)
- [ ] 5. Method development
- [ ] 6. Sealed test + REPORT.md
