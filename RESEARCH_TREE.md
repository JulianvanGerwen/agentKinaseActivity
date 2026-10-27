# RESEARCH TREE — Kinase activity inference from phosphoproteomics

Legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding

## Status summary (2026-10-27, M0: kickoff)
Project initialised. Hardware: 64 CPU, 375 GB RAM, 4× RTX 4090; Python env `/home/plm_env/bin/python` (torch 2.6+cu124, pandas, sklearn, statsmodels; no lightgbm/decoupler yet). Internet available. Provided data in `data/` (compendium of log2FC: 1.09 GB long table; PSP KS/site/regulatory dumps; Hijazi inhibitor→kinase target matrix; Klaeger tables; papers Hijazi, KSEA, JUMPsem). The predicted-KS project lives at `/cluster/project/beltrao/JvG/Projects/hackathonKinaseSubs/kinaseSubMain` (NOT under /workspace).

Key early observation: kinaseSubMain's deployed library uses (F6b) ridge deconvolution of inhibitor conditions incl. Hijazi and (F6a) co-regulation over the compendium → leaks into our evaluation. They also built a "clean R2" library (benchmark PMIDs + CPTAC removed) and a benchmarKIN harness reproducing Müller-Dott 2025 (PSP+z-score perturbation AUROC 0.782). Phase 1 = literature + data profiling + audit of kinaseSubMain reuse/leakage, in parallel.

## Decision log
- D0 (M0): Work in /workspace (git). `data/`, `external/`, `cache/` git-ignored (PSP is CC BY-NC-SA; compendium large). Python env = /home/plm_env.

## Tree
- [~] 1. Phase 1 — groundwork
  - [x] 1.1 LITERATURE.md (594 lines; [V]/[M] verification tags)
    - [!] benchmarKIN (Müller-Dott 2025): z-score best on perturbation (AUROC ~0.79); library choice > algorithm; naive top-k predicted KS hurts (~0.65). JUMPsem = multi-sample one-factor SEM per kinase (not single-sample).
    - Hijazi targets profiled at 200 nM but dosed 1–10 µM → off-targets under-annotated; decryptM dose-resolved + Klaeger Kd → dose-aware GT.
  - [~] 1.2 Data profile: compendium, ground-truth candidates, PSP coverage → docs/data_profile.md, DATA.md
  - [x] 1.3 kinaseSubMain audit → docs/kinasesub_audit.md
    - [!] Their data/01-03 + kinhib file are byte-identical to ours. SAFE predicted-KS: R2 'nophospho' (DEV-trained, s0, full 350k×503 matrix) + R3 (no phospho feats) + Cantley atlas. UNSAFE: deployment release & 'clean_bench' (still uses 37 compendium studies incl. decryptM via F6a/F6b).
    - [!] benchmarKIN available locally (230 perturbation exps: Hernandez 94, Hijazi 118, Tyr 18; 7 CPTAC cohorts) + validated harness. Independent of compendium: Tyr set, 3 Hernandez studies, CPTAC.
    - Reusable: site universe mapping, 692 inhibitor conditions with occupancy-based targets (Klaeger/Reinecke/Davis), kinase families/paralogs, ESM kinase embeddings.
- [ ] 2. DESIGN.md + ≤3 review rounds
- [ ] 3. Evaluation harness + verification (toy data, controls)
- [ ] 4. Baselines (z-score, KSEA, JUMPsem, …)
- [ ] 5. Method development
- [ ] 6. Sealed test + REPORT.md
