# RESEARCH TREE — Kinase activity inference from phosphoproteomics

Legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding

## Status summary (2026-10-27, M1: G1 frozen)
Design closed after 3 review rounds (DESIGN v3 + C7/C8). G1 frozen (results/G1.json, commit 7de8eef): GT built from rule tables and independently audited (12 fixes), split seed 90 (TEST sealed in cache/sealed: 468 conds, 470 M1-sel pairs, 67 kinases, 15 unseen), harness verified independently (from-scratch reimplementation identical; benchmarKIN reference 0.775 reproduced on their data), bias-corrected bootstrap (C8). MDE on TEST ≈ 0.06 M1-sel. 96 tests pass. Libraries built (PSP variants, leakage-filtered curated union, Sugiyama, predicted top-500 for 503 kinases).

Key finding for strategy: under the pre-registered coverage-aware primary metric (unscored = ½), PSP z-score scores M1-sel 0.515 (scored-only 0.633, coverage 0.58) — barely above random (Δ +0.010, n.s.) — while the popularity prior-only control scores 0.687. Winning therefore requires accurate scores for MANY kinases (incl. understudied ones), i.e. extended/predicted libraries + a statistic that ranks well across kinases with very different substrate counts. Next: baselines finishing; then method branch 1 (library variants) and 2 (Bayesian core).

## Decision log
- D7: core-v1 ships z_vif (calibrated) not z_post: −0.007 M1-sel is within the −0.02 non-inferiority margin; calibration is a design desideratum; DEV-confirm look 1 pre-declared z_vif primary. Reference number for further branches = 0.663.
- D6: Adopted bundle 'core-v1' = curated_union + r2np@100 (α=1), Bayesian core, family kernel ρ_f=3, ρ_g nested, τ=0.125; no winsor (failed adoption rule).
- D5: Do NOT add any kinase-popularity prior (literature counts etc.) to close the gap to prior-only: it exploits annotation/drug-target bias the brief warns against. Prior-only remains a reported control; T-unseen and harmony carry the evidence.
- D4: G1 frozen at commit 7de8eef. MDE ≈0.06.
- D3: ½ rule confirmed (coverage is part of the objective; M2 ≥ 0 required in 'beats'). Bootstrap bias fixed by C8 before any method comparison.
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
- [x] 3. G1 build — FROZEN (results/G1.json)
  - [x] 3.1 GT + split + DEV-only cache: 1,047 GT conds; M1-sel 1,660 pairs / 140 kinases; split seed 79 (TEST 396 conds, 445 M1-sel pairs, 78 kinases, 18 unseen). Hijazi doses unobtainable → 10 µM default. 16 tests pass.
    - [x] 3.1b GT audit: 12 issues (4 high) fixed, sample pass 89.7% pre-fix; combination doses ambiguous → low; cellular-potency override (gefitinib/EGFR, dasatinib/ABL1). FINAL split seed 90: TEST 468 conds / 470 M1-sel pairs / 67 kinases (15 unseen); DEV-CV 877 pairs; devconfirm 119. GT_FREEZE sha 1e6f791c…
  - [x] 3.2 harness src/kinact/eval.py + tests (20 pass; paired bootstrap B=2000 0.3 s)
    - [x] 3.2b integration: loaders, smoke run. DEV: zscore M1-sel 0.521 (scored 0.666, cov 0.63); prior-only 0.706; random 0.476; MDE 0.036 (M1-sel groups) / 0.059 (all groups)
    - [x] 3.2c harness verification: no code bug; gap 0.78→0.52 = data processing −0.06, broader GT −0.09, per-exp/macro +0.04, ½ rule −0.15. Bootstrap composition bias → C8. Random unbiased (0.495±0.018). Prior-only strictly OOF.
    - [!] PSP z-score ≈ random on M1-sel (coverage-limited); prior-only 0.687 is the bar.
  - [x] 3.3 KS libraries: psp_all 401 kin (225 with ≥5 quantified), curated_union 439 kin (296 ≥5), Sugiyama in vitro 347 kin (342 ≥5), predicted top-500 r2np/r3 (503), cantley (396). 405 PMIDs excluded (compendium + HT-MS). 27 tests.
- [x] 4. Baselines on DEV (results/baselines/baselines_dev.md; refs: decoupler/KSEAapp Spearman ≈1.0)
  - best per method M1-sel: zscore curated min1 0.563 · KSEA 0.562 · PhosX-like 0.557 · KinLib enrich 0.550 · motif z 0.541 · viper/ulm/wmean/gsea ~0.53 · RoKAI 0.526 · mlm 0.522 · JUMPsem 0.507 (coverage 0.49; 6–12 s/cond).
  - [!] All published methods ≈ chance on coverage-aware M1-sel; library/min-substrate drive gains. Candidate 'best-other-baseline' (excl. z/KSEA/JUMPsem) for G2: PhosX-like (0.557).
- [~] 5. Method development (DEV only; adoption per DESIGN §7)
  - [x] 5.1 B1 library variants (408 configs; results/b1/b1_report.md)
    - [!] z-score + r2np top-200: M1-sel 0.664 (nested 0.638; optimism +0.026), coverage ≈1.0, M2 vs PSP-z 0.699 vs 0.633; Δ vs PSP-z +0.12 [0.07,0.16]. Still < prior-only 0.687 (n.s.). r2np > r3 > cantley > curated; Hijazi weak (0.57), decryptM 0.89 (suspicious).
    - [x] 5.1a r2np/r3/cantley leakage audit: SAFE (docs/audit/r2np_leakage.md). decryptM high = dasatinib/Tyr composition (residue-matched random library already 0.70 on D) → quote D only with caveat; residue-matched null added as descriptive control.
  - [x] 5.2 B2 Bayesian core (results/b2/b2_report.md): curated_union + r2np@100 (α=1), family kernel ρ_f=3, ρ_g=0.25, τ=0.125, no noise model. Nested M1-sel 0.655 (H .644 D .768 O .674), M2 +0.033 vs PSP-z, T-unseen 0.589 (prior-only 0.324). Ablations: −K 0.617, −predicted 0.607, neither 0.553; noise model hurts; posterior-mean ranking 0.672 (post hoc). Calibration 4.6% overall but grows with substrate count (1.2%→7.6%). 8.9 ms/cond CPU-1T.
    - [!] Complementary to prior-only (prior wins H; B2 wins D and unseen kinases).
  - [x] 5.3 B2.1+B3 (results/b21/b21_report.md)
    - [!] ADOPTED: ρ_f ≤ 3 → nested 0.670 with z_post ranking (Δ +0.015 [0.006,0.035], p=0.001; optimistic — cap decided post hoc). z_post is ~10× conservative (0.6% |z|>1.96).
    - Output = z_vif no winsor (D7): DEV nested M1-sel 0.663 (M1-all 0.615, M2 0.687 vs PSP-z 0.633, harmony .581/.688/.737/.645, T-unseen 0.647); calibration 6.4% DEV nulls / 5.4% site-permuted. (Earlier '5.0%, no M1 cost' figures were the non-adopted winsor variant — corrected.)
    - NOT adopted: winsorising (+0.004, threshold chosen post hoc), posterior-mean output (not robust), r3/cantley edges (flat), similarity kernel (flat). Bigger search spaces lose ~0.03 to selection noise → keep final search small.
    - [x] 5.3a audit (docs/audit/core_v1_audit.md): reproduced exactly, no leakage. DEV-confirm look #1 (179 M1-sel pairs): core-v1 0.803 vs PSP-z 0.545 (+0.26 [0.19,0.31]), curated-z 0.618 (+0.19), prior-only 0.856 (−0.05, n.s.). M2 gain vs PSP-z only +0.008 on devconfirm. 4 looks left.
  - [-] 5.4 B4 co-regulation expansion (results/b4/b4_report.md): on B2 core ΔM1-sel −0.010 [−0.021,−0.003] (expanded sites shared across kinases → noise); neighbour fill-in n.s.; on weighted z +0.026 (only reaches B1 level). Kept as stacking candidate (T-unseen 0.684). NOTE for G2: exclude DEV conditions sharing a control with TEST conditions from any co-regulation fit (metadata only).
  - [ ] 5.4b shared-function edge weights
  - [-] 5.5 B5 self-sites (results/b5/b5_report.md): alone 0.726 on pairs with measured self-site (35% of M1-sel pairs); paradoxical up for PKC/PKA/ERKi/TAK1; on core +0.003 [−0.001,0.008] n.s. (half of self-sites already library edges). GT-A helper mask_self_sites.
  - [ ] 5.6 B6 stacking
- [ ] 6. Sealed test + REPORT.md
