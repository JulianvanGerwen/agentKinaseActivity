# RESEARCH TREE — Kinase activity inference from phosphoproteomics

Legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding

## Status summary (2026-10-28, M2: G2 frozen, sealed TEST run done)
The sealed one-shot TEST run is complete (FREEZE.json + results/test/, commit be2ef15; freeze checks passed; independently reproduced bit-identically from commit 03ab343, runtime fields aside). On 468 TEST conditions (245 evaluable, 470 M1-sel pairs, 67 kinases, 15 unseen), core-v3 scores M1-sel 0.714 vs z-score 0.584, KSEA 0.582, PhosX-like 0.549, JUMPsem 0.505; all four "beats" tests pass (Holm p 0.003; Δ +0.13..+0.21; T-unseen Δ > 0; M2 ≥ 0). Harmony 0.725 (bins 0.72/0.66/0.70/0.82). Runtime 6.7 ms/cond (1 thread), 2.0 ms (16 threads), 1.7 ms (GPU); full compendium 2.9 s.

Not beaten: the prior-only control (0.712; Δ +0.002, p 0.50). It wins on well-studied kinases (T-seen 0.844 vs 0.759) but is far below chance on unseen kinases (0.253 vs core-v3 0.557). Calibration is 4.8% on site-permuted data but 8.1% on real null conditions. core-v3 vs core-v1 +0.021, carried by signature kinases. On the external Hernández set, core-v3 is below core-v1 (0.552 vs 0.613). Retro done and REPORT.md written. Remaining ideas need a fresh independent test set (TEST is spent); see REPORT §6.

## Decision log
- D14 (G2/TEST): FREEZE.json written at 03ab343 and TEST run exactly once; results committed with the freeze file (be2ef15) and reproduced by an independent agent. Claims are restricted to the pre-registered Holm family: core-v3 beats z-score, KSEA, JUMPsem and PhosX-like; there is no claim vs prior-only, and the M2 vs z-score test is not significant (p 0.15 after Holm). No further method changes after TEST.
- D13: Final TEST model = core-v3 with signature settings (n_top, κ)=(100,0) exactly as frozen for DEV-confirm look #2 (results/g2_prep/devconfirm_look2_frozen.json); build verifies ρ_sig/E1/counts reproduce. A1 secondary coreg config = majority of per-fold DEV picks. Comparators fixed: z-score curated min1 (+ PSP min3 reported), KSEA curated min1, JUMPsem best DEV, best-other = PhosX-like, prior-only.
- D12: Final box B8 (results/b8/b8_report.md): core-v3 = core-v2 + E1 (per-signature-kinase permutation-null scale of the fused z; signature-kinase |z|>1.96 on permuted data 7.0%→5.1%, ΔM1-sel −0.0007). E2 shared-factor removal (Δ −0.0016, p .80) and E3 top-3 bagging (Δ −0.0012, p .73) rejected. Development stopped; next: seeds 1–4 re-split, DEV-confirm look #2 on core-v3, G2.
- D11: core-v2 = core-v1 + late-fusion signatures (B7b, α=1, signatures fitted on DEV-CV conditions only). Conditional on DEV-confirm look #2; revert to core-v1 if Δ ≤ 0 there.
- D10: Near diminishing returns (fresh-reviewer verdict): remaining ideas ≤ +0.01 ≪ TEST MDE 0.06. Plan: B7b audit → one final idea box → seeds 1–4 re-split check → DEV-confirm look → G2 → TEST.
- D9: Near-misses (A1 p=0.057) are NOT adopted by relaxing the rule post hoc; instead pre-registered as secondary TEST analyses.
- D8: Amendment C10 (pre-registered) allows per-kinase perturbation signatures in B7 with strict safeguards; rationale: data-driven KS evidence is in scope (brief: 'other sources of KS relationships'), popularity priors remain banned.
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
  - [-] 5.4b B4b learned shared edge weights (results/b4b/b4b_report.md): nested picks θ=0 every fold; 1/promiscuity^β +0.007 [−0.007,0.019] p=0.16 (M1-all +0.019 sig.) but calibration worse (7.6%). Only stable learned signal: down-weight promiscuous sites.
  - [-] 5.5 B5 self-sites (results/b5/b5_report.md): alone 0.726 on pairs with measured self-site (35% of M1-sel pairs); paradoxical up for PKC/PKA/ERKi/TAK1; on core +0.003 [−0.001,0.008] n.s. (half of self-sites already library edges). GT-A helper mask_self_sites.
  - [-] 5.6 B6 stacking (results/b6/b6_report.md): only B4 coreg-z is complementary (rank corr 0.56; better in bin 0 and T-unseen). A1 = average(core z_vif, coreg z): +0.012 [−0.001,0.030] p=0.057 → NOT adopted (hard test); M1-all +0.020, harmony +0.015, T-unseen +0.036 sig., calibration 4.6%. Learned stacker −0.001. → A1 pre-registered as SECONDARY method for TEST (descriptive).
  - [-] 5.7 B7 signatures in joint posterior (results/b7/b7_report.md): +0.012 [−0.017,0.026] p=0.15; kinases with signature +0.076 (p=.003), Hijazi +0.036; T-unseen −0.026 (TYK2 via JAK-kernel borrowing) → fails C10.
    - [!] 5.7b B7b late fusion (pre-declared; results/b7b/b7b_report.md): M1-sel 0.663→0.680, Δ +0.018 [−0.000,0.023] p=0.0085; T-unseen +0.0006; all bins +; shuffled control −0.002..−0.008; calibration 6.2%/5.5%. Signature kinases +0.102. α at grid edge (1).
      - [x] audit (docs/audit/b7b_audit.md): reproduces byte-identical; no leakage (label randomisation invariant); α plateau 0.75–2; all 5 folds gain. Caveats: 71% of gain from 5 kinases (MAP3K7, ROCK1/2, MET, IKBKB); signature-kinase FP 8.1→9.8%; z over-dispersed for signature kinases (7.0%). → ADOPTED as core-v2, conditional on DEV-confirm.
  - [x] 5.8 B8 final box (results/b8/b8_report.md; pre-declared): E1 fused-z permutation-null recalibration ADOPTED (calibration fix; M1-sel 0.6796, sig-kinase calibration 5.1%); [-] E2 top-r factor removal −0.0016 p .80; [-] E3 top-3 core bagging −0.0012 p .73 (top configs r≥0.96). → core-v3. Development stopped.
- [x] 7. Deliverable core-v1: src/kinact/core.py + CLI `kinact score`; equivalence max|Δz| 8e-12; 4.6 ms/cond (CPU 1T, compendium avg), 2.1 ms (16T/GPU); DEV compendium 2.2–6.2 s. PSP-derived edges git-ignored, rebuilt by scripts/build_core_artifacts.py.
- [x] 9. Descriptive DEV analyses (results/descriptive/descriptive_dev.md; not used for selection)
  - M3 benchmarKIN-style (matched cells): core-v2 0.740 > core-v1 0.714 > curated-z 0.682 > PSP-z 0.668 > prior-only 0.598.
  - GT-A self-sites: cores most kinase-specific (80/145 kinases sig.). Kinase–kinase: AUROC 0.66–0.68 but 0.53–0.60 within substrate-overlap strata (mostly shared-substrate artefact) — honest limitation.
  - Residue-matched random library: M1-sel 0.529 (real 0.668).
- [x] 8. External test pipelines ready (docs/external_eval.md; run_external(..., confirm_test=True)); Tyr 18 exps/16 kin, Hernandez 15/7, CPTAC 7 cohorts/367 kin (gold = protein abundance, C9). No real method run.
- [x] 10. Pre-G2 checks (results/g2_prep/g2_prep.md)
  - Re-split seeds 1–4: core-v3 0.6800±0.0037 (seed 0: 0.6796) — no fold-selection inflation (caveat: partitions overlap seed 0, ARI 0.68–0.75). Δ vs core-v1 +0.011..+0.019; vs prior-only −0.006±0.009.
  - [!] DEV-confirm look #2 (179 pairs): core-v3 0.842 vs core-v1 0.803 (+0.039 [0.006,0.069], p .016) → D11 condition PASSED, core-v3 kept; vs PSP-z +0.298, KSEA +0.225, prior-only −0.014 (n.s.). 3 looks left.
- [x] 6. Sealed TEST (results/test/LEDGER.md; FREEZE.json; reproduced)
  - [!] core-v3 M1-sel 0.714: beats z-score (+0.130 [0.090,0.169]), KSEA (+0.132), JUMPsem (+0.209), PhosX-like (+0.165); all Holm p 0.003.
  - [!] vs prior-only +0.002 [−0.036,0.083] p 0.50 (no claim); T-unseen 0.557 vs 0.253; T-seen 0.759 vs 0.844.
  - core-v3 vs core-v1 +0.021 [0.001,0.031] (signature kinases +0.050); A1 secondary 0.686 (not better). Calibration: nulls 8.1%, site-permuted 4.8%. External: Tyr 0.737, Hernández 0.552 (core-v1 0.613), CPTAC 0.672.
  - [x] 6.1 LESSONS retro (docs/retro/m2_retro.md; LESSONS pruned to 12 rules) · [x] 6.2 REPORT.md
