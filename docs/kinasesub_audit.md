# Audit: kinaseSubMain (predicted kinase–substrate sites) for reuse

Sister repo (READ-ONLY): `/cluster/project/beltrao/JvG/Projects/hackathonKinaseSubs/kinaseSubMain` (below: `KSM/`).
Audited 2026-10-27 from REPORT.md v3, DATA.md, docs/{downstream,features_phospho,features_motif,model_r2,model_r3}.md,
run manifests and parquet metadata (read with `KSM/.venv/bin/python -I`). Nothing in KSM was modified.

**Key fact:** KSM's `data/01–03_*.csv` and `kinhib_inhib_comb__1.csv` are byte-identical to ours (`cmp`). Any KSM
feature that uses "the compendium" uses our evaluation data.

## 1. KS prediction outputs

Index: `kinase_idx` → `KSM/cache/index/kinases.parquet` (503 rows: uniprot_acc, gene, residue_class, Manning
group/family/subfamily, pseudokinase, KD bounds). `site_idx` → `KSM/cache/index/sites.parquet` (350,736 U_obs sites:
site_id `ACC_S123`, acc, res, pos, gene, window15, window31, in_U_ms; 347,354 in U_ms). Every full file has
**131,144,410 rows = all U_obs sites × the 503 kinases of matching residue class**, so full score matrices exist.

| Variant | File | Size | Columns | Training labels | Phospho features |
|---|---|---|---|---|---|
| Deployment final (R2 5 seeds + R3 5 seeds, fused) | `KSM/cache/release/r2r3_release__20261027-012616__s0/release_scores.parquet` (sha cceb32c7…) | 4.1 GB | kinase_idx, site_idx, score (fused z-logit), score_site, prob, pct, assignment_prob | DEV+TEST PSP (L2, 15,202 pairs) | **full F6 (incl. Hijazi, CPTAC)** |
| Deployment R2 only (5 seeds) | `KSM/cache/release/r2_deploy__20261027-010348__s0/r2_uobs_scores.parquet` | 1.2 GB | kinase_idx, site_idx, score | DEV+TEST | full F6 |
| Deployment R3 only (5 seeds) | `KSM/cache/release/r3_deploy_20261027-005816/uobs_predictions.parquet` (sha 2f62988c…; README.json) | 279 MB | kinase_idx, site_idx, score | DEV+TEST (14,747 pos) | **none** (ESM-2 650M windows + protein token + U_ms "phospho-channel" bit + Cantley distillation) |
| R2 clean_bench, DEV-only, seed 0 | `KSM/cache/runs/r2_final_cleanbench_uobs__20261026-204513__s0/uobs_predictions.parquet` (+ `_topk.parquet`, 2 MB) | 1.2 GB | kinase_idx, site_idx, score | DEV only | F6 minus 22 benchmarKIN PMIDs (incl. Hijazi) and minus all CPTAC |
| R2 nophospho, DEV-only, seed 0 | `KSM/cache/runs/r2_final_nophospho_uobs__20261026-213950__s0/uobs_predictions.parquet` (sha 1edbd269…; + `_topk.parquet` 129,302 rows, with kinase_acc/site_id) | 1.2 GB | kinase_idx, site_idx, score | DEV only | **none** (groups F1,F0,F2a,F3,F4,F5,F7) |
| Release top-k (from deployment final) | `KSM/results/release/kinase_top200.tsv.gz` (100,600 rows), `site_top5.tsv.gz` (1,753,680 rows) | 1.9 / 34.8 MB | kinase_acc, gene, rank, site_id, score, prob, pct, assignment_prob | DEV+TEST | full F6 |
| Cantley / Kinase Library atlas (B1) | `KSM/cache/features/motif/cantley_{unprimed,primed_soft,primed_hard}_{score,pct_published,pct_universe}.npy` (+ `cantley_imputed_*`, `cantley_meta.json`); load via `ksp.store.load_matrix("motif/<stem>")` | 12 GB dir | dense float32 503 × 350,736 | none (published atlas, kinase-library 1.8.0 reproduced within 1e-3) | none (but see §2: pct_universe/primed use U_ms) |

Score semantics (REPORT §4.5, model_r2 §7): `score` = raw logit (R2 LightGBM / R3) or fused z-scored logit.
`prob` = per-kinase cross-fitted tiered Platt calibration (fitted on DEV OOF; PU lower bound; not re-validated for the
DEV+TEST model). `pct` = within-kinase percentile. `assignment_prob` = softmax over kinases (T≈0.93) + ST/Y isotonic
map: no ≥0.8 tier; Y sites capped at 0.445. Only the deployment release has calibrated columns; the DEV-only runs have
raw `score` only. Unprimed Cantley has NaN for kinases without own matrix (nan_fraction 0.44 incl. wrong-class cells);
`cantley_imputed_*` fill 503 kinases by nearest atlas kinase.

Other caches: `KSM/cache/downstream/r2_{clean,nophospho}/{uobs_superset,measured_topk}.parquet` (top-k libraries);
TEST-only predictions under `KSM/cache/runs/r2r3_final__20261026-230138__s0/` (25.9 M rows; not useful to us).

Quality (REPORT §1, §7): final CV M1 (macro hard-negative AUROC) 0.786, TEST 0.794; Cantley 0.71. R2 nophospho CV M1
0.771 ≈ R2 v3 0.768 ≈ clean 0.771 — **phospho features add ~nothing to KS prediction**. Tyr kinases weak (TEST
M1 0.595 ≈ popularity floor 0.585).

## 2. Leakage map

How phospho data enter (features_phospho.md, REPORT §3):
- **F6a co-regulation** (label-derived, cross-fitted): correlation of each site with a kinase's training substrates
  across 66 compendium study blocks (all ~178 PMIDs, incl. Hijazi 31959955, decryptM 36926954) + 10 CPTAC cohorts.
- **F6b inhibitor deconvolution**: ridge of site log2FC on a condition × kinase occupancy matrix over **692 compendium
  conditions** (Hijazi 183 conditions using the *provided kinhib matrix = our ground-truth target annotation*;
  527 affinity-mapped via Klaeger/Reinecke/CurveCurator/Davis, mainly decryptM 290, PMID 32034161 54, 28188228 26,
  38047086 24 ...). This is direct leakage of our ground truth (drug→target→site response).
- **F6c**: CPTAC kinase abundance / A-loop correlation (CPTAC only).
- **Site existence U_ms** (all variants): sites with ≥1 MS observation (PSP MS, compendium, CPTAC, Ochoa). Restricts
  training positives/negatives and candidate universe; R3 "phospho-channel" bit marks flanking U_ms sites;
  `cantley_*_pct_universe` uses U_ms as reference; `primed_soft` weights flanks by # MS channels (incl. compendium).
- Ochoa 2020: only site existence (functional score excluded). benchmarKIN tumour data: never a feature (CPTAC F6 used
  LinkedOmics `external/cptac_pancan_v1.2`, 10 cohorts, superset of benchmarKIN's 7).

| Variant | Hijazi | Other compendium perturbations | CPTAC | Site existence | Safe for our eval? |
|---|---|---|---|---|---|
| Deployment final / R2 deploy / release top-k | F6a+F6b (target matrix itself) | F6a+F6b | F6a+F6c | yes | **NO** |
| R2 clean_bench (DEV) | removed | **still used** (37 PMIDs, 483 F6b conds; decryptM etc.) | removed | yes | only if eval is restricted to the 22 clean_bench PMIDs; NO for other compendium conditions |
| R2 nophospho (DEV) | none | none | none | yes | **YES** (residual: site existence only) |
| R3 deployment (DEV+TEST) | none | none | none | yes (phospho bit, U_ms rows) | **YES** (residual: site existence only) |
| Cantley `*_pct_published`, `*_unprimed_score` | none | none | none | no | **YES, fully** |
| Cantley `pct_universe`, `primed_*` | none | none | none | U_ms reference / flank weights | yes (trivial residual) |

Residual site-existence leakage: knowing that a site was quantified somewhere carries no sign/magnitude of the
perturbation response; it is negligible for activity inference but should be disclosed (KSM's own A31e disclosure).
Use our own "measured" universe anyway (only measured sites enter z-score/KSEA).

**PSP-label leakage.** All variants are trained on PSP (human + mouse/rat mapped by 15-mer). Our ground truth is
perturbation-based, so this is shared with every PSP baseline and not leakage per se. Two caveats: (i) PSP may curate
HTP evidence from benchmark studies (KSM downstream §8; undetectable since PSP KS rows carry no PMIDs); (ii) if we ever
evaluate on held-out PSP pairs, DEV+TEST models (deployment, R3 deploy) have seen all PSP; DEV-only runs (nophospho,
clean) held out TEST proteins (`KSM/results/splits/test_proteins.tsv`, 3,957 proteins, cluster-level).

**Recommendation.** Use **R2 nophospho (DEV)** and/or **R3 deploy** full matrices; optionally fuse them
(z-score per model, w≈0.6 R2 / 0.4 R3 as in KSM). Do not use the release files or deployment final.

**Re-running a no-phospho deployment variant (DEV+TEST, 5 seeds).** `scripts/models/r2_deploy.py` hardcodes the
frozen config (incl. F6a), so it needs a copied script with `groups=[F1,F0,F2a,F3,F4,F5,F7]` and `--families f2a_sep
f3nb`. Code resolves paths from `ROOT = parents[2]` of `src/ksp`, so work in a copy/symlink farm (it writes to
`cache/`, `results/runs/`). Steps: (1) `r2_deploy.py --stage xfit --families f2a_sep f3nb` ≈ 5 min parallel;
(2) `R2_THREADS=16 r2_deploy.py --stage train --seeds 0 1 2 3 4` ≈ 5 × (2.5 min train + 3–4 min inference)
≈ 30 min (seed-0 DEV nophospho took 152 s train + 251 s inference); (3) fuse with R3 deploy via a modified
`scripts/models/r2_release.py` ≈ 3 min. Peak RAM ~20 GB; needs `cache/features` (54 GB), `cache/index`,
`cache/labels`, `cache/splits`, `results/labels`. Total ≈ 45 min. Gain over the existing DEV-only seed-0 run is small;
prefer using the existing files first.

## 3. benchmarKIN

Downloaded: `KSM/external/benchmarkin/` (Zenodo 12566560, CC BY 4.0, 451 MB): `benchmark_data_{hernandez,hijazi}.csv`
+ metadata, CPTAC tumour `*_original_medcent_30plus.tsv` (BRCA, CCRCC, GBM, HNSCC, LSCC, LUAD, UCEC; 95–122 tumours
each), gold standard `benchmarking_gold_standard_top5per.Rds`, plus `networkin_human_predictions_3.1.tsv`.
Tyrosine set: `KSM/external/benchmarkin_github/tyrosine{Data,Meta}.rda` (18 experiments from 6 PMIDs, none in our compendium; GPL-3).
Mapped caches (python-ready): `KSM/cache/downstream/perturbation_data.parquet` (60,418 sites × 230 experiments + meta
cols), `perturbation_meta.parquet` (404 experiment–target rows: experiment, dataset, drug, PMID, Time_min, sign,
gene, kinase_acc), `tumour_<COHORT>.parquet`, `tumour_gs.parquet` (13,232 GS pairs, 400 kinases), mapping logs JSON.
Counts: 230 perturbation experiments (Hernandez 94, Hijazi 118 [HL60, MCF7], tyrosine 18), 80 target kinases; site
mapping 90.7% / 96.1% / 87.4% / 99.8%.

Harness (pure Python, validated vs paper: PSP+z-score 0.782 vs 0.79; tumour 0.664 vs 0.662; KL 0.650 vs 0.65;
shuffled 0.518): `KSM/src/ksp/downstream/{bench_data,libraries,inference,evaluate,paired,ldo,api}.py`, scripts
`KSM/scripts/downstream/{run_validation,run_predictions,pred_topk,run_final,selftest}.py`. Methods: zscore, ksea,
wmean, ulm (minsize 5); perturbation AUROC (RMS-scale, sign, mcauroc 1000×), scaled rank/PHit, tumour AUROC.
Depends on `ksp.index`/`ksp.kinases` canonical stores — copy those too. Uses `rdata` to read .Rds/.rda.

Overlap with our compendium (`KSM/results/features/phospho/benchmarkin_studies.tsv`): 22 of 24 Hernandez PMIDs are
in the compendium (276 conditions) and **Hijazi 31959955 is in both** (118 benchmarKIN exps vs 183 compendium
conditions). Not in compendium: PMIDs 23592840, 23612710, 24204825 (Hernandez) and all 18 tyrosine experiments.
→ benchmarKIN perturbation is largely *not* independent of our compendium; the tyrosine set, the 2 non-compendium
Hernandez studies and the CPTAC tumour benchmark are the independent parts.
`KSM/src/ksp/downstream/ldo.py` classifies experiments drug_seen/unseen/no_drug against the compendium.

## 4. Other reusable assets

| Asset | Path (under KSM/) | Note |
|---|---|---|
| Site universe + source flags | `cache/sites/universe.parquet` (350,740 × 33) | has `comp_uid` = our compendium peptide_UID → site_id (177,410 exact, 1,828 shifted); also PSP/CPTAC/Ochoa flags, windows |
| Site/protein/kinase index | `cache/index/{sites,proteins,kinases}.parquet`; `src/ksp/{index,sites,ids,uniprot}.py` | canonical UniProt accs, isoform→canonical by window; site_id `ACC_RES POS` |
| Source mapping log | `cache/sites/mapping_log.parquet` (684k rows) | per-source map_method/reason; reuse for isoform QC |
| Compendium long/wide | `cache/compendium_{long,wide}.parquet` | 17.0 M rows (peptide_UID, conditionID, log2FC); our data, just faster |
| Kinase list (503, incl./excl. reasons, residue class, Manning, pseudokinase) | `results/labels/kinase_list.tsv`, `cache/index/kinases.parquet` | ready-made kinase universe |
| Kinase name resolver | `results/labels/kinase_synonyms.tsv`; `src/ksp/kinases.py` (`resolve_kinase_names`) | maps PSP/Cantley/Manning/HGNC/benchmarKIN names → UniProt; key by accession |
| Manning / KinHub / KinCore / pkinfam | `external/{kinbase_manning2002,kinhub,kincore}`, `external/uniprot/pkinfam/` | raw sources; KinCore alignment of human PKs |
| Paralog map | `cache/splits/kinase_domain_identity.parquet` (15,404 pairs); `ksp.labels.NegativeMasks` | Par(k) = ≥60% KD identity or same Manning subfamily |
| Kinase LKO clusters | `results/splits/kinase_folds.tsv` | 120 kinase clusters (≥50% KD id ∪ family) — good for leave-family-out |
| PSP mouse/rat → human | `cache/psp_human_mapped_rows.parquet` (24,704 rows); labels `cache/labels/dev_pairs.parquet` | 15-mer mapping; DEV pairs only on disk (TEST sealed) |
| Substrate clusters / splits | `results/splits/{substrate_clusters,cv_folds,test_proteins}.tsv`; `src/ksp/splits.py`, `crossfit.py` | mmseqs ≥40% id clusters; cross-fit driver with LeakageError |
| ESM-2 kinase embeddings | `cache/features/esm2_{650m,3b}/kinase_{domain,full}.npy` (503 × 1280/2560) | kinase similarity for propagating to understudied kinases |
| ESM-2 site embeddings | `cache/features/esm2_*/site_pca64.npy`; `cache/r3/uobs/win_emb.npy` (9.4 GB) | per-site |
| PPI (systematic screens) | `cache/features/context/ppi/edges.parquet` + RWR; raw `external/{bioplex,huri,opencell,buljan2020,humancellmap}` | 279,677 edges, no curated KS; STRING/BioGRID shown leaky |
| Drug–kinase affinities | `data/Klaeger` (=ours), `external/{reinecke2024,curvecurator_zenodo_8399823,davis2011}`; loaders `src/ksp/features/phospho_inhib.py` (`build_affinity`, `occupancy`, `build_inhibition_matrix`) | **define inhibitor targets with occupancy dose/(dose+Kd)** for non-Hijazi compendium conditions |
| Per-condition drug/dose parse | `results/features/phospho/f6b_conditions.tsv` (1,937 conds: drugs, doses_nM, drug_ids, minutes, drop reason), `drug_aliases.tsv`, `conditions_meta.tsv` (ptype, minutes) | ready ground-truth candidates: 692 inhibitor conditions with occupancy; time ≤120 min |
| Non-separable kinases | `results/features/phospho/f6b_nonseparable_pairs.tsv`, `f6b_kinases.tsv` | 191 pairs / 109 kinases with near-identical inhibition profiles — needed for fair scoring of promiscuous drugs |
| A-loop sites | `results/features/phospho/f6c_aloop_sites.tsv` | activation-loop site per kinase (452 ranges) — activating-site ground truth |
| Kinase Library engine | `src/ksp/features/motif_cantley.py`, `external/kinase_library/` | clean-room atlas scoring of any 15-mer |
| Other KS sources | `external/{sugiyama2019,iptmnet,omnipath,phosformer,selphi2,ikipdb_*}`; NetworKIN in `external/benchmarkin` | Sugiyama in vitro (hurts as labels per Müller-Dott); iKiP contains 61% of PSP |
| Evaluation code | `src/ksp/eval.py` (paired cluster bootstrap), `src/ksp/downstream/paired.py` | matched-universe paired bootstrap |

## 5. Lessons most relevant to us (KSM LESSONS.md + pre-prune list + downstream findings)

1. Treat label-derived features and "unbiased" resources as leaky until overlap-tested (STRING has 34% of PSP KS pairs vs 4% controls; iKiP-DB 61% of PSP). (LESSONS 3)
2. Key everything by UniProt accession, never symbol (174/404 Cantley names differ). (pre-prune 2)
3. Run null controls (shuffled libraries, random scores) through the real harness before the first model; controls must match positives on label-free axes — a shuffled PSP library still scores 0.52 because |z| grows with √n. (LESSONS 2; downstream §5)
4. Pre-register selection rule, adoption threshold, # configs and comparator set; include any arm that beats you on any stratum. (LESSONS 4)
5. Keep an independent auditor recomputing headline numbers from saved predictions; agent claims untrusted (caught 4 majors). (LESSONS 5)
6. Compare on matched kinase/experiment sets with paired CIs; native (own-coverage) contrasts were misleading (0.692 vs 0.664 became +0.016 n.s.). (REPORT §8)
7. Predicted libraries reached only **parity** with PSP for activity inference, and trailed PSP by 0.04–0.06 on unseen drugs; their main gain is **coverage** (50–65 vs 41 perturbed kinases; 274 vs 73 tumour kinases) and "PSP ∪ top100" helped slightly on tumours. Adding Kinase Library top-15 lowers performance (0.65). (REPORT §8)
8. Z-score (RoKAI, no centering) ≥ KSEA ≈ ULM > wmean on both benchmarks; reproduce this before claiming gains. (downstream §5)
9. Many kinases are not separable by inhibitor profiles (109 kinases, cosine ≥0.95) — define drug targets by occupancy and score kinase groups, not single kinases, for promiscuous drugs. (features_phospho §2)
10. Budget memory against the job cgroup (128 GiB), ≤60 GB per agent; stream 131 M-row parquets with `iter_batches` (the full file OOM-killed the harness). Turn PROMPT asks into an owned checklist at T+0. (LESSONS 6; downstream §7)
