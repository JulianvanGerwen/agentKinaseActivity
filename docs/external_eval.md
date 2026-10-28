# External (TEST-only) evaluation sets: benchmarKIN

These sets are part of the final TEST (DESIGN §2/§4). Do NOT score them with any real method before FREEZE.
Everything here was validated with random / oracle / sign-flipped synthetic scores only.

## Files
- `scripts/prep_external.py` -> `cache/external/` (writes `prep_report.json` with all numbers below).
- `src/kinact/external_eval.py` (metrics, controls, TEST driver), `tests/test_external.py`.
- Source: `external/benchmarkin_py/*` (benchmarKIN Zenodo 12566560 + GitHub tyrosine .rda, position-mapped to
  `{UniProt}_{res}{pos}`; CPTAC by 15-mer window). prep reuses those mapped parquets (not re-mapped).

## Sets and sizes
| set | content | n |
|---|---|---|
| `tyrosine` | 18 experiments, 6 PMIDs (Dasatinib x2, Imatinib, Bosutinib, Nilotinib, EGF x8, HRG, TCR x4), 2306 sites | 45 (exp, target) rows, 16 kinases |
| `hernandez_new` | 15 exps from PMIDs 23592840 (VEGF, sorafenib), 23612710 (BI2536), 24204825 (anti-CD3); 6432 sites | 23 rows, 7 kinases |
| tumour | 7 CPTAC cohorts (BRCA 122, CCRCC 103, GBM 99, HNSCC 108, LSCC 108, LUAD 110, UCEC 95 samples), 21.7k-35.7k sites each | 12,166 GS pairs in U (6786 GS+, 5380 GS-), 367 kinases |

Hijazi (118 exps) and the other 22 Hernandez PMIDs are in our compendium -> not external. All targets are in U.

## Mapping rates (source rows -> site key)
tyrosine 87.4% (2450/2803; loss = residue mismatch 328), Hernandez (whole file) 90.7%, CPTAC 99.75-99.81% per cohort.
Duplicate source rows mapping to one key are nan-averaged. 80.9% (tyrosine) / 91.0% (hernandez_new) of the mapped
site keys also exist in our compendium site universe (`cache/sites.parquet`).

## Drug identity (`pert_experiments.tsv`)
Columns: `drug`, `agent_keys` (canonical via `drug_aliases.tsv`), `stimulus_keys`, `drug_in_dev`. DEV = `gt_conditions.tsv`
(DEV + DEV-confirm only; asserted to hold no `test` rows); agents from `condition_doses.tsv`, stimuli by text match on
DEV Condition/Control (only EGF occurs). Seen in DEV: tyrosine 10/18 (dasatinib x2, EGF x8), hernandez_new 5/15 (BI2536
x4, sorafenib). Unseen: imatinib, bosutinib, nilotinib, HRG, TCR/anti-CD3, VEGF. Report seen/unseen strata (done in
`m1_perturbation`: `M1-all_drug_seen/unseen`).

## Tumour gold standard (note)
`tumour_gs` = benchmarKIN `benchmarking_gold_standard_top5per.Rds`: per cohort, samples with the top (GS+) / bottom (GS-)
5% protein level of that kinase. It is protein-abundance based, NOT activating-site based (DESIGN wording to be
corrected). Long format `tumour_gs_long.tsv`: cohort, sample, kinase, kinase_acc, label, sign (+1/-1), in_U.

## Metrics
Score convention: DataFrame condition(or sample) x HGNC kinase, NaN = unscored.
1. **M1 (ours)**. Perturbation: `Evaluator` on synthetic GT (group = PMID x drug, source `X`, M1-sel = 1..10 targets);
   negatives = all other U kinases, unscored negative/target = 1/2; macro over kinases. Tumour: per (cohort, kinase)
   AUROC of the score over GS+ vs GS- samples (1/2 if unscored), mean over cohorts, macro over kinases (`m1_tumour`).
2. **benchmarKIN own**. Perturbation (`benchmarkin_perturbation`): per-experiment RMS scaling (R `scale(center=FALSE)`),
   times sign, pooled; positives = (experiment, target), negatives = all other scored cells; experiments without a scored
   target dropped; decoupler mcauroc (1000 draws, negatives resampled to n_pos with replacement, seed 42). Tumour
   (`benchmarkin_tumour`): kinases with >=30 scored samples, per-cohort z over samples, GS+/GS- pooled over cohorts,
   1000 x 80% subsamples without replacement, mean AUROC.
Reference implementation: `$KINASESUB_DIR/src/ksp/downstream/`
(re-implemented, not imported).

## Validation (synthetic only)
Random scores: M1 0.49-0.50, benchmarKIN pooled 0.50-0.52 (tumour 0.50). Oracle: 1.00 on all; sign-flipped oracle < 0.01.
All-NaN scores -> M1 = 0.5 exactly. AUROC equals Mann-Whitney; RMS scaling checked. Reproducing benchmarKIN's published
numbers (PSP+z-score 0.79 perturbation, 0.662 tumour) is NOT possible here: no activity matrices or example outputs ship with
the data, and producing them requires running a real method (forbidden before TEST). The KSM reference harness
documents its own reproduction (0.782 / 0.664) on the same data.

## Running at TEST time (PI-launched script only)
```python
from kinact.external_eval import run_external
res = run_external(score_fn, confirm_test=True)   # score_fn(site x condition df) -> condition x kinase df
```
Runs both metrics on `tyrosine`, `hernandez_new`, tumour; for every compared method use the same call (the score matrices
must keep NaN for unscored kinases). For paired comparisons across methods, restrict to the common scored kinases
(M2-style) and bootstrap over `group` (PMID x drug) / cohort-kinase; not implemented here. Prep: `python scripts/prep_external.py`
(~10 s); tests: `pytest tests/test_external.py`.
