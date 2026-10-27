# Interfaces (shared contract between GT/split, library, harness, methods) — v1

Site key: `{UniProt}_{residue}{pos}` e.g. `P28482_T185` (canonical, from cache/sites.parquet). Kinase key: HGNC symbol in U.

## Ground truth / splits (written by scripts/build_gt.py, scripts/make_split.py)
Agent-visible files (DEV folds fold0..fold4 + devconfirm only) in `results/ground_truth/`; TEST rows live ONLY in `cache/sealed/` (same file names, never opened by development code; `cache/sealed/full/` holds all-split copies; `cache/sealed/draft_candidate_v0/` is the retired pre-G1 draft, unused).
Read every GT TSV with `kinact.eval.read_gt_tsv` (keep_default_na=False: `role == "null"` and kinase names such as `NA` are literal strings; only empty fields are NaN).
- `kinase_universe.tsv`: kinase, uniprot, manning_group, manning_family, is_pseudokinase
- `gt_pairs.tsv`: conditionID, kinase, cls ∈ {direct, downstream, masked}, sign ∈ {-1,+1} (0 only for masked), confidence ∈ {high, medium, low, -}, rule_id, note. Kinases in U absent for a GT condition = negatives.
- `gt_conditions.tsv`: conditionID, PMID, source ∈ {H,D,O}, group_id, n_targets, in_M1sel (bool), role ∈ {gt, null, specificity, excluded}, split ∈ {fold0..fold4, devconfirm, test}. Only role=="gt" with n_targets>0 enter M1.
- `harmony_bins.tsv`: conditionID, kinase, n_psp_quantified, bin (per condition-kinase); `harmony_bins_group.tsv`: kinase, group_id, n_psp_quantified, bin (frozen per (kinase, group): median over the group's positive conditions; bins 0 | 1-4 | 5-19 | >=20). `Evaluator(bins=...)` accepts either.
- `condition_groups.tsv`: conditionID, PMID, group_id, split, qc_status, gt_role, drugs.
- `GT_FREEZE.json`: sha256 of the above + kinase_universe + split.

## DEV compendium cache (`cache/dev/`, built by scripts/make_dev_cache.py; loaders in `kinact.data`)
- `quant_dev.npy` float32 [peptides x 1476 conditions] (NaN = not quantified), rows = `sites_dev.txt` (peptide_UID; map to canonical site keys via `cache/sites.parquet`), columns = `conditions_dev.txt`; `conditions_dev.tsv` = metadata + group_id, split (fold0..fold4), qc_status. TEST and devconfirm conditions are absent.
- `kinact.data.load_dev()` -> sites (canonical key) x conditions DataFrame; `load_dev_meta()`; `library_matrix(name, site_index)` -> sparse kinase x site CSR aligned to the site index.

## Libraries (written by scripts/build_libraries.py)
- `cache/libraries/<name>.parquet`: kinase, site, weight (float>0), source, evidence. Only sites present in the compendium need be kept, but keep all for external data.

## Method outputs
- A method is a function `score(X: pd.DataFrame sites×conditions log2FC (NaN=missing), lib) -> pd.DataFrame conditions×kinases` (NaN = unscored), larger = more active. Optional companion n_substrates frame.
- Stored as `results/scores/<method>__<config_hash>.parquet` (DEV only).

## Harness (src/kinact/eval.py)
`evaluate(scores, gt_pairs, gt_conditions, split_filter, ...) -> dict of metrics + per-kinase table`, with bootstrap utilities; refuses split=='test' unless FREEZE.json exists.
