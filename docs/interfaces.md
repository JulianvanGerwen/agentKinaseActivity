# Interfaces (shared contract between GT/split, library, harness, methods) — v1

Site key: `{UniProt}_{residue}{pos}` e.g. `P28482_T185` (canonical, from cache/sites.parquet). Kinase key: HGNC symbol in U.

## Ground truth / splits (written by scripts/build_gt.py, scripts/make_split.py)
- `results/ground_truth/kinase_universe.tsv`: kinase, uniprot, manning_group, manning_family, is_pseudokinase
- `results/ground_truth/gt_pairs.tsv`: conditionID, kinase, cls ∈ {direct, downstream, masked}, sign ∈ {-1,+1} (0 allowed only for masked), confidence, rule_id, note. Kinases in U absent for a GT condition = negatives.
- `results/ground_truth/gt_conditions.tsv`: conditionID, PMID, source ∈ {H,D,O}, group_id, n_targets, in_M1sel (bool), role ∈ {gt, null, specificity}, split ∈ {fold0..fold4, devconfirm, test}
- `results/ground_truth/GT_FREEZE.json`: sha256 of the above + kinase_universe + split.

## Libraries (written by scripts/build_libraries.py)
- `cache/libraries/<name>.parquet`: kinase, site, weight (float>0), source, evidence. Only sites present in the compendium need be kept, but keep all for external data.

## Method outputs
- A method is a function `score(X: pd.DataFrame sites×conditions log2FC (NaN=missing), lib) -> pd.DataFrame conditions×kinases` (NaN = unscored), larger = more active. Optional companion n_substrates frame.
- Stored as `results/scores/<method>__<config_hash>.parquet` (DEV only).

## Harness (src/kinact/eval.py)
`evaluate(scores, gt_pairs, gt_conditions, split_filter, ...) -> dict of metrics + per-kinase table`, with bootstrap utilities; refuses split=='test' unless FREEZE.json exists.
