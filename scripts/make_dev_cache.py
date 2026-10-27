#!/usr/bin/env python
"""Build the DEV-only data cache for development agents (DESIGN §4 physical sealing, B10).

Keeps only conditions whose split is fold0..fold4 (TEST and DEV-confirm groups removed), taken from
results/ground_truth/condition_groups.tsv (which itself lists no TEST conditions).

Outputs (cache/dev/):
  quant_dev.npy          float32 [n_sites x n_dev_conditions], NaN = not quantified (same row order as cache/sites.parquet)
  sites_dev.txt          row order (peptide_UID; == cache/quant_wide_sites.txt)
  conditions_dev.txt     column order (conditionID)
  conditions_dev.tsv     metadata of the DEV conditions (01_* columns + group_id, split, qc_status)
  README
Usage: python scripts/make_dev_cache.py
"""
import numpy as np
import pandas as pd

from gtlib import CACHE, GT, sha256

FOLDS = {f"fold{i}" for i in range(5)}


def main():
    cg = pd.read_csv(GT / "condition_groups.tsv", sep="\t", dtype=str, keep_default_na=False)
    dev = cg[cg["split"].isin(FOLDS)]
    conds = pd.read_parquet(CACHE / "conditions.parquet")
    cols = [l.strip() for l in open(CACHE / "quant_wide_conditions.txt")]
    pos = {c: j for j, c in enumerate(cols)}
    keep = [c for c in cols if c in set(dev["conditionID"])]
    X = np.load(CACHE / "quant_wide_f32.npy", mmap_mode="r")
    Xd = np.ascontiguousarray(X[:, [pos[c] for c in keep]])
    out = CACHE / "dev"
    out.mkdir(exist_ok=True)
    np.save(out / "quant_dev.npy", Xd)
    (out / "conditions_dev.txt").write_text("\n".join(keep) + "\n")
    (out / "sites_dev.txt").write_text(open(CACHE / "quant_wide_sites.txt").read())
    meta = conds[conds["conditionID"].isin(set(keep))].drop(columns=["cond_idx"]).merge(
        dev[["conditionID", "group_id", "split", "qc_status"]], on="conditionID")
    meta = meta.set_index("conditionID").loc[keep].reset_index()
    meta.to_csv(out / "conditions_dev.tsv", sep="\t", index=False)
    (out / "README").write_text(
        "DEV-only compendium cache (built by scripts/make_dev_cache.py from results/ground_truth/condition_groups.tsv).\n"
        f"Contains {len(keep)} of {len(cols)} conditions: splits fold0..fold4 only. TEST groups and DEV-confirm groups\n"
        "are removed (DESIGN §4, B10). Development agents must read compendium values ONLY from here.\n"
        "quant_dev.npy: float32 [n_sites x n_conditions], NaN = not quantified; rows = sites_dev.txt (peptide_UID,\n"
        "same order as cache/sites.parquet / quant_wide_sites.txt), columns = conditions_dev.txt.\n"
        "conditions_dev.tsv: metadata + group_id, split (CV fold), qc_status (scale_error/duplicate = excluded from\n"
        "evaluation; usable or not for unsupervised fits at the method's discretion, document it).\n"
        "Unsupervised use inside CV: exclude the held-out fold's conditions (DESIGN §4).\n"
        f"sha256(quant_dev.npy) = {sha256(out / 'quant_dev.npy')}\n")
    print(f"DEV cache: {Xd.shape}")


if __name__ == "__main__":
    main()
