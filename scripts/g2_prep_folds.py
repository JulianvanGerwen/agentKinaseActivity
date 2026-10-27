#!/usr/bin/env python
"""G2 prep, part 1a: re-split the DEV-CV groups into 5 folds with seeds 1-4 (DESIGN §4 winner's-curse check).

Same logic as scripts/make_split.py step 3, re-implemented from agent-visible files only (results/ground_truth DEV +
devconfirm rows; cache/sealed is never opened). TEST and DEV-confirm membership are NOT changed.
For seed s: rng = default_rng(s); the DEV-confirm greedy step's rng.shuffle calls are replayed on the same per-source
group lists (results discarded; DEV-confirm membership is fixed to the frozen split) so that seed 0 reproduces the
frozen folds exactly; then the frozen fold rule: shuffle DEV-CV groups, stable-sort by (M1-sel pairs*1000 + group
size) desc, assign each group to the fold with the smallest load of its stratum (M1-sel groups: source stratum,
weight = M1-sel pairs; other GT groups / non-GT groups: weight = group size in conditions).
Writes results/g2_prep/folds_seed{s}.tsv (conditionID, group_id, split) for s in 0..4 and checks seed 0 == frozen.
"""
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/workspace")
GT = ROOT / "results/ground_truth"
OUT = ROOT / "results/g2_prep"
FOLDS = [f"fold{i}" for i in range(5)]
N_FOLDS = 5


def rd(p):
    return pd.read_csv(p, sep="\t", keep_default_na=False, na_values=[""], dtype={"PMID": str})


def inputs():
    cg = rd(GT / "condition_groups.tsv")            # all non-TEST conditions
    gc = rd(GT / "gt_conditions.tsv")
    gp = rd(GT / "gt_pairs.tsv")
    assert not (cg["split"] == "test").any() and not (gc["split"] == "test").any()
    gc["in_M1sel"] = gc["in_M1sel"].astype(str) == "True"
    sel_c = set(gc.loc[gc["in_M1sel"], "conditionID"])
    src_of = gc.set_index("conditionID")["source"]
    gid = cg.set_index("conditionID")["group_id"]
    m1 = gp[(gp["cls"] == "direct") & gp["confidence"].isin(["high", "medium"]) & gp["conditionID"].isin(sel_c)].copy()
    m1["group_id"] = m1["conditionID"].map(gid)
    m1["source"] = m1["conditionID"].map(src_of)
    g_m1 = m1.groupby("group_id").size()
    g_src = m1.groupby("group_id")["source"].agg(lambda s: s.value_counts().sort_index().idxmax())
    gsize = cg["group_id"].value_counts()
    gt_groups = set(gc["group_id"])
    split_of_group = cg.groupby("group_id")["split"].agg(lambda s: s.iloc[0])
    assert (cg.groupby("group_id")["split"].nunique() == 1).all()
    return cg, g_m1, g_src, gsize, gt_groups, split_of_group


def folds_for_seed(seed, cg, g_m1, g_src, gsize, gt_groups, split_of_group):
    rng = np.random.default_rng(seed)
    nontest = set(split_of_group.index)
    strata = defaultdict(list)
    for g in g_m1.index:
        strata[g_src[g]].append(g)
    for st in sorted(strata):                          # replay the DEV-confirm greedy shuffles (membership fixed)
        dev_st = [g for g in sorted(strata[st]) if g in nontest]
        rng.shuffle(dev_st)
    dev_groups = sorted(g for g in nontest if split_of_group[g] in FOLDS)
    order = list(dev_groups)
    rng.shuffle(order)
    order = sorted(order, key=lambda g: -(g_m1.get(g, 0) * 1000 + gsize[g]))
    load, fold = defaultdict(float), {}
    for g in order:
        if g in g_m1.index:
            key, w = ("m1", g_src[g]), g_m1[g]
        elif g in gt_groups:
            key, w = ("gt", "all"), gsize[g]
        else:
            key, w = ("nongt", "all"), gsize[g]
        f = min(range(N_FOLDS), key=lambda i: (load[(key, i)], i))
        fold[g] = f"fold{f}"
        load[(key, f)] += w
    d = cg[cg["split"].isin(FOLDS)][["conditionID", "group_id"]].copy()
    d["split"] = d["group_id"].map(fold)
    return d


def main():
    args = inputs()
    cg = args[0]
    summ = {}
    for s in range(5):
        d = folds_for_seed(s, *args)
        if s == 0:
            ref = cg[cg["split"].isin(FOLDS)].set_index("conditionID")["split"]
            same = (d.set_index("conditionID")["split"] == ref).all()
            assert same, "seed 0 does not reproduce the frozen folds"
        d.to_csv(OUT / f"folds_seed{s}.tsv", sep="\t", index=False)
        summ[s] = d["split"].value_counts().sort_index().to_dict()
    print(json.dumps(summ, indent=1))
    print("seed 0 reproduces frozen folds: OK")


if __name__ == "__main__":
    main()
