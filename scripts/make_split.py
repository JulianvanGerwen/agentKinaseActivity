#!/usr/bin/env python
"""Grouped TEST / DEV split of the compendium (DESIGN v2 §4 + B4, B9, B10; PI decisions r3 2026-10-27).

1. Grouping (union-find over ALL 1,937 conditions; group = connected component):
   - every canonical drug ("agent", results/ground_truth/drug_aliases.tsv) found in the Condition OR Control text links
     all conditions mentioning it (class tokens such as 'MEKi' count only when no named compound is in the text);
     Hijazi <-> Wilkes (26060313) and every other cross-study reuse of a drug are linked this way;
   - conditions without any drug are linked by study (PMID);
   - identical-vector duplicates (results/ground_truth/condition_qc.tsv duplicate_of) are linked;
   - decryptM "0 nM vs 0 nM" null controls carry their drug name and therefore follow their drug group (B9).
   Group source = majority GT source (H/D/O) of the group's M1-sel pairs (else of its direct pairs; else 'none').

2. TEST selection (pre-registered; uses labels only, never method outputs). Units: M1-sel pairs = direct positives with
   confidence high/medium in conditions with in_M1sel. Harmony bin of a pair = bin of the median number of quantified
   PSP substrates of that kinase over the group's conditions where it is an M1-sel positive (r3 S4), bins 0/1-4/5-19/>=20.
   For seed s in 0..99: rng = default_rng(s); (i) within each source stratum the groups holding M1-sel pairs are shuffled
   and added greedily (a group is added iff it moves the stratum's TEST pair share closer to 0.30); (ii) two passes over
   the M1-sel groups in shuffled order toggle a group's TEST membership iff that lowers d(s) below; (iii) GT groups without
   M1-sel pairs are then added the same way on their high/medium direct pairs (or condition counts if they have none),
   pooled over sources. Groups without GT conditions stay in DEV.
   Distance (smaller = better):
       d(s) = mean over cells (source x harmony bin) with >= 20 M1-sel pairs of |TEST share_cell - 0.30|
              + P,   P = 1 + (0.20 - u) if u < 0.20 else 0,
   where u = fraction of TEST M1-sel kinases that are never a high/medium direct positive in any DEV ('gt') condition.
   The split with minimal d is chosen; ties -> lowest seed.

3. DEV: DEV-confirm = about 15 % of DEV M1-sel pairs per source stratum, chosen with the same greedy rule (seed 0);
   remaining DEV groups -> 5 folds (seed 0: shuffle, sort by size desc, assign to the fold with the fewest M1-sel
   pairs of the group's stratum; groups without M1-sel pairs balanced on condition counts).

Outputs
  results/ground_truth/gt_pairs.tsv, gt_conditions.tsv, harmony_bins.tsv, harmony_bins_group.tsv  (DEV + devconfirm rows)
  results/ground_truth/condition_groups.tsv   (all non-TEST conditions: conditionID, PMID, group_id, split, qc_status)
  results/ground_truth/test_condition_hashes.txt (salted sha256 of TEST conditionIDs; salt sealed)
  results/ground_truth/split_candidates.tsv, split_log.json
  cache/sealed/{gt_pairs,gt_conditions,harmony_bins,harmony_bins_group,condition_groups}.tsv (TEST rows only), salt.txt
Usage: python scripts/make_split.py [--outdir DIR] [--fulldir DIR]
"""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from gtlib import CACHE, GT, HIJAZI_PMID, SEALED, DrugMatcher, UnionFind, write_tsv

SEEDS = range(100)
TEST_SHARE, CONFIRM_SHARE, N_FOLDS = 0.30, 0.15, 5
UNSEEN_MIN, MIN_CELL = 0.20, 20
BINS = [-1, 0, 4, 19, 10 ** 9]
BIN_LABELS = ["0", "1-4", "5-19", ">=20"]
PAIR_COLS = ["conditionID", "kinase", "cls", "sign", "confidence", "rule_id", "note"]
COND_COLS = ["conditionID", "PMID", "source", "group_id", "n_targets", "in_M1sel", "role", "split"]


def build_groups(conds, qc):
    m = DrugMatcher()
    uf = UnionFind(conds["conditionID"])
    drug_of = {}
    for r in conds.itertuples():
        if r.PMID == HIJAZI_PMID:
            ags = m.find(str(r.Condition).split("_", 1)[1], r.PMID)
        else:
            ags = sorted(set(m.find(r.Condition, r.PMID)) | set(m.find(r.Control, r.PMID)))
        drug_of[r.conditionID] = ags
        if ags:
            for a in ags:
                uf.union(r.conditionID, "drug:" + a)
        else:
            uf.union(r.conditionID, "study:" + r.PMID)
    for cid, dup in zip(qc["conditionID"], qc["duplicate_of"]):
        if isinstance(dup, str) and dup:
            uf.union(cid, dup)
    root = {c: uf.find(c) for c in conds["conditionID"]}
    # deterministic group ids: order components by their smallest conditionID
    comp = defaultdict(list)
    for c, r in root.items():
        comp[r].append(c)
    order = sorted(comp.values(), key=lambda v: min(v))
    gid = {}
    for i, members in enumerate(order):
        for c in members:
            gid[c] = f"G{i:04d}"
    return gid, drug_of


def greedy(groups, size, rng, share, fixed_total=None):
    """Add groups (shuffled) iff they move the selected share of `size` closer to `share`."""
    g = list(groups)
    rng.shuffle(g)
    total = fixed_total if fixed_total is not None else sum(size[x] for x in g)
    target, cur, sel = share * total, 0.0, []
    for x in g:
        if abs(cur + size[x] - target) < abs(cur - target):
            sel.append(x)
            cur += size[x]
    return sel


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--fulldir", default=None)
    args = ap.parse_args(argv)
    full = Path(args.fulldir) if args.fulldir else SEALED / "full"
    out_res = Path(args.outdir) / "results" if args.outdir else GT
    out_seal = Path(args.outdir) / "sealed" if args.outdir else SEALED
    res_qc = (Path(args.outdir) / "results" if args.outdir else GT) / "condition_qc.tsv"

    conds = pd.read_parquet(CACHE / "conditions.parquet")
    conds["PMID"] = conds["PMID"].astype(str)
    qc = pd.read_csv(res_qc, sep="\t", dtype=str, keep_default_na=False)
    gc = pd.read_csv(full / "gt_conditions_full.tsv", sep="\t", dtype={"PMID": str}, keep_default_na=False)
    gc["in_M1sel"] = gc["in_M1sel"].astype(str) == "True"
    gp = pd.read_csv(full / "gt_pairs_full.tsv", sep="\t", dtype=str, keep_default_na=False)
    gp["sign"] = gp["sign"].astype(int)
    hb = pd.read_csv(full / "harmony_full.tsv", sep="\t")
    uni = pd.read_csv(out_res / "kinase_universe.tsv", sep="\t", keep_default_na=False)

    gid, drug_of = build_groups(conds, qc)
    gc["group_id"] = gc["conditionID"].map(gid)
    gsize = pd.Series(gid).value_counts()
    src_of = gc.set_index("conditionID")["source"]

    # M1-sel units
    sel_c = set(gc.loc[gc["in_M1sel"], "conditionID"])
    gt_c = set(gc.loc[gc["role"] == "gt", "conditionID"])
    dpos = gp[(gp["cls"] == "direct") & gp["confidence"].isin(["high", "medium"])].copy()
    dpos["group_id"] = dpos["conditionID"].map(gid)
    dpos["source"] = dpos["conditionID"].map(src_of)
    m1 = dpos[dpos["conditionID"].isin(sel_c)].merge(hb[["conditionID", "kinase", "n_psp_quantified"]],
                                                     on=["conditionID", "kinase"], how="left")
    kg_med = m1.groupby(["kinase", "group_id"])["n_psp_quantified"].median()
    m1["bin"] = pd.cut(m1.set_index(["kinase", "group_id"]).index.map(kg_med), BINS, labels=BIN_LABELS).astype(str)
    dpos_gt = dpos[dpos["conditionID"].isin(gt_c)]

    g_m1 = m1.groupby("group_id").size()
    g_dir = dpos_gt.groupby("group_id").size()
    g_ncond = gc.groupby("group_id").size()
    g_src = m1.groupby("group_id")["source"].agg(lambda s: s.value_counts().sort_index().idxmax())
    g_src2 = dpos_gt.groupby("group_id")["source"].agg(lambda s: s.value_counts().sort_index().idxmax())
    gt_groups = sorted(set(gc["group_id"]))
    m1_groups = sorted(g_m1.index)
    other_groups = [g for g in gt_groups if g not in g_m1.index]
    strata = defaultdict(list)
    for g in m1_groups:
        strata[g_src[g]].append(g)
    cells = m1.groupby(["source", "bin"]).size()
    cells = cells[cells >= MIN_CELL]
    pos_by_group = dpos_gt.groupby("group_id")["kinase"].apply(set).to_dict()
    m1_kin_by_group = m1.groupby("group_id")["kinase"].apply(set).to_dict()

    # dense per-group matrices for fast evaluation
    G = sorted(set(gt_groups))
    gi = {g: i for i, g in enumerate(G)}
    cell_idx = {c: j for j, c in enumerate(cells.index)}
    Gcell = np.zeros((len(G), len(cells)))
    for (g, s_, b), n in m1.groupby(["group_id", "source", "bin"]).size().items():
        if (s_, b) in cell_idx:
            Gcell[gi[g], cell_idx[(s_, b)]] += n
    kin_all = sorted(set(dpos_gt["kinase"]) | set(m1["kinase"]))
    ki = {k: j for j, k in enumerate(kin_all)}
    Gpos = np.zeros((len(G), len(kin_all)), bool)
    Gm1k = np.zeros((len(G), len(kin_all)), bool)
    for g, ks in pos_by_group.items():
        Gpos[gi[g], [ki[k] for k in ks]] = True
    for g, ks in m1_kin_by_group.items():
        Gm1k[gi[g], [ki[k] for k in ks]] = True
    tot = cells.to_numpy(float)
    m1w = np.array([g_m1.get(g, 0) for g in G], float)

    def evaluate_mask(t):
        share = (t @ Gcell) / tot
        d = float(np.mean(np.abs(share - TEST_SHARE)))
        dev_pos = Gpos[~t].any(0)
        test_k = Gm1k[t].any(0)
        nk = int(test_k.sum())
        nu = int((test_k & ~dev_pos).sum())
        u = nu / max(1, nk)
        pen = 1 + (UNSEEN_MIN - u) if u < UNSEEN_MIN else 0.0
        return d + pen, d, u, float(t @ m1w / m1w.sum()), nk, nu

    def evaluate(test):
        return evaluate_mask(np.array([g in test for g in G]))

    cand = []
    tests = {}
    for s in SEEDS:
        rng = np.random.default_rng(s)
        # pass 1: per-source greedy on pair share (starting point)
        test = []
        for st in sorted(strata):
            test += greedy(sorted(strata[st]), g_m1.to_dict(), rng, TEST_SHARE)
        t = np.array([g in set(test) for g in G])
        # passes 2-3: shuffled single-group toggles that lower the full pre-registered distance
        for _ in range(2):
            order = [gi[g] for g in m1_groups]
            rng.shuffle(order)
            cur = evaluate_mask(t)[0]
            for i in order:
                t[i] = ~t[i]
                new_d = evaluate_mask(t)[0]
                if new_d < cur:
                    cur = new_d
                else:
                    t[i] = ~t[i]
        test = [g for g in G if t[gi[g]] and g in g_m1.index]
        sz = {g: int(g_dir.get(g, 0)) for g in other_groups}
        with_d = [g for g in other_groups if sz[g] > 0]
        test += greedy(with_d, sz, rng, TEST_SHARE)
        no_d = [g for g in other_groups if sz[g] == 0]
        test += greedy(no_d, g_ncond.to_dict(), rng, TEST_SHARE)
        test = set(test)
        tests[s] = test
        dist, d, u, sh, nk, nu = evaluate(test)
        cand.append(dict(seed=s, distance=round(dist, 6), cell_dev=round(d, 6), unseen_frac=round(u, 4),
                         test_m1sel_pair_share=round(sh, 4), test_m1sel_kinases=nk, test_unseen_kinases=nu,
                         test_conditions=int(gc["group_id"].isin(test).sum())))
    cand = pd.DataFrame(cand)
    best = cand.sort_values(["distance", "seed"]).iloc[0]
    seed = int(best["seed"])
    test = tests[seed]

    # DEV-confirm and folds (seed 0)
    rng = np.random.default_rng(0)
    confirm = []
    for st in sorted(strata):
        dev_st = [g for g in sorted(strata[st]) if g not in test]
        confirm += greedy(dev_st, g_m1.to_dict(), rng, CONFIRM_SHARE)
    confirm = set(confirm)
    all_groups = sorted(set(gid.values()))
    dev_groups = [g for g in all_groups if g not in test and g not in confirm]
    fold = {}
    load = defaultdict(float)
    order = list(dev_groups)
    rng.shuffle(order)
    order = sorted(order, key=lambda g: -(g_m1.get(g, 0) * 1000 + gsize[g]))   # stable sort keeps shuffled ties
    for g in order:
        if g in g_m1.index:
            key, w = ("m1", g_src[g]), g_m1[g]
        elif g in set(gc["group_id"]):
            key, w = ("gt", "all"), gsize[g]
        else:
            key, w = ("nongt", "all"), gsize[g]
        f = min(range(N_FOLDS), key=lambda i: (load[(key, i)], i))
        fold[g] = f"fold{f}"
        load[(key, f)] += w
    split_of_group = {g: ("test" if g in test else "devconfirm" if g in confirm else fold[g]) for g in all_groups}

    # ---------------------------------------------------------------- write
    gc["split"] = gc["group_id"].map(split_of_group)
    allc = conds[["conditionID", "PMID"]].copy()
    allc["group_id"] = allc["conditionID"].map(gid)
    allc["split"] = allc["group_id"].map(split_of_group)
    allc = allc.merge(qc[["conditionID", "qc_status"]], on="conditionID", how="left")
    allc["drugs"] = allc["conditionID"].map(lambda c: ";".join(drug_of[c]))
    allc["gt_role"] = allc["conditionID"].map(gc.set_index("conditionID")["role"]).fillna("")
    is_test_c = allc["split"] == "test"
    test_ids = set(allc.loc[is_test_c, "conditionID"])

    gp["split"] = gp["conditionID"].map(gc.set_index("conditionID")["split"])
    hb["group_id"] = hb["conditionID"].map(gid)
    hb["split"] = hb["conditionID"].map(gc.set_index("conditionID")["split"])
    hbg = hb[hb["cls"] == "direct"].groupby(["kinase", "group_id"], as_index=False).agg(
        n_psp_quantified=("n_psp_quantified", "median"), split=("split", "first"))
    hbg["bin"] = pd.cut(hbg["n_psp_quantified"], BINS, labels=BIN_LABELS).astype(str)

    def put(df, cols, name, sort):
        dv = df[df["split"] != "test"].sort_values(sort)
        ts = df[df["split"] == "test"].sort_values(sort)
        write_tsv(dv[cols], out_res / f"{name}.tsv")
        write_tsv(ts[cols], out_seal / f"{name}.tsv")

    gc_out = gc.copy()
    gc_out["in_M1sel"] = gc_out["in_M1sel"].astype(bool)
    put(gc_out, COND_COLS, "gt_conditions", ["conditionID"])
    put(gp, PAIR_COLS, "gt_pairs", ["conditionID", "kinase"])
    put(hb, ["conditionID", "kinase", "n_psp_quantified", "bin"], "harmony_bins", ["conditionID", "kinase"])
    put(hbg, ["kinase", "group_id", "n_psp_quantified", "bin"], "harmony_bins_group", ["kinase", "group_id"])
    cg_cols = ["conditionID", "PMID", "group_id", "split", "qc_status", "gt_role", "drugs"]
    put(allc, cg_cols, "condition_groups", ["conditionID"])
    salt = hashlib.sha256((out_seal / "gt_pairs.tsv").read_bytes()).hexdigest()
    (out_seal / "salt.txt").write_text(salt + "\n")
    hashes = sorted(hashlib.sha256((salt + c).encode()).hexdigest() for c in test_ids)
    (out_res / "test_condition_hashes.txt").write_text("\n".join(hashes) + "\n")
    write_tsv(cand, out_res / "split_candidates.tsv")

    # summary log (aggregates only)
    fam = dict(zip(uni["kinase"], [f if f else k for k, f in zip(uni["kinase"], uni["manning_family"])]))
    dev_pos = set(dpos_gt.loc[~dpos_gt["group_id"].isin(test), "kinase"])
    test_m1k = set(m1.loc[m1["group_id"].isin(test), "kinase"])
    dev_fam = {fam.get(k, k) for k in dev_pos}
    largest = gsize.idxmax()
    lg_pmids = sorted(allc.loc[allc["group_id"] == largest, "PMID"].unique())
    log = dict(chosen_seed=seed, chosen_distance=float(best["distance"]), cell_deviation=float(best["cell_dev"]),
               unseen_frac=float(best["unseen_frac"]), n_groups=len(all_groups), n_gt_groups=len(gt_groups),
               largest_group=dict(group_id=largest, n_conditions=int(gsize[largest]), n_pmids=len(lg_pmids),
                                  pmids=lg_pmids, split=split_of_group[largest]),
               n_conditions_by_split=allc["split"].value_counts().sort_index().to_dict(),
               n_gt_conditions_by_split=gc["split"].value_counts().sort_index().to_dict(),
               n_m1sel_pairs_by_split=m1.assign(split=m1["group_id"].map(split_of_group))["split"].value_counts().sort_index().to_dict(),
               test_m1sel_kinases=len(test_m1k), T_unseen=len(test_m1k - dev_pos),
               T_family_unseen=len({k for k in test_m1k if fam.get(k, k) not in dev_fam}),
               n_feasible_candidates=int((cand["unseen_frac"] >= UNSEEN_MIN).sum()))
    (out_res / "split_log.json").write_text(json.dumps(log, indent=1, sort_keys=True) + "\n")
    print(json.dumps(log, indent=1))


if __name__ == "__main__":
    main()
