#!/usr/bin/env python
"""Descriptive residue-matched null library (follow-up of docs/audit/r2np_leakage.md). Not used for selection.

For each kinase, every core-v1 library edge (curated_union + r2np top-100, models/core_v1/edges.tsv.gz) whose site is
measured in >= 1 DEV condition is replaced by a random DEV-measured site with the SAME residue (S/T/Y), sampled
without replacement within the kinase; weight and source are kept; unmeasured edges are dropped (they never enter
the likelihood). So per kinase the library size (in measured sites) and residue composition are matched, while site
identity is random. Scored with the core-v1 model (same config / kernel / rho_e) on the DEV GT conditions, 3 seeds.
Reported: M1-sel, M1-sel-scored, coverage and M1-sel by source (H/D/O), next to the real methods.
"""
import json

import numpy as np
import pandas as pd

from desc_common import METHODS, OUT, ROOT, Data, score_all
from kinact.core import KinaseActivityModel

SEEDS = [0, 1, 2]


def matched_edges(e, measured, seed):
    rng = np.random.default_rng(seed)
    e = e[e["site"].isin(measured)].copy()
    res_of = lambda s: pd.Series(s).str.extract(r"_([STY])\d+$")[0].to_numpy()
    e["res"] = res_of(e["site"])
    pool = {r: np.asarray(sorted(s for s in measured if s.rsplit("_", 1)[1][0] == r)) for r in "STY"}
    out = []
    for k, g in e.groupby("kinase", sort=False):
        g = g.drop_duplicates("site").copy()
        for r, gg in g.groupby("res"):
            g.loc[gg.index, "site"] = rng.choice(pool[r], len(gg), replace=False)
        out.append(g)
    return pd.concat(out)[["site", "kinase", "weight", "source"]]


def summ(d, S):
    P = d.ev.pairs(S)
    r = {"M1-sel": d.ev.metric(P, "M1-sel")["value"], "M1-sel-scored": d.ev.metric(P, "M1-sel-scored")["value"],
         "coverage": float(P[P["in_M1sel"]]["target_scored"].mean())}
    for s in "HDO":
        r[s] = d.ev.metric(P, f"M1-sel:source={s}")["value"]
    return r


def main():
    d = Data()
    X = d.X[d.gt_conds].astype(np.float64)
    measured = set(d.X.index[np.isfinite(d.X.to_numpy()).any(1)])
    base = KinaseActivityModel.from_config()
    e = pd.read_csv(ROOT / "models/core_v1/edges.tsv.gz", sep="\t")
    kin = pd.read_csv(ROOT / "models/core_v1/kinases.tsv", sep="\t", keep_default_na=False, na_values=[""])
    res = {}
    S = score_all(d, masked=False)
    for m in METHODS:
        res[m] = summ(d, S[m].reindex(index=d.gt_conds))
        print(m, res[m], flush=True)
    em = e[e["site"].isin(measured)]
    res["_lib"] = {"edges_measured": int(len(em)), "residue_share": em["site"].str.extract(r"_([STY])")[0]
                   .value_counts(normalize=True).round(3).to_dict()}
    nulls = []
    for sd in SEEDS:
        mdl = KinaseActivityModel(base.cfg, matched_edges(e, measured, sd), kin)
        Z = mdl.score(X).z.reindex(columns=d.kin)
        r = summ(d, Z); nulls.append(r)
        print("resnull seed", sd, r, flush=True)
    res["residue-matched null (core-v1 model)"] = {k: float(np.mean([n[k] for n in nulls])) for k in nulls[0]}
    res["residue-matched null (core-v1 model)"]["sd_M1-sel"] = float(np.std([n["M1-sel"] for n in nulls], ddof=1))
    res["residue-matched null per seed"] = nulls
    (OUT / "resnull.json").write_text(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main()
