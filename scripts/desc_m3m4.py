#!/usr/bin/env python
"""Descriptive M3 (benchmarKIN-compatible) and M4 (kinase-centric within study) on DEV GT conditions (folds 0-4).
DESIGN §5 'Descriptive' / §12 C4. Not used for model selection.

M3: kinaseSubMain reference `ksp.downstream.evaluate.perturbation_eval` (read-only import): experiments = DEV GT
    conditions with >= 1 direct high/medium positive; targets = those positives; experiment sign = majority sign
    (kinact.controls.condition_sign). Pooled AUROC = scores R-scale()d per experiment, times sign, all other scored
    kinases are negatives (benchmarKIN: no masking), background experiments without a scored target removed,
    decoupler-style MC subsampling (1000 x, negatives balanced to #positives, seed 42). Per-experiment ranking:
    rank of each target among scored kinases (sign-adjusted), scaled rank = rank / n_scored (lower = better),
    mean over (kinase) of per-kinase mean; P(hit top-k).
M4: per (kinase k, study = PMID) among DEV GT conditions where k is a direct high/medium positive in >= 1 and an
    evaluable negative (not listed in gt_pairs) in >= 1 condition, both scored: AUROC of s_k * score_k (s_k =
    majority sign of k's positives in that study) positives vs negatives; mean over studies per kinase, macro over
    kinases. Scored kinases only.
"""
import json
import os
import sys

import numpy as np
import pandas as pd

from desc_common import FOLDS, METHODS, OUT, Data, score_all
from kinact import controls as ctl
from kinact.external_eval import auroc

sys.path.insert(0, os.path.join(os.environ.get("KINASESUB_DIR", "kinaseSubMain"), "src"))
from ksp.downstream.evaluate import perturbation_eval  # noqa: E402


class Bench:
    def __init__(self, experiments, meta):
        self.experiments, self.meta = experiments, meta


def m3(d, S, conds):
    pos = d.P[(d.P["cls"] == "direct") & d.P["confidence"].isin(["high", "medium"]) & d.P["conditionID"].isin(conds)]
    exps = sorted(set(pos["conditionID"]))
    sign = ctl.condition_sign(d.P, exps)
    bench = Bench(pd.DataFrame({"experiment": exps, "sign": sign.reindex(exps).to_numpy()}),
                  pd.DataFrame({"experiment": pos["conditionID"], "kinase_acc": pos["kinase"], "gene": pos["kinase"]}))
    act = S.reindex(index=exps).T                       # kinase x experiment
    r, _ = perturbation_eval(act, bench, n_iter=1000, seed=42)
    return r


def m4(d, S):
    gc = d.C[d.C["split"].isin(FOLDS) & (d.C["role"] == "gt") & (d.C["n_targets"] > 0)]
    gc = gc[gc["conditionID"].isin(S.index)]
    pos = d.P[(d.P["cls"] == "direct") & d.P["confidence"].isin(["high", "medium"])
              & d.P["conditionID"].isin(gc["conditionID"])]
    listed = set(zip(d.P["conditionID"], d.P["kinase"]))
    pmid = gc.set_index("conditionID")["PMID"]
    rows = []
    for (st, k), g in pos.assign(PMID=pos["conditionID"].map(pmid)).groupby(["PMID", "kinase"]):
        if k not in S.columns:
            continue
        sk = 1.0 if g["sign"].astype(float).sum() >= 0 else -1.0
        study = list(pmid.index[pmid == st])
        negc = [c for c in study if (c, k) not in listed]
        tp = (sk * S.loc[list(g["conditionID"]), k]).dropna()
        tn = (sk * S.loc[negc, k]).dropna() if negc else pd.Series(dtype=float)
        if len(tp) and len(tn):
            rows.append(dict(kinase=k, PMID=st, auc=auroc(tp.to_numpy(), tn.to_numpy()), n_pos=len(tp), n_neg=len(tn)))
    t = pd.DataFrame(rows)
    if t.empty:
        return {"M4": np.nan, "n_kinases": 0, "n_kinase_studies": 0}, t
    pk = t.groupby("kinase")["auc"].mean()
    return {"M4": float(pk.mean()), "n_kinases": int(len(pk)), "n_kinase_studies": int(len(t)),
            "frac_kinases_gt_0.5": float((pk > 0.5).mean())}, t


def main():
    d = Data()
    S = score_all(d, masked=False)
    res, m4t = {}, []
    allc = list(d.gt_conds)
    sel = list(d.C.loc[d.C["conditionID"].isin(allc) & d.C["in_M1sel"].astype(str).eq("True"), "conditionID"])
    for m in METHODS:
        Sm = S[m].reindex(index=allc)
        r = {"M1-sel_check": d.ev.metric(d.ev.pairs(Sm), "M1-sel")["value"]}
        r["M3_all"] = m3(d, Sm, allc)
        r["M3_M1sel"] = m3(d, Sm, sel)
        r["M4"], t = m4(d, Sm)
        m4t.append(t.assign(method=m))
        res[m] = r
        print(m, round(r["M1-sel_check"], 4), {k: round(v, 3) for k, v in r["M3_all"].items() if isinstance(v, float)},
              r["M4"], flush=True)
    # M3 on the matched universe: (condition, kinase) cells scored by every phospho method (others -> NaN)
    ph = [m for m in METHODS if m != "prior-only"]
    common_mask = np.logical_and.reduce([S[m].reindex(index=allc).notna().to_numpy() for m in ph])
    for m in METHODS:
        Sm = S[m].reindex(index=allc).where(common_mask)
        res[m]["M3_matched"] = m3(d, Sm, allc)
    # M4 on kinase-study units evaluable for every phospho method (common set)
    T = pd.concat(m4t)
    cnt = T[T.method.isin(ph)].groupby(["kinase", "PMID"]).size()
    common = set(cnt.index[cnt == len(ph)])
    for m in METHODS:
        tt = T[(T.method == m).to_numpy() & np.array([x in common for x in zip(T.kinase, T.PMID)], bool)]
        res[m]["M4_common"] = {"M4": float(tt.groupby("kinase")["auc"].mean().mean()),
                               "n_kinases": int(tt.kinase.nunique()), "n_kinase_studies": int(len(tt))}
    T.to_csv(OUT / "m4_units.tsv", sep="\t", index=False)
    (OUT / "m3m4.json").write_text(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main()
