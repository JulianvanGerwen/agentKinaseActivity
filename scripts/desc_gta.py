#!/usr/bin/env python
"""Descriptive GT-A (DESIGN §3 / §12 B13): per-kinase Spearman between the kinase's score and its own activating-site
log2FC across non-inhibitor DEV conditions. Not used for model selection.

- Activating self-sites: selfsites.psp_self_sites(U) with sign +1 (PSP Regulatory_sites '(enzymatic) activity,
  induced'; conflicting annotations dropped). Per condition: mean raw log2FC over the kinase's measured activating sites.
- Scores: computed on input with ALL self-sites masked (selfsites.self_site_keys: psp+aloop + every PSP 'activity'
  regulatory site on a U kinase), every method (desc_common.score_all(masked=True)).
- Conditions: desc_common.non_inhibitor_conditions (qc ok, no kinase-inhibitor agent, not Hijazi): stimuli, genetic,
  other drugs, time courses. Pairs with both values finite; kinase evaluable if >= 20 such conditions in >= 2 studies.
- Null: 2000 permutations of the self-site vector WITHIN study (PMID) (keeps study structure); one-sided
  p = (1 + #{rho* >= rho}) / 2001; BH over each method's evaluable kinases; significant = q <= 0.10.
- prior-only has no within-kinase variation outside GT conditions -> not applicable.
"""
import json

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from desc_common import OUT, Data, non_inhibitor_conditions, score_all
from kinact.explore import selfsites as SS

N_PERM = 2000
MIN_N = 20
METHODS = ["core-v1", "core-v2", "psp-z-min3", "cur-z-min1"]


def bh(p):
    p = np.asarray(p, float); n = len(p)
    o = np.argsort(p); q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(q, 1); return out


def spearman_perm(x, y, strata, rng):
    rx = rankdata(x); ry = rankdata(y)
    rx = (rx - rx.mean()) / rx.std(); ry = (ry - ry.mean()) / ry.std()
    obs = float(np.mean(rx * ry))
    codes = np.unique(strata, return_inverse=True)[1].astype(float)
    base = np.argsort(codes, kind="stable")
    order = np.argsort(codes[None, :] * 2 + rng.random((N_PERM, len(ry))), axis=1)
    perm = np.empty_like(order); perm[:, base] = order
    P = ry[perm]
    null = P @ rx / len(rx)
    return obs, (1 + int((null >= obs - 1e-12).sum())) / (N_PERM + 1)


def main():
    d = Data()
    conds = non_inhibitor_conditions(d)
    S = score_all(d, masked=True)
    act = SS.psp_self_sites(d.Udf)
    act = act[(act["sign"] > 0) & act["site"].isin(d.X.index)]
    pmid = d.meta["PMID"].reindex(conds).astype(str).to_numpy()
    Xa = d.X.loc[act["site"].unique(), conds].astype(np.float64)
    selfv = {}
    for k, g in act.groupby("kinase"):
        selfv[k] = Xa.loc[g["site"]].mean(axis=0, skipna=True).to_numpy()   # NaN if none measured
    rows = []
    rng = np.random.default_rng(0)
    for m in METHODS:
        Sm = S[m].reindex(index=conds)
        for k, v in selfv.items():
            s = Sm[k].to_numpy() if k in Sm.columns else np.full(len(conds), np.nan)
            ok = np.isfinite(v) & np.isfinite(s)
            if ok.sum() < MIN_N or len(set(pmid[ok])) < 2 or np.ptp(s[ok]) == 0:
                continue
            rho, p = spearman_perm(s[ok], v[ok], pmid[ok], rng)
            rows.append(dict(method=m, kinase=k, n=int(ok.sum()), n_studies=len(set(pmid[ok])), rho=rho, p=p))
    # specificity control: same self-site vector vs the score of 5 random OTHER scored kinases (no permutation test)
    ctrl = {}
    for m in METHODS:
        Sm = S[m].reindex(index=conds)
        scored = [k for k in Sm.columns if Sm[k].notna().sum() >= MIN_N]
        vals = []
        for k, v in selfv.items():
            if not any((r["method"] == m) and (r["kinase"] == k) for r in rows):
                continue
            for j in rng.choice([x for x in scored if x != k], 5, replace=False):
                s = Sm[j].to_numpy(); ok = np.isfinite(v) & np.isfinite(s)
                if ok.sum() >= MIN_N and np.ptp(s[ok]) > 0:
                    vals.append(pd.Series(s[ok]).corr(pd.Series(v[ok]), method="spearman"))
        ctrl[m] = float(np.mean(vals))
    T = pd.DataFrame(rows)
    T["q"] = T.groupby("method")["p"].transform(bh)
    T.to_csv(OUT / "gta_per_kinase.tsv", sep="\t", index=False)
    common = set.intersection(*[set(T.loc[T.method == m, "kinase"]) for m in METHODS])
    c12 = set.intersection(*[set(T.loc[T.method == m, "kinase"]) for m in ["core-v1", "core-v2"]])
    res = {"n_conditions": len(conds), "n_kinases_with_activating_site": len(selfv),
           "n_common": len(common), "n_perm": N_PERM}
    for m in METHODS:
        t = T[T.method == m]; tc = t[t.kinase.isin(common)]
        res[m] = {"n_kinases": int(len(t)), "macro_rho": float(t.rho.mean()), "median_rho": float(t.rho.median()),
                  "frac_rho_pos": float((t.rho > 0).mean()), "n_sig_q10": int((t.q <= 0.10).sum()),
                  "macro_rho_common": float(tc.rho.mean()),
                  "n_sig_common_q10": int(((tc.q <= 0.10)).sum()),
                  "ctrl_random_other_kinase_macro_rho": ctrl[m]}
        print(m, res[m], flush=True)
    # paired v2 - v1 per kinase
    a = T[T.method == "core-v2"].set_index("kinase").rho; b = T[T.method == "core-v1"].set_index("kinase").rho
    dd = (a - b).dropna()
    res["v2_minus_v1"] = {"mean": float(dd.mean()), "n_up": int((dd > 0).sum()), "n_down": int((dd < 0).sum()),
                          "n_equal": int((dd == 0).sum())}
    (OUT / "gta.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
