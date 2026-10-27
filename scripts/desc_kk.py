#!/usr/bin/env python
"""Descriptive kinase-kinase relationship analysis (DESIGN §5, §12 B15). Not used for model selection.

- Activity profiles: scores of all qc-ok DEV conditions (core-v1 / core-v2 from desc_common; PSP-z and curated-z as
  references). Kinases with >= 100 scored conditions. Spearman correlation over pairwise-complete conditions
  (>= 100 shared).
- Gold (unordered pairs): PSP kinase->substrate edges (psp_ks_all, any organism, human-mapped site) A -> site on
  kinase B (A != B, both in U) where the site is a PSP regulatory site of B whose ON_FUNCTION mentions
  '(enzymatic) activity' (induced / inhibited / altered).
- Non-gold: every other profiled pair not linked by ANY PSP KS edge between the two kinases.
- AUROC of |rho| gold vs non-gold, overall and stratified by substrate-library Jaccard of the two kinases (core-v1
  library: curated_union + r2np top-100; strata 0 | (0,.05] | (.05,.2] | >.2), and within non-same-family pairs.
- Top 20 novel pairs: |rho| ranking under core-v2 among pairs with Jaccard = 0, different Manning family, not gold and
  no PSP KS edge in either direction.
"""
import json

import numpy as np
import pandas as pd

from desc_common import OUT, ROOT, Data, qc_ok_conditions, score_all
from kinact.external_eval import auroc

MIN_C = 100
METHODS = ["core-v1", "core-v2", "psp-z-min3", "cur-z-min1"]
STRATA = [(-1, 0, "J=0"), (0, 0.05, "0<J<=.05"), (0.05, 0.2, ".05<J<=.2"), (0.2, 1.01, "J>.2")]


def spearman_matrix(S):
    """Pairwise-complete Spearman (ranks recomputed per pair would be exact; here ranks per kinase over its own scored
    conditions, then Pearson over shared conditions - the standard pairwise approximation)."""
    R = S.rank(axis=0)
    A = R.to_numpy(np.float64); M = np.isfinite(A); A0 = np.where(M, A, 0.0); Mf = M.astype(np.float64)
    n = Mf.T @ Mf
    sx = A0.T @ Mf; sxx = (A0 ** 2).T @ Mf; sxy = A0.T @ A0
    with np.errstate(all="ignore"):
        cov = sxy - sx * sx.T / n
        vx = sxx - sx ** 2 / n
        r = cov / np.sqrt(vx * vx.T)
    r[n < MIN_C] = np.nan
    return pd.DataFrame(r, index=S.columns, columns=S.columns), n


def gold_pairs(d):
    ks = pd.read_parquet(ROOT / "cache/psp_ks_all.parquet", columns=["kinase", "site_key_psp", "hs_acc"])
    reg = pd.read_parquet(ROOT / "cache/psp_regulatory.parquet", columns=["site_key_psp", "hs_acc", "ON_FUNCTION"])
    u2k = d.Udf.dropna(subset=["uniprot"]).set_index("uniprot")["kinase"]
    U = set(d.kin)
    ks = ks[ks["kinase"].isin(U) & ks["hs_acc"].isin(u2k.index)].assign(sub=lambda t: t["hs_acc"].map(u2k))
    ks = ks[ks["kinase"] != ks["sub"]]
    anyedge = {frozenset(p) for p in zip(ks["kinase"], ks["sub"])}
    regact = set(reg.loc[reg["ON_FUNCTION"].fillna("").str.contains("activity") & reg["hs_acc"].isin(u2k.index),
                         "site_key_psp"].dropna())
    g = ks[ks["site_key_psp"].isin(regact)]
    return {frozenset(p) for p in zip(g["kinase"], g["sub"])}, anyedge


def jaccard(d, kins):
    e = pd.read_csv(ROOT / "models/core_v1/edges.tsv.gz", sep="\t")
    e = e[e["kinase"].isin(kins)]
    ki = pd.Index(kins); si = pd.Index(e["site"].unique())
    import scipy.sparse as sp
    B = sp.csr_matrix((np.ones(len(e)), (si.get_indexer(e["site"]), ki.get_indexer(e["kinase"]))),
                      shape=(len(si), len(ki)))
    B.data[:] = 1; B.sum_duplicates(); B.data[:] = 1
    inter = (B.T @ B).toarray(); sz = np.diag(inter)
    with np.errstate(all="ignore"):
        J = inter / (sz[:, None] + sz[None, :] - inter)
    return pd.DataFrame(np.nan_to_num(J), index=kins, columns=kins)


def main():
    d = Data()
    conds = qc_ok_conditions(d)
    S = score_all(d, masked=False)
    gold, anyedge = gold_pairs(d)
    fam = d.Udf.set_index("kinase")["manning_family"].fillna("").to_dict()
    res = {"n_conditions": len(conds), "n_gold_pairs_U": len(gold)}
    tables = {}
    for m in METHODS:
        Sm = S[m].reindex(index=conds)
        Sm = Sm.loc[:, Sm.notna().sum() >= MIN_C]
        R, n = spearman_matrix(Sm)
        kins = list(R.index)
        J = jaccard(d, kins).to_numpy()
        iu = np.triu_indices(len(kins), 1)
        t = pd.DataFrame({"a": np.asarray(kins)[iu[0]], "b": np.asarray(kins)[iu[1]], "rho": R.to_numpy()[iu],
                          "n": n[iu], "J": J[iu]}).dropna(subset=["rho"])
        t["gold"] = [frozenset(p) in gold for p in zip(t.a, t.b)]
        t["anyedge"] = [frozenset(p) in anyedge for p in zip(t.a, t.b)]
        t["same_fam"] = [fam.get(a, "") != "" and fam.get(a) == fam.get(b) for a, b in zip(t.a, t.b)]
        t = t[t.gold | ~t.anyedge]
        tables[m] = t
        r = {"n_kinases": len(kins), "n_pairs": int(len(t)), "n_gold": int(t.gold.sum()),
             "auroc_abs": auroc(t.rho[t.gold].abs(), t.rho[~t.gold].abs()),
             "auroc_signed": auroc(t.rho[t.gold], t.rho[~t.gold]),
             "mean_abs_rho_gold": float(t.rho[t.gold].abs().mean()), "mean_abs_rho_other": float(t.rho[~t.gold].abs().mean())}
        nf = t[~t.same_fam]
        r["auroc_abs_diff_family"] = auroc(nf.rho[nf.gold].abs(), nf.rho[~nf.gold].abs()); r["n_gold_diff_family"] = int(nf.gold.sum())
        for lo, hi, lab in STRATA:
            s = t[(t.J > lo) & (t.J <= hi)] if lo >= 0 else t[t.J == 0]
            r[f"stratum {lab}"] = {"n_gold": int(s.gold.sum()), "n_other": int((~s.gold).sum()),
                                   "auroc_abs": auroc(s.rho[s.gold].abs(), s.rho[~s.gold].abs()) if s.gold.any() else None}
        res[m] = r
        print(m, json.dumps(r), flush=True)
    # common pair set for core-v1 vs core-v2 (identical kinase sets expected)
    a, b = tables["core-v1"], tables["core-v2"]
    m2 = a.merge(b, on=["a", "b"], suffixes=("_v1", "_v2"))
    res["core v1-v2 rho agreement (pearson over pairs)"] = float(m2.rho_v1.corr(m2.rho_v2))
    # top-20 novel
    t = m2[(m2.J_v2 == 0) & ~m2.same_fam_v2 & ~m2.gold_v2 & ~m2.anyedge_v2].copy()
    t["absr"] = t.rho_v2.abs()
    top = t.sort_values("absr", ascending=False).head(20)
    top = top.assign(fam_a=top.a.map(fam), fam_b=top.b.map(fam))[["a", "fam_a", "b", "fam_b", "rho_v2", "rho_v1", "n_v2"]]
    top.to_csv(OUT / "kk_top20_novel.tsv", sep="\t", index=False)
    res["n_novel_candidates"] = int(len(t))
    pd.concat([x.assign(method=k) for k, x in tables.items()]).to_parquet(ROOT / "cache/descriptive/kk_pairs.parquet")
    (OUT / "kk.json").write_text(json.dumps(res, indent=1, default=float))
    print(top.round(3).to_string())


if __name__ == "__main__":
    main()
