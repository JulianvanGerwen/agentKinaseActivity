#!/usr/bin/env python
"""Independent audit of core-v1 (D6) on DEV folds 0-4 (docs/audit/core_v1_audit.md, part 1).

1. Re-scores the core-v1 search space (curated_union + r2np@100 alpha=1, no winsor, wsd scale, borrow; tau x rho_f
   (<=3) x rho_g) with kinact.explore.bayes2, recomputes inner-fold M1-sel with an independent implementation of the
   nested selection (and checks it against results/b21/grid_out.tsv), assembles OOF, evaluates with the frozen
   harness.
2. Leakage checks: (a) the fold-f pick is invariant to arbitrary corruption of fold-f scores; (b) rho_e uses only
   DEV role=null conditions, and a leave-outer-fold-out rho_e gives the same OOF; (c) z_vif contains no label input.
3. B2 reference (cache/b2 archives) re-evaluated; paired bootstrap of M1-sel and T-unseen deltas.
4. Independent minimal numpy reimplementation of the scorer on 50 GT conditions (Spearman vs bayes2).
No ledger rows are written by this script (it re-evaluates already-ledgered configurations; the audit row is
appended once at the end, see main()).
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path("/workspace")
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))
from exp_b21_common import Ctx, FOLDS, B2, b2_reference, calib_bins  # noqa: E402
import exp_b21_grid as G  # noqa: E402
from kinact import bootstrap as bs, ledger  # noqa: E402
from kinact.eval import to_jsonable  # noqa: E402
from kinact.library import load_library  # noqa: E402
from kinact.methods.linear import zscore  # noqa: E402

OUT = ROOT / "results/audit_core_v1"; OUT.mkdir(parents=True, exist_ok=True)
TAUS = [0.06, 0.125, 0.25, 0.5, 1.0]; RFS = [0.0, 1.0, 3.0]; RGS = [0.0, 0.25, 1.0]
OUTS = ["z_post", "z_vif", "mean", "z_freq"]


def cfg_of(tau, rf, rg):
    return dict(winsor=None, alphas={"r2np@100": 1.0}, tau=tau, rho_f=rf, rho_g=rg, scale="wsd")


def my_inner(ev, S, fold_of):
    """Independent inner metric: drop every pair whose condition is in fold g, M1-sel on the rest."""
    P = ev.pairs(S)
    f = P["conditionID"].map(fold_of)
    return {g: ev.metric(P[f != g], "M1-sel")["value"] for g in FOLDS}, P


def t_unseen_mask(ctx, P):
    C, Pg = ctx.C, ctx.P
    gc = C[C["split"].isin(FOLDS) & (C["role"] == "gt")].set_index("conditionID")["split"]
    pos = Pg[(Pg["cls"] == "direct") & Pg["confidence"].isin(["high", "medium"]) & Pg["conditionID"].isin(gc.index)]
    pos = pos.assign(fold=pos["conditionID"].map(gc))
    seen = {f: set(pos.loc[pos["fold"] != f, "kinase"]) for f in FOLDS}
    f = P["conditionID"].map(ctx.fold_of)
    return np.array([k not in seen[ff] for k, ff in zip(P["kinase"], f)])


def summary(ctx, S, psp=None):
    r = ctx.ev.evaluate(S, n_perm=2000)
    v = lambda x: x["value"] if isinstance(x, dict) and "value" in x else x
    P = r["pairs"]
    out = {"M1-sel": v(r["M1-sel"]), "M1-all": v(r["M1-all"]), "M1-sel-scored": v(r["M1-sel-scored"]),
           "coverage_M1-sel": r["coverage_M1-sel"], "harmony": v(r["harmony"]), "harmony_per_bin": r["harmony"]["per_bin"],
           "by_source": {k: v(x) for k, x in r["by_source"].items()},
           "by_sign": {str(k): v(x) for k, x in r["by_sign"].items()},
           "n_sig": r["coverage_significant"]["n_significant"],
           "T-unseen": ctx.ev.metric(P[t_unseen_mask(ctx, P)], "M1-sel")["value"],
           "T-unseen_n_pairs": int((t_unseen_mask(ctx, P) & P["in_M1sel"]).sum())}
    if psp is not None:
        out["M2_vs_pspz"] = ctx.ev.m2(S, psp)
    return out


# ------------------------------------------------------------------ independent minimal scorer (numpy, CPU)
def indep_scores(X, conds, kin, Udf, tau, rf, rg, rho_e):
    """From scratch: library -> W; Manning kernel; per condition dense posterior; z_post and z_vif."""
    cur = load_library("curated_union"); r2 = load_library("r2np", k=100)
    sidx = pd.Index(X.index); kidx = pd.Index(kin)
    W = np.zeros((len(sidx), len(kidx)))
    ci, ki = sidx.get_indexer(cur["site"]), kidx.get_indexer(cur["kinase"])
    ok = (ci >= 0) & (ki >= 0)
    for s, k, w in zip(ci[ok], ki[ok], cur["weight"].to_numpy(float)[ok]):
        W[s, k] = max(W[s, k], w)
    curmask = W != 0
    ci, ki = sidx.get_indexer(r2["site"]), kidx.get_indexer(r2["kinase"])
    ok = (ci >= 0) & (ki >= 0)
    for s, k in zip(ci[ok], ki[ok]):
        if not curmask[s, k]:
            W[s, k] = 1.0
    rows = np.flatnonzero((W != 0).any(1)); W = W[rows]
    U = Udf.set_index("kinase").reindex(kin)
    fam = U["manning_family"].fillna("").astype(str).values; grp = U["manning_group"].fillna("").astype(str).values
    K = np.eye(len(kin))
    for i in range(len(kin)):
        for j in range(len(kin)):
            if fam[i] and fam[i] == fam[j]:
                K[i, j] += rf
            if grp[i] and grp[i] == grp[j] and grp[i] not in ("Other", "Atypical"):
                K[i, j] += rg
    d = np.sqrt(np.diag(K)); K = K / np.outer(d, d)
    P0 = np.linalg.inv(K) / tau ** 2
    rel = (K > 1e-12) & ~np.eye(len(kin), dtype=bool)
    zp, zv = [], []
    for c in conds:
        x = X[c].to_numpy(float)
        m_all = ~np.isnan(x)
        r = x[m_all] - np.median(x[m_all])
        y = np.full(len(x), np.nan); y[m_all] = r / np.sqrt(np.mean(r ** 2))   # wsd without winsor = RMS scale
        yl = y[rows]; ml = ~np.isnan(yl)
        Wm = W[ml]
        G_ = Wm.T @ Wm
        L = G_ + P0
        Li = np.linalg.inv(L)
        mu = Li @ (Wm.T @ yl[ml])
        A = Li @ Wm.T                                   # p x m
        q = (A * A).sum(1)
        suma = A.sum(1)
        me = np.where(q > 0, suma ** 2 / q, 0)
        has = (Wm != 0).sum(0) >= 1
        rep = has | ((has.astype(float) @ rel.astype(float)) > 0)
        zp.append(np.where(rep, mu / np.sqrt(np.diag(Li)), np.nan))
        zv.append(np.where(rep, mu / np.sqrt(q * (1 + rho_e * np.maximum(me - 1, 0))), np.nan))
    return np.array(zp), np.array(zv)


def main():
    ctx = Ctx(ks={"r2np": [100]})
    rep = {"nulls": ctx.nulls, "null_roles": ctx.C.set_index("conditionID").loc[ctx.nulls, ["role", "split"]].to_dict("index")}
    assert all(ctx.C.set_index("conditionID").loc[c, "role"] == "null" for c in ctx.nulls)
    assert len(set(ctx.nulls) & set(ctx.conds)) == 0
    psp = zscore(ctx.X[ctx.conds], "psp_all", min_substrates=3).reindex(columns=ctx.kin)

    # ---- 1. re-score the core-v1 space
    store, rows = {}, []
    rho_cache = {}
    for tau in TAUS:
        for rf in RFS:
            for rg in RGS:
                cfg = cfg_of(tau, rf, rg)
                r = G.compute(ctx, cfg, "gt")
                rho = G.fit_rho(ctx, cfg, None); rho_cache[(tau, rf, rg)] = rho
                for o in OUTS:
                    arr = G.outputs_from(r, rho, [(o, None)])[o]
                    S = pd.DataFrame(arr, index=ctx.conds, columns=ctx.kin)
                    inner, P = my_inner(ctx.ev, S, ctx.fold_of)
                    store[(tau, rf, rg, o)] = arr
                    rows.append(dict(tau=tau, rho_f=rf, rho_g=rg, output=o, rho_e=rho,
                                     M1sel=ctx.ev.metric(P, "M1-sel")["value"], **{f"in_{g}": v for g, v in inner.items()}))
    g = pd.DataFrame(rows)
    g.to_csv(OUT / "dev_grid_rescored.tsv", sep="\t", index=False)
    # compare with producer grid
    pg = pd.read_csv(ROOT / "results/b21/grid_out.tsv", sep="\t", keep_default_na=False, na_values=[""])
    pg = pg[(pg["winsor"].astype(str).isin(["", "nan", "None"])) & (pg["alphas"] == "{'r2np@100': 1.0}") & pg["output"].isin(OUTS)]
    mg = g.merge(pg, on=["tau", "rho_f", "rho_g", "output"], suffixes=("", "_prod"))
    diffs = {"n_matched": len(mg), "max_abs_M1sel": float((mg["M1sel"] - mg["M1-sel"]).abs().max()),
             "max_abs_inner": float(max((mg[f"in_{f}"] - mg[f"inner_{f}"]).abs().max() for f in FOLDS)),
             "max_abs_rho_e": float((mg["rho_e"] - mg["rho_e_prod"]).abs().max()) if "rho_e_prod" in mg else None}
    rep["grid_vs_producer"] = diffs
    print("grid vs producer", diffs, flush=True)

    # ---- nested (independent)
    def nested(mask, store_=store):
        sub = g[mask]
        parts, chosen = [], {}
        for f in FOLDS:
            b = sub.loc[sub[f"in_{f}"].idxmax()]
            key = (b.tau, b.rho_f, b.rho_g, b.output)
            chosen[f] = dict(tau=b.tau, rho_f=b.rho_f, rho_g=b.rho_g, output=b.output, inner=b[f"in_{f}"])
            idx = ctx.cf == f
            parts.append(pd.DataFrame(store_[key][idx], index=np.asarray(ctx.conds)[idx], columns=ctx.kin))
        return pd.concat(parts).loc[ctx.conds], chosen

    V = {
        "core-v1 rank z_post (tau,rho_f<=3,rho_g nested)": g.output == "z_post",
        "core-v1 z_vif (tau,rho_f<=3,rho_g nested)": g.output == "z_vif",
        "core-v1 z_post fixed tau=.125 rho_f=3 (rho_g nested)": (g.output == "z_post") & (g.tau == 0.125) & (g.rho_f == 3),
        "core-v1 z_vif fixed tau=.125 rho_f=3 (rho_g nested)": (g.output == "z_vif") & (g.tau == 0.125) & (g.rho_f == 3),
        "core-v1 z_vif fixed tau=.125 rho_f=3 rho_g=1": (g.output == "z_vif") & (g.tau == 0.125) & (g.rho_f == 3) & (g.rho_g == 1),
        "core-v1 z_post fixed tau=.125 rho_f=3 rho_g=1": (g.output == "z_post") & (g.tau == 0.125) & (g.rho_f == 3) & (g.rho_g == 1),
        "core-v1 mean (tau,rho_f<=3,rho_g nested)": g.output == "mean",
    }
    S_v, ch_v = {}, {}
    for n, m in V.items():
        S_v[n], ch_v[n] = nested(m)
    S_ref, ch_ref = b2_reference(ctx)
    plan = bs.BootstrapPlan.from_evaluator(ctx.ev, B=2000, seed=0)
    res = {"B2 reference": summary(ctx, S_ref, psp), "PSP z-score min3": summary(ctx, psp)}
    res["B2 reference"]["chosen"] = ch_ref
    Pref = ctx.ev.pairs(S_ref); uref = t_unseen_mask(ctx, Pref)
    for n, S in S_v.items():
        r = summary(ctx, S, psp); r["chosen"] = ch_v[n]
        d = bs.paired_delta(ctx.ev, S, S_ref, "M1-sel", plan)
        r["dM1sel_vs_B2ref"] = {k: d[k] for k in ["estimate", "ci_low", "ci_high", "p_one_sided"]}
        for ep in ["M1-all", "harmony", "M2"]:
            dd = bs.paired_delta(ctx.ev, S, S_ref, ep, plan)
            r[f"d{ep}_vs_B2ref"] = {k: dd[k] for k in ["estimate", "ci_low", "ci_high", "p_one_sided"]}
        # T-unseen paired bootstrap (unseen mask identical for both: depends on kinase & fold only)
        Pa = ctx.ev.pairs(S); ua = t_unseen_mask(ctx, Pa)
        assert (ua == uref).all()
        ra = bs.metric_replicates(ctx.ev, Pa[ua], "M1-sel", plan); rb = bs.metric_replicates(ctx.ev, Pref[uref], "M1-sel", plan)
        pt = ctx.ev.metric(Pa[ua], "M1-sel")["value"] - ctx.ev.metric(Pref[uref], "M1-sel")["value"]
        s_ = bs._summ(pt, ra - rb)
        r["dTunseen_vs_B2ref"] = {"estimate": pt, "ci_low": s_["ci_low"], "ci_high": s_["ci_high"],
                                  "p_one_sided": bs.shifted_null_p(pt, ra - rb)}
        res[n] = r
        print(n, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items() if k in
                  ["M1-sel", "M1-all", "T-unseen", "dM1sel_vs_B2ref", "dTunseen_vs_B2ref"]}, flush=True)

    # decomposition: which folds changed between B2 ref and core-v1 z_post
    S_c = S_v["core-v1 rank z_post (tau,rho_f<=3,rho_g nested)"]
    Pc = ctx.ev.pairs(S_c)
    Pc["fold"] = Pc["conditionID"].map(ctx.fold_of); Pref["fold"] = Pref["conditionID"].map(ctx.fold_of)
    per_fold = {}
    for f in FOLDS:
        a = Pc[Pc.fold == f]; b = Pref[Pref.fold == f]
        per_fold[f] = {"core": ctx.ev.metric(a, "M1-sel")["value"], "b2ref": ctx.ev.metric(b, "M1-sel")["value"],
                       "core_unseen": ctx.ev.metric(a[t_unseen_mask(ctx, a)], "M1-sel")["value"],
                       "b2ref_unseen": ctx.ev.metric(b[t_unseen_mask(ctx, b)], "M1-sel")["value"],
                       "b2ref_rho_f": ch_ref[f]["rho_f"], "core_rho_f": ch_v["core-v1 rank z_post (tau,rho_f<=3,rho_g nested)"][f]["rho_f"]}
    rep["per_fold_core_vs_b2ref"] = per_fold
    # unseen kinases contributing most to the T-unseen change
    ua = t_unseen_mask(ctx, Pc)
    ka = Pc[ua & Pc.in_M1sel].groupby("kinase")["auroc"].mean(); kb = Pref[ua & Pref.in_M1sel].groupby("kinase")["auroc"].mean()
    dk = (ka - kb).sort_values()
    rep["Tunseen_kinase_deltas"] = {"n_kinases": int(len(dk)), "top_gain": dk.tail(6).round(3).to_dict(),
                                    "top_loss": dk.head(4).round(3).to_dict(),
                                    "n_kinases_changed_>0.01": int((dk.abs() > 0.01).sum())}

    # ---- 2a. leakage: corrupt fold-f scores, the pick for f must not change
    leak = {}
    rng = np.random.default_rng(1)
    for f in FOLDS:
        idx = ctx.cf == f
        rows2 = []
        for key, arr in store.items():
            if key[3] != "z_post":
                continue
            a2 = arr.copy(); a2[idx] = rng.standard_normal(a2[idx].shape)
            P2 = ctx.ev.pairs(pd.DataFrame(a2, index=ctx.conds, columns=ctx.kin))
            fo = P2["conditionID"].map(ctx.fold_of)
            rows2.append((ctx.ev.metric(P2[fo != f], "M1-sel")["value"], key))
        best = max(rows2)[1]
        orig = ch_v["core-v1 rank z_post (tau,rho_f<=3,rho_g nested)"][f]
        leak[f] = {"pick_after_corruption": best, "pick_orig": (orig["tau"], orig["rho_f"], orig["rho_g"], "z_post"),
                   "same": best == (orig["tau"], orig["rho_f"], orig["rho_g"], "z_post")}
    rep["leak_corruption_test"] = leak
    print("corruption test", {f: v["same"] for f, v in leak.items()}, flush=True)

    # ---- 2b. rho_e leave-outer-fold-out (nulls in fold f excluded when scoring fold f)
    nf = pd.Series(ctx.nulls).map(ctx.fold_of).to_numpy()
    rep["null_folds"] = dict(zip(ctx.nulls, nf))
    chz = ch_v["core-v1 z_vif (tau,rho_f<=3,rho_g nested)"]
    parts, rlofo = [], {}
    for f in FOLDS:
        c = chz[f]; cfg = cfg_of(c["tau"], c["rho_f"], c["rho_g"])
        rn = G.compute(ctx, cfg, "null")
        keep = nf != f
        zf = B2.output_stat(rn["post"], "z_freq")
        rho_f_, _ = B2.fit_rho_null(zf[keep], rn["me"][keep], rn["rep"][keep])
        rlofo[f] = {"rho_e_all_nulls": rho_cache[(c["tau"], c["rho_f"], c["rho_g"])], "rho_e_lofo": rho_f_}
        r = G.compute(ctx, cfg, "gt")
        arr = G.outputs_from(r, rho_f_, [("z_vif", None)])["z_vif"]
        idx = ctx.cf == f
        parts.append(pd.DataFrame(arr[idx], index=np.asarray(ctx.conds)[idx], columns=ctx.kin))
    S_lofo = pd.concat(parts).loc[ctx.conds]
    rep["rho_e_lofo"] = {"per_fold": rlofo, "M1-sel_lofo": ctx.ev.metric(ctx.ev.pairs(S_lofo), "M1-sel")["value"],
                         "M1-sel_allnulls": res["core-v1 z_vif (tau,rho_f<=3,rho_g nested)"]["M1-sel"]}
    print("rho_e lofo", rep["rho_e_lofo"], flush=True)

    # ---- calibration per bin (fold-wise chosen config on DEV nulls + site-permuted GT conds)
    def calib(ch, out_name):
        cal = {}
        for which, conds in [("null", ctx.nulls), ("perm", ctx.conds)]:
            fo = pd.Series(conds).map(ctx.fold_of).to_numpy()
            Zs, Ns = [], []
            for f in FOLDS:
                idx = np.flatnonzero(fo == f)
                if not len(idx):
                    continue
                c = ch[f]; cfg = cfg_of(c["tau"], c["rho_f"], c["rho_g"])
                r = G.compute(ctx, cfg, which)
                z = G.outputs_from(r, rho_cache[(c["tau"], c["rho_f"], c["rho_g"])], [(out_name, None)])[out_name]
                Zs.append(z[idx]); Ns.append(r["n"][idx])
            cal[which] = calib_bins(np.vstack(Zs), np.vstack(Ns))
            cal[which + "_ncells"] = int(np.isfinite(np.vstack(Zs)).sum())
        return cal
    for n, o in [("core-v1 z_vif (tau,rho_f<=3,rho_g nested)", "z_vif"),
                 ("core-v1 rank z_post (tau,rho_f<=3,rho_g nested)", "z_post"),
                 ("core-v1 z_vif fixed tau=.125 rho_f=3 (rho_g nested)", "z_vif")]:
        res[n]["calibration"] = calib(ch_v[n], o)
        print(n, res[n]["calibration"], flush=True)

    # ---- 4. independent reimplementation, 50 random GT conditions
    rs = np.random.default_rng(7)
    c50 = list(rs.choice(ctx.conds, 50, replace=False))
    tau, rf, rg = 0.125, 3.0, 0.25
    rho = rho_cache[(tau, rf, rg)]
    zp_i, zv_i = indep_scores(ctx.X, c50, ctx.kin, ctx.Udf, tau, rf, rg, rho)
    pos = [ctx.conds.index(c) for c in c50]
    zp_b = store[(tau, rf, rg, "z_post")][pos]; zv_b = store[(tau, rf, rg, "z_vif")][pos]
    def cmp(a, b):
        same_nan = bool((np.isnan(a) == np.isnan(b)).all())
        ok = ~np.isnan(a) & ~np.isnan(b)
        per = [stats.spearmanr(a[i, ok[i]], b[i, ok[i]])[0] for i in range(a.shape[0])]
        return {"same_scored_set": same_nan, "spearman_pooled": float(stats.spearmanr(a[ok], b[ok])[0]),
                "spearman_min_per_cond": float(np.min(per)), "max_abs_diff": float(np.nanmax(np.abs(a - b)))}
    rep["reimpl_50"] = {"z_post": cmp(zp_i, zp_b), "z_vif": cmp(zv_i, zv_b), "cfg": [tau, rf, rg, rho]}
    print("reimpl", rep["reimpl_50"], flush=True)

    rep["results"] = res
    (OUT / "dev_audit.json").write_text(json.dumps(to_jsonable(rep), indent=1, default=str))
    # pooled (all-DEV) grid for the frozen devconfirm config
    g.to_csv(OUT / "dev_grid_rescored.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
