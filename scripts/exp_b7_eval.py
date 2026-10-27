#!/usr/bin/env python
"""B7 evaluation: nested selection per variant, OOF metrics, paired group bootstrap vs core-v1 (alpha=0 of the same
pipeline == core-v1 z_vif nested, checked), T-seen / T-unseen / no-signature splits, shuffled- and generic-signature
controls, calibration (DEV nulls + site-permuted GT conditions), signature examples, drug-disjointness check, ledger.
Writes results/b7/b7_eval.json and results/b7/b7_tables.md.
"""
import itertools
import json
import pickle

import numpy as np
import pandas as pd

from exp_b7_common import (B7Ctx, CACHE, FOLDS, OUT, ALPHAS, NTOPS, KAPPAS, VARIANTS, PRIMARY, BRANCH, frame,
                           config_dict)
from exp_b21_common import calib_bins
from audit_core_v1_dev import summary, t_unseen_mask
from kinact import bootstrap as bs, ledger
from kinact.eval import to_jsonable, read_gt_tsv
from kinact.methods.linear import zscore

N_SHUFFLE = 3


def m1sel(ev, S, keep_folds, fold_of):
    P = ev.pairs(S)
    f = P["conditionID"].map(fold_of)
    return ev.metric(P[f.isin(keep_folds)], "M1-sel")["value"]


def grid_table(ctx, variant, store):
    c = ctx.c
    cfgs = [(0.0, 0, 0.0)] + [(a, nt, kp) for kp, nt, a in itertools.product(KAPPAS, NTOPS, ALPHAS)]
    rows, outer = [], {}
    for a, nt, kp in cfgs:
        get = (lambda fo, g: store[("core", fo, g)][0]) if a == 0 else \
              (lambda fo, g: store[(variant, fo, g, kp, nt, a)][0])
        r = {"alpha": a, "n_top": nt, "kappa": kp}
        for f in FOLDS:
            parts = {g: get(f, g) for g in FOLDS if g != f}
            r[f"inner_{f}"] = m1sel(c.ev, frame(ctx, parts), [g for g in FOLDS if g != f], c.fold_of)
        S = frame(ctx, {f: get(f, f) for f in FOLDS})
        outer[(a, nt, kp)] = S
        P = c.ev.pairs(S)
        r["oof_M1-sel"] = c.ev.metric(P, "M1-sel")["value"]
        r["oof_M1-all"] = c.ev.metric(P, "M1-all")["value"]
        rows.append(r)
        if a:
            ledger.append("bayes2+sig", config_dict(variant, a, nt, kp), "dev",
                          {"M1-sel": r["oof_M1-sel"], "M1-all": r["oof_M1-all"],
                           **{k: v for k, v in r.items() if k.startswith("inner_")},
                           "note": "per-config OOF; signatures fitted on the 4 other folds"},
                          branch=BRANCH, folds=range(5))
    return pd.DataFrame(rows), outer


def nested(ctx, g, outer):
    parts, chosen = {}, {}
    for f in FOLDS:
        b = g.loc[g[f"inner_{f}"].idxmax()]
        key = (float(b.alpha), int(b.n_top), float(b.kappa))
        chosen[f] = dict(alpha=key[0], n_top=key[1], kappa=key[2], inner=float(b[f"inner_{f}"]))
        parts[f] = outer[key].to_numpy()[ctx.fold_pos[f]]
    return frame(ctx, parts), chosen


def delta(ev, A, B, ep, plan, mask=None):
    if mask is None:
        d = bs.paired_delta(ev, A, B, ep, plan)
        return {k: d[k] for k in ["estimate", "ci_low", "ci_high", "p_one_sided"]}
    Pa, Pb = ev.pairs(A), ev.pairs(B)
    ma, mb = mask(Pa), mask(Pb)
    assert (ma == mb).all()
    ra = bs.metric_replicates(ev, Pa[ma], ep, plan); rb = bs.metric_replicates(ev, Pb[mb], ep, plan)
    pt = ev.metric(Pa[ma], ep)["value"] - ev.metric(Pb[mb], ep)["value"]
    s = bs._summ(pt, ra - rb)
    return {"estimate": pt, "ci_low": s["ci_low"], "ci_high": s["ci_high"], "p_one_sided": bs.shifted_null_p(pt, ra - rb),
            "n_pairs": int((ma & Pa["in_M1sel"]).sum())}


def sig_kinases(variant):
    d = json.loads((OUT / f"fits_{variant}.json").read_text())
    return {f: set(d[f]["kinases"]) for f in FOLDS}


def control_scores(ctx, variant, chosen, kind, seed=0, which="gt"):
    """Fold-wise chosen config; signatures replaced by a control (shuffled / generic / none=real)."""
    c = ctx.c
    parts, rhos, ns = {}, {}, {}
    for f in FOLDS:
        ch = chosen[f]
        idx = ctx.fold_pos[f]
        if ch["alpha"] == 0:
            res = ctx.score(idx, [("x", None)], ctx.picks[f], [], which=which, return_n=True)
            Z, rho, n = res[("x", 0.0)]
        else:
            fits = pickle.load(open(CACHE / f"sig_{variant}_{f}.pkl", "rb"))
            s = fits[ch["kappa"]]
            s = {"real": s, "shuffled": s.shuffled(seed) if kind == "shuffled" else None,
                 "generic": s.generic_control() if kind == "generic" else None}[kind]
            Ws = s.matrix(len(c.X.index), c.kin, ch["n_top"])
            Z, rho, n = ctx.score(idx, [("x", Ws)], ctx.picks[f], [ch["alpha"]], which=which, return_n=True)[("x", ch["alpha"])]
        parts[f], rhos[f], ns[f] = Z, rho, n
    return frame(ctx, parts), rhos, ns


def calibration(ctx, variant, chosen):
    """|z|>1.96: DEV nulls (null in fold f scored with the fold-f outer design) and site-permuted GT conds."""
    c = ctx.c
    Zs, Ns = [], []
    for f in FOLDS:
        ni = np.flatnonzero(ctx.null_fold == f)
        if not len(ni):
            continue
        ch = chosen[f]
        if ch["alpha"] == 0:
            Wl, al = [("x", None)], []
        else:
            fits = pickle.load(open(CACHE / f"sig_{variant}_{f}.pkl", "rb"))
            Wl, al = [("x", fits[ch["kappa"]].matrix(len(c.X.index), c.kin, ch["n_top"]))], [ch["alpha"]]
        # score the fold-f nulls as the "evaluated" set: reuse score() with which='null'
        res = ctx.score(ni, Wl, ctx.picks[f], al, which="null", return_n=True)
        Z, _, n = res[("x", ch["alpha"] if ch["alpha"] else 0.0)]
        Zs.append(Z); Ns.append(n)
    out = {"null": calib_bins(np.vstack(Zs), np.vstack(Ns)), "null_ncells": int(np.isfinite(np.vstack(Zs)).sum())}
    S, _, ns = control_scores(ctx, variant, chosen, "real", which="perm")
    N = np.full(S.shape, np.nan)
    for f in FOLDS:
        N[ctx.fold_pos[f]] = ns[f]
    out["perm"] = calib_bins(S.to_numpy(), N)
    return out


def gene_map():
    h = pd.read_parquet("/workspace/cache/hgnc_symbols.parquet")
    h = h[h["kind"] == "symbol"].dropna(subset=["uniprot_ids"])
    m = {}
    for s, u in zip(h["symbol"], h["uniprot_ids"]):
        for x in str(u).replace(";", ",").split(","):
            m.setdefault(x.strip(), s)
    return m


def examples(ctx, variant, chosen, kinases=("MTOR", "MAPK1", "MAP2K1", "EGFR", "GSK3B", "PRKACA", "CDK7", "MET"), n=8):
    c = ctx.c
    gm = gene_map()
    psp = pd.read_parquet("/workspace/cache/psp_ks_human.parquet", columns=["kinase", "site_key_psp"]).dropna()
    psp_k = psp.groupby("site_key_psp")["kinase"].apply(set).to_dict()
    f = "fold0"
    kp = chosen[f]["kappa"] if chosen[f]["alpha"] else 1.0
    fits = pickle.load(open(CACHE / f"sig_{variant}_{f}.pkl", "rb"))[kp]
    ex = {}
    for k in kinases:
        if k not in fits.sites:
            continue
        rows = []
        for si, w in zip(fits.sites[k][:n], fits.weights[k][:n]):
            sk = c.X.index[si]; acc, res = sk.split("_", 1)
            ks = psp_k.get(sk, set())
            rows.append(f"{gm.get(acc, acc)} {res} ({w:+.2f}{'; PSP:' + ','.join(sorted(ks)[:3]) if ks else ''})")
        ex[k] = {"n_groups": fits.n_groups[k], "top": rows}
    # how many top-50 signature sites are PSP substrates of k / of any kinase; overlap between kinases
    stats = {}
    for k, s in fits.sites.items():
        keys = c.X.index[s[:50]]
        stats[k] = {"psp_of_k": int(sum(k in psp_k.get(x, set()) for x in keys)),
                    "psp_any": int(sum(x in psp_k for x in keys))}
    def jac(a, b, n=50):
        A, B = set(fits.sites[a][:n]), set(fits.sites[b][:n]); return len(A & B) / len(A | B)
    pairs = {f"{a}-{b}": jac(a, b) for a, b in [("MAPK1", "MAP2K1"), ("MAPK1", "MAPK3"), ("MAP2K1", "MAP2K2"),
                                                ("MTOR", "RPS6KB1"), ("MTOR", "MAPK1"), ("EGFR", "MAPK1")]
             if a in fits.sites and b in fits.sites}
    ks_ = sorted(fits.sites)
    allj = [jac(a, b) for a, b in itertools.combinations(ks_, 2)]
    psp_rand = float(np.mean([x in psp_k for x in c.X.index[fits.eligible_sites]]))
    return {"fold": f, "kappa": kp, "examples": ex, "psp_stats_top50": stats, "jaccard_top50": pairs,
            "jaccard_all_pairs_mean": float(np.mean(allj)) if allj else None,
            "psp_any_rate_eligible_sites": psp_rand, "n_kinases": len(ks_)}


def drug_disjointness(ctx):
    cg = read_gt_tsv("/workspace/results/ground_truth/condition_groups.tsv").set_index("conditionID")
    m = ctx.c.meta
    out = {}
    for f in FOLDS:
        ev_c = [x for x in m.index[m["split"] == f]]
        tr_c = [x for x in m.index[m["split"].isin([g for g in FOLDS if g != f])]]
        dr = lambda cs: set(d for x in cs for d in str(cg.loc[x, "drugs"]).split(";") if d and d != "nan") if len(cs) else set()
        de, dt = dr(ev_c), dr(tr_c)
        out[f] = {"n_drugs_eval": len(de), "shared_drugs": sorted(de & dt),
                  "shared_groups": len(set(m.loc[ev_c, "group_id"]) & set(m.loc[tr_c, "group_id"]))}
    return out


def main():
    ctx = B7Ctx()
    c = ctx.c
    psp = zscore(c.X[c.conds], "psp_all", min_substrates=3).reindex(columns=c.kin)
    plan = bs.BootstrapPlan.from_evaluator(c.ev, B=2000, seed=0)
    rep = {"gt_hash": ctx.gt_hash, "core_picks": ctx.picks, "drug_disjointness": drug_disjointness(ctx)}
    S_core = None
    res, grids, chosen_all, S_nested = {}, {}, {}, {}
    for v in VARIANTS:
        store = pickle.load(open(CACHE / f"grid_{v}.pkl", "rb"))
        g, outer = grid_table(ctx, v, store)
        g.to_csv(OUT / f"grid_{v}.tsv", sep="\t", index=False)
        grids[v] = g
        if S_core is None:
            S_core = outer[(0.0, 0, 0.0)]
        S, ch = nested(ctx, g, outer)
        S_nested[v], chosen_all[v] = S, ch
        best = g.loc[g["oof_M1-sel"].idxmax()]
        print(v, "nested", round(m1sel(c.ev, S, FOLDS, c.fold_of), 4), ch, "pooled-best", dict(best), flush=True)
        del store
    # ---- reference
    ref = summary(c, S_core, psp)
    rep["core_v1_reproduced_M1-sel"] = ref["M1-sel"]
    res["core-v1 (alpha=0)"] = ref
    print("core-v1", ref["M1-sel"], ref["T-unseen"], flush=True)

    def masks_for(variant):
        sk = sig_kinases(variant)
        has = lambda P: np.array([k in sk[c.fold_of[cid]] for k, cid in zip(P["kinase"], P["conditionID"])])
        uns = lambda P: t_unseen_mask(c, P)
        return {"T-unseen": uns, "T-seen": lambda P: ~uns(P), "has-signature": has, "no-signature": lambda P: ~has(P)}

    for v in VARIANTS:
        S = S_nested[v]
        r = summary(c, S, psp)
        r["chosen"] = chosen_all[v]
        r["M2_vs_core"] = c.ev.m2(S, S_core)
        r["delta_vs_core"] = {ep: delta(c.ev, S, S_core, ep, plan) for ep in ["M1-sel", "M1-all", "harmony", "M2"]}
        for name, mk in masks_for(v).items():
            r["delta_vs_core"][name] = delta(c.ev, S, S_core, "M1-sel", plan, mk)
            r[name] = c.ev.metric(c.ev.pairs(S)[mk(c.ev.pairs(S))], "M1-sel")["value"]
        for src in ["H", "D", "O"]:
            r["delta_vs_core"][f"source={src}"] = delta(c.ev, S, S_core, f"M1-sel:source={src}", plan)
        for sg in ["-1", "1"]:
            r["delta_vs_core"][f"sign={sg}"] = delta(c.ev, S, S_core, f"M1-sel:sign={sg}", plan)
        for b in ["0", "1-4", "5-19", ">=20"]:
            try:
                r["delta_vs_core"][f"bin={b}"] = delta(c.ev, S, S_core, f"M1-sel:bin={b}", plan)
            except Exception as e:  # noqa
                r["delta_vs_core"][f"bin={b}"] = str(e)
        res[f"B7 {v} (nested)"] = r
        ledger.append("bayes2+sig nested", {"variant": v, "grid": "alpha x n_top x kappa + core", "chosen": chosen_all[v]},
                      "dev", {k: r[k] for k in ["M1-sel", "M1-all", "M1-sel-scored", "coverage_M1-sel", "T-unseen"]},
                      branch=BRANCH, folds=range(5))
        print(v, {k: r[k] for k in ["M1-sel", "M1-all", "T-unseen", "harmony"]}, r["delta_vs_core"]["M1-sel"], flush=True)

    # ---- controls + calibration (primary and the best other variant)
    best_other = max([v for v in VARIANTS if v != PRIMARY], key=lambda v: res[f"B7 {v} (nested)"]["M1-sel"])
    for v in [PRIMARY, best_other]:
        ch = chosen_all[v]
        ctl = {}
        Sreal, _, _ = control_scores(ctx, v, ch, "real")
        ctl["real_rescored_M1-sel"] = m1sel(c.ev, Sreal, FOLDS, c.fold_of)
        sh = []
        for seed in range(N_SHUFFLE):
            Ssh, _, _ = control_scores(ctx, v, ch, "shuffled", seed)
            sh.append(m1sel(c.ev, Ssh, FOLDS, c.fold_of))
            if seed == 0:
                ctl["shuffled_seed0_vs_core"] = delta(c.ev, Ssh, S_core, "M1-sel", plan)
                ctl["shuffled_seed0_T-unseen"] = summary(c, Ssh)["T-unseen"]
        ctl["shuffled_M1-sel"] = sh
        Sg, _, _ = control_scores(ctx, v, ch, "generic")
        ctl["generic_M1-sel"] = m1sel(c.ev, Sg, FOLDS, c.fold_of)
        ctl["generic_vs_core"] = delta(c.ev, Sg, S_core, "M1-sel", plan)
        res[f"B7 {v} (nested)"]["controls"] = ctl
        res[f"B7 {v} (nested)"]["calibration"] = calibration(ctx, v, ch)
        print(v, "controls", ctl, res[f"B7 {v} (nested)"]["calibration"], flush=True)
    # core calibration via the same code path
    ch0 = {f: {"alpha": 0.0, "n_top": 0, "kappa": 0.0} for f in FOLDS}
    res["core-v1 (alpha=0)"]["calibration"] = calibration(ctx, PRIMARY, ch0)

    # ---- fixed-config sensitivity (primary): pooled grid of OOF M1-sel, top rows
    rep["grid_top"] = {v: grids[v].sort_values("oof_M1-sel", ascending=False).head(5).to_dict("records") for v in VARIANTS}
    rep["fits"] = {v: {k: {"n_kinases": d["n_kinases"], "n_train_conds": d["n_train_conds"]}
                       for k, d in json.loads((OUT / f"fits_{v}.json").read_text()).items()} for v in VARIANTS}
    rep["examples"] = {v: examples(ctx, v, chosen_all[v]) for v in [PRIMARY, "all+contrast"]}
    rep["results"] = res
    (OUT / "b7_eval.json").write_text(json.dumps(to_jsonable(rep), indent=1, default=str))
    print("done", flush=True)


if __name__ == "__main__":
    main()
