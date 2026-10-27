#!/usr/bin/env python
"""G2 prep, part 1b: winner's-curse / fold-sensitivity check (DESIGN §4). DEV-CV only (cache/dev).

For one fold seed s (folds from results/g2_prep/folds_seed{s}.tsv, made by g2_prep_folds.py; seed 0 = frozen folds)
the FULL nested core-v3 pipeline is re-run with the s-folds, exactly as in the producing scripts:
  1. core-v1: z_vif grid tau{.06,.125,.25,.5,1} x rho_f{0,1,3} x rho_g{0,.25,1} (scores fold-independent; rho_e on the
     7 DEV nulls), per-outer-fold argmax of inner M1-sel  (scripts/audit_core_v1_dev.py)
  2. B7 joint-posterior grid, PRIMARY variant clean+contrast, (alpha_sig, n_top, kappa) + core row; signatures for
     (outer f, eval g) fitted on folds not in {f, g}; per-fold nested pick -> (n_top, kappa)_f   (exp_b7_grid/eval)
  3. B7b late fusion, alpha in {.25,.5,1} nested with (n_top, kappa)_f inherited -> core-v2      (exp_b7b_eval)
  4. E1 permutation-null scale v_k from site-permuted GT conditions of folds != f -> core-v3   (exp_b8_eval run_e1)
Reference baselines (fixed configs, OOF where they have a training part): PSP z-score min3, curated_union z-score
min1, prior-only (freq from the other 4 folds, oracle sign); plus z-score with (library, min) nested over
{psp_all, psp_human_only, curated_union} x {1,3,5} (selection-symmetric baseline, descriptive).
Writes results/g2_prep/resplit_seed{s}.json; ledger branch G2PREP (one row per seed for core-v3 and core-v1).
Usage: g2_prep_resplit.py --seed S
"""
import argparse
import itertools
import json
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/workspace")
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))
import exp_b21_common as E21  # noqa: E402
from exp_b7_common import B7Ctx, FOLDS, frame, ALPHAS as B7_ALPHAS, NTOPS, KAPPAS  # noqa: E402
from exp_b7_eval import delta, m1sel  # noqa: E402
from exp_b7b_eval import Late, VARIANT, ALPHAS as B7B_ALPHAS  # noqa: E402
from audit_core_v1_dev import summary, t_unseen_mask  # noqa: E402
from kinact import bootstrap as bs, controls as ctl, ledger  # noqa: E402
from kinact.eval import to_jsonable  # noqa: E402
from kinact.explore import signatures_late as SL, final_box as FB  # noqa: E402
from kinact.methods.linear import zscore  # noqa: E402

OUT = ROOT / "results/g2_prep"
TAUS = [0.06, 0.125, 0.25, 0.5, 1.0]; RFS = [0.0, 1.0, 3.0]; RGS = [0.0, 0.25, 1.0]
ZLIBS = ["psp_all", "psp_human_only", "curated_union"]; ZMINS = [1, 3, 5]


def set_folds(ctx, split):
    """Re-assign DEV-CV fold labels everywhere the pipeline reads them (TEST/devconfirm untouched)."""
    c = ctx.c
    m = c.C["split"].isin(FOLDS)
    c.C.loc[m, "split"] = c.C.loc[m, "conditionID"].map(split).to_numpy()
    assert c.C.loc[m, "split"].isin(FOLDS).all()
    mm = c.meta["split"].isin(FOLDS)
    c.meta.loc[mm, "split"] = c.meta.loc[mm, "conditionID"].map(split).to_numpy()
    assert c.meta.loc[mm, "split"].isin(FOLDS).all()
    c.fold_of = c.meta["split"].copy()
    c.cf = pd.Series(c.conds).map(c.fold_of).to_numpy()
    ctx.fold_pos = {f: np.flatnonzero(c.cf == f) for f in FOLDS}
    ctx.null_fold = pd.Series(c.nulls).map(c.fold_of).to_numpy()


def nested_from_full(c, S_by_key, keys):
    """Fold-independent scores: inner_f = M1-sel over folds != f; pick argmax (first in key order)."""
    rows = []
    for k in keys:
        P = c.ev.pairs(S_by_key[k]); f = P["conditionID"].map(c.fold_of)
        rows.append({"key": k, **{f"inner_{g}": c.ev.metric(P[f != g], "M1-sel")["value"] for g in FOLDS},
                     "oof": c.ev.metric(P, "M1-sel")["value"]})
    g = pd.DataFrame(rows)
    parts, chosen = {}, {}
    for f in FOLDS:
        b = g.loc[g[f"inner_{f}"].idxmax()]
        chosen[f] = b["key"]
        parts[f] = S_by_key[b["key"]]
    A = np.full((len(c.conds), len(c.kin)), np.nan)
    for f in FOLDS:
        pos = c.cf == f
        A[pos] = parts[f].to_numpy()[pos]
    return pd.DataFrame(A, index=c.conds, columns=parts[FOLDS[0]].columns), chosen, g


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--seed", type=int, required=True)
    args = ap.parse_args()
    s = args.seed
    t0 = time.time()
    ctx = B7Ctx(); c = ctx.c
    split = pd.read_csv(OUT / f"folds_seed{s}.tsv", sep="\t").set_index("conditionID")["split"]
    set_folds(ctx, split)
    rep = {"seed": s, "gt_hash": ctx.gt_hash,
           "fold_sizes_gt_conds": pd.Series(c.cf).value_counts().sort_index().to_dict()}
    P_all = c.ev.pairs(pd.DataFrame(0.0, index=c.conds, columns=c.kin))
    rep["fold_m1sel_pairs"] = P_all[P_all["in_M1sel"]]["conditionID"].map(c.fold_of).value_counts().sort_index().to_dict()
    idx_all = np.arange(len(c.conds))

    # ---------------- 1. core-v1 grid (fold-independent scores) + nested
    core_S = {}
    for tau, rf, rg in itertools.product(TAUS, RFS, RGS):
        Z, _ = ctx.score(idx_all, [("x", None)], (tau, rf, rg), [])[("x", 0.0)]
        core_S[(tau, rf, rg)] = pd.DataFrame(Z, index=c.conds, columns=c.kin)
    S_v1, ch_core, _ = nested_from_full(c, core_S, list(core_S))
    ctx.picks = {f: tuple(float(x) for x in ch_core[f]) for f in FOLDS}
    print(f"[s{s}] core picks {ctx.picks}; core-v1 M1-sel {m1sel(c.ev, S_v1, FOLDS, c.fold_of):.4f} "
          f"({time.time()-t0:.0f}s)", flush=True)
    zcore = {(f, g): core_S[ctx.picks[f]].to_numpy()[ctx.fold_pos[g]] for f in FOLDS for g in FOLDS}

    # ---------------- 2. B7 joint grid (primary variant) + nested (n_top, kappa)
    fits, store = {}, {}
    jobs = [([x for x in FOLDS if x not in (f, g)], [f, g], [(f, g), (g, f)]) for f, g in itertools.combinations(FOLDS, 2)]
    jobs += [([x for x in FOLDS if x != f], [f], [(f, f)]) for f in FOLDS]
    for tr, held, uses in jobs:
        evc = list(c.meta.index[c.meta["split"].isin(held)])
        ft = ctx.fit(tr, evc, VARIANT)
        for f_out, g_ev in uses:
            fits[(f_out, g_ev)] = ft
            Wl = [((kp, nt), ft[kp].matrix(len(c.X.index), c.kin, nt)) for kp, nt in itertools.product(KAPPAS, NTOPS)]
            res = ctx.score(ctx.fold_pos[g_ev], Wl, ctx.picks[f_out], B7_ALPHAS)
            for (lab, a), (Z, _) in res.items():
                store[(f_out, g_ev, lab[0], lab[1], a)] = Z
    print(f"[s{s}] B7 grid done ({time.time()-t0:.0f}s)", flush=True)
    cfgs = [(0.0, 0, 0.0)] + [(a, nt, kp) for kp, nt, a in itertools.product(KAPPAS, NTOPS, B7_ALPHAS)]
    get = lambda cf, fo, g: zcore[(fo, g)] if cf[0] == 0 else store[(fo, g, cf[2], cf[1], cf[0])]
    rows = []
    for cf in cfgs:
        r = {"cfg": cf}
        for f in FOLDS:
            r[f"inner_{f}"] = m1sel(c.ev, frame(ctx, {g: get(cf, f, g) for g in FOLDS if g != f}), [g for g in FOLDS if g != f], c.fold_of)
        r["oof"] = m1sel(c.ev, frame(ctx, {f: get(cf, f, f) for f in FOLDS}), FOLDS, c.fold_of)
        rows.append(r)
    g7 = pd.DataFrame(rows)
    b7_chosen, picks = {}, {}
    for f in FOLDS:
        b = g7.loc[g7[f"inner_{f}"].idxmax()]
        b7_chosen[f] = b["cfg"]
        if b["cfg"][0] == 0:     # core row won: inherit (n_top, kappa) from the best signature row (documented)
            b = g7.iloc[1:].loc[g7.iloc[1:][f"inner_{f}"].idxmax()]
        picks[f] = (int(b["cfg"][1]), float(b["cfg"][2]))
    S_b7 = frame(ctx, {f: get(b7_chosen[f], f, f) for f in FOLDS})
    del store
    print(f"[s{s}] B7 picks {b7_chosen}", flush=True)

    # ---------------- 3. B7b late fusion, alpha nested -> core-v2
    L = Late(ctx)
    zs = {k: L.zsig(fits[k][picks[k[0]][1]], picks[k[0]][0], "gt", ctx.fold_pos[k[1]])[0] for k in fits}
    rows, outer = [], {}
    for a in B7B_ALPHAS:
        r = {"alpha": a}
        for f in FOLDS:
            r[f"inner_{f}"] = m1sel(c.ev, frame(ctx, {g: SL.fuse(zcore[(f, g)], zs[(f, g)], a) for g in FOLDS if g != f}),
                                    [g for g in FOLDS if g != f], c.fold_of)
        outer[a] = {f: SL.fuse(zcore[(f, f)], zs[(f, f)], a) for f in FOLDS}
        r["oof"] = m1sel(c.ev, frame(ctx, outer[a]), FOLDS, c.fold_of)
        rows.append(r)
    gb = pd.DataFrame(rows)
    alpha = {f: float(gb.loc[gb[f"inner_{f}"].idxmax(), "alpha"]) for f in FOLDS}
    S_v2 = frame(ctx, {f: outer[alpha[f]][f] for f in FOLDS})

    # ---------------- 4. E1 -> core-v3
    parts = {}
    sigmask = {f: np.isin(np.asarray(c.kin), list(fits[(f, f)][picks[f][1]].sites)) for f in FOLDS}
    for f in FOLDS:
        tp = np.concatenate([ctx.fold_pos[g] for g in FOLDS if g != f])
        zc_p = ctx.score(tp, [("x", None)], ctx.picks[f], [], which="perm")[("x", 0.0)][0]
        zs_p = L.zsig(fits[(f, f)][picks[f][1]], picks[f][0], "perm", tp)[0]
        v, vs = FB.perm_scales(zc_p.astype(np.float64), zs_p, alpha[f], sigmask[f])
        parts[f] = FB.fuse_calibrated(zcore[(f, f)], zs[(f, f)], alpha[f], v, vs)
    S_v3 = frame(ctx, parts)
    print(f"[s{s}] alpha {alpha}; v1 {m1sel(c.ev, S_v1, FOLDS, c.fold_of):.4f} v2 {m1sel(c.ev, S_v2, FOLDS, c.fold_of):.4f} "
          f"v3 {m1sel(c.ev, S_v3, FOLDS, c.fold_of):.4f} ({time.time()-t0:.0f}s)", flush=True)

    # ---------------- baselines
    X = c.X[c.conds]
    psp = zscore(X, "psp_all", min_substrates=3).reindex(columns=c.kin)
    cur = zscore(X, "curated_union", min_substrates=1).reindex(columns=c.kin)
    po = []
    for f in FOLDS:
        ec = [x for x in c.conds if c.fold_of[x] == f]
        tc = list(c.C.loc[c.C["split"].isin([g for g in FOLDS if g != f]), "conditionID"])
        po.append(ctl.prior_only(c.P, c.C, tc, ec, c.kin))
    prior = pd.concat(po).loc[c.conds]
    zgrid = {(lib, mn): zscore(X, lib, min_substrates=mn).reindex(columns=c.kin) for lib in ZLIBS for mn in ZMINS}
    S_zn, ch_z, _ = nested_from_full(c, zgrid, list(zgrid))

    meth = {"core-v3": S_v3, "core-v2": S_v2, "core-v1": S_v1, "B7 joint (nested)": S_b7,
            "PSP z-score min3": psp, "curated z-score min1": cur, "prior-only": prior, "z-score nested lib/min": S_zn}
    res = {}
    for n, S in meth.items():
        r = summary(c, S, psp if n != "PSP z-score min3" else None)
        r.pop("harmony_per_bin", None)
        res[n] = r
    plan = bs.BootstrapPlan.from_evaluator(c.ev, B=2000, seed=0)
    deltas = {}
    for comp in ["core-v1", "PSP z-score min3", "curated z-score min1", "prior-only", "z-score nested lib/min"]:
        deltas[comp] = {ep: delta(c.ev, S_v3, meth[comp], ep, plan) for ep in ["M1-sel", "M1-all", "harmony"]}
        deltas[comp]["T-unseen"] = delta(c.ev, S_v3, meth[comp], "M1-sel", plan, lambda P: t_unseen_mask(c, P))
    deltas["core-v2 vs core-v1"] = {"M1-sel": delta(c.ev, S_v2, S_v1, "M1-sel", plan)}
    rep.update({"core_picks": ctx.picks, "b7_chosen": {f: list(v) for f, v in b7_chosen.items()},
                "ntop_kappa": picks, "alpha": alpha, "zscore_nested_chosen": {f: list(v) for f, v in ch_z.items()},
                "b7_grid_top": g7.sort_values("oof", ascending=False).head(5).astype(str).to_dict("records"),
                "b7b_grid": gb.to_dict("records"),
                "n_sig_kinases_outer": {f: int(sigmask[f].sum()) for f in FOLDS},
                "results": res, "delta_core_v3": deltas, "wall_s": time.time() - t0})
    (OUT / f"resplit_seed{s}.json").write_text(json.dumps(to_jsonable(rep), indent=1, default=str))
    pickle.dump({n: S for n, S in meth.items()}, open(ROOT / f"cache/g2_prep/resplit_seed{s}_scores.pkl", "wb"))
    for n in ["core-v3", "core-v1"]:
        r = res[n]
        ledger.append(f"{n} nested (fold seed {s})", {"bundle": n, "fold_seed": s, "pipeline": "g2_prep_resplit.py"},
                      "dev", {k: r[k] for k in ["M1-sel", "M1-all", "M1-sel-scored", "coverage_M1-sel", "harmony", "T-unseen"]},
                      branch="G2PREP", folds=range(5))
    print(f"[s{s}] done {time.time()-t0:.0f}s", json.dumps({n: round(r["M1-sel"], 4) for n, r in res.items()}), flush=True)


if __name__ == "__main__":
    main()
