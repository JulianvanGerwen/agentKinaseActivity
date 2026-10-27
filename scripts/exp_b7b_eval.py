#!/usr/bin/env python
"""B7b: late fusion of B7 signatures (pre-declared in results/b7b/b7b_report.md). DEV only.
Nested alpha in {0.25, 0.5, 1}; (n_top, kappa) fixed per outer fold to the B7 PRIMARY nested picks. One primary
comparison vs core-v1 nested + C10 checks. Writes results/b7b/b7b_eval.json.
"""
import itertools
import json
import pickle
import time

import numpy as np
import pandas as pd

from exp_b7_common import B7Ctx, CACHE, FOLDS, frame
from exp_b7_eval import delta, drug_disjointness, m1sel
from exp_b21_common import calib_bins
from audit_core_v1_dev import summary, t_unseen_mask
from kinact import bootstrap as bs, ledger
from kinact.eval import to_jsonable
from kinact.explore import signatures_late as SL
from kinact.methods.linear import zscore
from pathlib import Path

ROOT = Path("/workspace")
OUT = ROOT / "results/b7b"; OUT.mkdir(parents=True, exist_ok=True)
BRANCH = "B7b"
VARIANT = "clean+contrast"
ALPHAS = [0.25, 0.5, 1.0]
N_SHUFFLE = 3


def b7_picks():
    d = json.loads((ROOT / "results/b7/b7_eval.json").read_text())
    ch = d["results"][f"B7 {VARIANT} (nested)"]["chosen"]
    return {f: (int(c["n_top"]), float(c["kappa"])) for f, c in ch.items()}


def cfg(alpha, picks):
    return dict(model="core-v1 z_vif + late signature z", variant=VARIANT, alpha=alpha, fusion="(zc+a zs)/sqrt(1+a^2)",
                ntop_kappa_per_fold={f: list(v) for f, v in picks.items()}, clip=SL.CLIP, min_sites=SL.MIN_SITES,
                rho_sig="RMS over DEV nulls")


class Late:
    def __init__(self, ctx):
        self.ctx = ctx
        c = ctx.c
        self.gt_cols = ctx.allconds.get_indexer(c.conds)
        self.null_cols = ctx.allconds.get_indexer(c.nulls)
        self.nsites = len(c.X.index)

    def zsig(self, sigs, n_top, which, idx):
        """Standardised z_sig for 'gt' (positions in c.conds), 'null' (positions in c.nulls) or 'perm'."""
        ctx = self.ctx
        Ws = sigs.matrix(self.nsites, ctx.c.kin, n_top)
        zn = SL.sig_project(Ws, ctx.Y[:, self.null_cols], ctx.M[:, self.null_cols])
        rho = SL.rho_sig(zn)
        if which == "gt":
            cols = self.gt_cols[idx]; Z = SL.sig_project(Ws, ctx.Y[:, cols], ctx.M[:, cols])
        elif which == "null":
            Z = zn[idx]
        else:
            Z = SL.sig_project(Ws, ctx.Yp[:, idx], ctx.Mp[:, idx])
        return Z / rho, rho


def main():
    t0 = time.time()
    ctx = B7Ctx()
    c = ctx.c
    L = Late(ctx)
    picks = b7_picks()
    print("ntop/kappa picks", picks, flush=True)
    store = pickle.load(open(CACHE / f"grid_{VARIANT}.pkl", "rb"))
    zcore = {(f, g): store[("core", f, g)][0] for f in FOLDS for g in FOLDS}
    del store
    S_core = frame(ctx, {f: zcore[(f, f)] for f in FOLDS})
    # ---------------- signature fits + z_sig (inner: folds not in {f,g}; outer: folds != f)
    zs, rhos, outer_fits = {}, {}, {}
    for f, g in itertools.combinations(FOLDS, 2):
        held = [f, g]
        tr = [x for x in FOLDS if x not in held]
        evc = list(c.meta.index[c.meta["split"].isin(held)])
        kps = sorted({picks[f][1], picks[g][1]})
        fits = ctx.fit(tr, evc, VARIANT, kappas=kps)
        for fo, ge in [(f, g), (g, f)]:
            nt, kp = picks[fo]
            zs[(fo, ge)], rhos[(fo, ge)] = L.zsig(fits[kp], nt, "gt", ctx.fold_pos[ge])
    for f in FOLDS:
        nt, kp = picks[f]
        outer_fits[f] = pickle.load(open(CACHE / f"sig_{VARIANT}_{f}.pkl", "rb"))[kp]
        zs[(f, f)], rhos[(f, f)] = L.zsig(outer_fits[f], nt, "gt", ctx.fold_pos[f])
    print(f"fits done {time.time()-t0:.0f}s; rho_sig range {min(rhos.values()):.2f}-{max(rhos.values()):.2f}", flush=True)
    # ---------------- per-alpha inner / OOF, nested pick
    rows, outer = [], {}
    for a in ALPHAS:
        r = {"alpha": a}
        for f in FOLDS:
            parts = {g: SL.fuse(zcore[(f, g)], zs[(f, g)], a) for g in FOLDS if g != f}
            r[f"inner_{f}"] = m1sel(c.ev, frame(ctx, parts), [g for g in FOLDS if g != f], c.fold_of)
        S = frame(ctx, {f: SL.fuse(zcore[(f, f)], zs[(f, f)], a) for f in FOLDS})
        outer[a] = S
        P = c.ev.pairs(S)
        r["oof_M1-sel"] = c.ev.metric(P, "M1-sel")["value"]; r["oof_M1-all"] = c.ev.metric(P, "M1-all")["value"]
        rows.append(r)
        ledger.append("core-v1+late-sig", cfg(a, picks), "dev",
                      {"M1-sel": r["oof_M1-sel"], "M1-all": r["oof_M1-all"],
                       **{k: v for k, v in r.items() if k.startswith("inner_")}, "note": "per-alpha OOF"},
                      branch=BRANCH, folds=range(5))
    g = pd.DataFrame(rows)
    g.to_csv(OUT / "grid.tsv", sep="\t", index=False)
    print(g.round(4).to_string(), flush=True)
    chosen, parts = {}, {}
    for f in FOLDS:
        b = g.loc[g[f"inner_{f}"].idxmax()]
        chosen[f] = float(b.alpha)
        parts[f] = outer[chosen[f]].to_numpy()[ctx.fold_pos[f]]
    S = frame(ctx, parts)
    print("chosen", chosen, flush=True)
    # sanity: kinases without a signature are unchanged
    has_k = {f: set(outer_fits[f].sites) for f in FOLDS}
    for f in FOLDS:
        no = ~np.isin(np.asarray(c.kin), list(has_k[f]))
        A, B0 = S.to_numpy()[ctx.fold_pos[f]][:, no], S_core.to_numpy()[ctx.fold_pos[f]][:, no]
        assert np.allclose(np.nan_to_num(A, nan=9e9), np.nan_to_num(B0, nan=9e9), atol=1e-6)
    # ---------------- primary comparison + C10 checks
    plan = bs.BootstrapPlan.from_evaluator(c.ev, B=2000, seed=0)
    psp = zscore(c.X[c.conds], "psp_all", min_substrates=3).reindex(columns=c.kin)
    res = {"core-v1": summary(c, S_core, psp)}
    r = summary(c, S, psp)
    r["chosen_alpha"] = chosen
    r["M2_vs_core"] = c.ev.m2(S, S_core)
    d = {ep: delta(c.ev, S, S_core, ep, plan) for ep in ["M1-sel", "M1-all", "harmony", "M2"]}
    has = lambda P: np.array([k in has_k[c.fold_of[cid]] for k, cid in zip(P["kinase"], P["conditionID"])])
    uns = lambda P: t_unseen_mask(c, P)
    masks = {"T-unseen": uns, "T-seen": lambda P: ~uns(P), "has-signature": has, "no-signature": lambda P: ~has(P)}
    for name, mk in masks.items():
        d[name] = delta(c.ev, S, S_core, "M1-sel", plan, mk)
        Pp = c.ev.pairs(S); r[name] = c.ev.metric(Pp[mk(Pp)], "M1-sel")["value"]
        Pc = c.ev.pairs(S_core); res["core-v1"][name] = c.ev.metric(Pc[mk(Pc)], "M1-sel")["value"]
    for src in ["H", "D", "O"]:
        d[f"source={src}"] = delta(c.ev, S, S_core, f"M1-sel:source={src}", plan)
    for sg in ["-1", "1"]:
        d[f"sign={sg}"] = delta(c.ev, S, S_core, f"M1-sel:sign={sg}", plan)
    for b in ["0", "1-4", "5-19", ">=20"]:
        d[f"bin={b}"] = delta(c.ev, S, S_core, f"M1-sel:bin={b}", plan)
    r["delta_vs_core"] = d
    print("B7b", {k: r[k] for k in ["M1-sel", "M1-all", "harmony", "T-unseen"]}, d["M1-sel"], flush=True)
    # shuffled-signature control (outer fits, fold-wise chosen alpha)
    sh = []
    for seed in range(N_SHUFFLE):
        p2 = {}
        for f in FOLDS:
            z, _ = L.zsig(outer_fits[f].shuffled(seed), picks[f][0], "gt", ctx.fold_pos[f])
            p2[f] = SL.fuse(zcore[(f, f)], z, chosen[f])
        Ssh = frame(ctx, p2)
        sh.append(delta(c.ev, Ssh, S_core, "M1-sel", plan))
        sh[-1]["M1-sel"] = m1sel(c.ev, Ssh, FOLDS, c.fold_of)
    ctl = {"shuffled": sh, "shuffled_mean_delta": float(np.mean([x["estimate"] for x in sh]))}
    ctl["real_minus_shuffled"] = d["M1-sel"]["estimate"] - ctl["shuffled_mean_delta"]
    r["controls"] = ctl
    print("shuffled", ctl, flush=True)
    # calibration: fold-f DEV nulls with the fold-f design; site-permuted GT conditions
    cal = {}
    for lab, use_sig in [("core-v1", False), ("B7b", True)]:
        Zs, Ns, Zp, Np = [], [], [], []
        for f in FOLDS:
            ni = np.flatnonzero(ctx.null_fold == f)
            if len(ni):
                Zc, _, n = ctx.score(ni, [("x", None)], ctx.picks[f], [], which="null", return_n=True)[("x", 0.0)]
                if use_sig:
                    Zc = SL.fuse(Zc, L.zsig(outer_fits[f], picks[f][0], "null", ni)[0], chosen[f])
                Zs.append(Zc); Ns.append(n)
            idx = ctx.fold_pos[f]
            Zc, _, n = ctx.score(idx, [("x", None)], ctx.picks[f], [], which="perm", return_n=True)[("x", 0.0)]
            if use_sig:
                Zc = SL.fuse(Zc, L.zsig(outer_fits[f], picks[f][0], "perm", idx)[0], chosen[f])
            Zp.append(Zc); Np.append(n)
        cal[lab] = {"null": calib_bins(np.vstack(Zs), np.vstack(Ns)), "null_ncells": int(np.isfinite(np.vstack(Zs)).sum()),
                    "perm": calib_bins(np.vstack(Zp), np.vstack(Np))}
    # z_sig alone calibration (outer fits; perm conditions) for diagnosis
    zp = np.vstack([L.zsig(outer_fits[f], picks[f][0], "perm", ctx.fold_pos[f])[0] for f in FOLDS])
    v = zp[np.isfinite(zp)]
    cal["z_sig_perm"] = {"rate": float((np.abs(v) > 1.96).mean()), "sd": float(v.std()), "n": int(len(v))}
    r["calibration"] = cal
    print("calibration", cal, flush=True)
    res["B7b (nested)"] = r
    ledger.append("core-v1+late-sig nested", {**cfg("nested", picks), "grid": ALPHAS, "chosen": chosen}, "dev",
                  {k: r[k] for k in ["M1-sel", "M1-all", "M1-sel-scored", "coverage_M1-sel", "T-unseen"]},
                  branch=BRANCH, folds=range(5))
    rep = {"gt_hash": ctx.gt_hash, "ntop_kappa_picks": picks, "rho_sig": {f"{a}|{b}": v for (a, b), v in rhos.items()},
           "drug_disjointness": drug_disjointness(ctx), "grid": rows, "results": res,
           "n_sig_kinases_outer": {f: len(has_k[f]) for f in FOLDS}, "wall_s": time.time() - t0}
    (OUT / "b7b_eval.json").write_text(json.dumps(to_jsonable(rep), indent=1, default=str))
    print(f"done {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
