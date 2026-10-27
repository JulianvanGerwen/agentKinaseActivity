#!/usr/bin/env python
"""B8 final box (pre-declared in results/b8/b8_report.md). DEV only.
Reproduces core-v2 (B7b nested), then E1 (permutation-null fused calibration), E2 (top-r shared-factor removal,
r nested in {1,2,3}), E3 (top-3 core-config bagging). One comparison each vs core-v2.
Writes results/b8/b8_eval.json and cache/b8/*.pkl (score frames + components for the bundle script).
Usage: exp_b8_eval.py [--only E1,E2,E3]
"""
import argparse
import itertools
import json
import pickle
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from exp_b7_common import B7Ctx, CACHE as B7CACHE, FOLDS, frame
from exp_b7_eval import delta, m1sel
from exp_b7b_eval import Late, b7_picks, VARIANT
from exp_b21_common import calib_bins
from audit_core_v1_dev import summary, t_unseen_mask
from kinact import bootstrap as bs, ledger
from kinact.eval import to_jsonable, psp_bin
from kinact.explore import bayes2 as B2, signatures_late as SL, final_box as FB, coreg as CR
from kinact.explore.bayes import report_mask
from kinact.methods.linear import zscore

ROOT = Path("/workspace")
OUT = ROOT / "results/b8"; OUT.mkdir(parents=True, exist_ok=True)
CACHE = ROOT / "cache/b8"; CACHE.mkdir(parents=True, exist_ok=True)
BRANCH = "B8"
RS = [1, 2, 3]
TAUS = [0.06, 0.125, 0.25, 0.5, 1.0]; RGS = [0.0, 0.25, 1.0]


def b7b_alpha():
    d = json.loads((ROOT / "results/b7b/b7b_eval.json").read_text())
    return {f: float(a) for f, a in d["results"]["B7b (nested)"]["chosen_alpha"].items()}


class B8:
    def __init__(self):
        self.t0 = time.time()
        self.ctx = ctx = B7Ctx()
        self.c = c = ctx.c
        self.L = Late(ctx)
        self.picks = b7_picks()
        self.alpha = b7b_alpha()
        store = pickle.load(open(B7CACHE / f"grid_{VARIANT}.pkl", "rb"))
        self.zcore = {(f, g): store[("core", f, g)][0] for f in FOLDS for g in FOLDS}
        del store
        self.fits = {}            # (f, g) -> Signatures used for outer f scoring eval fold g
        for f, g in itertools.combinations(FOLDS, 2):
            tr = [x for x in FOLDS if x not in (f, g)]
            evc = list(c.meta.index[c.meta["split"].isin([f, g])])
            kps = sorted({self.picks[f][1], self.picks[g][1]})
            fits = ctx.fit(tr, evc, VARIANT, kappas=kps)
            self.fits[(f, g)] = fits[self.picks[f][1]]; self.fits[(g, f)] = fits[self.picks[g][1]]
        for f in FOLDS:
            self.fits[(f, f)] = pickle.load(open(B7CACHE / f"sig_{VARIANT}_{f}.pkl", "rb"))[self.picks[f][1]]
        self.zs = {k: self.L.zsig(fit, self.picks[k[0]][0], "gt", ctx.fold_pos[k[1]])[0] for k, fit in self.fits.items()}
        self.sigmask = {f: np.isin(np.asarray(c.kin), list(self.fits[(f, f)].sites)) for f in FOLDS}
        self.Wb_csc = ctx.Wb.T.tocsr()                         # p x sites
        self.nsites = len(c.X.index)
        self.plan = bs.BootstrapPlan.from_evaluator(c.ev, B=2000, seed=0)
        self.psp = zscore(c.X[c.conds], "psp_all", min_substrates=3).reindex(columns=c.kin)
        print(f"[b8] init {time.time()-self.t0:.0f}s", flush=True)

    # ------------------------------------------------------------------ core-v2
    def core_v2(self):
        parts = {f: SL.fuse(self.zcore[(f, f)], self.zs[(f, f)], self.alpha[f]) for f in FOLDS}
        return frame(self.ctx, parts)

    def train_pos(self, f):
        return np.concatenate([self.ctx.fold_pos[g] for g in FOLDS if g != f])

    # ------------------------------------------------------------------ scoring with an arbitrary input (projection)
    def h_of(self, Yp, Mp):
        return torch.as_tensor(np.asarray(self.Wb_csc @ (Mp * Yp)).T, dtype=torch.float32, device="cuda")

    def zcore_custom(self, which, idx, core, YE, ME, YN, MN):
        """z_vif of conditions (which, idx) with input YE (sites x len(idx)) and nulls YN, rho_e refit on YN."""
        ctx = self.ctx
        tau, rf, rg = core
        K = ctx.kernel(rf, rg)
        res = []
        for wh, ix, Y, M in (("null", np.arange(len(self.c.nulls)), YN, MN), (which, idx, YE, ME)):
            b = ctx._sub(wh, ix)
            post = B2.posterior2(b["G"], self.h_of(Y, M), b["g0"], K, tau)
            res.append((post, report_mask(b["n"], K, 1, True), b["n"]))
        (pN, rN, _), (pE, rE, nE) = res
        rho, _ = B2.fit_rho_null(B2.output_stat(pN, "z_freq"), B2.m_eff(pN["q"], pN["suma"]), rN)
        return np.where(rE, B2.output_stat(pE, "z_vif", rho_e=rho), np.nan), nE

    def zsig_custom(self, fit, n_top, YE, ME, YN, MN):
        Ws = fit.matrix(self.nsites, self.c.kin, n_top)
        rho = SL.rho_sig(SL.sig_project(Ws, YN, MN))
        return SL.sig_project(Ws, YE, ME) / rho

    def YM(self, which, idx):
        return self.ctx._YM(which, idx)

    # ------------------------------------------------------------------ evaluation helpers
    def compare(self, S, S_ref, label):
        c, plan = self.c, self.plan
        d = {ep: delta(c.ev, S, S_ref, ep, plan) for ep in ["M1-sel", "M1-all", "harmony", "M2"]}
        uns = lambda P: t_unseen_mask(c, P)
        for name, mk in {"T-unseen": uns, "T-seen": lambda P: ~uns(P)}.items():
            d[name] = delta(c.ev, S, S_ref, "M1-sel", plan, mk)
        for src in ["H", "D", "O"]:
            d[f"source={src}"] = delta(c.ev, S, S_ref, f"M1-sel:source={src}", plan)
        for sg in ["-1", "1"]:
            d[f"sign={sg}"] = delta(c.ev, S, S_ref, f"M1-sel:sign={sg}", plan)
        for b in ["0", "1-4", "5-19", ">=20"]:
            d[f"bin={b}"] = delta(c.ev, S, S_ref, f"M1-sel:bin={b}", plan)
        ni = {"M1-all": d["M1-all"]["estimate"] >= -0.02, "M2": d["M2"]["estimate"] >= -0.02,
              "T-unseen": d["T-unseen"]["estimate"] >= -0.01}
        for k in ["source=H", "source=D", "source=O", "bin=0", "bin=1-4", "bin=5-19", "bin=>=20"]:
            ni[k] = d[k]["estimate"] >= -0.02
        d["noninferiority"] = ni
        d["hard_test_pass"] = bool(d["M1-sel"]["p_one_sided"] < 0.05 and d["M1-sel"]["estimate"] > 0)
        print(label, "dM1-sel", {k: round(v, 4) for k, v in d["M1-sel"].items()}, "NI", all(ni.values()), flush=True)
        return d

    def full(self, S):
        r = summary(self.c, S, self.psp)
        Pp = self.c.ev.pairs(S); m = t_unseen_mask(self.c, Pp)
        r["T-seen"] = self.c.ev.metric(Pp[~m], "M1-sel")["value"]
        return r

    def calib(self, Zn, Nn, Sn, Zp, Np_, Sp):
        """Z: stacked z; N: n_substrates; S: signature-kinase mask (same shape)."""
        out = {}
        for lab, Z, N, Sm in (("null", Zn, Nn, Sn), ("perm", Zp, Np_, Sp)):
            out[lab] = calib_bins(Z, N)
            out[lab + "_sigkin"] = calib_bins(Z, N, ok=Sm)
            out[lab + "_sigkin_ncells"] = int((Sm & np.isfinite(Z)).sum())
        # per signature kinase (perm, >= 50 cells)
        ok = Sp & np.isfinite(Zp)
        rates = [(np.abs(Zp[ok[:, j], j]) > 1.96).mean() for j in range(Zp.shape[1]) if ok[:, j].sum() >= 50]
        if rates:
            out["perm_sigkin_per_kinase"] = {"n": len(rates), "median": float(np.median(rates)),
                                             "p90": float(np.quantile(rates, .9)), "max": float(np.max(rates))}
        return out


def ledger_row(method, cfg, r, extra=None):
    ledger.append(method, cfg, "dev", {**{k: r[k] for k in ["M1-sel", "M1-all", "M1-sel-scored", "coverage_M1-sel",
                                                             "T-unseen"] if k in r}, **(extra or {})},
                  branch=BRANCH, folds=range(5))


# ====================================================================== E1
def run_e1(b: B8, S_ref, res):
    ctx, c = b.ctx, b.c
    parts, scales, cal_parts = {}, {}, {"v2": [], "e1": []}
    for f in FOLDS:
        a, nt = b.alpha[f], b.picks[f][0]
        tp = b.train_pos(f)
        zc_p = ctx.score(tp, [("x", None)], ctx.picks[f], [], which="perm")[("x", 0.0)][0]
        zs_p = b.L.zsig(b.fits[(f, f)], nt, "perm", tp)[0]
        v, vs = FB.perm_scales(zc_p.astype(np.float64), zs_p, a, b.sigmask[f])
        scales[f] = {"v_median": float(np.nanmedian(v)), "v_min": float(np.nanmin(v)), "v_max": float(np.nanmax(v)),
                     "fixed_1+a2": 1 + a * a}
        parts[f] = FB.fuse_calibrated(b.zcore[(f, f)], b.zs[(f, f)], a, v, vs)
        # calibration conditions: fold-f nulls and fold-f permuted GT conditions
        ni = np.flatnonzero(ctx.null_fold == f)
        idx = ctx.fold_pos[f]
        for lab in ("v2", "e1"):
            fu = (lambda zc, zs: SL.fuse(zc, zs, a)) if lab == "v2" else (lambda zc, zs: FB.fuse_calibrated(zc, zs, a, v, vs))
            rec = {}
            if len(ni):
                Zc, _, n = ctx.score(ni, [("x", None)], ctx.picks[f], [], which="null", return_n=True)[("x", 0.0)]
                rec["null"] = (fu(Zc, b.L.zsig(b.fits[(f, f)], nt, "null", ni)[0]), n)
            Zc, _, n = ctx.score(idx, [("x", None)], ctx.picks[f], [], which="perm", return_n=True)[("x", 0.0)]
            rec["perm"] = (fu(Zc, b.L.zsig(b.fits[(f, f)], nt, "perm", idx)[0]), n)
            rec["sig"] = b.sigmask[f]
            cal_parts[lab].append(rec)
    S = frame(ctx, parts)
    cal = {}
    for lab, recs in cal_parts.items():
        Zn = np.vstack([r["null"][0] for r in recs if "null" in r]); Nn = np.vstack([r["null"][1] for r in recs if "null" in r])
        Sn = np.vstack([np.broadcast_to(r["sig"], r["null"][0].shape) for r in recs if "null" in r])
        Zp = np.vstack([r["perm"][0] for r in recs]); Np_ = np.vstack([r["perm"][1] for r in recs])
        Sp = np.vstack([np.broadcast_to(r["sig"], r["perm"][0].shape) for r in recs])
        cal[lab] = b.calib(Zn, Nn, Sn, Zp, Np_, Sp)
    r = b.full(S)
    r["delta_vs_core_v2"] = b.compare(S, S_ref, "E1")
    r["scales"] = scales
    r["calibration"] = cal
    e_v2 = abs(cal["v2"]["perm_sigkin"]["all"] - 0.05); e_e1 = abs(cal["e1"]["perm_sigkin"]["all"] - 0.05)
    nl = cal["e1"]["null_sigkin"]["all"]
    r["adoption"] = {"a_perm_sigkin_closer_to_5pct": bool(e_e1 < e_v2), "b_null_sigkin_in_[2,8]%": bool(0.02 <= nl <= 0.08),
                     "c_dM1sel_ge_-0.01": bool(r["delta_vs_core_v2"]["M1-sel"]["estimate"] >= -0.01)}
    r["adopt"] = all(r["adoption"].values())
    print("E1 calibration", json.dumps(to_jsonable(cal))[:1500], "\nE1 adoption", r["adoption"], flush=True)
    ledger_row("core-v2+E1 perm-null fused calibration", dict(model="core-v2", e1="F/sqrt(mean F^2 on permuted train conds)",
               alpha=b.alpha, min_perm_cells=FB.MIN_PERM_CELLS), r)
    res["E1"] = r
    return S, cal


# ====================================================================== E2
def run_e2(b: B8, S_ref, res):
    ctx, c = b.ctx, b.c
    allc = np.asarray(ctx.allconds)
    gt_role = set(c.C.loc[c.C["role"] == "gt", "conditionID"])
    YN, MN = b.YM("null", np.arange(len(c.nulls)))
    factors, sv = {}, {}

    def factors_for(held, nongt=False):
        key = (tuple(sorted(held)), nongt)
        if key not in factors:
            ho = list(c.meta.index[c.meta["split"].isin(held)])
            ok = CR.allowed_mask(c.meta, allc, ho)
            if nongt:
                ok &= ~np.isin(allc, list(gt_role))
            rows, V, s = FB.fit_factors(ctx.Y, ctx.M, np.flatnonzero(ok))
            factors[key] = (rows, V); sv["|".join(key[0]) + ("|nonGT" if nongt else "")] = {
                "n_conds": int(ok.sum()), "n_sites": int(len(rows)), "sv": s.round(1).tolist()}
        return factors[key]

    def score_use(fo, ge, r, nongt=False, which="gt"):
        held = [fo] if fo == ge else [fo, ge]
        rows, V = factors_for(held, nongt)
        idx = ctx.fold_pos[ge]
        YE, ME = b.YM(which, idx)
        YEp = FB.project(YE, ME, rows, V, r); YNp = FB.project(YN, MN, rows, V, r)
        zc, n = b.zcore_custom(which, idx, ctx.picks[fo], YEp, ME, YNp, MN)
        zs = b.zsig_custom(b.fits[(fo, ge)], b.picks[fo][0], YEp, ME, YNp, MN)
        return SL.fuse(zc, zs, b.alpha[fo]), n

    # sanity: h from Wb equals the cached block h (unprojected)
    YE, ME = b.YM("gt", ctx.fold_pos["fold0"][:20])
    h1 = b.h_of(YE.astype(np.float64), ME).cpu().numpy(); h0 = ctx.base["gt"]["h"][ctx.fold_pos["fold0"][:20]].cpu().numpy()
    assert np.allclose(h1, h0, rtol=1e-3, atol=1e-3), np.abs(h1 - h0).max()
    rows_, outer = [], {}
    inner_z = {}
    for f in FOLDS:
        for g in FOLDS:
            if g != f:
                for r in RS:
                    inner_z[(f, g, r)] = score_use(f, g, r)[0]
        print(f"E2 inner {f} done {time.time()-b.t0:.0f}s", flush=True)
    for r in RS:
        rec = {"r": r}
        for f in FOLDS:
            rec[f"inner_{f}"] = m1sel(c.ev, frame(ctx, {g: inner_z[(f, g, r)] for g in FOLDS if g != f}),
                                      [g for g in FOLDS if g != f], c.fold_of)
        outer[r] = {f: score_use(f, f, r)[0] for f in FOLDS}
        S_r = frame(ctx, outer[r]); P = c.ev.pairs(S_r)
        rec["oof_M1-sel"] = c.ev.metric(P, "M1-sel")["value"]; rec["oof_M1-all"] = c.ev.metric(P, "M1-all")["value"]
        rows_.append(rec)
        ledger.append("core-v2+E2 factor removal", dict(model="core-v2", r=r, site_frac=FB.SITE_FRAC, clip=FB.CLIP,
                      min_factor_sites=FB.MIN_FACTOR_SITES), "dev",
                      {"M1-sel": rec["oof_M1-sel"], "M1-all": rec["oof_M1-all"],
                       **{k: v for k, v in rec.items() if k.startswith("inner_")}, "note": "per-r OOF"},
                      branch=BRANCH, folds=range(5))
    g = pd.DataFrame(rows_)
    g.to_csv(OUT / "e2_grid.tsv", sep="\t", index=False)
    print(g.round(4).to_string(), flush=True)
    chosen = {f: int(g.loc[g[f"inner_{f}"].idxmax(), "r"]) for f in FOLDS}
    S = frame(ctx, {f: outer[chosen[f]][f] for f in FOLDS})
    r = b.full(S)
    r["chosen_r"] = chosen; r["grid"] = rows_; r["svd"] = sv
    r["delta_vs_core_v2"] = b.compare(S, S_ref, "E2")
    # pre-declared check: factors from non-GT allowed conditions only, same r*
    S_ng = frame(ctx, {f: score_use(f, f, chosen[f], nongt=True)[0] for f in FOLDS})
    r["check_nonGT_factors"] = {"M1-sel": b.c.ev.metric(c.ev.pairs(S_ng), "M1-sel")["value"],
                                "delta_vs_core_v2": delta(c.ev, S_ng, S_ref, "M1-sel", b.plan)}
    # calibration (fold-f nulls / permuted fold-f GT conditions, projected with the outer factors)
    Zn, Nn, Sn, Zp, Np_, Sp = [], [], [], [], [], []
    for f in FOLDS:
        z, n = score_use(f, f, chosen[f], which="perm")
        Zp.append(z); Np_.append(n); Sp.append(np.broadcast_to(b.sigmask[f], z.shape))
        ni = np.flatnonzero(ctx.null_fold == f)
        if len(ni):
            rows, V = factors_for([f])
            YNp = FB.project(YN, MN, rows, V, chosen[f])
            zc, n2 = b.zcore_custom("null", np.arange(len(c.nulls)), ctx.picks[f], YNp, MN, YNp, MN)
            zs = b.zsig_custom(b.fits[(f, f)], b.picks[f][0], YNp, MN, YNp, MN)
            z2 = SL.fuse(zc, zs, b.alpha[f])[ni]
            Zn.append(z2); Nn.append(n2[ni]); Sn.append(np.broadcast_to(b.sigmask[f], z2.shape))
    r["calibration"] = b.calib(np.vstack(Zn), np.vstack(Nn), np.vstack(Sn), np.vstack(Zp), np.vstack(Np_), np.vstack(Sp))
    d = r["delta_vs_core_v2"]
    r["adopt"] = bool(d["hard_test_pass"] and all(d["noninferiority"].values()))
    ledger_row("core-v2+E2 factor removal nested", dict(model="core-v2", r_grid=RS, chosen=chosen,
               site_frac=FB.SITE_FRAC), r)
    res["E2"] = r
    pickle.dump({"chosen": chosen, "factors": {k: v for k, v in factors.items() if not k[1]}}, open(CACHE / "e2.pkl", "wb"))
    return S


# ====================================================================== E3
def run_e3(b: B8, S_ref, res):
    ctx, c = b.ctx, b.c
    allidx = np.arange(len(c.conds))
    Z, Zp, rows_ = {}, {}, []
    for tau in TAUS:
        for rg in RGS:
            core = (tau, 3.0, rg)
            Z[core] = ctx.score(allidx, [("x", None)], core, [], which="gt")[("x", 0.0)][0]
            Zp[core] = ctx.score(allidx, [("x", None)], core, [], which="perm")[("x", 0.0)][0]
            S = pd.DataFrame(Z[core], index=c.conds, columns=c.kin)
            P = c.ev.pairs(S); fo = P["conditionID"].map(c.fold_of)
            rows_.append({"tau": tau, "rho_g": rg, **{f"inner_{f}": c.ev.metric(P[fo != f], "M1-sel")["value"] for f in FOLDS}})
    g = pd.DataFrame(rows_)
    g.to_csv(OUT / "e3_grid.tsv", sep="\t", index=False)
    parts, top, corr = {}, {}, {}
    for f in FOLDS:
        t3 = g.sort_values(f"inner_{f}", ascending=False).head(3)
        cfgs = [(float(x.tau), 3.0, float(x.rho_g)) for x in t3.itertuples()]
        assert cfgs[0] == tuple(ctx.picks[f]), (f, cfgs[0], ctx.picks[f])
        # top-1 reproduces the cached core-v1 z of fold f
        assert np.allclose(np.nan_to_num(Z[cfgs[0]][ctx.fold_pos[f]], nan=9e9), np.nan_to_num(b.zcore[(f, f)], nan=9e9), atol=1e-4)
        tp = b.train_pos(f)
        zb, cm = FB.bag_z([Z[k][ctx.fold_pos[f]].astype(np.float64) for k in cfgs], [Zp[k][tp].astype(np.float64) for k in cfgs])
        top[f] = [list(k) for k in cfgs]; corr[f] = cm.round(4).tolist()
        parts[f] = SL.fuse(zb, b.zs[(f, f)], b.alpha[f])
    S = frame(ctx, parts)
    r = b.full(S)
    r["top3"] = top; r["corr"] = corr
    r["delta_vs_core_v2"] = b.compare(S, S_ref, "E3")
    # calibration on fold-f perm GT conditions + nulls
    Zn, Nn, Sn, Zpp, Np_, Sp = [], [], [], [], [], []
    for f in FOLDS:
        cfgs = [tuple(x) for x in top[f]]; tp = b.train_pos(f); idx = ctx.fold_pos[f]; a, nt = b.alpha[f], b.picks[f][0]
        _, _, n = ctx.score(idx, [("x", None)], ctx.picks[f], [], which="perm", return_n=True)[("x", 0.0)]
        zb, _ = FB.bag_z([Zp[k][idx].astype(np.float64) for k in cfgs], [Zp[k][tp].astype(np.float64) for k in cfgs])
        z = SL.fuse(zb, b.L.zsig(b.fits[(f, f)], nt, "perm", idx)[0], a)
        Zpp.append(z); Np_.append(n); Sp.append(np.broadcast_to(b.sigmask[f], z.shape))
        ni = np.flatnonzero(ctx.null_fold == f)
        if len(ni):
            zl, nn = [], None
            for k in cfgs:
                zz, _, nn = ctx.score(ni, [("x", None)], k, [], which="null", return_n=True)[("x", 0.0)]; zl.append(zz.astype(np.float64))
            zb, _ = FB.bag_z(zl, [Zp[k][tp].astype(np.float64) for k in cfgs])
            z = SL.fuse(zb, b.L.zsig(b.fits[(f, f)], nt, "null", ni)[0], a)
            Zn.append(z); Nn.append(nn); Sn.append(np.broadcast_to(b.sigmask[f], z.shape))
    r["calibration"] = b.calib(np.vstack(Zn), np.vstack(Nn), np.vstack(Sn), np.vstack(Zpp), np.vstack(Np_), np.vstack(Sp))
    d = r["delta_vs_core_v2"]
    r["adopt"] = bool(d["hard_test_pass"] and all(d["noninferiority"].values()))
    ledger_row("core-v2+E3 top-3 core bagging", dict(model="core-v2", n_bag=3, rho_f=3.0, taus=TAUS, rgs=RGS,
               scale="sum z / sqrt(sum corr_perm)"), r)
    res["E3"] = r
    pickle.dump({"top": top, "Z": Z, "Zp": Zp}, open(CACHE / "e3.pkl", "wb"))
    return S


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--only", default="E1,E2,E3")
    which = ap.parse_args().only.split(",")
    b = B8()
    S_ref = b.core_v2()
    m = b.c.ev.metric(b.c.ev.pairs(S_ref), "M1-sel")["value"]
    assert abs(m - 0.6803082554884371) < 1e-6, m
    print("core-v2 reproduced", m, flush=True)
    res = {"core-v2": b.full(S_ref)}
    frames = {"core-v2": S_ref}
    fn = OUT / "b8_eval.json"
    old = json.loads(fn.read_text()) if fn.exists() else {}
    for e, run in [("E1", run_e1), ("E2", run_e2), ("E3", run_e3)]:
        if e in which:
            out = run(b, S_ref, res)
            frames[e] = out[0] if isinstance(out, tuple) else out
            old.setdefault("results", {})[e] = to_jsonable(res[e])
    old.setdefault("results", {})["core-v2"] = to_jsonable(res["core-v2"])
    old["gt_hash"] = b.ctx.gt_hash; old["alpha"] = b.alpha; old["ntop_kappa"] = b.picks
    old.setdefault("wall_s", {})["+".join(which)] = time.time() - b.t0
    fn.write_text(json.dumps(to_jsonable(old), indent=1, default=str))
    for k, S in frames.items():
        S.to_pickle(CACHE / f"S_{k}.pkl")
    print(f"done {time.time()-b.t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
