#!/usr/bin/env python
"""Branch B6 step 2/3: complementarity of components and nested-CV stacking vs core-v1 z_vif (D7).

  --stage comp   : per-component OOF metrics + per-condition rank correlations  -> results/b6/complementarity.json
  --stage stack  : pre-declared stackers (results/b6/predeclared.json, written before this stage is run)
                   -> results/b6/stack_eval.json, ledger branch B6.

Normalisation / recalibration / calibration are cross-fitted: for outer fold f, sigma_j and the null correlation R
come from site-permuted DEV GT conditions NOT in fold f (label-free), using fold f's version of each fold-specific
component. Calibration (|z|>1.96) is reported on fold-f permuted rows and on the DEV null conditions of fold f.
"""
import argparse
import json
import time

import numpy as np
import pandas as pd

from exp_b21_common import FOLDS, ROOT, Ctx, bs, calib_bins, full_eval, gt_freeze_hash  # noqa: F401
from exp_b4b_eval import unseen_mask
from kinact import ledger
from kinact.bootstrap import metric_replicates, shifted_null_p
from kinact.eval import to_jsonable
from kinact.explore import stacking as ST

OUT = ROOT / "results/b6"
COMP = ["zvif", "mean", "prom", "wzexp", "phosx", "kinlib", "curz", "self"]
PERFOLD = {"zvif", "mean", "prom", "wzexp"}


class Data:
    def __init__(self):
        z = np.load(ROOT / "cache/b6/components.npz", allow_pickle=True)
        self.z = z
        self.ctx = Ctx(ks={"r2np": [100]})
        c = self.ctx
        assert list(z["conds"]) == c.conds and list(z["kin"]) == c.kin and list(c.ev.kinases) == c.kin
        self.cf = c.cf
        self.null_fold = z["null_fold"]
        self.nulls = list(z["nulls"])
        # PSP-substrate counts for calibration bins (n of core-v1, from b4b fits)
        f0 = np.load(ROOT / "results/b4b/fit_fold0.npz")
        self.n_gt, self.n_null, self.n_perm = f0["n_gt"], f0["n_null"], f0["n_perm"]

    def raw(self, comp, which, f):
        a = self.z[f"{comp}__{which}"]
        return a[FOLDS.index(f)] if comp in PERFOLD else a

    def normalised(self, comps, f, transform="z"):
        """Per outer fold f: dict which -> list of normalised component matrices; sigma from perm rows not in f."""
        out = {w: [] for w in ["gt", "null", "perm"]}
        for cm in comps:
            P = self.raw(cm, "perm", f)
            if transform == "z":
                sig = ST.null_scale(P[self.cf != f])
                for w in out:
                    out[w].append(self.raw(cm, w, f) / sig)
            else:
                for w in out:
                    out[w].append(ST.rank_normal(self.raw(cm, w, f)))
        return out

    def oof(self, per_fold):
        A = np.full((len(self.ctx.conds), len(self.ctx.kin)), np.nan)
        for f, S in per_fold.items():
            A[self.cf == f] = S[self.cf == f]
        return pd.DataFrame(A, index=self.ctx.conds, columns=self.ctx.kin)


def metrics(D, S, S_ref=None, plan=None, um=None, rep_ref_u=None):
    c, ev = D.ctx, D.ctx.ev
    r = full_eval(c, S)
    Pp = ev.pairs(S)
    r["by_sign"] = {s: ev.metric(Pp, f"M1-sel:sign={s}")["value"] for s in ["-1", "1"]}
    if S_ref is not None:
        r["M2_vs_ref"] = ev.m2(S, S_ref)
        for ep in ["M1-sel", "M1-all", "harmony"]:
            d = bs.paired_delta(ev, S, S_ref, ep, plan)
            r[f"d{ep}"] = {k: d[k] for k in ["estimate", "ci_low", "ci_high", "p_one_sided"]}
        Pref = ev.pairs(S_ref)
        ru = metric_replicates(ev, Pp[um], "M1-sel", plan) - rep_ref_u
        pt = ev.metric(Pp[um], "M1-sel")["value"] - ev.metric(Pref[um], "M1-sel")["value"]
        q = np.quantile(ru, [0.025, 0.975])
        r["dT-unseen"] = {"estimate": pt, "ci_low": 2 * pt - q[1], "ci_high": 2 * pt - q[0],
                          "p_one_sided": shifted_null_p(pt, ru)}
    return r


def stage_comp(D):
    c, ev = D.ctx, D.ctx.ev
    res = {}
    oofs = {cm: D.oof({f: D.raw(cm, "gt", f) for f in FOLDS}) for cm in COMP}
    ref = oofs["zvif"]
    for cm in COMP:
        r = metrics(D, oofs[cm])
        r["M2_vs_zvif"] = ev.m2(oofs[cm], ref) if cm != "zvif" else None
        res[cm] = r
        print(cm, round(r["M1-sel"], 4), round(r["coverage"], 3), round(r["T-unseen"], 4), flush=True)
    # per-condition Spearman between components (cells scored by both), median over GT conditions
    rk = {cm: oofs[cm].rank(axis=1).to_numpy() for cm in COMP}
    rho = pd.DataFrame(np.eye(len(COMP)), index=COMP, columns=COMP)
    for i, a in enumerate(COMP):
        for b in COMP[i + 1:]:
            A, Bm = oofs[a].to_numpy(), oofs[b].to_numpy()
            vals = []
            for k in range(A.shape[0]):
                ok = np.isfinite(A[k]) & np.isfinite(Bm[k])
                if ok.sum() >= 10:
                    vals.append(pd.Series(A[k, ok]).corr(pd.Series(Bm[k, ok]), method="spearman"))
            rho.loc[a, b] = rho.loc[b, a] = float(np.nanmedian(vals))
    del rk
    # per-pair AUROC (scored-only) correlation with core-v1 on pairs where both score the target
    Pr = ev.pairs(ref)
    sel = Pr["in_M1sel"].to_numpy()
    pair_corr, wins = {}, {}
    for cm in COMP[1:]:
        Pc = ev.pairs(oofs[cm])
        ok = sel & Pr["target_scored"].to_numpy() & Pc["target_scored"].to_numpy()
        a, b = Pr["auroc"].to_numpy()[ok], Pc["auroc"].to_numpy()[ok]
        pair_corr[cm] = float(np.corrcoef(a, b)[0, 1])
        wins[cm] = {"n": int(ok.sum()), "comp_better_by_0.1": float(np.mean(b > a + 0.1)),
                    "core_better_by_0.1": float(np.mean(a > b + 0.1))}
    # null SDs (raw scale) per component (fold0 version)
    sig = {cm: ST.null_scale(D.raw(cm, "perm", "fold0")) for cm in COMP}
    out = {"gt_freeze_sha256": gt_freeze_hash(), "components": res, "spearman_median": rho.to_dict(),
           "pair_auroc_corr_with_zvif": pair_corr, "pair_wins_vs_zvif": wins, "perm_sd_fold0": sig}
    (OUT / "complementarity.json").write_text(json.dumps(to_jsonable(out), indent=1, default=str))
    print(rho.round(2).to_string())
    print(json.dumps(pair_corr, indent=0), json.dumps(wins, indent=0))


def stacked_fold(D, spec, f, w):
    comps = spec["components"]
    Nz = D.normalised(comps, f, spec.get("transform", "z"))
    R = ST.null_corr([P[D.cf != f] for P in Nz["perm"]])
    return {wh: ST.combine(Nz[wh], np.asarray(w, float), R) for wh in ["gt", "null", "perm"]}


def run_stacker(D, spec):
    c, ev = D.ctx, D.ctx.ev
    comps = spec["components"]
    J = len(comps)
    per, info = {}, {}
    if spec["kind"] == "average":
        for f in FOLDS:
            per[f] = stacked_fold(D, spec, f, np.ones(J))
            info[f] = {"w": [1.0] * J}
    else:
        w0 = np.array([1.0 if cm == "zvif" else 0.0 for cm in comps])
        lams = spec["lambdas"]
        for f in FOLDS:
            Nz = D.normalised(comps, f)
            Zgt = Nz["gt"]
            R = ST.null_corr([P[D.cf != f] for P in Nz["perm"]])
            T = [g for g in FOLDS if g != f]
            inner = {}
            for lam in lams:
                if lam == "inf":
                    preds = {g: ST.combine(Zgt, w0, R) for g in T}
                else:
                    preds = {}
                    for g in T:
                        tr = [cc for cc, ff in zip(c.conds, D.cf) if ff in T and ff != g]
                        wg = ST.fit_weights(ST.pair_data(ev, c.conds, tr, Zgt), w0, float(lam), R)
                        preds[g] = ST.combine(Zgt, wg, R)
                S = D.oof({g: preds[g] for g in T})
                Pp = ev.pairs(S)
                Pp = Pp[Pp["conditionID"].map(c.fold_of).isin(T)]
                inner[lam] = ev.metric(Pp, "M1-sel")["value"]
            best = max(lams, key=lambda l: inner[l])
            if best == "inf":
                wf = w0
            else:
                tr = [cc for cc, ff in zip(c.conds, D.cf) if ff != f]
                wf = ST.fit_weights(ST.pair_data(ev, c.conds, tr, Zgt), w0, float(best), R)
            # also record weights for every finite lambda (descriptive)
            wall = {}
            for lam in lams:
                if lam != "inf":
                    tr = [cc for cc, ff in zip(c.conds, D.cf) if ff != f]
                    wall[lam] = ST.fit_weights(ST.pair_data(ev, c.conds, tr, Zgt), w0, float(lam), R).tolist()
            per[f] = stacked_fold(D, spec, f, wf)
            info[f] = {"inner": inner, "lambda": best, "w": np.asarray(wf).tolist(), "w_all": wall}
            print(f"  {spec['name']} {f} inner {({k: round(v, 4) for k, v in inner.items()})} -> {best} "
                  f"w {np.round(wf, 3)}", flush=True)
    S = D.oof({f: per[f]["gt"] for f in FOLDS})
    Zn = np.vstack([per[f]["null"][D.null_fold == f] for f in FOLDS])
    Nn = np.vstack([D.n_null[D.null_fold == f] for f in FOLDS])
    Zp = np.vstack([per[f]["perm"][D.cf == f] for f in FOLDS])
    Np = np.vstack([D.n_perm[D.cf == f] for f in FOLDS])
    return S, {"calib_null": calib_bins(Zn, Nn), "calib_perm": calib_bins(Zp, Np), "folds": info}


def stage_stack(D):
    spec_all = json.loads((OUT / "predeclared.json").read_text())
    c, ev = D.ctx, D.ctx.ev
    plan = bs.BootstrapPlan.from_evaluator(ev, B=2000, seed=0)
    S_ref = D.oof({f: D.raw("zvif", "gt", f) for f in FOLDS})
    Pref = ev.pairs(S_ref)
    um = unseen_mask(c, Pref)
    rep_ref_u = metric_replicates(ev, Pref[um], "M1-sel", plan)
    res = {"core-v1 z_vif (ref)": metrics(D, S_ref)}
    Zn = np.vstack([D.raw("zvif", "null", f)[D.null_fold == f] for f in FOLDS])
    Zp = np.vstack([D.raw("zvif", "perm", f)[D.cf == f] for f in FOLDS])
    res["core-v1 z_vif (ref)"].update(
        calib_null=calib_bins(Zn, np.vstack([D.n_null[D.null_fold == f] for f in FOLDS])),
        calib_perm=calib_bins(Zp, np.vstack([D.n_perm[D.cf == f] for f in FOLDS])))
    print("ref", round(res["core-v1 z_vif (ref)"]["M1-sel"], 4), flush=True)
    for spec in spec_all["stackers"]:
        t0 = time.time()
        S, extra = run_stacker(D, spec)
        r = metrics(D, S, S_ref, plan, um, rep_ref_u)
        r.update(extra)
        r["spec"] = spec
        res[spec["name"]] = r
        ledger.append("b6_stack", spec, "dev", {k: r[k] for k in ["M1-sel", "M1-all", "M1-sel-scored", "coverage"]},
                      branch="B6", folds=range(5))
        S.to_parquet(OUT / f"oof_{spec['name']}.parquet")
        print(spec["name"], round(r["M1-sel"], 4), r["dM1-sel"], r["calib_null"]["all"], r["calib_perm"]["all"],
              f"{time.time()-t0:.0f}s", flush=True)
    out = {"gt_freeze_sha256": gt_freeze_hash(), "results": res}
    (OUT / "stack_eval.json").write_text(json.dumps(to_jsonable(out), indent=1, default=str))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    D = Data()
    {"comp": stage_comp, "stack": stage_stack}[a.stage](D)
