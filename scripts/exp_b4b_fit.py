#!/usr/bin/env python
"""Branch B4b: per outer fold f, fit learned edge-weight parameters theta (nested inner CV for lambda), score the
hand-specified promiscuity variant, and score DEV nulls / site-permuted conditions for calibration.
Writes results/b4b/fit_<fold>.npz (+ .json). One process per fold/GPU:
    CUDA_VISIBLE_DEVICES=i python scripts/exp_b4b_fit.py --fold fold<i>

PRE-SPECIFIED (written before the nested run; one exploratory look at fold1 with z_post chose the loss — see report):
  output     = z_vif (core-v1, D7). Training uses rho_e of the fold's core-v1 config; after fitting, rho_e is refitted
               on the 7 DEV null conditions for the weighted design (label-free) and the fold is scored with it.
  objective  = smooth-AUROC sigmoid pairwise loss, T = 0.5 (logistic loss reported as a fixed-lambda side variant)
  lambda     in {inf (theta = 0 = core-v1), 0.1, 0.01, 0.001}, chosen per outer fold by pooled inner-CV M1-sel
             (train on 3 folds, score the 4th; site statistics exclude outer and inner-validation fold)
  promiscuity variant: w = prom^-beta, beta in {0, .25, .5, .75, 1}; chosen per outer fold by pooled M1-sel on the
             other 4 folds (no fitted parameters, like the B2.1 grids)
"""
import argparse
import json
import time

import numpy as np
import pandas as pd
import torch

from exp_b4b_common import CORE_V1, EW, FOLDS, OUT, BCtx

LAMS = [0.1, 0.01, 0.001]
BETAS = [0.0, 0.25, 0.5, 0.75, 1.0]
MAX_ITER = 20


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fold", required=True)
    a = ap.parse_args()
    f = a.fold
    t0 = time.time()
    b = BCtx()
    c = b.ctx
    core, rep = b.cores[f], b.rep[f]
    rho0 = CORE_V1[f]["rho_e"]
    others = [g for g in FOLDS if g != f]
    out, info = {}, {"fold": f, "core": CORE_V1[f], "theta": {}, "rho_e": {}, "inner_m1": {}}

    def fit(train_folds, excl, lam, kind="sigmoid"):
        Phi = b.Phi(excl)
        th, h = EW.fit_theta(core, Phi, b.Yg, b.Mg, rep, c.ev, c.conds, b.conds_of(train_folds), lam=lam,
                             max_iter=MAX_ITER, kind=kind, out="z_vif", rho_e=rho0)
        return th, h, Phi

    def inner_m1(S):
        Pp = c.ev.pairs(pd.DataFrame(S, index=c.conds, columns=c.kin))
        return c.ev.metric(Pp[Pp["conditionID"].map(c.fold_of) != f], "M1-sel")["value"]

    def finish(name, th, Phi):
        """refit rho_e on nulls; store GT, null and perm scores."""
        r = b.rho_e(f, th, Phi)
        info["rho_e"][name] = r
        out[name] = b.score(f, th, Phi, rho_e=r)
        out[f"{name}__null"] = b.score(f, th, Phi, rho_e=r, which="null")
        out[f"{name}__perm"] = b.score(f, th, Phi, rho_e=r, which="perm")
        out[f"{name}__zpost"] = b.score(f, th, Phi, out="z_post")

    # core-v1 (theta = 0); its own rho_e (refit check vs grid value)
    Phi_f = b.Phi([f])
    finish("core", None, Phi_f)
    info["inner_m1"]["inf"] = inner_m1(out["core"])
    # inner CV: pooled OOF over the 4 other folds, per lambda
    for lam in LAMS:
        inner = np.full((len(c.conds), len(c.kin)), np.nan)
        for g in others:
            th, h, Phi = fit([x for x in others if x != g], [f, g], lam)
            r = b.rho_e(f, th, Phi)
            S = b.score(f, th, Phi, rho_e=r)
            inner[c.cf == g] = S[c.cf == g]
            info["theta"][f"inner_{g}_lam{lam}"] = th.tolist()
            print(f"[{f}] inner {g} lam {lam} {time.time()-t0:.0f}s loss {h[0]:.4f}->{h[-1]:.4f}", flush=True)
        info["inner_m1"][str(lam)] = inner_m1(inner)
    # outer fits (train on all 4 other folds; site stats exclude f)
    for lam in LAMS:
        th, h, Phi = fit(others, [f], lam)
        finish(f"outer_lam{lam}", th, Phi)
        info["theta"][f"outer_lam{lam}"] = th.tolist()
        print(f"[{f}] outer lam {lam} {time.time()-t0:.0f}s theta {np.round(th, 3)}", flush=True)
    th, h, Phi = fit(others, [f], 0.01, kind="logistic")
    finish("outer_logistic_lam0.01", th, Phi)
    info["theta"]["outer_logistic_lam0.01"] = th.tolist()
    # promiscuity variant (raw log promiscuity, no centring -> weight = prom^-beta exactly)
    Praw = torch.as_tensor(b.E[["log_prom"]].to_numpy(np.float32), device=b.device)
    for beta in BETAS:
        finish(f"prom_beta{beta}", np.array([-beta], np.float32), Praw)
    out["n_gt"], out["n_null"], out["n_perm"] = b.n, b.n_null, b.n_perm
    out["null_fold"] = b.null_fold.astype(str)
    np.savez_compressed(OUT / f"fit_{f}.npz", **out)
    info["wall_s"] = time.time() - t0
    info["features"] = EW.FEATURES
    (OUT / f"fit_{f}.json").write_text(json.dumps(info, indent=1))
    print(f"[{f}] done {time.time()-t0:.0f}s inner {info['inner_m1']}", flush=True)


if __name__ == "__main__":
    main()
