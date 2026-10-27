#!/usr/bin/env python
"""B7 grid: fold-dependent signature fits and z_vif scores for every (variant, outer fold, eval fold, kappa, n_top,
alpha). Writes cache/b7/grid_<variant>.pkl (scores) and results/b7/fits_<variant>.json (per-fit kinase counts).
Signature for (outer f, inner eval fold g) is fitted on folds not in {f, g}; for outer scoring of f on folds != f.
"""
import argparse
import itertools
import json
import pickle
import time

import numpy as np

from exp_b7_common import B7Ctx, CACHE, FOLDS, OUT, ALPHAS, NTOPS, KAPPAS, VARIANTS


def fold_conds(ctx, folds):
    m = ctx.c.meta
    return list(m.index[m["split"].isin(folds)])


def run_variant(ctx, variant, store, fits_info):
    t0 = time.time()
    c = ctx.c
    jobs = []      # (train_folds, [(outer_f, eval_g)])
    for f, g in itertools.combinations(FOLDS, 2):
        jobs.append(([x for x in FOLDS if x not in (f, g)], [f, g], [(f, g), (g, f)]))
    for f in FOLDS:
        jobs.append(([x for x in FOLDS if x != f], [f], [(f, f)]))
    for tr, held, uses in jobs:
        fits = ctx.fit(tr, fold_conds(ctx, held), variant)
        any_fit = fits[KAPPAS[0]]
        fits_info[f"{'+'.join(held)}"] = {"train_folds": tr, "n_kinases": len(any_fit.sites),
                                          "n_train_conds": len(any_fit.train_conds),
                                          "kinases": sorted(any_fit.sites), "n_groups": any_fit.n_groups}
        if len(held) == 1:
            with open(CACHE / f"sig_{variant}_{held[0]}.pkl", "wb") as fh:
                pickle.dump(fits, fh)
        for f_out, g_ev in uses:
            core = ctx.picks[f_out]
            idx = ctx.fold_pos[g_ev]
            Wl = [(("core",), None)]
            for kp, nt in itertools.product(KAPPAS, NTOPS):
                Wl.append(((kp, nt), fits[kp].matrix(len(c.X.index), c.kin, nt)))
            res = ctx.score(idx, Wl, core, ALPHAS)
            for (lab, a), (Z, rho) in res.items():
                if lab == ("core",):
                    store[("core", f_out, g_ev)] = (Z, rho)
                else:
                    store[(variant, f_out, g_ev, lab[0], lab[1], a)] = (Z, rho)
        print(f"[{variant}] held {held}: {len(any_fit.sites)} kinases, {len(any_fit.train_conds)} train conds; "
              f"{time.time()-t0:.0f}s", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variants", default=",".join(VARIANTS))
    args = ap.parse_args()
    ctx = B7Ctx()
    for v in args.variants.split(","):
        store, info = {}, {}
        run_variant(ctx, v, store, info)
        with open(CACHE / f"grid_{v}.pkl", "wb") as fh:
            pickle.dump(store, fh)
        (OUT / f"fits_{v}.json").write_text(json.dumps(info, indent=1, default=str))


if __name__ == "__main__":
    main()
