#!/usr/bin/env python
"""Branch B6 step 1: assemble component score matrices on DEV GT conditions (gt), the 7 DEV null conditions (null)
and site-permuted DEV GT conditions (perm, seed 0, within condition; the same permuted input for every component).

Components (all label-free given the fold; fold-specific ones stored per outer fold f = config chosen for f by the
producing branch's nested CV on folds != f, applied to every condition):
  zvif   core-v1 z_vif (results/b4b/fit_<f>.npz 'core'; tau/rho_g per fold as D7)          [per fold]
  mean   core-v1 posterior mean (same model, recomputed with exp_b4b_common.BCtx)          [per fold]
  prom   promiscuity-weighted core, beta = 0.25 fixed (fit_<f>.npz 'prom_beta0.25')        [per fold]
  wzexp  B4 weighted z + co-regulation expansion; config = B4 nested pick for fold f (cache/b4/grid_scores.npz);
         rows of fold g always use the expansion fitted without fold g (allowed conditions only)  [per fold]
  phosx  PhosX-like, cantley-site-top15 min3 (cache/scores)                                 [fixed]
  kinlib KinLib enrichment, cantley-site-top15 thr0.5 min3 (cache/scores)                   [fixed]
  curz   z-score curated_union min1 (cache/scores)                                          [fixed]
  self   self-site signal, PSP regulatory self-sites (B5), RMS-standardised per condition    [fixed]
Writes cache/b6/components.npz.
"""
import json
import os
import sys
import time

import numpy as np
import pandas as pd

os.environ.setdefault("CUDA_VISIBLE_DEVICES", "1")
from exp_b4b_common import FOLDS, ROOT, BCtx  # noqa: E402
from kinact.explore import bayes as B  # noqa: E402
from kinact.explore import coreg as CR  # noqa: E402
from kinact.explore import selfsites as SS  # noqa: E402
from kinact.library import load_library  # noqa: E402
from kinact.methods import METHODS  # noqa: E402

CACHE = ROOT / "cache/b6"
FIX = {"phosx": ("phosx_like", 15, dict(min_substrates=3), "phosx_like__cantley-site-top15__min3"),
       "kinlib": ("kl_enrich", 15, dict(min_substrates=3, thr=0.5), "kl_enrich__cantley-site-top15__thr0.5_min3"),
       "curz": ("zscore", "curated_union", dict(min_substrates=1), "zscore__curated_union__min1")}


def main():
    t0 = time.time()
    CACHE.mkdir(parents=True, exist_ok=True)
    b = BCtx(device="cuda")
    c = b.ctx
    conds, kin, nulls = c.conds, c.kin, c.nulls
    assert list(c.ev.kinases if hasattr(c.ev, "kinases") else kin) == kin
    out = {"conds": np.array(conds), "kin": np.array(kin), "nulls": np.array(nulls),
           "cf": c.cf, "null_fold": pd.Series(nulls).map(c.fold_of).to_numpy()}
    fits = {f: np.load(ROOT / f"results/b4b/fit_{f}.npz") for f in FOLDS}
    # ---- Bayesian-core components (per fold)
    for name, key in [("zvif", "core"), ("prom", "prom_beta0.25")]:
        out[f"{name}__gt"] = np.stack([fits[f][key] for f in FOLDS])
        out[f"{name}__null"] = np.stack([fits[f][key + "__null"] for f in FOLDS])
        out[f"{name}__perm"] = np.stack([fits[f][key + "__perm"] for f in FOLDS])
    for which in ["gt", "null", "perm"]:
        out[f"mean__{which}"] = np.stack([b.score(f, None, None, out="mean", which=which) for f in FOLDS])
    # sanity: recomputed z_vif equals stored
    chk = b.score("fold1", None, None, out="z_vif")
    print("zvif recompute max|d|", np.nanmax(np.abs(chk - fits["fold1"]["core"])), flush=True)
    print(f"[core] {time.time()-t0:.0f}s", flush=True)
    # ---- inputs shared by the remaining components
    X = c.X
    Xg, Xn, Xp = X[conds], X[nulls], c.Xp
    assert list(Xp.columns) == conds
    # ---- B4 wz + co-regulation expansion
    import exp_b4_grid as G4
    g4 = G4.Ctx()
    assert g4.conds == conds and g4.kin == kin
    ev4 = json.load(open(ROOT / "results/b4/b4_eval.json"))["chosen"]["B4exp_wz"]
    grid = pd.read_csv(ROOT / "results/b4/grid.tsv", sep="\t")
    store = np.load(ROOT / "cache/b4/grid_scores.npz", allow_pickle=True)
    sites = g4.sites
    fold_of = c.fold_of

    def wz(cfg, V, colfold):
        """weighted z (exp_b4_grid.score_fold 'wz') for columns V (sites x n), each column with its own fold's fit."""
        R = np.full((V.shape[1], len(kin)), np.nan)
        for g in FOLDS:
            j = np.flatnonzero(colfold == g)
            if not len(j):
                continue
            exp = CR.expansion_edges(g4.top(cfg["seed"], int(cfg["q"]), g), kin, sites, int(cfg["m"]), 0.0,
                                     cfg["sign"], "rs")
            d = CR.make_design(sites, kin, g4.base, exp)
            W = d.W(float(cfg["beta"]))
            Y, D, _, _ = B.standardise(V[:, j], np.ones(V.shape[0]))
            Dl, Yl = D[d.lib_rows], Y[d.lib_rows]
            W2 = W.multiply(W).tocsr()
            with np.errstate(all="ignore"):
                z = (W.T @ (Dl * Yl)).T / np.sqrt((W2.T @ Dl).T)
            Wb = (W != 0).astype(float)
            nm = (Wb.T @ (Dl > 0).astype(float)).T
            R[j] = np.where((nm >= 1) & np.isfinite(z), z, np.nan)
        return R

    for which, V, cols in [("gt", Xg, conds), ("null", Xn, nulls), ("perm", Xp, conds)]:
        Vs = V.reindex(sites).to_numpy(np.float64)
        colfold = pd.Series(cols).map(fold_of).to_numpy()
        arrs = []
        for f in FOLDS:
            cfg = ev4[f]
            A = wz(cfg, Vs, colfold)
            if which == "gt":
                sub = grid[(grid["kind"] == "expand") & (grid["output"] == "wz") & (grid["seed"] == cfg["seed"])
                           & (grid["q"] == cfg["q"]) & (grid["m"] == cfg["m"]) & (grid["sign"] == cfg["sign"])
                           & (grid["beta"] == cfg["beta"])]
                assert len(sub) == 1
                ref = store[sub["config_hash"].iloc[0]]
                print(f"wzexp {f} reproduce max|d| {np.nanmax(np.abs(ref - A)):.2e} nan-mismatch "
                      f"{int((np.isnan(ref) != np.isnan(A)).sum())}", flush=True)
            arrs.append(A)
        out[f"wzexp__{which}"] = np.stack(arrs)
    print(f"[wzexp] {time.time()-t0:.0f}s", flush=True)
    # ---- fixed baselines (gt/null from cache/scores; perm recomputed)
    for name, (m, lib, kw, fn) in FIX.items():
        S = pd.read_parquet(ROOT / f"cache/scores/{fn}.parquet")
        out[f"{name}__gt"] = S.reindex(index=conds, columns=kin).to_numpy(np.float64)
        out[f"{name}__null"] = S.reindex(index=nulls, columns=kin).to_numpy(np.float64)
        chk = METHODS[m](Xg.iloc[:, :20], lib, **kw).reindex(columns=kin).to_numpy(np.float64)
        print(f"{name} recompute (20 conds) max|d| {np.nanmax(np.abs(chk - out[f'{name}__gt'][:20])):.2e}", flush=True)
        out[f"{name}__perm"] = METHODS[m](Xp, lib, **kw).reindex(index=conds, columns=kin).to_numpy(np.float64)
        print(f"[{name}] {time.time()-t0:.0f}s", flush=True)
    # ---- self-site signal (B5, psp set)
    tab = SS.self_site_table(c.Udf, "psp")
    Wself = SS.self_matrix(tab, X.index, kin)
    for which, V in [("gt", Xg), ("null", Xn), ("perm", Xp)]:
        Y, D, _, _ = B.standardise(V.reindex(X.index).to_numpy(np.float64), np.ones(len(X.index)))
        s, _ = SS.self_signal(Y, D, Wself, True)
        out[f"self__{which}"] = s
    np.savez_compressed(CACHE / "components.npz", **out)
    print(f"done {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
