#!/usr/bin/env python
"""Branch 2 (Bayesian linear-Gaussian core): hyperparameter grid + nested CV on DEV folds 0-4.

Every condition c is scored with site noise variances v estimated on allowed DEV conditions EXCLUDING c's own
fold (qc_status == ok), so a config's score of c does not depend on the outer split; nested CV then selects, for
each outer fold f, the config maximising pooled M1-sel over the other four folds, and the pooled out-of-fold (OOF)
scores are evaluated once. Every (config, output, noise) evaluation on DEV is appended to results/dev_ledger.tsv
(branch B2). Outputs: results/b2/grid.tsv, results/b2/oof_<variant>.parquet, results/b2/nested.json.
"""
import argparse
import hashlib
import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from kinact import ledger  # noqa: E402
from kinact.data import load_dev, load_dev_meta  # noqa: E402
from kinact.eval import Evaluator, load_gt, read_gt_tsv  # noqa: E402
from kinact.explore import bayes as B  # noqa: E402
from kinact.library import load_library  # noqa: E402

FOLDS = [f"fold{i}" for i in range(5)]
OUT = ROOT / "results/b2"
CACHE = ROOT / "cache/b2"
LIB_CUR, LIB_PRED = "curated_union", ("r2np", 100)
GRID = dict(tau=[0.25, 0.5, 1.0, 2.0, 4.0], rho_f=[0.0, 0.5, 1.0, 3.0], rho_g=[0.0, 0.25, 1.0],
            alpha=[0.0, 0.25, 0.5, 1.0], intercept=["median", "joint"])


def gt_freeze_hash():
    return hashlib.sha256((ROOT / "results/ground_truth/GT_FREEZE.json").read_bytes()).hexdigest()


def setup():
    P, C, U, _ = load_gt()
    bins = read_gt_tsv(ROOT / "results/ground_truth/harmony_bins_group.tsv")
    Udf = read_gt_tsv(ROOT / "results/ground_truth/kinase_universe.tsv")
    ev = Evaluator(P, C, U, split_filter=FOLDS, bins=bins, zero_policy="ignore")
    meta = load_dev_meta()
    X = load_dev()
    return P, C, Udf, ev, meta, X


def noise_variances(X, meta, refresh=False):
    """v per held-out fold (allowed = DEV, qc ok, not that fold) + v_all (all DEV qc ok)."""
    CACHE.mkdir(parents=True, exist_ok=True)
    split = meta.loc[X.columns, "split"].to_numpy()
    ok = (meta.loc[X.columns, "qc_status"] == "ok").to_numpy()
    out = {}
    for f in FOLDS + ["all"]:
        p = CACHE / f"v_{f}.npy"
        if p.exists() and not refresh:
            out[f] = np.load(p); continue
        allowed = ok & (split != f)
        out[f] = B.site_noise_variance(X.to_numpy()[:, allowed])
        np.save(p, out[f])
    return out


def standardise_by_fold(X, conds, meta, vs, noise=True, v_override=None):
    V = X[conds].to_numpy(np.float64)
    Y = np.empty_like(V); D = np.empty_like(V)
    fold = meta.loc[conds, "split"].to_numpy()
    for f in FOLDS:
        j = np.flatnonzero(fold == f)
        if not len(j):
            continue
        v = np.ones(V.shape[0]) if not noise else (vs[f] if v_override is None else v_override)
        Y[:, j], D[:, j], _, _ = B.standardise(V[:, j], v)
    return Y, D


def inner_metrics(ev, P, fold_of):
    """M1-sel on the complement of each fold (selection criterion) + pooled metrics."""
    f = P["conditionID"].map(fold_of).to_numpy()
    res = {f"inner_{g}": ev.metric(P[f != g], "M1-sel")["value"] for g in FOLDS}
    for ep in ["M1-sel", "M1-all", "M1-sel-scored"]:
        res[ep] = ev.metric(P, ep)["value"]
    sel = P[P["in_M1sel"]]
    res["coverage"] = float(sel["target_scored"].mean())
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--noise", choices=["on", "off", "both"], default="both")
    ap.add_argument("--no-ledger", action="store_true")
    ap.add_argument("--taus", default=None, help="comma list overriding GRID['tau'] (grid extension)")
    ap.add_argument("--tag", default="")
    ap.add_argument("--rho_f", default=None, help="comma list overriding GRID['rho_f'] (grid extension)")
    args = ap.parse_args()
    if args.taus:
        GRID["tau"] = [float(t) for t in args.taus.split(",")]
    if args.rho_f:
        GRID["rho_f"] = [float(t) for t in args.rho_f.split(",")]
    t0 = time.time()
    P, C, Udf, ev, meta, X = setup()
    OUT.mkdir(parents=True, exist_ok=True)
    vs = noise_variances(X, meta)
    fold_of = meta["split"]
    conds = list(ev.cond_index)
    design = B.build_design(X.index, list(Udf["kinase"]), load_library(LIB_CUR), load_library(LIB_PRED[0], k=LIB_PRED[1]))
    print(f"setup {time.time()-t0:.0f}s; {len(conds)} conds; lib rows {len(design.lib_rows)}; GT_FREEZE {gt_freeze_hash()[:16]}",
          flush=True)
    rows, store = [], {}
    noises = {"on": [True], "off": [False], "both": [True, False]}[args.noise]
    for noise in noises:
        Y, D = standardise_by_fold(X, conds, meta, vs, noise=noise)
        for alpha in GRID["alpha"]:
            ts = time.time()
            G, h, n = B.suff_stats(design, alpha, Y, D)
            n = n.cpu().numpy()
            for tau, rf, rg, ic in itertools.product(GRID["tau"], GRID["rho_f"], GRID["rho_g"], GRID["intercept"]):
                K = B.kinase_kernel(Udf, design.kinases, rf, rg)
                tp = time.time()
                mu, sd = B.posterior(G, h, K, tau, ic)
                t_post = time.time() - tp
                for borrow in ([False, True] if (rf or rg) else [False]):
                    rep = B.report_mask(n, K, 1, borrow)
                    for output in ["z", "mean"]:
                        arr = np.where(rep, mu / sd if output == "z" else mu, np.nan).astype(np.float32)
                        S = pd.DataFrame(arr, index=conds, columns=design.kinases)
                        Pp = ev.pairs(S)
                        m = inner_metrics(ev, Pp, fold_of)
                        cfg = dict(model="bayes_core_v1", lib_cur=LIB_CUR, lib_pred=f"{LIB_PRED[0]}@{LIB_PRED[1]}",
                                   tau=tau, rho_f=rf, rho_g=rg, alpha=alpha, intercept=ic, borrow=borrow,
                                   output=output, noise=noise, min_sub=1, v="leave-own-fold-out")
                        key = ledger.config_hash(cfg)
                        rows.append({**cfg, "config_hash": key, **m, "t_post_ms_per_cond": 1e3 * t_post / len(conds)})
                        store[key] = arr
                        if not args.no_ledger:
                            ledger.append("bayes_core", cfg, "dev", {k: m[k] for k in ["M1-sel", "M1-all", "M1-sel-scored", "coverage"]},
                                          branch="B2", folds=range(5))
            print(f"noise={noise} alpha={alpha}: {time.time()-ts:.0f}s, best so far M1-sel "
                  f"{max(r['M1-sel'] for r in rows):.3f}", flush=True)
    grid = pd.DataFrame(rows)
    grid.to_csv(OUT / f"grid_noise-{args.noise}{args.tag}.tsv", sep="\t", index=False)
    np.savez_compressed(CACHE / f"grid_scores_noise-{args.noise}{args.tag}.npz", conds=np.array(conds),
                        kinases=np.array(design.kinases), **store)
    print(f"done {time.time()-t0:.0f}s, {len(grid)} evaluations", flush=True)


if __name__ == "__main__":
    main()
