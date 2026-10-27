#!/usr/bin/env python
"""Branch B4 grid: co-regulation expansion (idea 1/2) and co-regulation imputation (idea 3) on top of the B2
core (fixed B2 recommended config: tau .125, rho_f 3, rho_g .25, alpha 1 [curated_union + r2np@100], median,
borrow, no noise model) and with plain weighted z.

Every DEV GT condition of fold f is scored with co-regulation statistics fitted on the allowed conditions of fold f
only (cache/b4/*_<fold>.npz from scripts/exp_b4_precompute.py; coreg.allowed_mask). Outputs per config: z (B2
posterior z), mean (B2 posterior mean), wz. Every evaluation is appended to results/dev_ledger.tsv (branch B4).
Writes results/b4/grid.tsv and cache/b4/grid_scores.npz.
"""
import argparse
import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))
import exp_b2_grid as E2  # noqa: E402
from kinact import ledger  # noqa: E402
from kinact.explore import bayes as B  # noqa: E402
from kinact.explore import coreg as CR  # noqa: E402
from kinact.library import load_library  # noqa: E402

FOLDS = E2.FOLDS
CACHE = ROOT / "cache/b4"
OUT = ROOT / "results/b4"
CORE = dict(tau=0.125, rho_f=3.0, rho_g=0.25, intercept="median", borrow=True)
OUTPUTS = ["z", "mean", "wz"]


class Ctx:
    def __init__(self):
        self.P, self.C, self.Udf, self.ev, self.meta, X = E2.setup()
        self.sites = np.load(CACHE / "sites.npy", allow_pickle=True)
        self.X = X.loc[self.sites]
        self.kin = list(self.Udf["kinase"])
        self.conds = list(self.ev.cond_index)
        self.fold_of = self.meta["split"]
        cf = pd.Series(self.conds).map(self.fold_of).to_numpy()
        self.fold_conds = {f: [c for c, g in zip(self.conds, cf) if g == f] for f in FOLDS}
        cur = load_library("curated_union")[["kinase", "site", "weight"]]
        r2 = load_library("r2np", k=100)[["kinase", "site"]].assign(weight=1.0)
        r2 = r2.merge(cur[["kinase", "site"]], how="left", indicator=True)
        r2 = r2[r2["_merge"] == "left_only"].drop(columns="_merge")
        self.base = pd.concat([cur, r2], ignore_index=True)            # = B2 library at alpha = 1
        self.K = B.kinase_kernel(self.Udf, self.kin, CORE["rho_f"], CORE["rho_g"])
        self.YD = {}
        for f in FOLDS:
            V = self.X[self.fold_conds[f]].to_numpy(np.float64)
            Y, D, _, _ = B.standardise(V, np.ones(V.shape[0]))
            Yn = np.where(D > 0, Y, np.nan)
            self.YD[f] = (Y, D, Yn)
        self._top, self._nbr = {}, {}

    def top(self, seed, q, f):
        k = (seed, q, f)
        if k not in self._top:
            z = np.load(CACHE / f"top_{seed}_q{q}_{f}.npz")
            assert list(z["kinases"]) == self.kin
            self._top[k] = {x: z[x] for x in ["idx", "r", "n", "rs"]}
        return self._top[k]

    def nbr(self, q, f):
        k = (q, f)
        if k not in self._nbr:
            z = np.load(CACHE / f"nbr_q{q}_{f}.npz")
            self._nbr[k] = {x: z[x] for x in ["rows", "nbr", "r", "n"]}
        return self._nbr[k]


def score_fold(ctx, design, beta, Y, D, conds):
    W = design.W(beta)
    S = W.shape[0]
    chunk = max(1, int(24 * 27000 / max(S, 1)))
    G, h, n = B.suff_stats(design, beta, Y, D, chunk=chunk)
    mu, sd = B.posterior(G, h, ctx.K, CORE["tau"], CORE["intercept"])
    n = n.cpu().numpy()
    rep = B.report_mask(n, ctx.K, 1, CORE["borrow"])
    rows = design.lib_rows
    Dl, Yl = D[rows], Y[rows]
    W2 = W.multiply(W).tocsr()
    with np.errstate(all="ignore"):
        wz = (W.T @ (Dl * Yl)).T / np.sqrt((W2.T @ Dl).T)
    Wb = (W != 0).astype(float)
    nm = (Wb.T @ (Dl > 0).astype(float)).T
    out = {"z": np.where(rep, mu / sd, np.nan), "mean": np.where(rep, mu, np.nan),
           "wz": np.where((nm >= 1) & np.isfinite(wz), wz, np.nan)}
    return {k: pd.DataFrame(v, index=conds, columns=ctx.kin) for k, v in out.items()}, n


def score_config(ctx, cfg):
    """cfg: kind 'none' | 'expand' (seed, q, m, sign, beta, rmin, weight) | 'impute' (q, rmin, min_nbr, gamma)."""
    parts = {o: [] for o in OUTPUTS}
    nexp = []
    for f in FOLDS:
        Y, D, Yn = ctx.YD[f]
        conds = ctx.fold_conds[f]
        beta = 0.0
        exp = None
        if cfg["kind"] == "expand":
            exp = CR.expansion_edges(ctx.top(cfg["seed"], cfg["q"], f), ctx.kin, ctx.sites, cfg["m"], cfg["rmin"],
                                     cfg["sign"], cfg.get("weight", "rs"))
            beta = cfg["beta"]
            nexp.append(len(exp))
        design = CR.make_design(ctx.sites, ctx.kin, ctx.base, exp)
        if cfg["kind"] == "impute":
            nb = ctx.nbr(cfg["q"], f)
            yhat, qrel = CR.impute_from_neighbors(Yn, nb["rows"], nb["nbr"], nb["r"], cfg["min_nbr"], cfg["rmin"])
            miss = (D[nb["rows"]] == 0) & ~np.isnan(yhat)
            Y = Y.copy(); D = D.copy()
            Y[nb["rows"]] = np.where(miss, yhat, Y[nb["rows"]])
            D[nb["rows"]] = np.where(miss, cfg["gamma"], D[nb["rows"]])
            nexp.append(int(miss.sum()) / len(conds))
        res, _ = score_fold(ctx, design, beta if beta else 1.0, Y, D, conds) if exp is not None else \
            score_fold(ctx, design, 0.0, Y, D, conds)
        for o in OUTPUTS:
            parts[o].append(res[o])
    return {o: pd.concat(parts[o]).loc[ctx.conds] for o in OUTPUTS}, float(np.mean(nexp)) if nexp else 0.0


def grid_configs(which):
    cfgs = []
    if which in ("all", "ref"):
        cfgs.append(dict(kind="none"))
    if which in ("all", "expand"):
        for seed, q, m, sign, beta in itertools.product(["cur", "cur+r2np"], [0, 1, 3], [10, 25, 50, 100, 200],
                                                        ["signed", "pos"], [0.25, 0.5, 1.0]):
            cfgs.append(dict(kind="expand", seed=seed, q=q, m=m, sign=sign, beta=beta, rmin=0.0, weight="rs"))
    if which in ("all", "impute"):
        for q, rmin, mn, gamma in itertools.product([0, 3], [0.3, 0.5], [1, 3], [0.25, 0.5, 1.0]):
            cfgs.append(dict(kind="impute", q=q, rmin=rmin, min_nbr=mn, gamma=gamma))
    return cfgs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--which", default="all")
    ap.add_argument("--tag", default="")
    ap.add_argument("--no-ledger", action="store_true")
    args = ap.parse_args()
    t0 = time.time()
    ctx = Ctx()
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"setup {time.time()-t0:.0f}s GT_FREEZE {E2.gt_freeze_hash()[:16]}", flush=True)
    rows, store = [], {}
    for i, cfg in enumerate(grid_configs(args.which)):
        ts = time.time()
        S, nexp = score_config(ctx, cfg)
        for o in OUTPUTS:
            full = dict(model="b4_coreg_v1", core=CORE, base="curated_union+r2np@100", output=o, **cfg)
            key = ledger.config_hash(full)
            m = E2.inner_metrics(ctx.ev, ctx.ev.pairs(S[o]), ctx.fold_of)
            rows.append({**{k: (json.dumps(v) if isinstance(v, dict) else v) for k, v in full.items()},
                         "config_hash": key, **m, "n_exp_per_fold": nexp})
            store[key] = S[o].to_numpy(np.float32)
            if not args.no_ledger:
                ledger.append("b4_coreg", full, "dev", {k: m[k] for k in ["M1-sel", "M1-all", "M1-sel-scored", "coverage"]},
                              branch="B4", folds=range(5))
        if i % 10 == 0:
            best = max(rows, key=lambda r: r["M1-sel"])
            print(f"{i} {cfg} {time.time()-ts:.1f}s | M1-sel z/mean/wz "
                  f"{rows[-3]['M1-sel']:.3f}/{rows[-2]['M1-sel']:.3f}/{rows[-1]['M1-sel']:.3f} | best {best['M1-sel']:.3f}",
                  flush=True)
    g = pd.DataFrame(rows)
    g.to_csv(OUT / f"grid{args.tag}.tsv", sep="\t", index=False)
    np.savez_compressed(CACHE / f"grid_scores{args.tag}.npz", conds=np.array(ctx.conds), kinases=np.array(ctx.kin),
                        **store)
    print(f"done {time.time()-t0:.0f}s, {len(g)} evaluations", flush=True)


if __name__ == "__main__":
    main()
