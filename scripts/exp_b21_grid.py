#!/usr/bin/env python
"""Branch B2.1 (+B3) grids on DEV folds 0-4 (one experiment per --exp). Every evaluation -> results/dev_ledger.tsv
(branch B2.1); grid rows -> results/b21/grid_<exp>.tsv. Scores are recomputed for nested-CV choices in
exp_b21_eval.py (no score archives).

Experiments (base = B2 core: curated_union + r2np@100 alpha=1, median-centred, borrow, no noise model):
  out    items 1/2/5: output statistic x winsor x tau x rho_f x rho_g x alpha_r2np
  src    item 3: + r3 / cantley global top-k edge sources (alpha_r3, alpha_cantley, k per source)
  kern   item 4: + label-free substrate-profile similarity kernel S (r2np or cantley top-500 cosine^2), rho_s
Output statistics: see kinact.explore.bayes2 (z_post = B2; mean; z_freq; hybrid l; z_vif with rho_e fitted on the 7
DEV role=null conditions per config; empnull = per-condition empirical null).
"""
import argparse
import itertools
import time

import numpy as np
import pandas as pd

from exp_b21_common import B2, Ctx, OUT, inner_metrics, ledger_row
from kinact.explore.bayes import report_mask

OUTPUTS = [("z_post", None), ("mean", None), ("z_freq", None), ("hybrid", 0.25), ("hybrid", 0.5), ("hybrid", 0.75),
           ("z_vif", None), ("empnull", None)]


def out_name(o, l):
    return f"hybrid{l:g}" if o == "hybrid" else o


def compute(ctx, cfg, which="gt", S_mats=None, need_null_rho=False):
    """cfg keys: tau rho_f rho_g rho_s S alphas(dict name@k -> a) winsor scale. Returns dict(post, rep, n, me, K)."""
    names = ["cur"] + [n for n, a in cfg["alphas"].items() if a]
    bl = ctx.blocks(which, names, cfg.get("scale", "wsd"), cfg.get("winsor"))
    G, h, g0, n = B2.combine(bl, cfg["alphas"])
    S = S_mats.get(cfg.get("S")) if (S_mats and cfg.get("S")) else None
    K = B2.kernel2(ctx.Udf, ctx.kin, cfg["rho_f"], cfg["rho_g"], cfg.get("rho_s", 0.0), S)
    post = B2.posterior2(G, h, g0, K, cfg["tau"])
    rep = report_mask(n, K, 1, True)
    return dict(post=post, rep=rep, n=n, me=B2.m_eff(post["q"], post["suma"]), K=K)


def outputs_from(r, rho_e=0.0, which=OUTPUTS):
    out = {}
    for o, l in which:
        if o == "empnull":
            zf = B2.output_stat(r["post"], "z_freq")
            z, _, _ = B2.empirical_null(zf, r["me"], r["rep"])
        else:
            z = B2.output_stat(r["post"], o, l, rho_e)
        out[out_name(o, l)] = np.where(r["rep"], z, np.nan)
    return out


def fit_rho(ctx, cfg, S_mats):
    rn = compute(ctx, cfg, "null", S_mats)
    rho, _ = B2.fit_rho_null(B2.output_stat(rn["post"], "z_freq"), rn["me"], rn["rep"])
    return rho


def score_cfg(ctx, cfg, output, S_mats=None, which="gt"):
    """Score matrix (numpy, conds x kinases) for one (cfg, output name); z_vif refits rho_e on DEV nulls."""
    r = compute(ctx, cfg, which, S_mats)
    rho = fit_rho(ctx, cfg, S_mats) if output == "z_vif" else 0.0
    sel = [x for x in OUTPUTS if out_name(*x) == output]
    return outputs_from(r, rho, sel)[output], r


def similarity_mats(ctx):
    from kinact.library import load_library
    S = {n: B2.topk_similarity(load_library(n, k=500), ctx.kin, 500, power=2) for n in ["r2np", "cantley"]}
    for n, M in S.items():
        off = M[np.triu_indices_from(M, 1)]
        print(f"[S {n}] offdiag mean {off.mean():.4f}; q50/90/99/max {np.quantile(off, [.5, .9, .99, 1]).round(3)}", flush=True)
    return S


def run(ctx, exp, cfgs, outputs=OUTPUTS, S_mats=None, ledger=True):
    rows = []
    t0 = time.time()
    last_blocks = None
    for i, cfg in enumerate(cfgs):
        bkey = (tuple(sorted(k for k, a in cfg["alphas"].items() if a)), cfg.get("winsor"))
        if bkey != last_blocks:
            ctx.drop_blocks(); last_blocks = bkey
        r = compute(ctx, cfg, "gt", S_mats)
        rho = fit_rho(ctx, cfg, S_mats) if any(o == "z_vif" for o, _ in outputs) else 0.0
        for oname, arr in outputs_from(r, rho, outputs).items():
            S = pd.DataFrame(arr, index=ctx.conds, columns=ctx.kin)
            m = inner_metrics(ctx.ev, ctx.ev.pairs(S), ctx.fold_of)
            full = dict(model="bayes2", exp=exp, lib_cur="curated_union", scale=cfg.get("scale", "wsd"),
                        winsor=cfg.get("winsor"), tau=cfg["tau"], rho_f=cfg["rho_f"], rho_g=cfg["rho_g"],
                        rho_s=cfg.get("rho_s", 0.0), S=cfg.get("S"), alphas=cfg["alphas"], output=oname,
                        borrow=True, min_sub=1)
            from kinact.ledger import config_hash
            key = config_hash(full)
            rows.append({**full, "alphas": repr(cfg["alphas"]), "key": key, "rho_e": rho, **m})
            if ledger:
                ledger_row("bayes2", full, m)
        if i % 20 == 0:
            best = max(rows, key=lambda x: x["M1-sel"])
            print(f"[{exp}] {i+1}/{len(cfgs)} {time.time()-t0:.0f}s best {best['M1-sel']:.4f} {best['output']} "
                  f"tau={best['tau']} rf={best['rho_f']} rg={best['rho_g']} {best['alphas']} w={best['winsor']}",
                  flush=True)
    g = pd.DataFrame(rows)
    g.to_csv(OUT / f"grid_{exp}.tsv", sep="\t", index=False)
    return g


def cfgs_out():
    for w, a, tau, rf, rg in itertools.product([None, 5.0, 3.0, 2.0], [1.0, 2.0], [0.06, 0.125, 0.25, 0.5, 1.0],
                                               [0.0, 1.0, 3.0], [0.0, 0.25, 1.0]):
        yield dict(winsor=w, alphas={"r2np@100": a}, tau=tau, rho_f=rf, rho_g=rg)


def cfgs_src(winsor):
    for k2, k3, kc in itertools.product([100, 200], [100, 200], [100, 200]):
        for a3, ac in itertools.product([0.0, 0.25, 0.5, 1.0], [0.0, 0.25, 0.5, 1.0]):
            if (a3 == 0 and k3 == 200) or (ac == 0 and kc == 200):
                continue
            for tau, rf, rg in itertools.product([0.06, 0.125, 0.25], [1.0, 3.0], [0.25, 1.0]):
                yield dict(winsor=winsor, alphas={f"r2np@{k2}": 1.0, f"r3@{k3}": a3, f"cantley@{kc}": ac},
                           tau=tau, rho_f=rf, rho_g=rg)


def cfgs_kern(winsor):
    for Sn, rs, tau, rf, rg in itertools.product(["r2np", "cantley"], [0.0, 0.5, 1.0, 3.0], [0.06, 0.125, 0.25],
                                                 [0.0, 1.0, 3.0], [0.0, 0.25, 1.0]):
        if rs == 0 and Sn == "cantley":
            continue
        yield dict(winsor=winsor, alphas={"r2np@100": 1.0}, tau=tau, rho_f=rf, rho_g=rg, rho_s=rs,
                   S=Sn if rs else None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exp", required=True, choices=["out", "src", "kern"])
    ap.add_argument("--winsor", type=float, default=None)
    ap.add_argument("--outputs", default=None, help="comma list of output names to restrict to")
    args = ap.parse_args()
    outs = OUTPUTS if not args.outputs else [x for x in OUTPUTS if out_name(*x) in args.outputs.split(",")]
    if args.exp == "out":
        ctx = Ctx(ks={"r2np": [100]})
        run(ctx, "out", list(cfgs_out()), outs)
    elif args.exp == "src":
        ctx = Ctx(ks={"r2np": [100, 200], "r3": [100, 200], "cantley": [100, 200]})
        cf = sorted(cfgs_src(args.winsor), key=lambda c: tuple(sorted(k for k, a in c["alphas"].items() if a)))
        run(ctx, "src", cf, outs)
    else:
        ctx = Ctx(ks={"r2np": [100]})
        S = similarity_mats(ctx)
        run(ctx, "kern", list(cfgs_kern(args.winsor)), outs, S_mats=S)


if __name__ == "__main__":
    main()
