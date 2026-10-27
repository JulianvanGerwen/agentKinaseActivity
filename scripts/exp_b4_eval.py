#!/usr/bin/env python
"""Branch B4 evaluation: nested CV over the B4 grid (scripts/exp_b4_grid.py), comparison with the B2 fixed core
(same pipeline, no expansion), B2 nested OOF (re-assembled from cache/b2 grid scores), B1 best (wz r2np@200 global)
and prior-only; paired group bootstrap (B=2000, seed 0), M2, harmony bins, T-unseen, runtime.
Writes results/b4/b4_eval.json and results/b4/b4_tables.md."""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))
import exp_b2_eval as E2E  # noqa: E402
import exp_b4_grid as G4  # noqa: E402
from kinact import bootstrap as bs, controls as ctl, ledger  # noqa: E402
from kinact.explore import bayes as B  # noqa: E402
from kinact.explore import coreg as CR  # noqa: E402

FOLDS = G4.FOLDS
OUT, CACHE = G4.OUT, G4.CACHE


def nested_b4(g, zs, conds, fold_of, kin, mask):
    sub = g[mask]
    cf = pd.Series(conds).map(fold_of).to_numpy()
    parts, chosen = [], {}
    for f in FOLDS:
        best = sub.loc[sub[f"inner_{f}"].idxmax()]
        chosen[f] = {k: best[k] for k in ["kind", "seed", "q", "m", "sign", "beta", "rmin", "min_nbr", "gamma", "output"]
                     if k in best and pd.notna(best[k])} | {"inner": float(best[f"inner_{f}"])}
        arr = zs[best["config_hash"]]
        parts.append(pd.DataFrame(arr[cf == f], index=np.asarray(conds)[cf == f], columns=kin))
    return pd.concat(parts).loc[conds], chosen


def full_metrics(ctx, S, name):
    r = E2E.summarise(ctx.ev.evaluate(S))
    r["T-unseen"] = E2E.t_unseen(ctx.ev, S, ctx.P, ctx.C, ctx.fold_of)
    return r


def runtime(ctx, cfg_best):
    """CPU single-sample path (SingleSample.score_one, 1 thread) for the base design vs the chosen expansion
    (fold0 expansion), 200 random DEV GT conditions, median of 3."""
    from threadpoolctl import threadpool_limits
    rng = np.random.default_rng(0)
    cs = rng.choice(ctx.conds, min(200, len(ctx.conds)), replace=False)
    Xs = ctx.X[list(cs)].to_numpy(np.float64)
    out = {}
    exp = None
    if cfg_best and cfg_best.get("kind") == "expand":
        exp = CR.expansion_edges(ctx.top(cfg_best["seed"], int(cfg_best["q"]), "fold0"), ctx.kin, ctx.sites,
                                 int(cfg_best["m"]), float(cfg_best["rmin"]), cfg_best["sign"])
    for name, e, beta in [("base", None, 0.0), ("b4", exp, float(cfg_best.get("beta", 0) or 0))]:
        d = CR.make_design(ctx.sites, ctx.kin, ctx.base, e)
        cfg = B.BayesConfig(tau=G4.CORE["tau"], rho_f=G4.CORE["rho_f"], rho_g=G4.CORE["rho_g"], alpha=beta if beta else 1.0,
                            borrow=True, noise_model=False)
        ss = B.SingleSample(d, ctx.Udf, cfg, None)
        if not beta:
            ss.W = d.Wc.tocsr()
        ts = []
        with threadpool_limits(1):
            for rep in range(3):
                t0 = time.perf_counter()
                for j in range(Xs.shape[1]):
                    x = Xs[:, j]
                    med = np.nanmedian(x); sc = 1.4826 * np.nanmedian(np.abs(x - med))
                    ss.score_one(x[d.lib_rows], med, sc)
                ts.append((time.perf_counter() - t0) / Xs.shape[1] * 1e3)
        out[name] = {"ms_per_cond_1thread": float(np.median(ts)), "lib_rows": int(len(d.lib_rows)),
                     "edges": int(ss.W.nnz)}
    pre = json.load(open(CACHE / "allowed.json"))["timing"]
    out["precompute_s_per_fold_all_variants"] = {f: round(v["total_s"], 1) for f, v in pre.items()}
    return out


def main():
    t0 = time.time()
    ctx = G4.Ctx()
    conds, kin, fold_of = ctx.conds, ctx.kin, ctx.fold_of
    g = pd.read_csv(OUT / "grid.tsv", sep="\t")
    zs = np.load(CACHE / "grid_scores.npz", allow_pickle=True)
    assert list(zs["conds"]) == conds
    S, chosen = {}, {}
    isnone = g["kind"] == "none"
    for o in ["z", "mean", "wz"]:
        ref = g[isnone & (g["output"] == o)].iloc[0]
        S[f"ref_{o}"] = pd.DataFrame(zs[ref["config_hash"]], index=conds, columns=kin)
        S[f"B4exp_{o}"], chosen[f"B4exp_{o}"] = nested_b4(g, zs, conds, fold_of, kin, (g["output"] == o) & (isnone | (g["kind"] == "expand")))
        S[f"B4imp_{o}"], chosen[f"B4imp_{o}"] = nested_b4(g, zs, conds, fold_of, kin, (g["output"] == o) & (isnone | (g["kind"] == "impute")))
    S["B4all_zmean"], chosen["B4all_zmean"] = nested_b4(g, zs, conds, fold_of, kin, g["output"].isin(["z", "mean"]))
    # forced-on expansion (no 'none' option) to see the cost/benefit even when selection would decline it
    S["B4exp_forced_z"], chosen["B4exp_forced_z"] = nested_b4(g, zs, conds, fold_of, kin, (g["output"] == "z") & (g["kind"] == "expand"))
    S["B4imp_forced_z"], chosen["B4imp_forced_z"] = nested_b4(g, zs, conds, fold_of, kin, (g["output"] == "z") & (g["kind"] == "impute"))
    # B2 nested reference (re-assembled from the B2 grid archives)
    g2, store2 = E2E.load_grid()
    assert list(store2["conds"]) == conds
    S["B2_nested"], _ = E2E.nested(g2, store2, conds, fold_of, g2["output"] == "z")
    S["B2_nested"] = S["B2_nested"].reindex(columns=kin)
    # B1 best (apparent; nested B1 = 0.638)
    from kinact.explore import lib_variants as LV
    dd = LV.DevData()
    eng = LV.Engine(dd)
    S["B1_best"] = eng.score({"stat": "wz", "lib": {"type": "pred", "name": "r2np", "k": 200, "mode": "global", "gamma": 0},
                              "min_sub": 1}).loc[conds].reindex(columns=kin)
    po = []
    Cs = ctx.C.set_index("conditionID")["split"]
    for f in FOLDS:
        ec = [c for c in conds if fold_of[c] == f]
        tc = [c for c in ctx.C["conditionID"] if Cs[c] in FOLDS and Cs[c] != f]
        po.append(ctl.prior_only(ctx.P, ctx.C, tc, ec, kin))
    S["prior_only"] = pd.concat(po).loc[conds]
    ctx.ev.zero_policy = "ignore"
    res = {k: full_metrics(ctx, v, k) for k, v in S.items()}
    for k in S:
        if k.startswith("B4"):
            ledger.append("b4_nested", {"variant": k, "chosen": chosen[k]}, "dev",
                          {x: res[k][x] for x in ["M1-sel", "M1-all", "M1-sel-scored", "coverage"]}, branch="B4", folds=range(5))
    plan = bs.BootstrapPlan.from_evaluator(ctx.ev, B=2000, seed=0)
    boot = {}
    pairs = [("B4exp_z", "ref_z"), ("B4exp_mean", "ref_mean"), ("B4exp_wz", "ref_wz"), ("B4imp_z", "ref_z"),
             ("B4imp_mean", "ref_mean"), ("B4exp_forced_z", "ref_z"), ("B4imp_forced_z", "ref_z"),
             ("B4all_zmean", "ref_z"), ("B4all_zmean", "ref_mean"),
             ("B4exp_z", "B2_nested"), ("B4exp_mean", "B2_nested"), ("B4all_zmean", "B2_nested"), ("ref_z", "B2_nested"),
             ("B4exp_z", "B1_best"), ("B4exp_wz", "B1_best"), ("B4all_zmean", "prior_only")]
    for a, b in pairs:
        d = bs.paired_delta(ctx.ev, S[a], S[b], "M1-sel", plan)
        boot[f"{a}_vs_{b}"] = {k: float(d[k]) for k in ["estimate", "se", "ci_low", "ci_high", "p_one_sided"]}
        boot[f"{a}_vs_{b}"]["M2"] = ctx.ev.m2(S[a], S[b])
    # forced expansion: per-grid summary of expansion effect (pooled, non-nested) vs ref
    rt = runtime(ctx, chosen["B4exp_forced_z"]["fold0"])
    out = {"gt_freeze_sha256": G4.E2.gt_freeze_hash(), "results": res, "chosen": chosen, "bootstrap": boot,
           "runtime": rt, "n_grid_evals": int(len(g))}
    json.dump(out, open(OUT / "b4_eval.json", "w"), indent=1, default=lambda x: x.item() if hasattr(x, "item") else str(x))
    # tables
    L = ["| variant | M1-sel | M1-all | cov | harm | b0 | b1-4 | b5-19 | b>=20 | H | D | O | s- | s+ | T-unseen | nsig |",
         "|" + "---|" * 16]
    for k, r in res.items():
        hb = r["harmony_per_bin"]
        gb = lambda b: hb.get(b, {}).get("value", np.nan) if isinstance(hb.get(b), dict) else hb.get(b, np.nan)
        L.append(f"| {k} | {r['M1-sel']:.3f} | {r['M1-all']:.3f} | {r['coverage']:.3f} | {r['harmony']:.3f} | "
                 f"{gb('0'):.3f} | {gb('1-4'):.3f} | {gb('5-19'):.3f} | {gb('>=20'):.3f} | {r['by_source'].get('H', np.nan):.3f} | "
                 f"{r['by_source'].get('D', np.nan):.3f} | {r['by_source'].get('O', np.nan):.3f} | {r['by_sign'].get('-1', np.nan):.3f} | "
                 f"{r['by_sign'].get('1', np.nan):.3f} | {r['T-unseen']['value']:.3f} | {r['n_sig']} |")
    L += ["", "| comparison | ΔM1-sel | 95% CI | p | M2 A vs B (n pairs) |", "|---|---|---|---|---|"]
    for k, d in boot.items():
        m2 = d["M2"]
        m2s = f"{m2.get('A', np.nan):.3f} vs {m2.get('B', np.nan):.3f} ({m2.get('n_pairs', '')})" if isinstance(m2, dict) else str(m2)
        L.append(f"| {k} | {d['estimate']:+.3f} | [{d['ci_low']:+.3f}, {d['ci_high']:+.3f}] | {d['p_one_sided']:.4f} | {m2s} |")
    L += ["", "Chosen per fold:"] + [f"- {k}: " + "; ".join(f"{f}: {json.dumps(c, default=str)}" for f, c in v.items())
                                     for k, v in chosen.items()]
    L += ["", f"Runtime: {json.dumps(rt)}"]
    (OUT / "b4_tables.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))
    print(f"done {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
