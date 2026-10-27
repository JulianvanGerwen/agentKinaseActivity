#!/usr/bin/env python
"""Branch B2.1 (+B3): nested-CV variants from results/b21/grid_*.tsv, full OOF evaluation, paired group bootstrap
vs the B2 reference (B2 nested OOF scores, cache/b2), M2 vs PSP z-score, T-unseen within CV, calibration of each
variant's z (fold-wise chosen config applied to DEV role=null conditions and to the site-permuted GT conditions of
the same fold), runtime of the recommended config. Writes results/b21/b21_eval.json and results/b21/b21_tables.md.
"""
import ast
import json
import time

import numpy as np
import pandas as pd

import exp_b21_grid as G
from exp_b21_common import (B2, FOLDS, OUT, Ctx, b2_reference, bs, calib_bins, dump, full_eval, ledger_row,
                            nested)
from kinact.methods.linear import zscore

Z_OUTPUTS = {"z_post", "z_freq", "z_vif", "empnull"}


def load_grids():
    gs = []
    for e in ["out", "src", "kern"]:
        p = OUT / f"grid_{e}.tsv"
        if p.exists():
            g = pd.read_csv(p, sep="\t", keep_default_na=False, na_values=[""])
            gs.append(g)
    g = pd.concat(gs, ignore_index=True)
    g["winsor"] = pd.to_numeric(g["winsor"], errors="coerce")
    g["S"] = g["S"].where(g["S"].notna() & (g["S"] != ""), None)
    return g


def row_cfg(r):
    return dict(winsor=None if pd.isna(r["winsor"]) else float(r["winsor"]), alphas=ast.literal_eval(r["alphas"]),
                tau=float(r["tau"]), rho_f=float(r["rho_f"]), rho_g=float(r["rho_g"]),
                rho_s=float(r["rho_s"]), S=r["S"], scale=r["scale"])


class Scorer:
    def __init__(self, ctx, grid, S_mats):
        self.ctx, self.grid, self.S = ctx, grid.set_index("key"), S_mats
        self.cache = {}

    def __getitem__(self, key):
        if key not in self.cache:
            r = self.grid.loc[key]
            arr, _ = G.score_cfg(self.ctx, row_cfg(r), r["output"], self.S)
            self.cache[key] = arr
        return self.cache[key]


def calib_variant(ctx, grid, chosen, S_mats):
    """Fold-wise: config chosen for fold f scores the DEV nulls and site-permuted GT conditions of fold f.
    z = the variant's output if it is a z statistic, else the empirical-null z of the same config."""
    gi = grid.set_index("key")
    res = {}
    for which, conds in [("null", ctx.nulls), ("perm", ctx.conds)]:
        fo = pd.Series(conds).map(ctx.fold_of).to_numpy()
        Zs, Ns = [], []
        for f in FOLDS:
            idx = np.flatnonzero(fo == f)
            if not len(idx):
                continue
            r = gi.loc[chosen[f]["key"]]
            out = r["output"] if r["output"] in Z_OUTPUTS else "empnull"
            z, rr = G.score_cfg(ctx, row_cfg(r), out, S_mats, which=which)
            Zs.append(z[idx]); Ns.append(rr["n"][idx])
        res[which] = calib_bins(np.vstack(Zs), np.vstack(Ns))
    return res


def main():
    t0 = time.time()
    ctx = Ctx(ks={"r2np": [100, 200], "r3": [100, 200], "cantley": [100, 200]})
    grid = load_grids()
    S_mats = G.similarity_mats(ctx) if (grid["rho_s"] > 0).any() else None
    sc = Scorer(ctx, grid, S_mats)
    o = grid["exp"] == "out"
    w0 = grid["winsor"].isna()
    a1 = grid["alphas"] == "{'r2np@100': 1.0}"
    hyb = grid["output"].str.startswith("hybrid")
    variants = {
        "B2 config, new code (z_post, no winsor)": o & w0 & (grid["output"] == "z_post") & a1,
        "out: z_post": o & w0 & (grid["output"] == "z_post"),
        "out: mean": o & w0 & (grid["output"] == "mean"),
        "out: z_freq": o & w0 & (grid["output"] == "z_freq"),
        "out: hybrid (l free)": o & w0 & hyb,
        "out: z_vif (rho_e on DEV nulls)": o & w0 & (grid["output"] == "z_vif"),
        "out: empnull": o & w0 & (grid["output"] == "empnull"),
        "out: output free (pre-specified)": o & w0,
        "out+winsor: z_freq": o & (grid["output"] == "z_freq"),
        "out+winsor: mean": o & (grid["output"] == "mean"),
        "out+winsor: empnull": o & (grid["output"] == "empnull"),
        "out+winsor: output free": o,
    }
    for e in ["src", "kern"]:
        if (grid["exp"] == e).any():
            ge = grid["exp"] == e
            for out in sorted(grid.loc[ge, "output"].unique()):
                variants[f"{e}: {out}"] = ge & (grid["output"] == out)
            variants[f"{e}: output free"] = ge
            if e == "src":
                r2only = ge & grid["alphas"].map(lambda s: all(v == 0 for k, v in ast.literal_eval(s).items()
                                                               if not k.startswith("r2np")))
                variants["src: r2np-only subset, output free"] = r2only
    variants["ALL grids: output free (final)"] = pd.Series(True, index=grid.index)
    variants["ALL grids: empnull only (final-z)"] = grid["output"] == "empnull"

    S, chosen = {}, {}
    for name, m in variants.items():
        if m.sum() == 0:
            continue
        S[name], chosen[name] = nested(grid, sc, ctx, m, hp=["output", "winsor", "tau", "rho_f", "rho_g", "rho_s",
                                                                "S", "alphas"])
    S_ref, ch_ref = b2_reference(ctx)
    psp = zscore(ctx.X[ctx.conds], "psp_all", min_substrates=3).reindex(columns=ctx.kin)
    plan = bs.BootstrapPlan.from_evaluator(ctx.ev, B=2000, seed=0)
    res = {"B2 reference (nested OOF, cache/b2)": full_eval(ctx, S_ref)}
    res["B2 reference (nested OOF, cache/b2)"]["M2_vs_pspz"] = ctx.ev.m2(S_ref, psp)
    for name, s in S.items():
        r = full_eval(ctx, s)
        r["M2_vs_pspz"] = ctx.ev.m2(s, psp)
        d = bs.paired_delta(ctx.ev, s, S_ref, "M1-sel", plan)
        r["dM1sel_vs_B2"] = {k: d[k] for k in ["estimate", "ci_low", "ci_high", "p_one_sided"]}
        for ep in ["M1-all", "harmony"]:
            r[f"d{ep}_vs_B2"] = bs.paired_delta(ctx.ev, s, S_ref, ep, plan)["estimate"]
        r["chosen"] = chosen[name]
        r["calibration"] = calib_variant(ctx, grid, chosen[name], S_mats)
        res[name] = r
        ledger_row("bayes2_nested", {"variant": name, "chosen": {f: c["key"] for f, c in chosen[name].items()}}, r)
        print(f"{name}: M1-sel {r['M1-sel']:.4f} d {d['estimate']:+.4f} [{d['ci_low']:+.3f},{d['ci_high']:+.3f}] "
              f"p={d['p_one_sided']:.3f} calib {r['calibration']}", flush=True)
    # pooled-best (non-nested) per variant family, for reference
    best = grid.loc[grid["M1-sel"].idxmax()]
    rt = runtime(ctx, grid, chosen.get("ALL grids: output free (final)"), S_mats)
    out = {"gt_freeze_sha256": ctx.gt_hash, "n_grid_evaluations": int(len(grid)),
           "n_by_exp": grid["exp"].value_counts().to_dict(), "results": res,
           "pooled_best": best.to_dict(), "runtime": rt, "wall_s": time.time() - t0, "b2_ref_chosen": ch_ref}
    dump("b21_eval.json", out)
    tables(out)


def runtime(ctx, grid, chosen, S_mats, n_cond=200):
    from threadpoolctl import threadpool_limits
    if chosen is None:
        return {}
    # use the most frequently chosen config across folds
    keys = pd.Series([c["key"] for c in chosen.values()])
    r = grid.set_index("key").loc[keys.mode()[0]]
    c = row_cfg(r)
    cfg = B2.Bayes2Config(tau=c["tau"], rho_f=c["rho_f"], rho_g=c["rho_g"], rho_s=c["rho_s"], alphas=c["alphas"],
                          output=r["output"] if not r["output"].startswith("hybrid") else "hybrid",
                          lam=float(r["output"][6:]) if r["output"].startswith("hybrid") else 0.5,
                          scale=c["scale"], winsor=c["winsor"])
    ss = B2.SingleSample2(ctx.md, ctx.Udf, cfg, S_mats.get(c["S"]) if (S_mats and c["S"]) else None)
    rng = np.random.default_rng(0)
    cols = list(rng.choice(ctx.X.columns, n_cond, replace=False))
    V = ctx.X[cols].to_numpy(np.float64)
    out = {"cfg": {**c, "output": r["output"]}}
    for th in (1, 16):
        with threadpool_limits(th):
            ss.score_one(V[:, 0])
            ts = []
            for rep in range(5):
                t = time.time()
                for j in range(n_cond):
                    ss.score_one(V[:, j])
                ts.append(1e3 * (time.time() - t) / n_cond)
        out[f"cpu_{th}t_ms_per_cond"] = float(np.median(ts))
    return out


def tables(out):
    f3 = lambda x: "–" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{x:.3f}"
    L = ["| variant | M1-sel | M1-all | M2 vs PSP-z | harm 0 / 1-4 / 5-19 / ≥20 | T-uns | H / D / O | ΔM1-sel vs B2 [95% CI] p | "
         "calib null all (0/1-4/5-19/≥20) | calib perm all (0/1-4/5-19/≥20) |", "|" + "---|" * 10]
    for name, r in out["results"].items():
        hb = r["harmony_per_bin"]
        harm = " / ".join(f3(hb.get(b)) for b in ["0", "1-4", "5-19", ">=20"])
        src = " / ".join(f3(r["by_source"].get(s)) for s in ["H", "D", "O"])
        m2 = r["M2_vs_pspz"]
        m2s = f"{m2['A']:.3f} vs {m2['B']:.3f}" 
        d = r.get("dM1sel_vs_B2")
        ds = f"{d['estimate']:+.3f} [{d['ci_low']:+.3f}, {d['ci_high']:+.3f}] {d['p_one_sided']:.3f}" if d else "ref"
        cal = r.get("calibration", {})
        cs = lambda c: (f"{c['all']:.3f} ({'/'.join(f3(c.get(b)) for b in ['0', '1-4', '5-19', '>=20'])})" if c else "–")
        L.append(f"| {name} | {r['M1-sel']:.3f} | {r['M1-all']:.3f} | {m2s} | {harm} | {f3(r['T-unseen'])} | {src} | {ds} | "
                 f"{cs(cal.get('null'))} | {cs(cal.get('perm'))} |")
    L.append("")
    L.append("Chosen per fold:")
    for name, r in out["results"].items():
        if "chosen" in r:
            L.append(f"- {name}: " + "; ".join(
                f"{f}: {c['output']} w={c['winsor']} τ={c['tau']} ρf={c['rho_f']} ρg={c['rho_g']} ρs={c['rho_s']} "
                f"S={c['S']} {c['alphas']} (inner {c['inner']:.3f})" for f, c in r["chosen"].items()))
    (OUT / "b21_tables.md").write_text("\n".join(L))
    print("\n".join(L))


if __name__ == "__main__":
    main()
