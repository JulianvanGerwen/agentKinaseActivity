#!/usr/bin/env python
"""Branch 2: nested-CV selection on the grid (scripts/exp_b2_grid.py), full OOF evaluation, ablations,
calibration (DEV null conditions + site-permuted data), v-leakage check, runtime, same-family diagnostic.
Writes results/b2/b2_eval.json and results/b2/b2_tables.md (pasted into b2_report.md)."""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))
import exp_b2_grid as E  # noqa: E402
from kinact import bootstrap as bs, controls as ctl, ledger  # noqa: E402
from kinact.eval import calibration, psp_bin, to_jsonable  # noqa: E402
from kinact.explore import bayes as B  # noqa: E402
from kinact.library import load_library  # noqa: E402
from kinact.methods.linear import zscore  # noqa: E402

FOLDS = E.FOLDS
OUT, CACHE = E.OUT, E.CACHE
HP = ["tau", "rho_f", "rho_g", "alpha", "intercept", "borrow", "output", "noise"]


class _Stores:
    """Lazy union of several grid score archives."""
    def __init__(self, zs):
        self.zs = zs

    def __getitem__(self, k):
        if k in ("conds", "kinases"):
            return self.zs[0][k]
        for z in self.zs:
            if k in z.files:
                return z[k]
        raise KeyError(k)


def load_grid(tags=("", "_lowtau", "_rhof10")):
    gs, zs = [], []
    for t in tags:
        gs.append(pd.read_csv(OUT / f"grid_noise-both{t}.tsv", sep="\t"))
        zs.append(np.load(CACHE / f"grid_scores_noise-both{t}.npz", allow_pickle=True))
    for z in zs[1:]:
        assert list(z["conds"]) == list(zs[0]["conds"]) and list(z["kinases"]) == list(zs[0]["kinases"])
    g = pd.concat(gs, ignore_index=True).drop_duplicates("config_hash")
    return g, _Stores(zs)


def nested(g, store, conds, fold_of, mask=None):
    """Per outer fold f: argmax over (masked) configs of inner M1-sel on the other folds; OOF assembly."""
    sub = g if mask is None else g[mask]
    parts, chosen = [], {}
    kin = store["kinases"]
    cf = pd.Series(conds).map(fold_of).to_numpy()
    for f in FOLDS:
        best = sub.loc[sub[f"inner_{f}"].idxmax()]
        chosen[f] = {k: best[k] for k in HP} | {"config_hash": best["config_hash"], "inner": best[f"inner_{f}"]}
        arr = store[best["config_hash"]]
        parts.append(pd.DataFrame(arr[cf == f], index=np.asarray(conds)[cf == f], columns=kin))
    return pd.concat(parts).loc[conds], chosen


def summarise(r):
    v = lambda x: x["value"] if isinstance(x, dict) and "value" in x else x
    return {"M1-sel": v(r["M1-sel"]), "M1-all": v(r["M1-all"]), "M1-sel-scored": v(r["M1-sel-scored"]),
            "coverage": r["coverage_M1-sel"], "harmony": v(r["harmony"]), "harmony_per_bin": r["harmony"]["per_bin"],
            "by_source": {k: v(x) for k, x in r["by_source"].items()},
            "by_sign": {str(k): v(x) for k, x in r["by_sign"].items()},
            "n_sig": r["coverage_significant"]["n_significant"]}


def same_family_auroc(ev, S, Udf):
    """Diagnostic (r1_algo item 3): scored-only AUROC of each target vs evaluable negatives of its own Manning
    family; pair -> kinase mean -> macro (M1-sel pairs). Pairs without a scored same-family negative skipped."""
    fam = Udf.set_index("kinase")["manning_family"].reindex(ev.kinases).to_numpy()
    Sm = ev._matrix(S)
    vals = {}
    t = ev.targets[ev.targets["in_M1sel"]]
    for ci, ki, sg, k in zip(t["ci"], t["ki"], t["sign"], t["kinase"]):
        negs = ev.neg_mask[ci] & (fam == fam[ki])
        sk = Sm[ci, ki]
        v = Sm[ci, negs]; v = v[~np.isnan(v)]
        if np.isnan(sk) or len(v) == 0:
            continue
        a = (np.sum(sg * v < sg * sk) + 0.5 * np.sum(v == sk)) / len(v)
        vals.setdefault(k, []).append(a)
    km = [np.mean(x) for x in vals.values()]
    return {"value": float(np.mean(km)) if km else np.nan, "n_kinases": len(km)}


def t_unseen(ev, S, P, C, fold_of):
    """Pooled OOF T-unseen: M1-sel over pairs whose kinase is never a direct high/medium positive in the other
    DEV folds' GT conditions."""
    Pp = ev.pairs(S)
    gc = C[C["split"].isin(FOLDS) & (C["role"] == "gt")].set_index("conditionID")["split"]
    pos = P[(P["cls"] == "direct") & P["confidence"].isin(["high", "medium"]) & P["conditionID"].isin(gc.index)]
    pos = pos.assign(fold=pos["conditionID"].map(gc))
    seen = {f: set(pos.loc[pos["fold"] != f, "kinase"]) for f in FOLDS}
    f = Pp["conditionID"].map(fold_of)
    uns = np.array([k not in seen[ff] for k, ff in zip(Pp["kinase"], f)])
    r = ev.metric(Pp[uns], "M1-sel")
    return {"value": r["value"], "n_kinases": r["n_kinases"], "n_pairs": r["n_pairs"]}


def calib_bins(Z, N):
    z = Z.to_numpy(float); n = N.reindex(index=Z.index, columns=Z.columns).to_numpy(float)
    ok = ~np.isnan(z); lab = psp_bin(n.ravel()).reshape(n.shape)
    out = {"all": float((np.abs(z[ok]) > 1.96).mean()), "n": int(ok.sum())}
    for b in ["0", "1-4", "5-19", ">=20"]:
        m = ok & (lab == b)
        if m.any():
            out[b] = float((np.abs(z[m]) > 1.96).mean())
    out["sd_z"] = float(np.nanstd(z))
    return out


def score_cfg(X, conds, meta, vs, design, Udf, cfg, noise=None, v_override=None, return_n=False):
    noise = bool(cfg["noise"]) if noise is None else noise
    Y, D = E.standardise_by_fold(X, conds, meta, vs, noise=noise, v_override=v_override)
    G, h, n = B.suff_stats(design, cfg["alpha"], Y, D)
    K = B.kinase_kernel(Udf, design.kinases, cfg["rho_f"], cfg["rho_g"])
    mu, sd = B.posterior(G, h, K, cfg["tau"], cfg["intercept"])
    n = n.cpu().numpy()
    rep = B.report_mask(n, K, 1, bool(cfg["borrow"]))
    arr = np.where(rep, mu / sd if cfg["output"] == "z" else mu, np.nan)
    S = pd.DataFrame(arr, index=conds, columns=design.kinases)
    return (S, pd.DataFrame(n, index=conds, columns=design.kinases)) if return_n else S


def main():
    t0 = time.time()
    P, C, Udf, ev, meta, X = E.setup()
    vs = E.noise_variances(X, meta)
    g, store = load_grid()
    conds = list(store["conds"])
    assert conds == list(ev.cond_index)
    fold_of = meta["split"]
    design = B.build_design(X.index, list(Udf["kinase"]), load_library(E.LIB_CUR),
                            load_library(E.LIB_PRED[0], k=E.LIB_PRED[1]))
    kin = list(Udf["kinase"])
    # ---------------- references
    ref = {"zscore_psp_min3": zscore(X[conds], "psp_all", min_substrates=3).reindex(columns=kin),
           "zscore_cur_min1": zscore(X[conds], "curated_union", min_substrates=1).reindex(columns=kin)}
    po = []
    for f in FOLDS:
        ec = [c for c in conds if fold_of[c] == f]
        tc = [c for c in C["conditionID"] if C.set_index("conditionID")["split"][c] in FOLDS and
              C.set_index("conditionID")["split"][c] != f]
        po.append(ctl.prior_only(P, C, tc, ec, kin))
    ref["prior_only"] = pd.concat(po).loc[conds]
    # ---------------- nested variants
    core = g["output"] == "z"            # noise model on/off is selected by nested CV like a hyperparameter
    variants = {
        "B2_free": core,
        "B2_noborrow": core & ~g["borrow"],
        "B2_borrow": core & (g["borrow"] | ((g["rho_f"] == 0) & (g["rho_g"] == 0))),
        "abl_noK": core & (g["rho_f"] == 0) & (g["rho_g"] == 0),
        "abl_noise_on": core & g["noise"],
        "abl_noise_off": core & ~g["noise"],
        "abl_mean": g["output"] == "mean",
        "abl_alpha0": core & (g["alpha"] == 0),
        "abl_noK_alpha0": core & (g["alpha"] == 0) & (g["rho_f"] == 0) & (g["rho_g"] == 0),
        "abl_tau_ge1_noK_alpha0": core & (g["alpha"] == 0) & (g["rho_f"] == 0) & (g["rho_g"] == 0) & (g["tau"] >= 1),
    }
    S, chosen = {}, {}
    for name, m in variants.items():
        S[name], chosen[name] = nested(g, store, conds, fold_of, m)
    S.update(ref)
    res = {}
    for name, sc in S.items():
        ev.zero_policy = "ignore"
        res[name] = summarise(ev.evaluate(sc))
        res[name]["same_family"] = same_family_auroc(ev, sc, Udf)
        res[name]["T-unseen"] = t_unseen(ev, sc, P, C, fold_of)
        if name in ref:
            ledger.append(f"b2_ref_{name}", {"ref": name}, "dev", {k: res[name][k] for k in ["M1-sel", "M1-all", "M1-sel-scored", "coverage"]},
                          branch="B2", folds=range(5))
        else:
            ledger.append("bayes_core_nested", {"variant": name, "chosen": {f: c["config_hash"] for f, c in chosen[name].items()}},
                          "dev", {k: res[name][k] for k in ["M1-sel", "M1-all", "M1-sel-scored", "coverage"]},
                          branch="B2", folds=range(5))
    # M2 vs z-score, bootstrap
    plan = bs.BootstrapPlan.from_evaluator(ev, B=2000, seed=0)
    boot = {}
    for a in ["B2_free", "B2_noborrow", "B2_borrow"]:
        for b in ["zscore_psp_min3", "zscore_cur_min1", "prior_only", "abl_noK", "abl_noise_on", "abl_noise_off", "abl_mean", "abl_alpha0", "B2_noborrow"]:
            d = bs.paired_delta(ev, S[a], S[b], "M1-sel", plan)
            boot[f"{a}_vs_{b}"] = {k: d[k] for k in ["estimate", "se", "ci_low", "ci_high", "p_one_sided"]}
        for b in ["zscore_psp_min3", "zscore_cur_min1"]:
            boot[f"{a}_vs_{b}_M2"] = {"A": ev.m2(S[a], S[b])}
    # ---------------- pooled-best (optimistic, non-nested) config of the core grid
    gcore = g[core]
    best = gcore.loc[gcore["M1-sel"].idxmax()]
    best_cfg = {k: best[k] for k in HP}
    # ---------------- v leakage check: v from all DEV vs leave-own-fold-out, at the pooled-best config
    gn = g[core & g["noise"]]
    bn = gn.loc[gn["M1-sel"].idxmax()]
    S_vall = score_cfg(X, conds, meta, vs, design, Udf, {k: bn[k] for k in HP}, v_override=vs["all"])
    S_vlfo = pd.DataFrame(store[bn["config_hash"]], index=conds, columns=store["kinases"])
    vcheck = {"M1-sel_v_lfo": ev.metric(ev.pairs(S_vlfo), "M1-sel")["value"],
              "M1-sel_v_all": ev.metric(ev.pairs(S_vall), "M1-sel")["value"],
              "corr_logv_fold_vs_all": {f: float(np.corrcoef(np.log(vs[f]), np.log(vs["all"]))[0, 1]) for f in FOLDS},
              "spearman_scores": float(pd.Series(S_vlfo.to_numpy().ravel()).corr(pd.Series(S_vall.to_numpy().ravel()),
                                                                                   method="spearman"))}
    # ---------------- calibration: DEV null conditions + site-permuted eval conditions
    nulls = [c for c in C.loc[(C["role"] == "null") & C["split"].isin(FOLDS), "conditionID"] if c in X.columns]
    calib = {}
    cal_cfgs = {"pooled_best": best_cfg, **{f"nested_{k}_fold0": {kk: chosen[k]["fold0"][kk] for kk in HP}
                                            for k in ["B2_noborrow", "B2_free", "abl_noise_on"]}}
    Xp = ctl.site_permutation(X[conds], seed=0, mode="within_condition")
    for cname, cfg in cal_cfgs.items():
        Zn, Nn = score_cfg(X, nulls, meta, vs, design, Udf, {**cfg, "output": "z"}, noise=bool(cfg["noise"]), return_n=True)
        Zp, Np = score_cfg(Xp, conds, meta, vs, design, Udf, {**cfg, "output": "z"}, noise=bool(cfg["noise"]), return_n=True)
        calib[cname] = {"cfg": cfg, "null_conds": calib_bins(Zn, Nn), "site_perm": calib_bins(Zp, Np)}
    for zn, (lib, ms) in {"zscore_psp_min3": ("psp_all", 3), "zscore_cur_min1": ("curated_union", 1)}.items():
        Zn, Nn = zscore(X[nulls], lib, min_substrates=ms, return_n=True)
        Zp, Np = zscore(Xp, lib, min_substrates=ms, return_n=True)
        calib[zn] = {"null_conds": calib_bins(Zn, Nn), "site_perm": calib_bins(Zp, Np)}
    # ---------------- runtime
    rt = runtime(X, meta, vs, design, Udf, best_cfg)
    out = {"gt_freeze_sha256": E.gt_freeze_hash(), "n_grid_evaluations": int(len(g)),
           "n_core_configs": int(core.sum()), "results": res, "chosen": chosen, "bootstrap": boot,
           "pooled_best_cfg": best_cfg, "pooled_best_M1sel": float(best["M1-sel"]), "v_check": vcheck,
           "calibration": calib, "null_conditions": nulls, "runtime": rt, "wall_s": time.time() - t0}
    (OUT / "b2_eval.json").write_text(json.dumps(to_jsonable(out), indent=1, default=str))
    print(json.dumps(to_jsonable({k: out[k] for k in ["pooled_best_cfg", "pooled_best_M1sel", "v_check", "runtime"]}),
                     indent=1, default=str))
    for n_, r in res.items():
        print(n_, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items() if k not in ("harmony_per_bin",)})


def runtime(X, meta, vs, design, Udf, cfg, n_cond=200):
    from threadpoolctl import threadpool_limits
    import torch
    rng = np.random.default_rng(0)
    cols = list(rng.choice(X.columns, n_cond, replace=False))
    rt = {}
    # GPU batched (incl. standardisation, sufficient statistics, posterior)
    for rep in range(2):
        torch.cuda.synchronize(); t = time.time()
        score_cfg(X, cols, meta, vs, design, Udf, {**cfg}, noise=bool(cfg["noise"]))
        torch.cuda.synchronize(); rt["gpu_batched_ms_per_cond"] = 1e3 * (time.time() - t) / n_cond
    # CPU single-sample path, single thread and 16 threads
    v = vs["all"]
    ss = B.SingleSample(design, Udf, B.BayesConfig(
        tau=cfg["tau"], rho_f=cfg["rho_f"], rho_g=cfg["rho_g"], alpha=cfg["alpha"], intercept="median",
        borrow=bool(cfg["borrow"]), noise_model=bool(cfg["noise"])), v[design.lib_rows])
    V = X[cols].to_numpy(np.float64)
    w = 1 / np.sqrt(v)
    for th in (1, 16):
        with threadpool_limits(th):
            ts = []
            for rep in range(3):
                t = time.time()
                for j in range(n_cond):
                    x = V[:, j]; m = ~np.isnan(x)
                    med = np.median(x[m]); sig = 1.4826 * np.median(np.abs((x[m] - med) * w[m]))
                    ss.score_one(x[design.lib_rows], med, sig)
                ts.append(1e3 * (time.time() - t) / n_cond)
            rt[f"cpu_single_{th}t_ms_per_cond"] = float(np.median(ts))
    t = time.time(); B.site_noise_variance(X.to_numpy()[:, (meta.loc[X.columns, "qc_status"] == "ok").to_numpy()])
    rt["noise_model_fit_s_one_off"] = time.time() - t
    rt["note"] = "CPU path uses intercept='median' (joint intercept only in the GPU path)"
    return rt


if __name__ == "__main__":
    main()
