#!/usr/bin/env python
"""Smoke run on pooled DEV folds 0-4 (not devconfirm, not test): z-score baseline + controls.
Writes results/smoke/smoke_dev.{md,json}, results/smoke/mde.json, appends to results/dev_ledger.tsv."""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from kinact import bootstrap as bs, controls as ctl, ledger  # noqa: E402
from kinact.data import library_matrix, load_dev, load_dev_meta  # noqa: E402
from kinact.eval import Evaluator, load_gt, to_jsonable  # noqa: E402

FOLDS = [f"fold{i}" for i in range(5)]
OUT = ROOT / "results/smoke"
LIB, MIN_SUB = "psp_all", 3
G_TEST_ALL = 87  # GT_REPORT split table: groups in split 'test' (all groups, not only M1-sel groups)


def zscore(X: pd.DataFrame, lib: str, kinases, min_sub=MIN_SUB):
    """z_k = mean(x_s over measured substrates)*sqrt(n)/sd(x_c); NaN if n < min_sub."""
    W, kidx, _ = library_matrix(lib, X.index, kinases=kinases, binary=True)
    cols = np.asarray(W.sum(0)).ravel() > 0
    V = X.to_numpy(np.float64)
    sd = np.nanstd(V, axis=0, ddof=1)
    Vs = V[cols]; ok = ~np.isnan(Vs)
    Ws = W[:, np.flatnonzero(cols)]
    n = Ws @ ok.astype(np.float64)
    s = Ws @ np.where(ok, Vs, 0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        z = (s / n) * np.sqrt(n) / sd[None, :]
    z[n < min_sub] = np.nan
    return pd.DataFrame(z.T, index=X.columns, columns=list(kinases)), pd.DataFrame(n.T, index=X.columns, columns=list(kinases))


def main():
    t0 = time.time()
    P, C, U, bins_c = load_gt()
    from kinact.eval import read_gt_tsv
    bins = read_gt_tsv(ROOT / "results/ground_truth/harmony_bins_group.tsv")
    kin = list(U["kinase"])
    meta = load_dev_meta()
    ev = Evaluator(P, C, U, split_filter=FOLDS, bins=bins, zero_policy="ignore")
    conds = [c for c in ev.cond_index if c in meta.index]
    assert len(conds) == len(ev.cond_index), "eligible GT conditions missing from DEV cache"
    assert set(meta.loc[conds, "split"]) <= set(FOLDS)
    X = load_dev(conditions=conds)
    Z, N = zscore(X, LIB, kin)
    S = {"zscore": Z}
    S["random"] = ctl.random_scores(conds, kin, seed=0)
    # prior-only / shift-prior: frequencies from OTHER folds only (out-of-fold)
    fold_of = C.set_index("conditionID")["split"]
    po, sp_ = [], []
    for f in FOLDS:
        ec = [c for c in conds if fold_of[c] == f]
        tc = [c for c in C["conditionID"] if fold_of[c] in FOLDS and fold_of[c] != f]
        po.append(ctl.prior_only(P, C, tc, ec, kin))
        fr = ctl.training_frequency(P, C, tc, kin)
        sp_.append(ctl.shift_prior(X[ec], fr))
    S["prior_only"] = pd.concat(po).loc[conds]
    S["shift_prior"] = pd.concat(sp_).loc[conds]
    S["oracle"] = ctl.oracle(P, conds, kin)
    res, pol = {}, {"oracle": "ignore", "prior_only": "ignore", "shift_prior": "ignore"}
    for m, sc in S.items():
        ev.zero_policy = pol.get(m, "raise" if m == "zscore" else "ignore")
        r = ev.evaluate(sc)
        res[m] = r
    ev.zero_policy = "ignore"
    plan = bs.BootstrapPlan.from_evaluator(ev, B=2000, seed=0)
    deltas = {}
    for name, b in [("zscore_vs_prior_only", "prior_only"), ("zscore_vs_random", "random")]:
        d = bs.paired_delta(ev, S["zscore"], S[b], "M1-sel", plan)
        deltas[name] = {k: v for k, v in d.items() if k != "replicates"}
    # MDE
    G_dev = int(ev.targets[ev.targets["in_M1sel"]]["group"].nunique())
    se = deltas["zscore_vs_prior_only"]["se"]
    n_all_dev = int(C[C["split"].isin(FOLDS)]["group_id"].nunique())
    cg = read_gt_tsv(ROOT / "results/ground_truth/condition_groups.tsv")
    g_all_dev = int(cg[cg["split"].isin(FOLDS)]["group_id"].nunique())
    mde = {"SE_dev_dM1sel_zscore_vs_prior_only": se, "G_DEV_M1sel_groups": G_dev, "G_TEST": G_TEST_ALL,
           "G_TEST_note": "GT_REPORT split table: all groups in split 'test' (upper bound on M1-sel groups; "
                          "no label-derived count used)",
           "MDE_2.5": 2.5 * se * np.sqrt(G_dev / G_TEST_ALL),
           "MDE_3.5_alt": 3.5 * se * np.sqrt(G_dev / G_TEST_ALL),
           "alt_all_groups_ratio": {"G_DEV_all_groups": g_all_dev,
                                    "MDE_2.5": 2.5 * se * np.sqrt(g_all_dev / G_TEST_ALL)},
           "delta_estimate": deltas["zscore_vs_prior_only"]["estimate"], "B": 2000, "seed": 0,
           "folds": FOLDS}
    (OUT / "mde.json").write_text(json.dumps(to_jsonable(mde), indent=1))
    # outputs
    v = lambda x: x["value"] if isinstance(x, dict) and "value" in x else x
    summ = {}
    for m, r in res.items():
        summ[m] = {"M1-sel": v(r["M1-sel"]), "M1-all": v(r["M1-all"]), "M1-sel-scored": v(r["M1-sel-scored"]),
                   "M1-all-scored": v(r["M1-all-scored"]), "coverage_M1-sel": r["coverage_M1-sel"],
                   "coverage_M1-all": r["coverage_M1-all"], "harmony": v(r["harmony"]),
                   "harmony_per_bin": r["harmony"]["per_bin"],
                   "by_source": {k: v(x) for k, x in r["by_source"].items()},
                   "by_sign": {str(k): v(x) for k, x in r["by_sign"].items()},
                   "tie_share_M1-sel": r["tie_share_M1-sel"], "coverage_significant": r["coverage_significant"]}
    out = {"library": LIB, "min_substrates": MIN_SUB, "folds": FOLDS, "n_conditions": len(conds),
           "n_kinases": len(kin), "G_DEV_M1sel": G_dev, "metrics": summ, "deltas": deltas,
           "T_unseen": "not computed (pooled out-of-fold run; no single training-positive set)"}
    (OUT / "smoke_dev.json").write_text(json.dumps(to_jsonable(out), indent=1, default=str))
    L = [f"# Smoke run on DEV (folds 0-4 pooled; no devconfirm/test)\n",
         f"Library `{LIB}`, z-score = mean(x_s)*sqrt(n)/sd(x_c), min {MIN_SUB} substrates. "
         f"{len(conds)} eligible GT conditions, {G_dev} M1-sel groups. Controls: random (seed 0); prior-only and "
         "shift-prior use out-of-fold frequencies (other folds only); oracle. Runtime "
         f"{time.time()-t0:.0f} s. Script: `scripts/smoke_dev.py`.\n",
         "| method | M1-sel | M1-all | M1-sel-scored | pair coverage | harmony |", "|---|---|---|---|---|---|"]
    f3 = lambda v: "NA" if v is None or (isinstance(v, float) and np.isnan(v)) else f"{v:.3f}"
    for m, r in summ.items():
        L.append(f"| {m} | {f3(r['M1-sel'])} | {f3(r['M1-all'])} | {f3(r['M1-sel-scored'])} | "
                 f"{f3(r['coverage_M1-sel'])} | {f3(r['harmony'])} |")
    L += ["", "## By source / sign (M1-sel)", "", "| method | H | D | O | sign -1 | sign +1 |", "|---|---|---|---|---|---|"]
    for m, r in summ.items():
        bsr, bsg = r["by_source"] or {}, r["by_sign"] or {}
        g = lambda d, *ks: next((d[k] for k in ks if k in d), None)
        L.append(f"| {m} | {f3(g(bsr,'H'))} | {f3(g(bsr,'D'))} | {f3(g(bsr,'O'))} | "
                 f"{f3(g(bsg,'-1'))} | {f3(g(bsg,'1'))} |")
    L += ["", "## Harmony per bin (M1-sel) and coverage-significant kinases", "", "| method | 0 | 1-4 | 5-19 | >=20 | n sig (q<=0.10) |", "|---|---|---|---|---|---|"]
    for m, r in summ.items():
        hb = r["harmony_per_bin"]
        L.append(f"| {m} | " + " | ".join(f3(hb.get(k)) for k in ["0", "1-4", "5-19", ">=20"]) +
                 f" | {r['coverage_significant']['n_significant']} |")
    L += ["", "## Paired bootstrap (groups, B=2000, seed 0), M1-sel", ""]
    for k, d in deltas.items():
        L.append(f"- {k}: Δ={d['estimate']:+.3f}, SE={d['se']:.3f}, 95% CI [{d['ci_low']:+.3f}, {d['ci_high']:+.3f}], "
                 f"one-sided p={d['p_one_sided']:.4f}")
    L += ["", f"## MDE\n", f"SE(ΔM1-sel, z-score vs prior-only) = {se:.4f}; G_DEV = {G_dev}, G_TEST = {G_TEST_ALL}; "
          f"MDE = 2.5·SE·sqrt(G_DEV/G_TEST) = {mde['MDE_2.5']:.4f} (see results/smoke/mde.json).",
          "", "Raw metric dicts (incl. coverage-significant): `smoke_dev.json`. T-unseen not computed."]
    (OUT / "smoke_dev.md").write_text("\n".join(L) + "\n")
    for m, r in summ.items():
        lm = {k: r[k] for k in ["M1-sel", "M1-all", "M1-sel-scored", "coverage_M1-sel"]}
        ledger.append(f"smoke_{m}", {"library": LIB, "min_substrates": MIN_SUB, "method": m}, "dev", to_jsonable(lm),
                      branch="smoke", folds=range(5))
    print("\n".join(L))


if __name__ == "__main__":
    main()
