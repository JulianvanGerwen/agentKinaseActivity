#!/usr/bin/env python
"""DEV-confirm look #1 for core-v1 (DESIGN §7 DEV-confirm semantics; docs/audit/core_v1_audit.md part 2).

Frozen BEFORE loading any DEV-confirm value (pre-declared here):
  library   curated_union (w=1) + r2np global top-100 binary, alpha=1; borrow on; median-centred; wsd scale; NO winsor
  prior     tau=0.125, rho_f=3 (D6 fixed); rho_g = argmax pooled M1-sel over ALL 5 DEV folds of the z_vif output
            among {0, 0.25, 1} (results/audit_core_v1/dev_grid_rescored.tsv)
  output    PRIMARY  z_vif with rho_e fitted (fit_rho_null) on the 7 DEV role=null conditions only
            SECONDARY z_post (ranking statistic of the D6 adoption evidence), same posterior
Comparators: PSP z-score (psp_all, min 3), curated_union z-score (min 1), prior-only (freq from all DEV folds, oracle sign).
Data: only split=='devconfirm' columns of cache/quant_wide_f32.npy (memmap column subset) + DEV nulls from cache/dev.
"""
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path("/workspace")
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))
from kinact import bootstrap as bs, controls as ctl, ledger  # noqa: E402
from kinact.data import load_dev, SITES_PQ  # noqa: E402
from kinact.eval import Evaluator, load_gt, read_gt_tsv, to_jsonable, calibration  # noqa: E402
from kinact.explore import bayes2 as B2  # noqa: E402
from kinact.explore.bayes import report_mask  # noqa: E402
from kinact.library import load_library  # noqa: E402
from kinact.methods.linear import zscore  # noqa: E402
from exp_b21_common import calib_bins  # noqa: E402

OUT = ROOT / "results/audit_core_v1"
FOLDS = [f"fold{i}" for i in range(5)]


def load_devconfirm(conds):
    q = np.load(ROOT / "cache/quant_wide_f32.npy", mmap_mode="r")
    peps = (ROOT / "cache/quant_wide_sites.txt").read_text().split()
    allc = (ROOT / "cache/quant_wide_conditions.txt").read_text().split()
    assert peps == (ROOT / "cache/dev/sites_dev.txt").read_text().split()
    ci = pd.Index(allc).get_indexer(conds)
    assert (ci >= 0).all()
    M = np.asarray(q[:, ci], dtype=np.float32)          # only these columns are read
    sites = pd.read_parquet(SITES_PQ, columns=["peptide_UID", "site_key"]).set_index("peptide_UID")["site_key"].reindex(peps).to_numpy()
    codes, uniq = pd.factorize(sites)
    ok = ~np.isnan(M)
    S = sp.csr_matrix((np.ones(len(codes), np.float32), (codes, np.arange(len(codes)))), shape=(len(uniq), len(codes)))
    num = S @ np.where(ok, M, 0.0); cnt = S @ ok.astype(np.float32)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.where(cnt > 0, num / cnt, np.nan)
    return pd.DataFrame(out.astype(np.float32), index=uniq, columns=list(conds))


def main():
    P, C, U, _ = load_gt()
    Udf = read_gt_tsv(ROOT / "results/ground_truth/kinase_universe.tsv"); kin = list(Udf["kinase"])
    bins = read_gt_tsv(ROOT / "results/ground_truth/harmony_bins_group.tsv")
    assert not (C["split"].astype(str) == "test").any()

    # ---------------- freeze (DEV only)
    g = pd.read_csv(OUT / "dev_grid_rescored.tsv", sep="\t")
    sub = g[(g.output == "z_vif") & (g.tau == 0.125) & (g.rho_f == 3.0)]
    best = sub.loc[sub["M1sel"].idxmax()]
    frozen = dict(name="core-v1", lib_cur="curated_union", alphas={"r2np@100": 1.0}, scale="wsd", winsor=None,
                  borrow=True, tau=0.125, rho_f=3.0, rho_g=float(best.rho_g), output_primary="z_vif",
                  output_secondary="z_post", rho_e_source="7 DEV role=null conditions (folds 0-4)",
                  rho_g_selection="pooled all-DEV M1-sel, z_vif, rho_g in {0,0.25,1}",
                  dev_pooled_M1sel=float(best.M1sel))
    nulls = [c for c in C.loc[(C.role == "null") & C.split.isin(FOLDS), "conditionID"]]
    Xn = load_dev(conditions=nulls)
    lib_cur = load_library("curated_union"); r2 = load_library("r2np", k=100)
    md = B2.build_multi_design(Xn.index, kin, lib_cur, {"r2np@100": r2})
    K = B2.kernel2(Udf, kin, frozen["rho_f"], frozen["rho_g"])

    def post_of(X):
        Y, Mk, _ = B2.standardise(X.to_numpy(np.float64), "wsd", None)
        bl = B2.block_stats(md, Y, Mk, names=["cur", "r2np@100"])
        Gm, h, g0, n = B2.combine(bl, frozen["alphas"])
        post = B2.posterior2(Gm, h, g0, K, frozen["tau"])
        rep = report_mask(n, K, 1, True)
        return post, rep, n

    pn, rn, _ = post_of(Xn)
    rho_e, _ = B2.fit_rho_null(B2.output_stat(pn, "z_freq"), B2.m_eff(pn["q"], pn["suma"]), rn)
    frozen["rho_e"] = float(rho_e)
    chash = ledger.config_hash(frozen)
    frozen["config_hash"] = chash
    (OUT / "devconfirm_frozen_config.json").write_text(json.dumps(frozen, indent=1))
    print("FROZEN", frozen, flush=True)

    # ---------------- DEV-confirm data (only now)
    dc_all = list(C.loc[C.split == "devconfirm", "conditionID"])
    X = load_devconfirm(dc_all)
    assert X.index.equals(Xn.index)
    post, rep, n = post_of(X)
    Zv = np.where(rep, B2.output_stat(post, "z_vif", rho_e=rho_e), np.nan)
    Zp = np.where(rep, B2.output_stat(post, "z_post"), np.nan)
    core_v = pd.DataFrame(Zv, index=dc_all, columns=kin); core_p = pd.DataFrame(Zp, index=dc_all, columns=kin)
    psp = zscore(X, "psp_all", min_substrates=3).reindex(columns=kin)
    cur = zscore(X, "curated_union", min_substrates=1).reindex(columns=kin)
    dev_train = list(C.loc[C.split.isin(FOLDS), "conditionID"])
    gt_dc = list(C.loc[(C.split == "devconfirm") & (C.role == "gt"), "conditionID"])
    prior = ctl.prior_only(P, C, dev_train, gt_dc, kin)
    dpos = P[(P.cls == "direct") & P.confidence.isin(["high", "medium"]) &
             P.conditionID.isin(set(C.loc[C.split.isin(FOLDS) & (C.role == "gt"), "conditionID"]))]
    train_pos = set(dpos["kinase"])
    ev = Evaluator(P, C, U, split_filter=["devconfirm"], bins=bins, train_pos_kinases=train_pos, zero_policy="ignore")
    plan = bs.BootstrapPlan.from_evaluator(ev, B=2000, seed=0)
    meth = {"core-v1 (z_vif, primary)": core_v, "core-v1 (z_post, secondary)": core_p,
            "PSP z-score min3": psp, "curated_union z-score min1": cur, "prior-only (DEV freq)": prior}
    v = lambda x: x["value"] if isinstance(x, dict) and "value" in x else x
    res = {}
    for name, S in meth.items():
        r = ev.evaluate(S, n_perm=2000)
        res[name] = {"M1-sel": v(r["M1-sel"]), "M1-all": v(r["M1-all"]), "M1-sel-scored": v(r["M1-sel-scored"]),
                     "coverage_M1-sel": r["coverage_M1-sel"], "harmony": v(r.get("harmony")),
                     "harmony_per_bin": r["harmony"]["per_bin"] if "harmony" in r else None,
                     "by_source": {k: v(x) for k, x in r["by_source"].items()},
                     "by_sign": {str(k): v(x) for k, x in r["by_sign"].items()},
                     "T-unseen": v(r["T-unseen"]), "T-unseen_n_pairs": r["T-unseen"]["n_pairs"],
                     "n_sig": r["coverage_significant"]["n_significant"],
                     "n_M1sel_pairs": r["M1-sel"]["n_pairs"], "n_kinases": r["M1-sel"]["n_kinases"]}
    for prim in ["core-v1 (z_vif, primary)", "core-v1 (z_post, secondary)"]:
        for comp in ["PSP z-score min3", "curated_union z-score min1", "prior-only (DEV freq)"]:
            dd = {}
            for ep in ["M1-sel", "M1-all", "harmony", "M1-sel:unseen"]:
                d = bs.paired_delta(ev, meth[prim], meth[comp], ep, plan)
                dd[ep] = {k: d[k] for k in ["estimate", "ci_low", "ci_high", "p_one_sided"]}
            if not comp.startswith("prior"):
                m2 = ev.m2(meth[prim], meth[comp]); dd["M2"] = m2
            res[prim][f"vs {comp}"] = dd
    # calibration on devconfirm nulls (N-eval-like) + site-permuted devconfirm GT conditions
    dc_null = list(C.loc[(C.split == "devconfirm") & (C.role == "null"), "conditionID"])
    pos_n = [dc_all.index(c) for c in dc_null]; pos_g = [dc_all.index(c) for c in gt_dc]
    Xp = ctl.site_permutation(X[gt_dc], seed=0)
    pp, rp, n_p = post_of(Xp)
    cal = {"dc_nulls": calib_bins(Zv[pos_n], n[pos_n]), "dc_nulls_ncells": int(np.isfinite(Zv[pos_n]).sum()),
           "dc_perm": calib_bins(np.where(rp, B2.output_stat(pp, "z_vif", rho_e=rho_e), np.nan), n_p),
           "dc_nulls_zpost": calib_bins(Zp[pos_n], n[pos_n])}
    res["calibration_core_v1_zvif"] = cal
    out = {"frozen": frozen, "n_devconfirm_conditions_loaded": len(dc_all), "results": res,
           "date": datetime.now(timezone.utc).isoformat(timespec="seconds")}
    (OUT / "devconfirm_look1.json").write_text(json.dumps(to_jsonable(out), indent=1, default=str))
    for k, r in res.items():
        print(k, json.dumps(to_jsonable(r), default=str)[:1500], flush=True)

    # ---------------- logging
    for name in ["core-v1 (z_vif, primary)", "core-v1 (z_post, secondary)"]:
        r = res[name]
        m = {kk: r[kk] for kk in ["M1-sel", "M1-all", "M1-sel-scored", "coverage_M1-sel", "harmony", "T-unseen"]}
        m["dM1sel_vs_pspz"] = r["vs PSP z-score min3"]["M1-sel"]; m["dM1sel_vs_curz"] = r["vs curated_union z-score min1"]["M1-sel"]
        m["dM1sel_vs_prior"] = r["vs prior-only (DEV freq)"]["M1-sel"]; m["look"] = 1
        ledger.append(f"core-v1 {name.split('(')[1].rstrip(')')}", {**frozen, "look": 1, "statistic": name},
                      "devconfirm", m, branch="DEVCONFIRM", folds="devconfirm")
    lt = OUT.parent / "devconfirm_looks.tsv"
    rows = []
    for name in meth:
        r = res[name]
        rows.append({"look": 1, "date": out["date"][:10], "bundle": "core-v1", "config_hash": chash, "method": name,
                     "M1-sel": round(r["M1-sel"], 4), "M1-all": round(r["M1-all"], 4),
                     "M1-sel-scored": round(r["M1-sel-scored"], 4), "coverage": round(r["coverage_M1-sel"], 4),
                     "harmony": round(r["harmony"], 4), "T-unseen": round(r["T-unseen"], 4),
                     "n_M1sel_pairs": r["n_M1sel_pairs"],
                     "dM1sel_vs_pspz": json.dumps({k: round(x, 4) for k, x in r["vs PSP z-score min3"]["M1-sel"].items()}) if "vs PSP z-score min3" in r else "",
                     "note": "primary" if "primary" in name else ("secondary" if "secondary" in name else "comparator")})
    df = pd.DataFrame(rows)
    df.to_csv(lt, sep="\t", index=False, mode="a" if lt.exists() else "w", header=not lt.exists())


if __name__ == "__main__":
    main()
