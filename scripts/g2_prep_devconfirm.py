#!/usr/bin/env python
"""DEV-confirm look #2 (DESIGN §7; D11/D12). Logged; after this look 3 remain.

FROZEN BEFORE ANY DEV-CONFIRM VALUE IS LOADED (pre-declared here; written to results/g2_prep/devconfirm_look2_frozen.json
before the devconfirm columns are read):
  core (core-v1)  = DEV-confirm look-1 frozen config (results/audit_core_v1/devconfirm_frozen_config.json, hash
                    04b69b6513e3927d): curated_union + r2np@100 (a=1), wsd, no winsor, borrow, tau .125, rho_f 3,
                    rho_g .25, z_vif with rho_e refitted on the 7 DEV-CV nulls (must equal 0.01).
  signatures      = B7 clean+contrast (PRIMARY) fitted on ALL DEV-CV GT conditions (folds 0-4; C10 i: devconfirm groups
                    are disjoint from DEV-CV groups by construction, checked).
  (n_top, kappa)  = all-DEV rule: argmax of pooled 5-fold OOF M1-sel in the seed-0 B7 primary grid
                    (results/b7/grid_clean+contrast.tsv, signature rows only).
  alpha (fusion)  = all-DEV rule: argmax over {.25,.5,1} of pooled seed-0 OOF M1-sel of the late fusion with that
                    (n_top, kappa) in every fold (outer signature fits on the 4 other folds, per-fold core-v1 picks).
  rho_sig         = RMS of z_sig over the 7 DEV-CV nulls (final fit).
  core-v2         = (z_core + a z_sig)/sqrt(1+a^2).
  core-v3 (E1)    = F_k/sqrt(v_k), v_k = mean F_k^2 over site-permuted (seed 0) DEV-CV GT conditions, final fit.
  comparators     = PSP z-score min3, curated_union z-score min1, KSEA curated_union min1, prior-only (freq from DEV-CV
                    GT conditions, oracle sign).
  inference       = paired cluster bootstrap B=2000 seed 0, C8 (basic CI, shifted-null one-sided p).
  D11 condition   = point Delta M1-sel (devconfirm) of core-v2 AND core-v3 vs core-v1 > 0; else revert to core-v1.
Data: only split=='devconfirm' columns of cache/quant_wide_f32.npy (as look 1). Never cache/sealed.
"""
import json
import pickle
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/workspace")
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))
from exp_b7_common import B7Ctx, FOLDS, frame, CACHE as B7CACHE  # noqa: E402
from exp_b7_eval import m1sel  # noqa: E402
from exp_b7b_eval import Late, VARIANT  # noqa: E402
from exp_b21_common import calib_bins  # noqa: E402
from audit_core_v1_devconfirm import load_devconfirm  # noqa: E402
from kinact import bootstrap as bs, controls as ctl, ledger  # noqa: E402
from kinact.eval import Evaluator, to_jsonable  # noqa: E402
from kinact.explore import bayes2 as B2, signatures_late as SL, final_box as FB, signatures as SG  # noqa: E402
from kinact.explore.bayes import report_mask  # noqa: E402
from kinact.methods.linear import zscore, ksea  # noqa: E402

OUT = ROOT / "results/g2_prep"
LOOK = 2
FUSE_ALPHAS = [0.25, 0.5, 1.0]


def main():
    import argparse
    ap = argparse.ArgumentParser(); ap.add_argument("--dryrun", action="store_true",
                                                    help="code-path test: DEV fold0 stands in for devconfirm; nothing logged")
    DRY = ap.parse_args().dryrun
    global OUT
    if DRY:
        OUT = OUT / "dryrun"; OUT.mkdir(exist_ok=True)
    ctx = B7Ctx(); c = ctx.c
    L = Late(ctx)
    look1 = json.loads((ROOT / "results/audit_core_v1/devconfirm_frozen_config.json").read_text())
    core = (look1["tau"], look1["rho_f"], look1["rho_g"])
    K = ctx.kernel(core[1], core[2])

    # ---------------- all-DEV HP rules (DEV-CV only)
    g7 = pd.read_csv(ROOT / "results/b7/grid_clean+contrast.tsv", sep="\t")
    g7 = g7[g7["alpha"] > 0]
    b = g7.loc[g7["oof_M1-sel"].idxmax()]
    n_top, kappa = int(b["n_top"]), float(b["kappa"])
    zcore_f = {f: ctx.score(ctx.fold_pos[f], [("x", None)], ctx.picks[f], [])[("x", 0.0)][0] for f in FOLDS}
    outer_fits = {f: pickle.load(open(B7CACHE / f"sig_{VARIANT}_{f}.pkl", "rb"))[kappa] for f in FOLDS}
    zs_f = {f: L.zsig(outer_fits[f], n_top, "gt", ctx.fold_pos[f])[0] for f in FOLDS}
    a_oof = {a: m1sel(c.ev, frame(ctx, {f: SL.fuse(zcore_f[f], zs_f[f], a) for f in FOLDS}), FOLDS, c.fold_of)
             for a in FUSE_ALPHAS}
    alpha = max(FUSE_ALPHAS, key=lambda a: (a_oof[a], -FUSE_ALPHAS.index(a)))

    # ---------------- final fits on ALL DEV-CV
    C, P = c.C, c.P
    DC = "fold0" if DRY else "devconfirm"
    dc_all = [x for x in C.loc[C.split == DC, "conditionID"] if (not DRY) or x in c.X.columns]
    gt_dc = [x for x in C.loc[(C.split == DC) & (C.role == "gt"), "conditionID"] if x in dc_all]
    dc_null = [x for x in C.loc[(C.split == DC) & (C.role == "null"), "conditionID"] if x in dc_all]
    if not DRY:
        assert not set(C.loc[C.split == "devconfirm", "group_id"]) & set(C.loc[C.split.isin(FOLDS), "group_id"])
    fit = ctx.fit(FOLDS, dc_all, VARIANT, kappas=[kappa])[kappa]
    Ws = fit.matrix(len(c.X.index), c.kin, n_top)
    nidx = np.arange(len(c.nulls))
    rho_sig = SL.rho_sig(SL.sig_project(Ws, *ctx._YM("null", nidx)))
    # core on nulls -> rho_e; core + z_sig on permuted DEV-CV GT conditions -> v_k
    res0 = ctx.score(np.arange(len(c.conds)), [("x", None)], core, [], which="perm")
    zc_perm, rho_e = res0[("x", 0.0)]
    assert abs(rho_e - look1["rho_e"]) < 1e-9, (rho_e, look1["rho_e"])
    zs_perm = SL.sig_project(Ws, ctx.Yp, ctx.Mp) / rho_sig
    sigmask = np.isin(np.asarray(c.kin), list(fit.sites))
    v, vs = FB.perm_scales(zc_perm.astype(np.float64), zs_perm, alpha, sigmask)
    frozen = {"name": "core-v3", "core": {k: look1[k] for k in ["lib_cur", "alphas", "scale", "winsor", "borrow", "tau",
                                                                  "rho_f", "rho_g", "rho_e", "config_hash"]},
              "signature_variant": VARIANT, "signature_fit": "all DEV-CV GT conditions (folds 0-4)",
              "n_train_conds": len(fit.train_conds), "n_signature_kinases": int(sigmask.sum()),
              "n_top": n_top, "kappa": kappa, "ntop_kappa_rule": "argmax pooled seed-0 OOF M1-sel, B7 primary grid",
              "b7_grid_best_row": {k: float(b[k]) for k in ["alpha", "n_top", "kappa", "oof_M1-sel"]},
              "alpha": alpha, "alpha_rule": "argmax pooled seed-0 OOF M1-sel of late fusion with fixed (n_top,kappa)",
              "alpha_oof": a_oof, "rho_sig": rho_sig, "E1": "v_k on site-permuted DEV-CV GT conditions (seed 0)",
              "v_k_median": float(np.nanmedian(v)), "v_k_range": [float(np.nanmin(v)), float(np.nanmax(v))],
              "look": LOOK}
    chash = ledger.config_hash({k: x for k, x in frozen.items() if k not in ("alpha_oof",)} | {"sig_kinases": sorted(fit.sites)})
    frozen["config_hash"] = chash
    (OUT / "devconfirm_look2_frozen.json").write_text(json.dumps(to_jsonable(frozen), indent=1))
    pickle.dump({"fit": fit, "v": v, "vs": vs}, open(ROOT / "cache/g2_prep/core_v3_alldev_fit.pkl", "wb"))
    print("FROZEN", json.dumps(to_jsonable(frozen)), flush=True)

    # ---------------- DEV-confirm data (only now)
    X = c.X[dc_all] if DRY else load_devconfirm(dc_all)
    assert X.index.equals(c.X.index)
    Y, M, _ = B2.standardise(X.to_numpy(np.float64), "wsd", None)
    Y, M = Y.astype(np.float32), M.astype(np.float32)

    def core_z(Yx, Mx):
        bl = B2.block_stats(c.md, Yx, Mx, names=["cur", "r2np@100"])
        G, h, g0, n = B2.combine(bl, {"r2np@100": 1.0})
        post = B2.posterior2(G, h, g0, K, core[0])
        rep = report_mask(n, K, 1, True)
        nn = n.cpu().numpy() if hasattr(n, "cpu") else n
        return np.where(rep, B2.output_stat(post, "z_vif", rho_e=rho_e), np.nan), nn

    zc, n_dc = core_z(Y, M)
    zsg = SL.sig_project(Ws, Y, M) / rho_sig
    kin = c.kin
    mk = lambda A, idx=dc_all: pd.DataFrame(np.asarray(A), index=idx, columns=kin)
    v1, v2, v3 = mk(zc), mk(SL.fuse(zc, zsg, alpha)), mk(FB.fuse_calibrated(zc, zsg, alpha, v, vs))
    meth = {"core-v3": v3, "core-v2": v2, "core-v1": v1,
            "PSP z-score min3": zscore(X, "psp_all", min_substrates=3).reindex(columns=kin),
            "curated_union z-score min1": zscore(X, "curated_union", min_substrates=1).reindex(columns=kin),
            "KSEA curated_union min1": ksea(X, "curated_union", min_substrates=1).reindex(columns=kin),
            "prior-only (DEV-CV freq)": ctl.prior_only(P, C, list(C.loc[C.split.isin(FOLDS), "conditionID"]), gt_dc, kin)}
    dpos = P[(P.cls == "direct") & P.confidence.isin(["high", "medium"]) &
             P.conditionID.isin(set(C.loc[C.split.isin(FOLDS) & (C.role == "gt"), "conditionID"]))]
    bins = pd.read_csv(ROOT / "results/ground_truth/harmony_bins_group.tsv", sep="\t", keep_default_na=False, na_values=[""])
    ev = Evaluator(P, C, c.Udf, split_filter=[DC], bins=bins, train_pos_kinases=set(dpos["kinase"]),
                   zero_policy="ignore")
    plan = bs.BootstrapPlan.from_evaluator(ev, B=2000, seed=0)
    val = lambda x: x["value"] if isinstance(x, dict) and "value" in x else x
    res = {}
    for name, S in meth.items():
        r = ev.evaluate(S, n_perm=2000)
        res[name] = {"M1-sel": val(r["M1-sel"]), "M1-all": val(r["M1-all"]), "M1-sel-scored": val(r["M1-sel-scored"]),
                     "coverage_M1-sel": r["coverage_M1-sel"], "harmony": val(r["harmony"]),
                     "harmony_per_bin": r["harmony"]["per_bin"], "by_source": {k: val(x) for k, x in r["by_source"].items()},
                     "T-unseen": val(r["T-unseen"]), "T-unseen_n_pairs": r["T-unseen"]["n_pairs"],
                     "n_sig": r["coverage_significant"]["n_significant"], "n_M1sel_pairs": r["M1-sel"]["n_pairs"],
                     "n_kinases": r["M1-sel"]["n_kinases"]}
    sigk = set(fit.sites)
    for name in ["core-v3", "core-v2", "core-v1"]:
        Pp = ev.pairs(meth[name]); m = Pp["kinase"].isin(sigk).to_numpy()
        res[name]["has-signature"] = ev.metric(Pp[m], "M1-sel")["value"]
        res[name]["no-signature"] = ev.metric(Pp[~m], "M1-sel")["value"]
        res[name]["n_pairs_has_signature"] = int((m & Pp["in_M1sel"].to_numpy()).sum())
    pick = lambda d: {k: d[k] for k in ["estimate", "ci_low", "ci_high", "p_one_sided"]}
    deltas = {}
    for a_, b_ in [("core-v3", "core-v1"), ("core-v2", "core-v1"), ("core-v3", "core-v2")] + \
                  [(y, x) for y in ["core-v3", "core-v2", "core-v1"] for x in list(meth)[3:]]:
        dd = {ep: pick(bs.paired_delta(ev, meth[a_], meth[b_], ep, plan)) for ep in ["M1-sel", "M1-all", "harmony"]}
        if not b_.startswith("prior"):
            dd["M2"] = ev.m2(meth[a_], meth[b_])
        deltas[f"{a_} vs {b_}"] = dd
    d11 = {"dM1sel_v2_vs_v1": deltas["core-v2 vs core-v1"]["M1-sel"]["estimate"],
           "dM1sel_v3_vs_v1": deltas["core-v3 vs core-v1"]["M1-sel"]["estimate"]}
    d11["pass"] = bool(d11["dM1sel_v2_vs_v1"] > 0 and d11["dM1sel_v3_vs_v1"] > 0)
    d11["outcome"] = "keep core-v3" if d11["pass"] else "REVERT to core-v1"
    # calibration: devconfirm nulls + site-permuted devconfirm GT conditions
    pn = [dc_all.index(x) for x in dc_null]
    Xp = ctl.site_permutation(X[gt_dc], seed=0)
    Yp, Mp, _ = B2.standardise(Xp.to_numpy(np.float64), "wsd", None)
    Yp, Mp = Yp.astype(np.float32), Mp.astype(np.float32)
    zcp, n_p = core_z(Yp, Mp); zsp = SL.sig_project(Ws, Yp, Mp) / rho_sig
    sm_n = np.broadcast_to(sigmask, (len(pn), len(kin))); sm_p = np.broadcast_to(sigmask, zcp.shape)
    cal = {}
    for lab, Zn, Zp in [("core-v1", zc[pn], zcp), ("core-v2", SL.fuse(zc, zsg, alpha)[pn], SL.fuse(zcp, zsp, alpha)),
                        ("core-v3", FB.fuse_calibrated(zc, zsg, alpha, v, vs)[pn], FB.fuse_calibrated(zcp, zsp, alpha, v, vs))]:
        cal[lab] = {"dc_nulls": calib_bins(Zn, n_dc[pn]), "dc_perm": calib_bins(Zp, n_p),
                    "dc_nulls_sigkin": calib_bins(Zn, n_dc[pn], ok=sm_n), "dc_perm_sigkin": calib_bins(Zp, n_p, ok=sm_p)}
    cal["n_null_conds"] = len(pn); cal["n_null_cells"] = int(np.isfinite(zc[pn]).sum())
    date = datetime.now(timezone.utc).isoformat(timespec="seconds")
    out = {"frozen": frozen, "look": LOOK, "date": date, "n_devconfirm_conditions_loaded": len(dc_all),
           "results": res, "deltas": deltas, "D11": d11, "calibration": cal}
    (OUT / "devconfirm_look2.json").write_text(json.dumps(to_jsonable(out), indent=1, default=str))
    print(json.dumps(to_jsonable({k: {kk: res[k][kk] for kk in ["M1-sel", "M1-all", "harmony", "T-unseen"]} for k in res})), flush=True)
    print("D11", d11, flush=True)

    # ---------------- logging (ledger + devconfirm_looks.tsv)
    if DRY:
        print("dry run: nothing logged"); return
    for name, cfgd in [("core-v3", frozen), ("core-v2", {**frozen, "name": "core-v2", "E1": None}),
                       ("core-v1", {**look1, "look": LOOK})]:
        r = res[name]
        m = {kk: r[kk] for kk in ["M1-sel", "M1-all", "M1-sel-scored", "coverage_M1-sel", "harmony", "T-unseen"]}
        m["dM1sel_vs_pspz"] = deltas[f"{name} vs PSP z-score min3"]["M1-sel"]
        if name != "core-v1":
            m["dM1sel_vs_core_v1"] = deltas[f"{name} vs core-v1"]["M1-sel"]
        m["look"] = LOOK
        ledger.append(f"{name} (DEV-confirm look {LOOK})", cfgd, "devconfirm", m, branch="DEVCONFIRM", folds="devconfirm")
    lt = ROOT / "results/devconfirm_looks.tsv"
    rows = []
    for name in meth:
        r = res[name]
        dkey = f"core-v3 vs {name}" if name != "core-v3" else None
        note = {"core-v3": "primary (bundle)", "core-v2": "D11 check", "core-v1": "D11 reference (look-1 config)"}.get(name, "comparator")
        rows.append({"look": LOOK, "date": date[:10], "bundle": "core-v3", "config_hash": chash, "method": name,
                     "M1-sel": round(r["M1-sel"], 4), "M1-all": round(r["M1-all"], 4),
                     "M1-sel-scored": round(r["M1-sel-scored"], 4), "coverage": round(r["coverage_M1-sel"], 4),
                     "harmony": round(r["harmony"], 4), "T-unseen": round(r["T-unseen"], 4), "n_M1sel_pairs": r["n_M1sel_pairs"],
                     "dM1sel_vs_pspz": json.dumps({k: round(x, 4) for k, x in deltas[f"{name} vs PSP z-score min3"]["M1-sel"].items()})
                     if f"{name} vs PSP z-score min3" in deltas else "",
                     "note": note + (f"; dM1sel core-v3 minus this {deltas[dkey]['M1-sel']['estimate']:+.4f}" if dkey in deltas else "")
                     + ("; D11 " + d11["outcome"] if name == "core-v3" else "")})
    pd.DataFrame(rows).to_csv(lt, sep="\t", index=False, mode="a", header=False)


if __name__ == "__main__":
    main()
