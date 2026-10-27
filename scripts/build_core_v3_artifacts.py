#!/usr/bin/env python
"""Build the core-v3 model artefacts (models/core_v3/) for kinact.core.KinaseActivityModel (RESEARCH_TREE D11/D12).

core-v3 = core-v1 (same library / kernel / rho_e; the core-v1 artefacts are rebuilt into models/core_v3 by the same
code) + late-fused perturbation-response signatures (B7b, alpha = 1) + the E1 per-signature-kinase permutation-null
scale (B8). Everything is DEV-CV only (folds 0-4; cache/dev + the DEV rows of results/ground_truth):

  1. (n_top, kappa, alpha) are TAKEN from the pre-declared DEV-confirm look #2 freeze
     (results/g2_prep/devconfirm_look2_frozen.json: n_top/kappa = argmax pooled seed-0 OOF M1-sel of the B7 primary
     grid; alpha = argmax OOF of the late fusion), so the deployed/TEST model is exactly the DEV-confirmed one. The
     build re-derives rho_sig and the E1 scale range and fails if they disagree with that freeze record.
     Descriptive only: models/core_v3/sig_selection.tsv = pooled DEV-CV OOF M1-sel of the B7b late fusion at alpha=1
     over n_top x kappa (cache/b7 outer fits); it is flat (0.676-0.683) and is NOT used for the choice.
  2. Signatures: kinact.explore.signatures.fit_signatures, variant clean+contrast (n_targets <= 10, generic-response
     contrast), on ALL DEV-CV GT training conditions (C10 i: "at TEST, fitted on DEV-CV conditions only"); DEV-confirm
     and TEST groups are not in cache/dev (drug-group disjoint by construction; run_test.py re-asserts it).
  3. rho_sig = RMS of the raw z_sig over the 7 DEV role=null conditions (B7b rule), computed with the deployed scorer.
  4. E1 scales: v_k = mean (z_core + alpha z_sig)^2 and v_sig_k = mean z_sig^2 over site-permuted DEV-CV GT
     conditions (kinact.controls.site_permutation seed 0, within condition, the exact permuted matrix of B7b/B8),
     >= 30 cells else pooled over signature kinases (kinact.explore.final_box.perm_scales). Label-free.
  5. rho_e refit on the DEV nulls (must reproduce the config, as for core-v1).
Nothing from cache/sealed/, no DEV-confirm or TEST values.

Outputs (models/core_v3/): edges.tsv.gz (PSP-derived -> git-ignored), kinases.tsv, signatures.tsv.gz (compendium
derived, small), e1_scales.tsv, sig_selection.tsv, signature_training_conditions.txt, manifest.json (sha256 of every
artefact + inputs; committed).
Usage: build_core_v3_artifacts.py [--config config/core_v3.json] [--write-config]   (--write-config fills the
selected n_top/kappa and rho_sig into the config on first build; afterwards the build only verifies them.)
"""
import argparse
import hashlib
import json
import pickle
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))
from build_core_artifacts import build_edges  # noqa: E402
from kinact import controls as ctl, core  # noqa: E402
from kinact.data import load_dev, load_dev_meta  # noqa: E402
from kinact.eval import Evaluator, load_gt, read_gt_tsv  # noqa: E402
from kinact.explore import bayes2 as B2, final_box as FB, signatures as SG, signatures_late as SL  # noqa: E402
from kinact.library import LIBDIR, sha256_file  # noqa: E402

FOLDS = [f"fold{i}" for i in range(5)]
VARIANT = "clean+contrast"
NTOPS, KAPPAS = [20, 50, 100, 200], [0.0, 1.0, 4.0]


def sha_bytes(b):
    return hashlib.sha256(b).hexdigest()


class DevData:
    """DEV-CV inputs exactly as the B7/B7b/B8 research context (exp_b21_common.Ctx / exp_b7_common.B7Ctx)."""

    def __init__(self):
        P, C, U, _ = load_gt()
        assert not (C["split"].astype(str) == "test").any()
        self.P, self.C = P, C
        self.Udf = read_gt_tsv(ROOT / "results/ground_truth/kinase_universe.tsv")
        self.kin = list(self.Udf["kinase"])
        bins = read_gt_tsv(ROOT / "results/ground_truth/harmony_bins_group.tsv")
        self.ev = Evaluator(P, C, U, split_filter=FOLDS, bins=bins, zero_policy="ignore")
        self.meta = load_dev_meta()
        assert set(self.meta["split"]) <= set(FOLDS)
        self.X = load_dev()                                        # float32, sites x DEV-cache conditions
        self.conds = list(self.ev.cond_index)
        self.fold_of = self.meta["split"]
        self.cf = pd.Series(self.conds).map(self.fold_of).to_numpy()
        self.nulls = [c for c in C.loc[(C["role"] == "null") & C["split"].isin(FOLDS), "conditionID"]
                      if c in self.X.columns]
        Y, M, _ = B2.standardise(self.X.to_numpy(np.float64), "wsd", None)
        self.Yfit = np.where(M > 0, Y, np.nan).astype(np.float32)  # signature fitting input (B7Ctx.Yfit)
        self.Y, self.M = Y.astype(np.float32), M.astype(np.float32)
        self.allconds = pd.Index(self.X.columns)
        self.sites = pd.Index(self.X.index)

    def gt_hash(self):
        return sha256_file(ROOT / "results/ground_truth/GT_FREEZE.json")


# ----------------------------------------------------------------------------------------------- 1. selection
def selection_table(d: DevData, alpha: float) -> pd.DataFrame:
    """Pooled DEV-CV OOF M1-sel of the B7b late fusion (alpha fixed) for every (n_top, kappa) of the B7 grid."""
    cache = ROOT / "cache/b7"
    store = pickle.load(open(cache / f"grid_{VARIANT}.pkl", "rb"))
    zcore = {f: store[("core", f, f)][0] for f in FOLDS}
    del store
    fits = {f: pickle.load(open(cache / f"sig_{VARIANT}_{f}.pkl", "rb")) for f in FOLDS}
    gt_cols = d.allconds.get_indexer(d.conds)
    null_cols = d.allconds.get_indexer(d.nulls)
    fold_pos = {f: np.flatnonzero(d.cf == f) for f in FOLDS}
    rows = []
    for nt in NTOPS:
        for kp in KAPPAS:
            A = np.full((len(d.conds), len(d.kin)), np.nan, np.float32)
            for f in FOLDS:
                Ws = fits[f][kp].matrix(len(d.sites), d.kin, nt)
                rho = SL.rho_sig(SL.sig_project(Ws, d.Y[:, null_cols], d.M[:, null_cols]))
                cols = gt_cols[fold_pos[f]]
                zs = SL.sig_project(Ws, d.Y[:, cols], d.M[:, cols]) / rho
                A[fold_pos[f]] = SL.fuse(zcore[f], zs, alpha)
            S = pd.DataFrame(A, index=d.conds, columns=d.kin)
            P = d.ev.pairs(S)
            rows.append({"n_top": nt, "kappa": kp, "alpha": alpha, "oof_M1-sel": d.ev.metric(P, "M1-sel")["value"],
                         "oof_M1-all": d.ev.metric(P, "M1-all")["value"]})
            print(f"  n_top {nt:3d} kappa {kp}: OOF M1-sel {rows[-1]['oof_M1-sel']:.4f}", flush=True)
    t = pd.DataFrame(rows)
    # reproduction anchor: the B7 per-fold picks at alpha give the B7b per-alpha OOF (results/b7b/grid.tsv)
    picks = {f: (int(c["n_top"]), float(c["kappa"])) for f, c in json.loads(
        (ROOT / "results/b7/b7_eval.json").read_text())["results"][f"B7 {VARIANT} (nested)"]["chosen"].items()}
    A = np.full((len(d.conds), len(d.kin)), np.nan, np.float32)
    for f in FOLDS:
        Ws = fits[f][picks[f][1]].matrix(len(d.sites), d.kin, picks[f][0])
        rho = SL.rho_sig(SL.sig_project(Ws, d.Y[:, null_cols], d.M[:, null_cols]))
        cols = gt_cols[fold_pos[f]]
        A[fold_pos[f]] = SL.fuse(zcore[f], SL.sig_project(Ws, d.Y[:, cols], d.M[:, cols]) / rho, alpha)
    ref = pd.read_csv(ROOT / "results/b7b/grid.tsv", sep="\t").set_index("alpha").loc[alpha, "oof_M1-sel"]
    got = d.ev.metric(d.ev.pairs(pd.DataFrame(A, index=d.conds, columns=d.kin)), "M1-sel")["value"]
    print(f"  anchor: B7b fold picks alpha={alpha}: {got:.6f} vs results/b7b/grid.tsv {ref:.6f}", flush=True)
    assert abs(got - ref) < 1e-6, (got, ref)
    return t


# ----------------------------------------------------------------------------------------------- 2. signatures
def fit_final_signatures(d: DevData, n_top: int, kappa: float):
    v = {"clean+contrast": dict(nt_max=10, contrast=True)}[VARIANT]
    tr = SG.training_conditions(d.C, d.meta, FOLDS, eval_conds=[])
    pos = SG.positives(d.P, d.C, tr, v["nt_max"])
    sig = SG.fit_signatures(d.Yfit, d.allconds, tr, pos, d.meta["group_id"], eval_conds=[], kappa=kappa,
                            n_max=n_top, contrast=v["contrast"])
    rows = []
    for k, s in sig.sites.items():
        w = sig.weights[k]
        for r, (si, wi) in enumerate(zip(s[:n_top], w[:n_top])):
            rows.append((k, d.sites[si], float(wi), r + 1, sig.n_groups[k]))
    tab = pd.DataFrame(rows, columns=["kinase", "site", "weight", "rank", "n_groups"])
    return sig, tab, tr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=str(ROOT / "config/core_v3.json"))
    ap.add_argument("--write-config", action="store_true", help="fill n_top/kappa/rho_sig into the config")
    a = ap.parse_args()
    t0 = time.time()
    cfg_path = Path(a.config)
    cfg = json.loads(cfg_path.read_text())
    sgc = cfg["signatures"]
    alpha = float(sgc["alpha"])
    out = (cfg_path.resolve().parent / cfg["artifacts"]).resolve()
    out.mkdir(parents=True, exist_ok=True)
    d = DevData()
    print(f"[data] {time.time()-t0:.0f}s; {len(d.conds)} DEV-CV GT conds, {len(d.nulls)} nulls, GT_FREEZE "
          f"{d.gt_hash()[:16]}", flush=True)

    # ---- core-v1 part (identical artefacts / rho_e check)
    kin = d.Udf[["kinase", "uniprot", "manning_group", "manning_family"]]
    kin.to_csv(out / "kinases.tsv", sep="\t", index=False)
    edges = build_edges(cfg, list(kin["kinase"]))
    edges.to_csv(out / "edges.tsv.gz", sep="\t", index=False)
    base = core.KinaseActivityModel(cfg, edges, kin)
    Xn = d.X[d.nulls].astype(np.float64)
    po = base.posterior(base.prepare(Xn))
    with np.errstate(all="ignore"):
        zf = po["mu"] / np.sqrt(po["q"]); me = np.where(po["q"] > 0, po["suma"] ** 2 / po["q"], 0.0)
    rho_e, rho_tab = core.fit_rho_e(zf, me, po["report"])
    if abs(rho_e - cfg["rho_e"]) > 1e-12:
        raise SystemExit(f"rho_e refit {rho_e} != config {cfg['rho_e']}")

    # ---- 1. (n_top, kappa) freeze rule
    sel = selection_table(d, alpha)
    sel.to_csv(out / "sig_selection.tsv", sep="\t", index=False, float_format="%.6f")
    frz = json.loads((ROOT / "results/g2_prep/devconfirm_look2_frozen.json").read_text())
    assert frz["name"] == "core-v3" and frz["signature_variant"] == VARIANT and float(frz["alpha"]) == alpha
    assert (frz["core"]["tau"], frz["core"]["rho_f"], frz["core"]["rho_g"], frz["core"]["rho_e"]) == \
        (cfg["tau"], cfg["rho_f"], cfg["rho_g"], cfg["rho_e"])
    n_top, kappa = int(frz["n_top"]), float(frz["kappa"])
    best = sel[(sel["n_top"] == n_top) & (sel["kappa"] == kappa)].iloc[0]
    print(f"[freeze] DEV-confirm look-2 freeze: n_top {n_top} kappa {kappa} alpha {alpha} "
          f"(descriptive late-fusion OOF {best['oof_M1-sel']:.4f})", flush=True)

    # ---- 2. final signatures
    sig, tab, tr = fit_final_signatures(d, n_top, kappa)
    tab.to_csv(out / "signatures.tsv.gz", sep="\t", index=False, float_format="%.9g")
    (out / "signature_training_conditions.txt").write_text("\n".join(tr) + "\n")
    print(f"[signatures] {tab['kinase'].nunique()} kinases, {len(tab)} edges, {len(tr)} training conditions", flush=True)

    # ---- 3. rho_sig on DEV nulls (deployed scorer; written float -> reloaded weights)
    tab_rt = pd.read_csv(out / "signatures.tsv.gz", sep="\t")
    m = core.KinaseActivityModel(cfg, edges, kin).attach_signatures(tab_rt, None, 1.0, alpha,
                                                                    float(sgc.get("clip", 4.0)), int(sgc.get("min_sites", 5)))
    zsn = m.sig_z(m.prepare(Xn))
    rho_sig = SL.rho_sig(zsn)
    # research cross-check (float32 standardised input as in B7b)
    Ws = sig.matrix(len(d.sites), d.kin, n_top)
    nc = d.allconds.get_indexer(d.nulls)
    rho_ref = SL.rho_sig(SL.sig_project(Ws, d.Y[:, nc], d.M[:, nc]))
    print(f"[rho_sig] {rho_sig:.6f} (research path {rho_ref:.6f})", flush=True)
    assert abs(rho_sig - rho_ref) < 1e-4 * rho_ref

    # ---- 4. E1 scales on site-permuted DEV-CV GT conditions
    Xp = ctl.site_permutation(d.X[d.conds], seed=0, mode="within_condition")
    m.attach_signatures(tab_rt, None, rho_sig, alpha, float(sgc.get("clip", 4.0)), int(sgc.get("min_sites", 5)))
    rp = m.score(Xp.astype(np.float64))
    zc_p, zs_p = rp.z_core.to_numpy(), rp.z_sig.to_numpy()
    v, vs = FB.perm_scales(zc_p, zs_p, alpha, m.sig["has"])
    sc = pd.DataFrame({"kinase": d.kin, "v_both": v, "v_sig": vs})
    sc = sc[m.sig["has"]].reset_index(drop=True)
    sc["n_cells_both"] = (np.isfinite(zc_p) & np.isfinite(zs_p)).sum(0)[m.sig["has"]]
    sc["n_cells_sig"] = np.isfinite(zs_p).sum(0)[m.sig["has"]]
    sc["n_groups"] = sc["kinase"].map(tab.groupby("kinase")["n_groups"].first())
    sc.to_csv(out / "e1_scales.tsv", sep="\t", index=False, float_format="%.9g")
    assert abs(rho_sig - frz["rho_sig"]) < 1e-6 * frz["rho_sig"], (rho_sig, frz["rho_sig"])
    assert abs(np.median(sc.v_both) - frz["v_k_median"]) < 1e-4 * frz["v_k_median"], (np.median(sc.v_both), frz["v_k_median"])
    assert len(tr) == frz["n_train_conds"] and tab["kinase"].nunique() == frz["n_signature_kinases"]
    print(f"[E1] v_k median {np.median(sc.v_both):.3f} range {sc.v_both.min():.3f}-{sc.v_both.max():.3f}", flush=True)

    # ---- config: fill or verify
    want = {"n_top": n_top, "kappa": kappa, "rho_sig": float(rho_sig)}
    if a.write_config:
        sgc.update(want)
        cfg_path.write_text(json.dumps(cfg, indent=1) + "\n")
    else:
        bad = {k: (sgc.get(k), v_) for k, v_ in want.items()
               if sgc.get(k) is None or abs(float(sgc[k]) - float(v_)) > 1e-9 * max(1.0, abs(float(v_)))}
        if bad:
            raise SystemExit(f"config/core_v3.json disagrees with the rebuilt artefacts: {bad} (use --write-config "
                             "only on the first build)")

    files = ["edges.tsv.gz", "kinases.tsv", "signatures.tsv.gz", "e1_scales.tsv", "sig_selection.tsv",
             "signature_training_conditions.txt"]
    man = {
        "name": cfg["name"], "built": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "config_sha256": sha_bytes(cfg_path.read_bytes()),
        "artifact_sha256": {f: sha256_file(out / f) for f in files},
        "git_ignored": ["edges.tsv.gz"],
        "inputs_sha256": {
            "GT_FREEZE.json": d.gt_hash(),
            "gt_pairs.tsv": sha256_file(ROOT / "results/ground_truth/gt_pairs.tsv"),
            "gt_conditions.tsv": sha256_file(ROOT / "results/ground_truth/gt_conditions.tsv"),
            "kinase_universe.tsv": sha256_file(ROOT / "results/ground_truth/kinase_universe.tsv"),
            "quant_dev.npy": sha256_file(ROOT / "cache/dev/quant_dev.npy"),
            **{f"lib:{n}": sha256_file(LIBDIR / f"{n}.parquet")
               for n in [cfg["library"]["curated"], *cfg["library"].get("predicted", {})]}},
        "n_edges": int(len(edges)), "n_edges_by_source": edges["source"].value_counts().to_dict(),
        "rho_e_refit": float(rho_e), "rho_e_fit_conditions": d.nulls, "rho_e_fit_table": rho_tab,
        "signatures": {"variant": VARIANT, "n_top": n_top, "kappa": kappa, "alpha": alpha,
                       "n_kinases": int(tab["kinase"].nunique()), "n_edges": int(len(tab)),
                       "n_training_conditions": len(tr), "training_splits": FOLDS,
                       "selection_rule": "taken from results/g2_prep/devconfirm_look2_frozen.json (look-2 freeze, "
                                         f"config_hash {frz['config_hash']}); sig_selection.tsv descriptive only",
                       "devconfirm_look2_frozen_sha256": sha256_file(ROOT / "results/g2_prep/devconfirm_look2_frozen.json"),
                       "late_fusion_oof_M1-sel_descriptive": float(best["oof_M1-sel"])},
        "rho_sig_refit": float(rho_sig), "rho_sig_research_path": float(rho_ref),
        "e1": {"n_perm_conditions": int(Xp.shape[1]), "permutation": "site_permutation seed 0 within_condition, "
               "DEV-CV GT conditions (Evaluator folds 0-4 order)", "min_perm_cells": FB.MIN_PERM_CELLS,
               "v_both_median": float(np.median(sc.v_both)), "v_both_min": float(sc.v_both.min()),
               "v_both_max": float(sc.v_both.max())},
        "data": "cache/dev (DEV-CV folds 0-4), results/ground_truth DEV rows, cache/libraries, cache/b7 (selection)",
    }
    (out / "manifest.json").write_text(json.dumps(man, indent=1, default=float) + "\n")
    print(f"wrote {out} ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
