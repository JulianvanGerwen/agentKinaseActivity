#!/usr/bin/env python
"""THE single sealed-TEST run (DESIGN §5, §8, §12 C7-C10). To be launched ONCE, by the PI, after G2.

    /home/plm_env/bin/python scripts/run_test.py 2>&1 | tee results/test/run.log      # real TEST run
    /home/plm_env/bin/python scripts/run_test.py --dry-run-dev                         # smoke test (no TEST data)
    /home/plm_env/bin/python scripts/run_test.py --print-freeze-template               # hashes the PI puts in FREEZE.json

Refuses to run unless /workspace/FREEZE.json exists and ALL of these match the working tree:
    gt_hash               == kinact.eval.sealed_gt_hash()          (also enforced by load_sealed_test)
    gt_freeze_sha256      == sha256(results/ground_truth/GT_FREEZE.json)
    config_sha256         == sha256(config/core_v3.json)
    model_manifest_sha256 == sha256(models/core_v3/manifest.json)  (+ every artefact sha256 in the manifest, incl. the
                                                                     git-ignored edges.tsv.gz, is re-checked on disk)
    code_commit           == `git rev-parse HEAD`, and no tracked file under src/ scripts/ config/ models/ tests/ is
                             modified (git status --porcelain --untracked-files=no)
FREEZE.json is written by the PI only; this script never writes it.

What is computed (all fixed at G2; nothing is chosen on TEST):
  final (PRIMARY)   core-v3 = KinaseActivityModel.from_config('config/core_v3.json')
  Holm comparators  z-score curated_union min1 (DESIGN B14: better DEV config than the PSP min3 default)
                    KSEA curated_union min1
                    JUMPsem curated_union min3, 'aligned' variant (best DEV config), native R 1.1.0 per study unit
                      (PMID x Sample x Control, pooled per PMID if < 4) over the TEST conditions only
                    best-other-baseline = PhosX-like Cantley site-top15 min3 (DEV-best excluding z/KSEA/JUMPsem)
                    prior-only = s_c * freq_k (freq over all DEV-CV GT conditions; oracle condition sign)
  controls          shift-prior (freq_k * mean_c x), random N(0,1) seed 0; PSP z-score min3 and JUMPsem native
                    (descriptive)
  secondary         core-v1 (config/core_v1.json); A1 = avg(core-v1 z_vif, B4 co-regulation wz) (kinact.secondary;
                    co-regulation fitted on DEV-CV conditions EXCLUDING any condition that shares a control key
                    (PMID, Control, Sample) with a TEST condition — metadata only)
  metrics           frozen core: M1-sel (primary), M1-all, M1-sel-scored + coverage, M2 (vs each comparator), harmony
                    (+ per bin), by source / sign, T-seen / T-unseen / T-family-unseen, coverage-significant kinases,
                    calibration (|z| > 1.96, overall / per bin / signature kinases) on the N-eval TEST nulls and on
                    site-permuted TEST GT conditions, runtime; has-/no-signature split (C10 iv).
  inference         paired cluster bootstrap over GT groups (B = 2000, seed 0, shared resamples, C8 basic CI +
                    shifted-null one-sided p). TEST Holm family (6): final vs z-score, KSEA, JUMPsem,
                    best-other-baseline, prior-only on M1-sel; final vs z-score on M2. "Beats" = Holm p < 0.05 AND
                    dM1-sel >= 0.01 AND dM1-sel(T-unseen) > 0 AND point dM2 >= 0. MDE reference = results/G1.json
                    (pre-registered, C7); the realised 2.5 SE on TEST is reported descriptively.
  external          kinact.external_eval.run_external(..., confirm_test=True) for core-v3 / core-v1 / z-score / KSEA /
                    PhosX-like / PSP z / random, + paired bootstrap per set (perturbation: GT groups = PMID x drug;
                    tumour: kinases), final vs each comparator.
Outputs: results/test/{LEDGER.md, test_results.json, per_kinase/*.tsv, scores/*.parquet, freeze_check.json}.

--dry-run-dev: substitutes the DEV-confirm split for TEST (labels from results/ground_truth, nulls = DEV-confirm
role=null conditions; A1 control exclusion against DEV-confirm conditions; T-unseen vs DEV-CV positives), skips the
FREEZE check (hashes are still logged) and the external sets, and writes to results/test_dryrun/. Its numbers are a
SMOKE TEST ONLY: they are not a DEV-confirm look and must not be logged or used for any decision.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from kinact import bootstrap as bs, controls as ctl  # noqa: E402
from kinact.core import KinaseActivityModel  # noqa: E402
from kinact.data import SITES_PQ, load_dev, load_dev_meta  # noqa: E402
from kinact.eval import (Evaluator, calibration, load_gt, psp_bin, read_gt_tsv, sealed_gt_hash,  # noqa: E402
                         to_jsonable)
from kinact.methods import ksea, phosx_like, run_jumpsem, study_units, zscore  # noqa: E402

FOLDS = [f"fold{i}" for i in range(5)]
FREEZE = ROOT / "FREEZE.json"
CFG3 = ROOT / "config/core_v3.json"
CFG1 = ROOT / "config/core_v1.json"
MANIFEST = ROOT / "models/core_v3/manifest.json"
B, SEED = 2000, 0
JOBS = 32

FINAL = "core-v3"
HOLM_COMPARATORS = [("z-score", "zscore curated_union min1"), ("KSEA", "KSEA curated_union min1"),
                    ("JUMPsem", "JUMPsem curated_union min3 (aligned)"),
                    ("best-other-baseline", "PhosX-like cantley site-top15 min3"), ("prior-only", "prior-only")]
DESCRIPTIVE = ["core-v1", "A1 (core-v1 + coreg wz)", "zscore psp_all min3", "JUMPsem curated_union min3 (native)",
               "shift-prior", "random"]


# =============================================================================================== provenance
def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while b := f.read(1 << 24):
            h.update(b)
    return h.hexdigest()


def git(*args) -> str:
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, check=True).stdout.strip()


def expected_hashes(include_sealed: bool) -> dict:
    out = {"gt_freeze_sha256": sha256(ROOT / "results/ground_truth/GT_FREEZE.json"),
           "config_sha256": sha256(CFG3), "model_manifest_sha256": sha256(MANIFEST),
           "code_commit": git("rev-parse", "HEAD")}
    if include_sealed:
        out = {"gt_hash": sealed_gt_hash(), **out}
    return out


def dirty_tracked() -> list:
    s = git("status", "--porcelain", "--untracked-files=no", "--", "src", "scripts", "config", "models", "tests")
    return [l for l in s.splitlines() if l.strip()]


def check_manifest_artifacts() -> dict:
    man = json.loads(MANIFEST.read_text())
    bad = {f: h for f, h in man["artifact_sha256"].items() if sha256(MANIFEST.parent / f) != h}
    return {"n_checked": len(man["artifact_sha256"]), "mismatch": bad}


def verify_freeze(F: dict, want: dict, dirty: list, art: dict) -> None:
    """Raise SystemExit unless every expected hash is present in FREEZE.json and equal, the tracked tree is clean
    and every model artefact matches its manifest hash."""
    mism = {k: {"FREEZE": F.get(k), "working_tree": v} for k, v in want.items() if F.get(k) != v}
    if mism or dirty or art["mismatch"]:
        raise SystemExit("REFUSED: FREEZE.json does not match the working tree:\n"
                         + json.dumps({"hash_mismatch": mism, "modified_tracked_files": dirty,
                                       "artifact_mismatch": art["mismatch"]}, indent=1))


def check_freeze() -> dict:
    env = os.environ.get("KINACT_FREEZE_PATH")
    if env and Path(env).resolve() != FREEZE.resolve():
        raise SystemExit(f"REFUSED: KINACT_FREEZE_PATH={env} points away from {FREEZE}")
    if not FREEZE.exists():
        raise SystemExit(f"REFUSED: {FREEZE} does not exist (written by the PI at G2 only)")
    F = json.loads(FREEZE.read_text())
    want = expected_hashes(include_sealed=True)
    art = check_manifest_artifacts()
    verify_freeze(F, want, dirty_tracked(), art)
    return {"FREEZE.json": F, "FREEZE_sha256": sha256(FREEZE), "verified": want, "artifacts": art}


def versions() -> dict:
    import scipy
    import threadpoolctl
    v = {"python": sys.version.split()[0], "platform": platform.platform(), "numpy": np.__version__,
         "pandas": pd.__version__, "scipy": scipy.__version__, "threadpoolctl": threadpoolctl.__version__}
    try:
        import torch
        v["torch"] = torch.__version__
        v["cuda"] = torch.cuda.get_device_name(0) if torch.cuda.is_available() else None
    except ImportError:                                                                  # pragma: no cover
        v["torch"] = None
    from kinact.methods.jumpsem import RSCRIPT
    try:
        r = subprocess.run([str(RSCRIPT), "-e", 'cat(R.version.string, as.character(packageVersion("JUMPsem")), '
                            'as.character(packageVersion("lavaan")), sep="|")'], capture_output=True, text=True,
                           timeout=120, env=dict(os.environ, LC_ALL="C.UTF-8"))
        v["R|JUMPsem|lavaan"] = r.stdout.strip()
    except Exception as e:                                                               # pragma: no cover
        v["R|JUMPsem|lavaan"] = f"unavailable: {e}"
    return v


# =============================================================================================== data
def load_columns(conds) -> pd.DataFrame:
    """Sites x conditions (log2FC float64) for the requested compendium columns only (memmap column subset);
    peptides collapsed to site keys exactly as kinact.data.load_dev (row set == the DEV cache)."""
    import scipy.sparse as sp
    q = np.load(ROOT / "cache/quant_wide_f32.npy", mmap_mode="r")
    peps = (ROOT / "cache/quant_wide_sites.txt").read_text().split()
    allc = (ROOT / "cache/quant_wide_conditions.txt").read_text().split()
    assert peps == (ROOT / "cache/dev/sites_dev.txt").read_text().split()
    ci = pd.Index(allc).get_indexer(list(conds))
    if (ci < 0).any():
        raise KeyError(f"{int((ci < 0).sum())} conditions missing from the compendium cache")
    M = np.asarray(q[:, np.sort(ci)], dtype=np.float32)[:, np.argsort(np.argsort(ci))]
    sites = pd.read_parquet(SITES_PQ, columns=["peptide_UID", "site_key"]).set_index("peptide_UID")["site_key"]
    codes, uniq = pd.factorize(sites.reindex(peps).to_numpy())
    ok = ~np.isnan(M)
    S = sp.csr_matrix((np.ones(len(codes), np.float32), (codes, np.arange(len(codes)))), shape=(len(uniq), len(codes)))
    num = S @ np.where(ok, M, 0.0); cnt = S @ ok.astype(np.float32)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.where(cnt > 0, num / cnt, np.nan)
    return pd.DataFrame(out.astype(np.float32), index=uniq, columns=list(conds)).astype(np.float64)


class Split:
    """Evaluated split: TEST (sealed labels) or, in --dry-run-dev, DEV-confirm (visible labels)."""

    def __init__(self, dry: bool):
        self.dry = dry
        P, C, U, _ = load_gt()
        self.Udf = read_gt_tsv(ROOT / "results/ground_truth/kinase_universe.tsv")
        self.kin = list(self.Udf["kinase"])
        self.P_vis, self.C_vis = P, C
        cg_vis = read_gt_tsv(ROOT / "results/ground_truth/condition_groups.tsv")
        if dry:
            self.name = "devconfirm"
            self.P = P[P["conditionID"].isin(set(C.loc[C["split"] == "devconfirm", "conditionID"]))]
            self.C = C[C["split"] == "devconfirm"]
            self.bins = read_gt_tsv(ROOT / "results/ground_truth/harmony_bins_group.tsv")
            self.cg = cg_vis[cg_vis["split"] == "devconfirm"]
            dev_splits = FOLDS
            self.freeze = None
        else:
            from kinact.eval import load_sealed_test
            self.name = "test"
            s = load_sealed_test(freeze=FREEZE)
            self.P, self.C, self.bins = s["gt_pairs"], s["gt_conditions"], s["harmony_bins_group"]
            self.cg = s["condition_groups"]
            self.cg = self.cg[self.cg["split"].astype(str) == "test"]
            self.C = self.C[self.C["split"].astype(str) == "test"]
            dev_splits = FOLDS + ["devconfirm"]                      # 'DEV' = every non-TEST split (make_split.py)
            self.freeze = FREEZE
        self.conds = list(self.cg["conditionID"])
        self.gt_conds = list(self.C.loc[self.C["role"].astype(str) == "gt", "conditionID"])
        self.nulls = list(self.C.loc[self.C["role"].astype(str) == "null", "conditionID"])
        assert set(self.gt_conds) | set(self.nulls) <= set(self.conds)
        dev_meta = load_dev_meta()
        assert not set(self.conds) & set(dev_meta.index), "evaluated conditions found in the DEV cache"
        assert not set(self.cg["group_id"]) & set(dev_meta["group_id"]), "evaluated groups overlap DEV-CV groups"
        # DEV positives (T-unseen / T-family-unseen) and DEV-CV prior frequencies
        dgt = C[C["split"].isin(dev_splits) & (C["role"] == "gt")]
        dpos = P[(P["cls"] == "direct") & P["confidence"].isin(["high", "medium"]) & P["conditionID"].isin(set(dgt["conditionID"]))]
        self.dev_pos = set(dpos["kinase"])
        fam = dict(zip(self.Udf["kinase"], [f if isinstance(f, str) and f else k
                                            for k, f in zip(self.Udf["kinase"], self.Udf["manning_family"])]))
        self.fam = fam
        self.dev_fam = {fam.get(k, k) for k in self.dev_pos}
        self.dev_cv_train = list(C.loc[C["split"].isin(FOLDS), "conditionID"])
        self.ev = Evaluator(self.P, self.C, self.Udf, split_filter=[self.name], bins=self.bins,
                            train_pos_kinases=self.dev_pos, zero_policy="ignore", freeze=self.freeze)
        allmeta = pd.read_parquet(ROOT / "cache/conditions.parquet").set_index("conditionID", drop=False)
        qc = read_gt_tsv(ROOT / "results/ground_truth/condition_qc.tsv").set_index("conditionID")["qc_status"]
        self.meta = allmeta.loc[self.conds].assign(qc_status=qc.reindex(self.conds).to_numpy())
        self.allmeta = allmeta


# =============================================================================================== scoring
def zframe(S, kin, index):
    return S.reindex(index=index, columns=kin)


def score_all(sp_: Split, X: pd.DataFrame, log) -> tuple[dict, dict]:
    kin, conds = sp_.kin, list(X.columns)
    m3 = KinaseActivityModel.from_config(CFG3)
    m1 = KinaseActivityModel.from_config(CFG1)
    S, aux = {}, {}
    t = time.time(); R3 = m3.score(X)
    S[FINAL] = R3.z; aux["core-v3"] = R3; log(f"core-v3 {time.time()-t:.1f}s")
    R1 = m1.score(X); S["core-v1"] = R1.z; aux["core-v1"] = R1
    assert np.array_equal(np.nan_to_num(R1.z.to_numpy(), nan=9e9), np.nan_to_num(R3.z_core.to_numpy(), nan=9e9))
    for name, fn in [("zscore curated_union min1", lambda: zscore(X, "curated_union", 1, return_n=True)),
                     ("zscore psp_all min3", lambda: zscore(X, "psp_all", 3, return_n=True)),
                     ("KSEA curated_union min1", lambda: ksea(X, "curated_union", 1, return_n=True)),
                     ("PhosX-like cantley site-top15 min3", lambda: phosx_like(X, 15, min_substrates=3, return_n=True,
                                                                               n_jobs=JOBS))]:
        t = time.time(); Z, N = fn()
        S[name] = zframe(Z, kin, conds); aux[name] = zframe(N, kin, conds); log(f"{name} {time.time()-t:.1f}s")
    # JUMPsem: native R per study unit over the evaluated split's qc-ok conditions only
    t = time.time()
    ok = list(sp_.meta.index[sp_.meta["qc_status"] == "ok"])
    units = study_units(sp_.meta, ok)
    r = run_jumpsem(X[ok], "curated_union", units=units, min_substrates=3, n_jobs=48, log=lambda *a: None)
    S["JUMPsem curated_union min3 (aligned)"] = zframe(r["aligned"], kin, conds)
    S["JUMPsem curated_union min3 (native)"] = zframe(r["native"], kin, conds)
    aux["jumpsem_units"] = r["units"]
    log(f"JUMPsem {time.time()-t:.1f}s; {len(r['units'])} units, {int((r['units']['error'] != '').sum())} failed, "
        f"{int(units.notna().sum())} conditions in units")
    # controls
    freq = ctl.training_frequency(sp_.P_vis, sp_.C_vis, sp_.dev_cv_train, kin)
    s_c = ctl.condition_sign(sp_.P, conds)
    S["prior-only"] = pd.DataFrame(np.outer(s_c.to_numpy(), freq.to_numpy()), index=conds, columns=kin)
    S["shift-prior"] = zframe(ctl.shift_prior(X, freq), kin, conds)
    S["random"] = ctl.random_scores(conds, kin, seed=SEED)
    aux["freq"] = freq
    # A1 (secondary): co-regulation fitted on DEV-CV minus conditions sharing a control with an evaluated condition
    t = time.time()
    S["A1 (core-v1 + coreg wz)"], aux["A1"] = score_a1(sp_, X, m1, log)
    log(f"A1 {time.time()-t:.1f}s")
    return S, aux


def dev_permuted():
    """Site-permuted DEV-CV GT conditions (seed 0) — the label-free calibration input of B6/B7b/B8."""
    P, C, U, _ = load_gt()
    bins = read_gt_tsv(ROOT / "results/ground_truth/harmony_bins_group.tsv")
    gconds = list(Evaluator(P, C, U, split_filter=FOLDS, bins=bins, zero_policy="ignore").cond_index)
    Xd = load_dev()
    return Xd, ctl.site_permutation(Xd[gconds], seed=0, mode="within_condition")


def score_a1(sp_: Split, X: pd.DataFrame, m1, log):
    from kinact import secondary as S2
    Xd, Xp = dev_permuted()
    dmeta = load_dev_meta()
    Xs = Xd[Xd.notna().any(axis=1)]
    m = dmeta.reindex(Xs.columns)
    eval_keys = set(S2.control_keys(sp_.allmeta.loc[sp_.conds]))
    shares = S2.control_keys(m).isin(eval_keys).to_numpy()
    allowed = (m["split"].isin(FOLDS) & (m["qc_status"] == "ok")).to_numpy() & ~shares
    cw = S2.CoregWZ(Xd, allowed, sp_.kin)
    a1 = S2.A1(m1.score(Xp.astype(np.float64)).z.to_numpy(), cw.score(Xp).to_numpy())
    info = {"n_allowed": int(allowed.sum()), "n_dev_sharing_control_with_eval": int(shares.sum()),
            "n_expansion_edges": int(len(cw.exp)), "sigma": list(a1.s), "R": a1.R.tolist(), "cfg": cw.cfg}
    log(f"A1 coreg: allowed {info['n_allowed']} DEV-CV conditions (excluded {info['n_dev_sharing_control_with_eval']} "
        "sharing a control with evaluated conditions)")
    zv = m1.score(X).z
    return a1.combine(zv, cw.score(X)), info


# =============================================================================================== evaluation
def masks(sp_: Split):
    uns = lambda P: ~P["kinase"].isin(sp_.dev_pos).to_numpy()
    fun = lambda P: ~P["kinase"].map(lambda k: sp_.fam.get(k, k)).isin(sp_.dev_fam).to_numpy()
    sig = set(pd.read_csv(ROOT / "models/core_v3/signatures.tsv.gz", sep="\t")["kinase"])
    has = lambda P: P["kinase"].isin(sig).to_numpy()
    return {"T-seen": lambda P: ~uns(P), "T-unseen": uns, "T-family-unseen": fun,
            "has-signature": has, "no-signature": lambda P: ~has(P)}


def delta(ev, A, Bm, ep, plan, mask=None) -> dict:
    keys = ["estimate", "se", "ci_low", "ci_high", "p_one_sided"]
    if mask is None:
        d = bs.paired_delta(ev, A, Bm, ep, plan)
        return {k: d[k] for k in keys}
    Pa, Pb = ev.pairs(A), ev.pairs(Bm)
    ma, mb = mask(Pa), mask(Pb)
    assert (ma == mb).all()
    ra = bs.metric_replicates(ev, Pa[ma], ep, plan); rb = bs.metric_replicates(ev, Pb[mb], ep, plan)
    pt = ev.metric(Pa[ma], ep)["value"] - ev.metric(Pb[mb], ep)["value"]
    s = bs._summ(pt, ra - rb)
    return {"estimate": pt, "se": s["se"], "ci_low": s["ci_low"], "ci_high": s["ci_high"],
            "p_one_sided": bs.shifted_null_p(pt, ra - rb), "n_pairs_M1sel": int((ma & Pa["in_M1sel"]).sum())}


def method_metrics(sp_: Split, S: pd.DataFrame, mk: dict, n_perm=2000) -> tuple[dict, pd.DataFrame]:
    ev = sp_.ev
    r = ev.evaluate(S, n_perm=n_perm)
    v = lambda x: x["value"] if isinstance(x, dict) and "value" in x else x
    P = r["pairs"]
    out = {"M1-sel": v(r["M1-sel"]), "M1-all": v(r["M1-all"]), "M1-sel-scored": v(r["M1-sel-scored"]),
           "coverage_M1-sel": r["coverage_M1-sel"], "coverage_M1-all": r["coverage_M1-all"],
           "harmony": v(r.get("harmony")), "harmony_per_bin": (r.get("harmony") or {}).get("per_bin"),
           "by_source": {k: v(x) for k, x in r["by_source"].items()},
           "by_sign": {str(k): v(x) for k, x in r["by_sign"].items()},
           "coverage_significant": r["coverage_significant"],
           "n_M1sel_pairs": r["M1-sel"]["n_pairs"], "n_M1sel_kinases": r["M1-sel"]["n_kinases"],
           "scored_cell_fraction": float(S.reindex(index=ev.cond_index).notna().to_numpy().mean())}
    for name, f in mk.items():
        m = f(P)
        out[name] = v(ev.metric(P[m], "M1-sel")); out[f"{name}_n_pairs"] = int((m & P["in_M1sel"]).sum())
    return out, r["per_kinase"]


def calib_block(Z: pd.DataFrame, N: pd.DataFrame | None, conds, sig_cols=None) -> dict:
    conds = [c for c in conds if c in Z.index]
    if not conds:
        return {"n_conditions": 0}
    c = calibration(Z, conds, N)
    out = {k: v for k, v in c.items() if k != "per_kinase"}
    if sig_cols is not None:
        z = Z.loc[conds, sig_cols].to_numpy(float)
        ok = np.isfinite(z)
        out["signature_kinases"] = {"frac_exceed": float((np.abs(z[ok]) > 1.96).mean()) if ok.any() else None,
                                    "n_cells": int(ok.sum())}
    return out


def evaluate_all(sp_: Split, S: dict, aux: dict, X: pd.DataFrame, log) -> dict:
    ev, kin = sp_.ev, sp_.kin
    plan = bs.BootstrapPlan.from_evaluator(ev, B=B, seed=SEED)
    mk = masks(sp_)
    res = {"n": {"conditions": len(sp_.conds), "gt_conditions_evaluable": len(ev.cond_index), "nulls": len(sp_.nulls),
                 "M1sel_pairs": int(ev.targets["in_M1sel"].sum()), "groups": len(ev.group_index),
                 "kinases_M1sel": int(ev.targets.loc[ev.targets["in_M1sel"], "kinase"].nunique()),
                 "T-unseen_kinases_M1sel": len(set(ev.targets.loc[ev.targets["in_M1sel"], "kinase"]) - sp_.dev_pos)},
           "methods": {}, "per_kinase": {}}
    for name, Sm in S.items():
        t = time.time()
        res["methods"][name], res["per_kinase"][name] = method_metrics(sp_, Sm, mk)
        log(f"metrics {name}: M1-sel {res['methods'][name]['M1-sel']:.4f} ({time.time()-t:.0f}s)")
    F = S[FINAL]
    # ---------------- Holm family (6) + beats
    fam = []
    for lab, comp in HOLM_COMPARATORS:
        fam.append({"test": f"{FINAL} vs {lab} (M1-sel)", "comparator": comp, "endpoint": "M1-sel",
                    **delta(ev, F, S[comp], "M1-sel", plan)})
    fam.append({"test": f"{FINAL} vs z-score (M2)", "comparator": HOLM_COMPARATORS[0][1], "endpoint": "M2",
                **delta(ev, F, S[HOLM_COMPARATORS[0][1]], "M2", plan)})
    hp = bs.holm([x["p_one_sided"] for x in fam])
    for x, h in zip(fam, hp):
        x["p_holm"] = float(h)
    beats = {}
    for (lab, comp), x in zip(HOLM_COMPARATORS, fam[:5]):
        m2 = ev.m2(F, S[comp])
        du = delta(ev, F, S[comp], "M1-sel", plan, mk["T-unseen"])
        beats[lab] = {"comparator": comp, "beats": bs.beats(x, x["p_holm"], m2["delta"], du["estimate"]),
                      "criteria": {"holm_p<0.05": x["p_holm"] < 0.05, "dM1-sel>=0.01": x["estimate"] >= 0.01,
                                   "dM1-sel(T-unseen)>0": du["estimate"] > 0, "point dM2>=0": m2["delta"] >= 0},
                      "dM1-sel": x["estimate"], "p_holm": x["p_holm"], "dM2_point": m2["delta"], "M2": m2,
                      "dT-unseen": du}
    res["holm_family"], res["beats"] = fam, beats
    # ---------------- descriptive comparisons (not in the family)
    desc = {}
    pairs_ = [(FINAL, c) for _, c in HOLM_COMPARATORS] + [(FINAL, d) for d in DESCRIPTIVE] + \
             [("A1 (core-v1 + coreg wz)", "core-v1"), ("core-v1", HOLM_COMPARATORS[0][1]),
              (HOLM_COMPARATORS[0][1], "prior-only")]
    for a, b_ in pairs_:
        d = {ep: delta(ev, S[a], S[b_], ep, plan) for ep in ["M1-sel", "M1-all", "harmony", "M2"]}
        d["M2_detail"] = ev.m2(S[a], S[b_])
        for tn in ["T-seen", "T-unseen", "T-family-unseen"]:
            d[tn] = delta(ev, S[a], S[b_], "M1-sel", plan, mk[tn])
        if (a, b_) == (FINAL, "core-v1"):
            for tn in ["has-signature", "no-signature"]:
                d[tn] = delta(ev, S[a], S[b_], "M1-sel", plan, mk[tn])
            for src in sorted(ev.targets["source"].unique()):
                d[f"source={src}"] = delta(ev, S[a], S[b_], f"M1-sel:source={src}", plan)
            for bn in ["0", "1-4", "5-19", ">=20"]:
                if (ev.targets["bin"] == bn).any():
                    d[f"bin={bn}"] = delta(ev, S[a], S[b_], f"M1-sel:bin={bn}", plan)
        desc[f"{a} vs {b_}"] = d
    res["comparisons"] = desc
    # ---------------- MDE reference (C7) + realised SE
    g1 = json.loads((ROOT / "results/G1.json").read_text())
    zp = desc[f"{HOLM_COMPARATORS[0][1]} vs prior-only"]["M1-sel"]
    res["mde"] = {"preregistered_mde_m1sel (results/G1.json, C7)": g1.get("mde_m1sel"),
                  "realised_2.5se_zscore_vs_prior_only (descriptive)": 2.5 * zp["se"],
                  "realised_2.5se_final_vs_zscore (descriptive)": 2.5 * fam[0]["se"]}
    # ---------------- calibration: N-eval nulls + site-permuted evaluated GT conditions
    sig_cols = [k for k in kin if k in set(pd.read_csv(ROOT / "models/core_v3/signatures.tsv.gz", sep="\t")["kinase"])]
    N3 = aux["core-v3"].n_substrates
    cal = {"nulls": sp_.nulls}
    cal["N-eval nulls"] = {FINAL: calib_block(S[FINAL], N3, sp_.nulls, sig_cols),
                           "core-v1": calib_block(S["core-v1"], N3, sp_.nulls, sig_cols),
                           "A1 (core-v1 + coreg wz)": calib_block(S["A1 (core-v1 + coreg wz)"], N3, sp_.nulls, sig_cols),
                           "zscore curated_union min1": calib_block(S["zscore curated_union min1"],
                                                                    aux["zscore curated_union min1"], sp_.nulls),
                           "KSEA curated_union min1": calib_block(S["KSEA curated_union min1"],
                                                                  aux["KSEA curated_union min1"], sp_.nulls)}
    gt = [c for c in sp_.gt_conds if c in X.columns]
    Xp = ctl.site_permutation(X[gt], seed=SEED, mode="within_condition")
    m3 = KinaseActivityModel.from_config(CFG3)
    Rp = m3.score(Xp)
    cal["site-permuted GT"] = {FINAL: calib_block(Rp.z, Rp.n_substrates, gt, sig_cols),
                               "core-v1": calib_block(Rp.z_core, Rp.n_substrates, gt, sig_cols)}
    res["calibration"] = cal
    return res


# =============================================================================================== runtime
def runtime_block(sp_: Split, log) -> dict:
    from threadpoolctl import threadpool_limits

    from kinact.runtime import bench_runtime
    m3 = KinaseActivityModel.from_config(CFG3)
    Xd = load_dev().astype(np.float64)
    rng = np.random.default_rng(0)
    bench = sorted(rng.choice(sorted(Xd.columns), 200, replace=False))

    class GPU:
        def prepare(self, X, lib=None):
            return m3.prepare(X)

        def score(self, st):
            return m3.score(st, device="cuda")
    out = {"bench": "kinact.runtime.bench_runtime, 200 fixed DEV conditions (seed 0), warm, median of 5",
           "cpu": bench_runtime(m3, Xd[bench], threads=(1, 16))}
    try:
        import torch
        if torch.cuda.is_available():
            out["gpu"] = bench_runtime(GPU(), Xd[bench], threads=(16,))
    except ImportError:                                                                  # pragma: no cover
        pass
    allc = (ROOT / "cache/quant_wide_conditions.txt").read_text().split()
    if sp_.dry:     # never touch TEST data in the smoke test: non-TEST compendium conditions only
        allc = [c for c in allc if c in set(read_gt_tsv(ROOT / "results/ground_truth/condition_groups.tsv")["conditionID"])]
    Xa = load_columns(allc)
    comp = {}
    for nt in (1, 16):
        with threadpool_limits(nt):
            ts = []
            for _ in range(3):
                t = time.perf_counter(); m3.score(m3.prepare(Xa)); ts.append(time.perf_counter() - t)
            comp[f"cpu_{nt}t_total_s"] = float(np.median(ts))
    comp["n_conditions"] = len(allc)
    out["compendium"] = comp
    log(f"runtime: compendium ({len(allc)} conds) 1t {comp['cpu_1t_total_s']:.1f}s, 16t {comp['cpu_16t_total_s']:.1f}s")
    return out


# =============================================================================================== external
def external_block(log) -> dict:
    from kinact import external_eval as XE
    m3 = KinaseActivityModel.from_config(CFG3)
    m1 = KinaseActivityModel.from_config(CFG1)
    U = XE.load_universe()
    fns = {FINAL: lambda d: m3.score(d.astype(np.float64)).z,
           "core-v1": lambda d: m1.score(d.astype(np.float64)).z,
           "zscore curated_union min1": lambda d: zscore(d.astype(np.float64), "curated_union", 1),
           "KSEA curated_union min1": lambda d: ksea(d.astype(np.float64), "curated_union", 1),
           "PhosX-like cantley site-top15 min3": lambda d: phosx_like(d.astype(np.float64), 15, min_substrates=3,
                                                                      n_jobs=JOBS),
           "zscore psp_all min3": lambda d: zscore(d.astype(np.float64), "psp_all", 3)}
    rec, res = {}, {}
    for name, fn in fns.items():
        store = []

        def wrapped(d, _fn=fn, _st=store):
            s = _fn(d).reindex(columns=U)
            _st.append(s)
            return s
        t = time.time()
        res[name] = XE.run_external(wrapped, confirm_test=True)
        rec[name] = store                       # order: PERT_SETS..., then COHORTS
        log(f"external {name} {time.time()-t:.0f}s")
    res["random"] = XE.run_external(lambda d: XE.random_scores(d.columns, U, seed=SEED), confirm_test=True)
    # ---------------- paired bootstrap per set: final vs each comparator
    paired = {}
    for i, sname in enumerate(XE.PERT_SETS):
        _, exps, tg = XE.load_pert(sname)
        gp, gc = XE.pert_gt_tables(exps, tg)
        ev = Evaluator(gp, gc, U, zero_policy="ignore")
        plan = bs.BootstrapPlan.from_evaluator(ev, B=B, seed=SEED)
        A = rec[FINAL][i]
        paired[sname] = {c: {ep: delta(ev, A, rec[c][i], ep, plan) for ep in ("M1-sel", "M1-all")}
                         for c in fns if c != FINAL}
    _, gs = XE.load_tumour()
    nP = len(XE.PERT_SETS)

    def km(stores):
        sc = {c: stores[nP + j] for j, c in enumerate(XE.COHORTS)}
        return XE.m1_tumour(sc, gs)["per_kinase"]
    kA = km(rec[FINAL])
    rng = np.random.default_rng(SEED)
    idx = rng.integers(0, len(kA), size=(B, len(kA)))
    paired["tumour"] = {}
    for c in fns:
        if c == FINAL:
            continue
        kB = km(rec[c]).reindex(kA.index)
        dk = (kA - kB).to_numpy()
        reps = dk[idx].mean(1)
        s = bs._summ(float(dk.mean()), reps)
        paired["tumour"][c] = {"M1-kinase": {"estimate": s["estimate"], "se": s["se"], "ci_low": s["ci_low"],
                                             "ci_high": s["ci_high"], "p_one_sided": bs.shifted_null_p(s["estimate"], reps),
                                             "n_kinases": int(len(dk)), "unit": "kinase (cluster bootstrap)"}}
    return {"results": res, "paired_vs_final": paired}


# =============================================================================================== report
def f3(x):
    return "NA" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.3f}"


def fd(d):
    return f"{d['estimate']:+.3f} [{d['ci_low']:+.3f}, {d['ci_high']:+.3f}]"


def write_ledger(out: Path, R: dict, dry: bool):
    L = [f"# {'DRY RUN (DEV-confirm substituted for TEST) — SMOKE TEST ONLY, NOT A RESULT' if dry else 'Sealed TEST run'}",
         "", f"Generated by `scripts/run_test.py` at {R['provenance']['finished']}. Wall time {R['provenance']['wall_s']:.0f} s.",
         f"Evaluated split: `{R['split']}`. Final (primary): **{FINAL}** (`config/core_v3.json`).", ""]
    if dry:
        L += ["> These numbers are a smoke test of the TEST pipeline on DEV-confirm. They are NOT a DEV-confirm look,",
              "> are not logged in results/devconfirm_looks.tsv or the ledger, and must not be used for any decision.", ""]
    pv = R["provenance"]
    L += ["## Provenance", "", "| item | value |", "|---|---|"]
    for k, v in {**pv["hashes"], **{f"version:{a}": b for a, b in pv["versions"].items()}}.items():
        L.append(f"| {k} | `{v}` |")
    L.append(f"| FREEZE check | {'skipped (dry run)' if dry else 'passed'} |")
    n = R["n"]
    L += ["", "## Data", "", f"{n['conditions']} conditions ({n['gt_conditions_evaluable']} evaluable GT, {n['nulls']} N-eval nulls); "
          f"{n['M1sel_pairs']} M1-sel pairs, {n['kinases_M1sel']} kinases ({n['T-unseen_kinases_M1sel']} T-unseen), "
          f"{n['groups']} GT groups.", ""]
    L += ["## Methods (frozen-core metrics)", "",
          "| method | M1-sel | M1-all | M1-sel-scored | cov | harmony | 0 / 1-4 / 5-19 / ≥20 | H / D / O | T-seen | T-unseen | T-fam-unseen | n sig |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for name, m in R["methods"].items():
        hb = m.get("harmony_per_bin") or {}
        bsrc = m["by_source"]
        L.append(f"| {name} | {f3(m['M1-sel'])} | {f3(m['M1-all'])} | {f3(m['M1-sel-scored'])} | {f3(m['coverage_M1-sel'])} | "
                 f"{f3(m['harmony'])} | {' / '.join(f3(hb.get(b)) for b in ['0', '1-4', '5-19', '>=20'])} | "
                 f"{' / '.join(f3(bsrc.get(s)) for s in ['H', 'D', 'O'])} | {f3(m['T-seen'])} | {f3(m['T-unseen'])} "
                 f"({m['T-unseen_n_pairs']}) | {f3(m['T-family-unseen'])} ({m['T-family-unseen_n_pairs']}) | "
                 f"{m['coverage_significant']['n_significant']} |")
    L += ["", "## TEST Holm family (6 tests; one-sided, B = 2000, seed 0, C8)", "",
          "| test | Δ [95% basic CI] | p | Holm p |", "|---|---|---|---|"]
    for x in R["holm_family"]:
        L.append(f"| {x['test']} | {fd(x)} | {x['p_one_sided']:.4f} | {x['p_holm']:.4f} |")
    L += ["", "## \"Beats\" rule (Holm p < 0.05 AND ΔM1-sel ≥ 0.01 AND ΔM1-sel(T-unseen) > 0 AND point ΔM2 ≥ 0)", "",
          "| comparator | beats | Holm p | ΔM1-sel | ΔT-unseen | ΔM2 (point) |", "|---|---|---|---|---|---|"]
    for lab, b_ in R["beats"].items():
        L.append(f"| {lab} ({b_['comparator']}) | **{'YES' if b_['beats'] else 'no'}** | {b_['p_holm']:.4f} | "
                 f"{b_['dM1-sel']:+.3f} | {b_['dT-unseen']['estimate']:+.3f} | {b_['dM2_point']:+.3f} |")
    L += ["", "If a test fails: estimate + CI only, no claim (DESIGN §5).", "", "## MDE reference", ""]
    for k, v in R["mde"].items():
        L.append(f"- {k}: {f3(v)}")
    L += ["", "## Descriptive comparisons (paired, not in the Holm family)", "",
          "| comparison | ΔM1-sel | ΔM1-all | Δharmony | ΔM2 (matched) | ΔT-seen | ΔT-unseen | ΔT-family-unseen |",
          "|---|---|---|---|---|---|---|---|"]
    for k, d in R["comparisons"].items():
        L.append(f"| {k} | {fd(d['M1-sel'])} | {fd(d['M1-all'])} | {fd(d['harmony'])} | {fd(d['M2'])} | "
                 f"{fd(d['T-seen'])} | {fd(d['T-unseen'])} | {fd(d['T-family-unseen'])} |")
    cv = R["comparisons"].get(f"{FINAL} vs core-v1", {})
    if "has-signature" in cv:
        L += ["", f"{FINAL} vs core-v1 by signature status (C10 iv): has-signature {fd(cv['has-signature'])} "
              f"({cv['has-signature'].get('n_pairs_M1sel')} pairs), no-signature {fd(cv['no-signature'])}."]
    L += ["", "## Calibration (|z| > 1.96; target 5%)", "", "| set | method | overall | per bin | signature kinases |",
          "|---|---|---|---|---|"]
    for setn in ["N-eval nulls", "site-permuted GT"]:
        for meth, c in R["calibration"][setn].items():
            pb = c.get("per_bin") or {}
            sk = c.get("signature_kinases") or {}
            L.append(f"| {setn} | {meth} | {f3(c.get('frac_exceed'))} ({c.get('n_scored', 0)} cells) | "
                     + " / ".join(f"{b}: {f3(pb.get(b))}" for b in ['0', '1-4', '5-19', '>=20'])
                     + f" | {f3(sk.get('frac_exceed'))} ({sk.get('n_cells', '-')}) |")
    rt = R["runtime"]
    L += ["", "## Runtime", "", f"- {rt['bench']}: CPU 1t {rt['cpu']['threads_1']['ms_per_condition']:.2f} ms/cond, "
          f"16t {rt['cpu']['threads_16']['ms_per_condition']:.2f} ms/cond"
          + (f", GPU {rt['gpu']['threads_16']['ms_per_condition']:.2f} ms/cond" if "gpu" in rt else ""),
          f"- compendium ({rt['compendium']['n_conditions']} conditions, prepare + score, median of 3): "
          f"1t {rt['compendium']['cpu_1t_total_s']:.1f} s, 16t {rt['compendium']['cpu_16t_total_s']:.1f} s"]
    if "external" in R:
        L += ["", "## External sets (secondary; benchmarKIN; run_external(confirm_test=True))", "",
              "| method | tyrosine M1-sel | tyr benchmarKIN | hernandez_new M1-sel | hern benchmarKIN | tumour M1-kinase | tumour benchmarKIN |",
              "|---|---|---|---|---|---|---|"]
        for name, r in R["external"]["results"].items():
            g = lambda s, k: r[s]["M1"][k]["value"] if isinstance(r[s]["M1"].get(k), dict) else r[s]["M1"].get(k)
            L.append(f"| {name} | {f3(g('tyrosine', 'M1-sel'))} | {f3(r['tyrosine']['benchmarkin']['auroc_mean'])} | "
                     f"{f3(g('hernandez_new', 'M1-sel'))} | {f3(r['hernandez_new']['benchmarkin']['auroc_mean'])} | "
                     f"{f3(r['tumour']['M1']['M1-kinase'])} | {f3(r['tumour']['benchmarkin']['auroc_mean'])} |")
        L += ["", "Paired bootstrap (final − comparator; B = 2000, seed 0, C8):", "", "| set | comparator | endpoint | Δ [95% CI] | p |",
              "|---|---|---|---|---|"]
        for sname, d in R["external"]["paired_vs_final"].items():
            for c, eps in d.items():
                for ep, x in eps.items():
                    L.append(f"| {sname} | {c} | {ep} | {fd(x)} | {x['p_one_sided']:.4f} |")
    L += ["", "Per-kinase tables: `per_kinase/<method>.tsv`; all numbers: `test_results.json`."]
    (out / "LEDGER.md").write_text("\n".join(L) + "\n")


# =============================================================================================== main
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dry-run-dev", action="store_true", help="smoke test on DEV-confirm; no TEST data, no externals")
    ap.add_argument("--print-freeze-template", action="store_true",
                    help="print the hashes FREEZE.json must contain (does not write anything)")
    ap.add_argument("--skip-runtime", action="store_true", help="dry run only: skip the runtime block")
    a = ap.parse_args()
    if a.print_freeze_template:
        print(json.dumps({**expected_hashes(include_sealed=True), "dirty_tracked_files": dirty_tracked()}, indent=1))
        return
    if a.skip_runtime and not a.dry_run_dev:
        raise SystemExit("--skip-runtime is only allowed with --dry-run-dev")
    dry = a.dry_run_dev
    t0 = time.time()
    np.random.seed(SEED)
    warnings.filterwarnings("ignore", category=FutureWarning)
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    out = ROOT / ("results/test_dryrun" if dry else "results/test")
    log = lambda msg: print(f"[{time.time()-t0:7.1f}s] {msg}", flush=True)
    # ---------------- FREEZE guard (before ANY sealed file is read)
    if dry:
        freeze = {"skipped": "dry run", "working_tree": expected_hashes(include_sealed=False),
                  "modified_tracked_files": dirty_tracked(), "artifacts": check_manifest_artifacts()}
    else:
        freeze = check_freeze()
        if (out / "test_results.json").exists():
            raise SystemExit(f"REFUSED: {out / 'test_results.json'} exists — the TEST run is done once")
    out.mkdir(parents=True, exist_ok=True)
    (out / "per_kinase").mkdir(exist_ok=True); (out / "scores").mkdir(exist_ok=True)
    (out / "freeze_check.json").write_text(json.dumps(freeze, indent=1, default=str))
    log(f"{'DRY RUN on devconfirm' if dry else 'TEST RUN'}; freeze check {'skipped' if dry else 'passed'}")
    vers = versions()
    # ---------------- split, data, leakage checks
    sp_ = Split(dry)
    sig_train = set((ROOT / "models/core_v3/signature_training_conditions.txt").read_text().split())
    dev_groups = read_gt_tsv(ROOT / "results/ground_truth/condition_groups.tsv").set_index("conditionID")["group_id"]
    assert not {dev_groups.get(c) for c in sig_train} & set(sp_.cg["group_id"]), \
        "C10: a signature training condition shares a drug group with an evaluated condition"
    assert not sig_train & set(sp_.conds)
    X = load_columns(sp_.conds)
    log(f"{sp_.name}: {len(sp_.conds)} conditions loaded ({X.shape[0]} sites); {len(sp_.ev.cond_index)} evaluable GT; "
        f"{len(sp_.nulls)} nulls")
    # ---------------- scoring + evaluation
    S, aux = score_all(sp_, X, log)
    for name, Sm in S.items():
        Sm.astype(np.float64).to_parquet(out / "scores" / (name.replace(" ", "_").replace("/", "-") + ".parquet"))
    R = evaluate_all(sp_, S, aux, X, log)
    for name, tab in R.pop("per_kinase").items():
        tab.to_csv(out / "per_kinase" / (name.replace(" ", "_").replace("/", "-") + ".tsv"), sep="\t", index=False,
                   float_format="%.6g")
    R["split"] = sp_.name
    R["a1_fit"] = aux["A1"]
    R["jumpsem_units"] = aux["jumpsem_units"].to_dict(orient="records")
    R["runtime"] = runtime_block(sp_, log) if not a.skip_runtime else {"skipped": True}
    if not dry:
        R["external"] = external_block(log)
    R["provenance"] = {"hashes": {**(freeze.get("verified") or freeze.get("working_tree")),
                                  "FREEZE_sha256": freeze.get("FREEZE_sha256", "n/a (dry run)")},
                       "versions": vers, "argv": sys.argv, "seed": SEED, "B": B,
                       "finished": datetime.now(timezone.utc).isoformat(timespec="seconds"), "wall_s": time.time() - t0,
                       "dry_run": dry}
    if a.skip_runtime:
        R["runtime"] = {"bench": "skipped", "cpu": {"threads_1": {"ms_per_condition": float("nan")},
                                                    "threads_16": {"ms_per_condition": float("nan")}},
                        "compendium": {"n_conditions": 0, "cpu_1t_total_s": float("nan"), "cpu_16t_total_s": float("nan")}}
    (out / "test_results.json").write_text(json.dumps(to_jsonable(R), indent=1, default=str))
    write_ledger(out, R, dry)
    log(f"wrote {out}")


if __name__ == "__main__":
    main()
