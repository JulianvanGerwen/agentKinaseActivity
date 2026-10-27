"""Shared setup / evaluation helpers for branch B2.1 (+B3) experiments (scripts/exp_b21_*.py).

DEV only: cache/dev + results/ground_truth DEV files. Scoring of a condition never depends on the outer split (no
site-noise model), so one score matrix per config serves every nested-CV fold.

Nested CV: for each outer fold f, choose the config maximising pooled M1-sel over the other 4 folds (inner_<f>),
take its scores on fold f, pool OOF, evaluate once.

PRE-SPECIFIED choice criterion for the output statistic (item 1; written before any B2.1 result was seen):
the output (mean / z_post / z_freq / hybrid_l / z_vif / empnull) is treated as a hyperparameter and chosen by the
same nested-CV inner M1-sel as tau/rho. The calibrated z (reported on the z scale) is the B3 statistic regardless
of which ranking output wins; if the winning ranking output is not a calibrated z, the deliverable reports both
(score for ranking, z for calibration) and this is stated in the report.
"""
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))
from kinact import bootstrap as bs, controls as ctl, ledger  # noqa: E402
from kinact.data import load_dev, load_dev_meta  # noqa: E402
from kinact.eval import Evaluator, load_gt, psp_bin, read_gt_tsv, to_jsonable  # noqa: E402
from kinact.explore import bayes2 as B2  # noqa: E402
from kinact.library import load_library  # noqa: E402

FOLDS = [f"fold{i}" for i in range(5)]
OUT = ROOT / "results/b21"
CACHE = ROOT / "cache/b21"
BRANCH = "B2.1"


def gt_freeze_hash():
    return hashlib.sha256((ROOT / "results/ground_truth/GT_FREEZE.json").read_bytes()).hexdigest()


class Ctx:
    """Data, evaluator, design and calibration inputs."""

    def __init__(self, ks=None):
        t0 = time.time()
        OUT.mkdir(parents=True, exist_ok=True); CACHE.mkdir(parents=True, exist_ok=True)
        P, C, U, _ = load_gt()
        bins = read_gt_tsv(ROOT / "results/ground_truth/harmony_bins_group.tsv")
        self.Udf = read_gt_tsv(ROOT / "results/ground_truth/kinase_universe.tsv")
        self.P, self.C = P, C
        self.ev = Evaluator(P, C, U, split_filter=FOLDS, bins=bins, zero_policy="ignore")
        self.meta = load_dev_meta()
        self.X = load_dev()
        self.conds = list(self.ev.cond_index)
        self.fold_of = self.meta["split"]
        self.cf = pd.Series(self.conds).map(self.fold_of).to_numpy()
        self.kin = list(self.Udf["kinase"])
        self.nulls = [c for c in C.loc[(C["role"] == "null") & C["split"].isin(FOLDS), "conditionID"]
                      if c in self.X.columns]
        self.Xp = ctl.site_permutation(self.X[self.conds], seed=0, mode="within_condition")
        ks = ks or {"r2np": [100, 200], "r3": [50, 100, 200], "cantley": [50, 100, 200]}
        preds = {}
        for n, kl in ks.items():
            lib = load_library(n, k=max(kl))
            for k in kl:
                preds[f"{n}@{k}"] = lib[lib["rank"] <= k]
        self.lib_cur = load_library("curated_union")
        self.md = B2.build_multi_design(self.X.index, self.kin, self.lib_cur, preds)
        self.gt_hash = gt_freeze_hash()
        print(f"[ctx] {time.time()-t0:.0f}s; {len(self.conds)} GT conds, {len(self.nulls)} DEV nulls, "
              f"lib rows {len(self.md.lib_rows)}; GT_FREEZE {self.gt_hash[:16]}", flush=True)
        self._blocks = {}

    def blocks(self, which: str, names, scale="mad", winsor=None):
        """which: gt | null | perm. Cached block statistics for the listed sources."""
        key = (which, tuple(names), scale, winsor)
        if key not in self._blocks:
            while len(self._blocks) >= 3:          # small LRU: GPU memory
                self._blocks.pop(next(iter(self._blocks)))
                import torch
                torch.cuda.empty_cache()
            X = {"gt": self.X[self.conds], "null": self.X[self.nulls], "perm": self.Xp}[which]
            Y, M, _ = B2.standardise(X.to_numpy(np.float64), scale, winsor)
            self._blocks[key] = B2.block_stats(self.md, Y, M, names=list(names))
        return self._blocks[key]

    def drop_blocks(self):
        self._blocks = {}
        import torch
        torch.cuda.empty_cache()


def inner_metrics(ev, Pp, fold_of):
    f = Pp["conditionID"].map(fold_of).to_numpy()
    res = {f"inner_{g}": ev.metric(Pp[f != g], "M1-sel")["value"] for g in FOLDS}
    for ep in ["M1-sel", "M1-all", "M1-sel-scored"]:
        res[ep] = ev.metric(Pp, ep)["value"]
    res["coverage"] = float(Pp[Pp["in_M1sel"]]["target_scored"].mean())
    return res


def nested(grid: pd.DataFrame, store: dict, ctx: Ctx, mask=None, hp=()):
    sub = grid if mask is None else grid[mask]
    parts, chosen = [], {}
    for f in FOLDS:
        best = sub.loc[sub[f"inner_{f}"].idxmax()]
        chosen[f] = {k: best[k] for k in hp if k in best} | {"key": best["key"], "inner": float(best[f"inner_{f}"])}
        arr = store[best["key"]]
        parts.append(pd.DataFrame(arr[ctx.cf == f], index=np.asarray(ctx.conds)[ctx.cf == f], columns=ctx.kin))
    return pd.concat(parts).loc[ctx.conds], chosen


def t_unseen(ctx, S):
    ev, P, C = ctx.ev, ctx.P, ctx.C
    Pp = ev.pairs(S)
    gc = C[C["split"].isin(FOLDS) & (C["role"] == "gt")].set_index("conditionID")["split"]
    pos = P[(P["cls"] == "direct") & P["confidence"].isin(["high", "medium"]) & P["conditionID"].isin(gc.index)]
    pos = pos.assign(fold=pos["conditionID"].map(gc))
    seen = {f: set(pos.loc[pos["fold"] != f, "kinase"]) for f in FOLDS}
    f = Pp["conditionID"].map(ctx.fold_of)
    uns = np.array([k not in seen[ff] for k, ff in zip(Pp["kinase"], f)])
    r = ev.metric(Pp[uns], "M1-sel")
    return r["value"]


def full_eval(ctx, S, n_perm=2000):
    r = ctx.ev.evaluate(S, n_perm=n_perm)
    v = lambda x: x["value"] if isinstance(x, dict) and "value" in x else x
    return {"M1-sel": v(r["M1-sel"]), "M1-all": v(r["M1-all"]), "M1-sel-scored": v(r["M1-sel-scored"]),
            "coverage": r["coverage_M1-sel"], "harmony": v(r["harmony"]), "harmony_per_bin": r["harmony"]["per_bin"],
            "by_source": {k: v(x) for k, x in r["by_source"].items()},
            "n_sig": r["coverage_significant"]["n_significant"], "T-unseen": t_unseen(ctx, S)}


def calib_bins(Z: np.ndarray, N: np.ndarray, ok=None):
    ok = (~np.isnan(Z)) if ok is None else ok & ~np.isnan(Z)
    lab = psp_bin(N.ravel()).reshape(N.shape)
    out = {"all": float((np.abs(Z[ok]) > 1.96).mean())}
    for b in ["0", "1-4", "5-19", ">=20"]:
        m = ok & (lab == b)
        if m.any():
            out[b] = float((np.abs(Z[m]) > 1.96).mean())
    return out


def ledger_row(method, cfg, m):
    ledger.append(method, cfg, "dev", {k: m[k] for k in ["M1-sel", "M1-all", "M1-sel-scored", "coverage"] if k in m},
                  branch=BRANCH, folds=range(5))


def b2_reference(ctx):
    """B2 nested OOF scores (B2_free variant) from the B2 grid archives (cache/b2)."""
    import exp_b2_eval as E2
    g, store = E2.load_grid()
    assert list(store["conds"]) == ctx.conds
    core = g["output"] == "z"
    S, chosen = E2.nested(g, store, ctx.conds, ctx.fold_of, core)
    return S.reindex(columns=ctx.kin), chosen


def dump(name, obj):
    (OUT / name).write_text(json.dumps(to_jsonable(obj), indent=1, default=str))
