"""Shared setup for branch B4b (learned shared-function edge weights on the core-v1 Bayesian core).

DEV only (cache/dev + results/ground_truth DEV files). Core-v1 (D7) = curated_union + r2np@100 (alpha 1), wsd scale,
no winsor, rho_f 3, output z_vif; tau and rho_g per fold as chosen by nested CV over the B2.1 out grid (z_vif, alpha 1,
no winsor; 45 configs; picks reproduced from results/b21/grid_out.tsv, identical to docs/audit/core_v1_audit.md);
rho_e per config fitted on the 7 DEV role=null conditions (label-free).
"""
import json
import time

import numpy as np
import pandas as pd

from exp_b21_common import B2, FOLDS, ROOT, Ctx, bs, full_eval, gt_freeze_hash, inner_metrics  # noqa: F401
from kinact import ledger
from kinact.explore import edgeweights as EW
from kinact.explore.bayes import report_mask
from kinact.library import load_library

OUT = ROOT / "results/b4b"
BRANCH = "B4b"
CORE_V1 = {f: dict(tau=0.125, rho_f=3.0, rho_g=0.25, rho_e=0.01) for f in FOLDS}
CORE_V1["fold0"]["rho_g"] = 0.0
CORE_V1["fold2"].update(tau=0.06, rho_e=0.002)


class BCtx:
    def __init__(self, device="cuda"):
        import torch
        t0 = time.time()
        OUT.mkdir(parents=True, exist_ok=True)
        self.ctx = c = Ctx(ks={"r2np": [100]})
        self.device = device
        md = c.md
        self.allowed = json.load(open(ROOT / "cache/b4/allowed.json"))["allowed"]
        # standardised values on lib rows, all DEV conditions (site stats) and GT conditions (scoring)
        Xall = c.X
        Y, M, _ = B2.standardise(Xall.to_numpy(np.float64), "wsd", None)
        self.all_conds = np.asarray(Xall.columns)
        self.Yall, self.Mall = Y[md.lib_rows].astype(np.float32), M[md.lib_rows].astype(np.float32)
        ci = pd.Index(self.all_conds).get_indexer(c.conds)
        assert (ci >= 0).all()
        self.Yg = torch.as_tensor(self.Yall[:, ci], device=device)
        self.Mg = torch.as_tensor(self.Mall[:, ci], device=device)
        pred500 = {n: load_library(n, k=500) for n in ["r2np", "r3", "cantley"]}
        self.E = EW.edge_table(md, Xall.index, c.Udf, c.lib_cur, load_library("psp_all"),
                               load_library("psp_invivo"), pred500)
        Wb = md.W({"r2np@100": 1.0}).tocsc()
        Wb.data[:] = 1.0
        Wb = Wb.toarray().astype(np.float32)
        self.n = (self.Mall[:, ci].T @ Wb)       # B x p measured edges
        # DEV null conditions (rho_e fit + calibration) and site-permuted GT conditions (calibration)
        ni = pd.Index(self.all_conds).get_indexer(c.nulls)
        self.Yn = torch.as_tensor(self.Yall[:, ni], device=device)
        self.Mn = torch.as_tensor(self.Mall[:, ni], device=device)
        self.n_null = self.Mall[:, ni].T @ Wb
        self.null_fold = pd.Series(c.nulls).map(c.fold_of).to_numpy()
        Yp, Mp, _ = B2.standardise(c.Xp.to_numpy(np.float64), "wsd", None)
        self.Yp = torch.as_tensor(Yp[md.lib_rows].astype(np.float32), device=device)
        self.Mp = torch.as_tensor(Mp[md.lib_rows].astype(np.float32), device=device)
        self.n_perm = Mp[md.lib_rows].T.astype(np.float32) @ Wb
        self.cores, self.rep, self.rep_null, self.rep_perm = {}, {}, {}, {}
        for f in FOLDS:
            K = B2.kernel2(c.Udf, c.kin, CORE_V1[f]["rho_f"], CORE_V1[f]["rho_g"])
            self.cores[f] = EW.DiffCore(self.E, len(md.lib_rows), len(c.kin), K, CORE_V1[f]["tau"], device)
            self.rep[f] = report_mask(self.n, K, 1, True)
            self.rep_null[f] = report_mask(self.n_null, K, 1, True)
            self.rep_perm[f] = report_mask(self.n_perm, K, 1, True)
        self.gt_hash = gt_freeze_hash()
        print(f"[b4b ctx] {time.time()-t0:.0f}s edges {len(self.E)} rows {len(md.lib_rows)}", flush=True)

    def site_cols(self, exclude_folds=()):
        """Columns of Yall usable for site stats: allowed for the first excluded (outer) fold, minus conditions of
        any further excluded (inner validation) fold."""
        ex = list(exclude_folds)
        ok = np.ones(len(self.all_conds), bool)
        if ex:
            ok &= np.isin(self.all_conds, self.allowed[ex[0]])
            fo = pd.Series(self.all_conds).map(self.ctx.fold_of).to_numpy()
            ok &= ~np.isin(fo, ex)
        return np.flatnonzero(ok)

    def Phi(self, exclude_folds, features=EW.FEATURES):
        import torch
        ls, mf = EW.site_stats(self.Yall, self.Mall, self.site_cols(exclude_folds))
        return torch.as_tensor(EW.feature_matrix(self.E, ls, mf, features), device=self.device)

    def conds_of(self, folds):
        return [c for c, f in zip(self.ctx.conds, self.ctx.cf) if f in folds]

    def score(self, f, theta, Phi, out="z_vif", rho_e=None, which="gt"):
        rho_e = CORE_V1[f]["rho_e"] if rho_e is None else rho_e
        Y, M, rep = {"gt": (self.Yg, self.Mg, self.rep[f]), "null": (self.Yn, self.Mn, self.rep_null[f]),
                     "perm": (self.Yp, self.Mp, self.rep_perm[f])}[which]
        return EW.score_theta(self.cores[f], theta, Phi, Y, M, rep, out=out, rho_e=rho_e)

    def rho_e(self, f, theta, Phi):
        return EW.fit_rho_e(self.cores[f], theta, Phi, self.Yn, self.Mn, self.rep_null[f])

    def oof(self, per_fold: dict):
        """per_fold: fold -> full score matrix (B x p); take fold rows -> DataFrame."""
        A = np.full((len(self.ctx.conds), len(self.ctx.kin)), np.nan)
        for f, S in per_fold.items():
            A[self.ctx.cf == f] = S[self.ctx.cf == f]
        return pd.DataFrame(A, index=self.ctx.conds, columns=self.ctx.kin)


def ledger_row(method, cfg, m):
    ledger.append(method, cfg, "dev", {k: m[k] for k in ["M1-sel", "M1-all", "M1-sel-scored", "coverage"] if k in m},
                  branch=BRANCH, folds=range(5))
