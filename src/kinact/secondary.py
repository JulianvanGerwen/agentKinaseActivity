"""Pre-registered SECONDARY TEST method A1 (RESEARCH_TREE D9; results/b6/b6_report.md, results/b6/predeclared.json).

A1 = 0-parameter average of two null-normalised components, recalibrated by their null correlation
(kinact.explore.stacking.combine):
  zvif  core-v1 z_vif (kinact.core, config/core_v1.json);
  wzexp B4 weighted z over the core library (curated_union + r2np@100, alpha 1) plus co-regulation expansion edges
        (kinact.explore.coreg; B4 config frozen for TEST = the per-fold nested pick shared by folds 0/2/4 and the
        majority of every field: seed profiles cur+r2np, q = 0 global factors, top m = 200 sites per kinase, sign 'pos',
        weight = shrunk r (n0 = 30, nmin = 30), beta = 1).
The co-regulation statistics are fitted label-free on ALLOWED conditions = DEV-CV (folds 0-4) conditions with
qc_status ok that do NOT share a control key (PMID, Control, Sample) with any evaluated condition (DESIGN §4; B4 note
for G2: metadata only). Normalisers sigma_j and the null correlation R come from site-permuted DEV-CV GT conditions
(seed 0, the B6 permuted input).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .explore import bayes as B
from .explore import coreg as CR
from .explore import stacking as ST

B4_FINAL = dict(seed="cur+r2np", q=0, m=200, sign="pos", beta=1.0, rmin=0.0, weight="rs", nmin=30, n0=30.0, mmax=200)


def kin_site_matrix(lib: pd.DataFrame, kinases, sites) -> sp.csr_matrix:
    ki = pd.Index(kinases).get_indexer(lib["kinase"]); si = pd.Index(sites).get_indexer(lib["site"])
    keep = (ki >= 0) & (si >= 0)
    e = pd.DataFrame({"k": ki[keep], "s": si[keep], "w": lib["weight"].to_numpy(float)[keep]})
    e = e.groupby(["k", "s"], as_index=False)["w"].max()
    return sp.csr_matrix((e["w"], (e["k"], e["s"])), shape=(len(kinases), len(sites)))


def base_libraries():
    """(seed library cur+r2np [kinase, site, weight max], base design library = B2 library at alpha 1)."""
    from .library import load_library
    cur = load_library("curated_union")[["kinase", "site", "weight"]]
    r2 = load_library("r2np", k=100)[["kinase", "site"]].assign(weight=1.0)
    seed = pd.concat([cur, r2]).groupby(["kinase", "site"], as_index=False)["weight"].max()
    r2b = r2.merge(cur[["kinase", "site"]], how="left", indicator=True)
    r2b = r2b[r2b["_merge"] == "left_only"].drop(columns="_merge")
    return seed, pd.concat([cur, r2b], ignore_index=True)


def control_keys(conditions: pd.DataFrame) -> pd.Series:
    return CR.control_key(conditions)


class CoregWZ:
    """B4 'wz + co-regulation expansion' fitted on the given allowed DEV columns (label-free)."""

    def __init__(self, X_dev: pd.DataFrame, allowed: np.ndarray, kinases, cfg=B4_FINAL, device="cuda"):
        X = X_dev[X_dev.notna().any(axis=1)]
        self.sites = pd.Index(X.index)
        self.kin = list(kinases)
        self.cfg = dict(cfg)
        assert cfg["seed"] == "cur+r2np" and cfg["q"] == 0
        seed, self.base = base_libraries()
        W = kin_site_matrix(seed, self.kin, self.sites)
        Y = CR.standardise(X.to_numpy(np.float64))
        Ya = Y[:, np.asarray(allowed, bool)]
        cand = np.flatnonzero((~np.isnan(Ya)).sum(1) >= cfg["nmin"])
        P = CR.seed_profiles(Ya, W)
        R, N = CR.profile_site_corr(P, Ya[cand], device=device)
        top = CR.topm_per_kinase(R, N, cfg["mmax"], cfg["nmin"], cfg["n0"], exclude=W[:, cand])
        del R, N
        top["idx"] = np.where(top["idx"] >= 0, cand[np.maximum(top["idx"], 0)], -1)
        self.top = top
        self.n_allowed = int(np.asarray(allowed, bool).sum())
        self.exp = CR.expansion_edges(top, self.kin, self.sites, cfg["m"], cfg["rmin"], cfg["sign"], cfg["weight"])
        self.design = CR.make_design(self.sites, self.kin, self.base, self.exp)
        self.W = self.design.W(float(cfg["beta"]))

    def score(self, V: pd.DataFrame) -> pd.DataFrame:
        """weighted z (exp_b4_grid.score_fold 'wz'); V: sites x conditions (any site index; reindexed)."""
        A = V.reindex(self.sites).to_numpy(np.float64)
        Y, D, _, _ = B.standardise(A, np.ones(A.shape[0]))
        rows = self.design.lib_rows
        Dl, Yl = D[rows], Y[rows]
        W = self.W
        W2 = W.multiply(W).tocsr()
        with np.errstate(all="ignore"):
            z = (W.T @ (Dl * Yl)).T / np.sqrt((W2.T @ Dl).T)
        nm = ((W != 0).astype(float).T @ (Dl > 0).astype(float)).T
        return pd.DataFrame(np.where((nm >= 1) & np.isfinite(z), z, np.nan), index=V.columns, columns=self.kin)


class A1:
    """A1 = combine([zvif/s1, wz/s2], w=(1,1), R) with s_j, R from permuted (label-free) DEV-CV GT conditions."""

    def __init__(self, zvif_perm: np.ndarray, wz_perm: np.ndarray):
        self.s = (ST.null_scale(zvif_perm), ST.null_scale(wz_perm))
        self.R = ST.null_corr([zvif_perm / self.s[0], wz_perm / self.s[1]])

    def combine(self, zvif: pd.DataFrame, wz: pd.DataFrame) -> pd.DataFrame:
        wz = wz.reindex(index=zvif.index, columns=zvif.columns)
        s = ST.combine([zvif.to_numpy(float) / self.s[0], wz.to_numpy(float) / self.s[1]], np.array([0.5, 0.5]), self.R)
        return pd.DataFrame(s, index=zvif.index, columns=zvif.columns)
