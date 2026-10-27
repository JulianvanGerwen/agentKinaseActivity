"""Motif-only (Kinase Library / PhosX-like) baselines on the Cantley atlas library (cache/libraries/cantley.parquet:
per kinase the top-500 candidate sites by the unprimed PSPA atlas score, weight = within-kinase percentile).

motif_zscore(k): RoKAI z-score over each kinase's top-k motif-predicted sites (lib 'cantley@k', binary).
phosx_like(n):   PhosX-style (Lussana & Petsalaki 2024): each site is assigned to its top-n kinases by percentile
                 (default n=10; restricted to the stored per-kinase top-500, i.e. percentile >= ~0.997), then the
                 GSEA running-sum NES (decoupler-style, `enrich.gsea`) per kinase.
kl_enrich(n, thr): Kinase Library enrichment (Johnson 2023): sites assigned to their top-n kinases (n=15); per
                 kinase one-sided hypergeometric (Fisher) tests up-regulated vs unchanged and down-regulated vs
                 unchanged (|log2FC| >= thr), score = log10(p_down) - log10(p_up) (positive = enriched among up sites).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

from ._common import DEFAULT_MIN_SUB, frame, ks_matrix
from .enrich import gsea
from .linear import zscore


def site_topn(n=10, name="cantley") -> pd.DataFrame:
    from ..library import load_library
    lib = load_library(name)
    lib = lib.sort_values(["site", "weight", "kinase"], ascending=[True, False, True])
    lib = lib.groupby("site", sort=False).head(n)
    return lib[["kinase", "site", "weight"]].assign(weight=1.0).reset_index(drop=True)


def motif_zscore(X, lib="cantley@100", min_substrates=DEFAULT_MIN_SUB, return_n=False):
    from ._common import load_lib_df
    df = load_lib_df(lib).assign(weight=1.0)
    return zscore(X, df, min_substrates, return_n)


def phosx_like(X, lib=10, min_substrates=DEFAULT_MIN_SUB, return_n=False, **kw):
    df = site_topn(int(lib)) if not isinstance(lib, pd.DataFrame) else lib
    return gsea(X, df, min_substrates, return_n, **kw)


def kl_enrich(X, lib=15, min_substrates=DEFAULT_MIN_SUB, return_n=False, thr=0.5):
    df = site_topn(int(lib)) if not isinstance(lib, pd.DataFrame) else lib
    K = ks_matrix(df, X.index, binary=True)
    Xa = X.to_numpy().astype(np.float64)
    up_all = (Xa >= thr).sum(0); dn_all = (Xa <= -thr).sum(0)
    meas = (~np.isnan(Xa)).sum(0)
    un_all = meas - up_all - dn_all
    V = Xa[K.site_rows]
    M = (~np.isnan(V)).astype(np.float64)
    U_ = (V >= thr).astype(np.float64); D_ = (V <= -thr).astype(np.float64)
    W = K.W
    n = W @ M; a_up = W @ U_; a_dn = W @ D_
    a_un = n - a_up - a_dn
    with np.errstate(all="ignore"):
        # population = up ∪ unchanged sites; successes = up sites; draws = kinase sites in the population
        p_up = hypergeom.sf(a_up - 1, (up_all + un_all)[None, :], up_all[None, :], a_up + a_un)
        p_dn = hypergeom.sf(a_dn - 1, (dn_all + un_all)[None, :], dn_all[None, :], a_dn + a_un)
        s = np.log10(np.clip(p_dn, 1e-300, 1)) - np.log10(np.clip(p_up, 1e-300, 1))
    return frame(s, X, K.kinases, n, min_substrates, return_n)
