"""Shared helpers for baseline methods (docs/interfaces.md "Method outputs").

Every method has the signature ``score(X, lib, min_substrates=3, return_n=False, **kw)``:
  X   sites x conditions DataFrame of log2FC (index = canonical site keys; NaN = not measured)
  lib library name (``psp_all``, ``cantley@100`` = predicted library top-k), a DataFrame
      (kinase, site, weight), or a prebuilt ``KSMatrix``.
Returns a conditions x kinases DataFrame (larger = more active; NaN = unscored, never 0); with
return_n=True also the conditions x kinases frame of measured-substrate counts n.
Each condition is scored on its own measured sites only (single-sample mode); a kinase is scored in a
condition iff it has >= min_substrates measured library substrates (benchmarKIN-like threshold).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
import scipy.sparse as sp

DEFAULT_MIN_SUB = 3


def load_lib_df(lib) -> pd.DataFrame:
    if isinstance(lib, pd.DataFrame):
        return lib
    if isinstance(lib, str):
        from ..library import load_library
        if "@" in lib:
            name, k = lib.split("@")
            return load_library(name, k=int(k))
        return load_library(lib)
    raise TypeError(f"unsupported library type {type(lib)}")


@dataclass
class KSMatrix:
    """Kinase x site sparse matrix restricted to the library sites present in X (columns = rows of X[site_rows])."""
    W: sp.csr_matrix          # kinases x len(site_rows), weights
    kinases: pd.Index
    site_rows: np.ndarray     # integer row positions into X


def ks_matrix(lib, site_index, binary=False) -> KSMatrix:
    if isinstance(lib, KSMatrix):
        return lib
    df = load_lib_df(lib)
    site_index = pd.Index(site_index)
    kinases = pd.Index(sorted(df["kinase"].unique()))
    ki = kinases.get_indexer(df["kinase"])
    si = site_index.get_indexer(df["site"])
    keep = si >= 0
    w = np.ones(keep.sum()) if binary else df["weight"].to_numpy(np.float64)[keep]
    e = pd.DataFrame({"k": ki[keep], "s": si[keep], "w": w}).groupby(["k", "s"], as_index=False)["w"].max()
    rows = np.unique(e["s"].to_numpy())
    col = np.searchsorted(rows, e["s"].to_numpy())
    W = sp.csr_matrix((e["w"].to_numpy(np.float64), (e["k"].to_numpy(), col)), shape=(len(kinases), len(rows)))
    return KSMatrix(W, kinases, rows)


def as_array(X: pd.DataFrame) -> np.ndarray:
    return X.to_numpy(np.float64, copy=False) if X.dtypes.iloc[0] == np.float64 else X.to_numpy(np.float64)


def condition_stats(X: pd.DataFrame):
    """Per-condition N (measured), mean, sd (ddof=1), sum of squares of centred values over ALL measured sites."""
    V = X.to_numpy()
    M = ~np.isnan(V)
    N = M.sum(0).astype(np.float64)
    s = np.nansum(V, axis=0, dtype=np.float64)
    mu = s / np.maximum(N, 1)
    ss = np.nansum((V.astype(np.float64) - mu) ** 2, axis=0) if V.shape[0] * V.shape[1] < 5e8 else \
        np.array([np.nansum((V[:, j].astype(np.float64) - mu[j]) ** 2) for j in range(V.shape[1])])
    sd = np.sqrt(ss / np.maximum(N - 1, 1))
    return N, mu, sd, ss


def frame(arr_kc: np.ndarray, X: pd.DataFrame, kinases, n_kc=None, min_sub=DEFAULT_MIN_SUB, return_n=False):
    """kinases x conditions array -> conditions x kinases frame, NaN where n < min_sub or non-finite."""
    a = np.array(arr_kc, dtype=np.float64)
    if n_kc is not None:
        a[n_kc < min_sub] = np.nan
    a[~np.isfinite(a)] = np.nan
    out = pd.DataFrame(a.T, index=X.columns, columns=list(kinases))
    if return_n:
        return out, pd.DataFrame(np.asarray(n_kc).T, index=X.columns, columns=list(kinases))
    return out
