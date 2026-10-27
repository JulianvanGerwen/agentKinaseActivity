"""Closed-form single-sample baselines, vectorised over conditions with per-condition measured-site masks.

Formulas (x = measured log2FC of condition c, N measured sites, mu/sd = mean/SD(ddof=1) over all N, w = library
weights (0 for non-targets), n = # measured substrates; identical to decoupler 2.3.1 run on the condition's
measured sites only, tmin = min_substrates; checked in scripts/check_baseline_refs.py):
  zscore (RoKAI / benchmarKIN / decoupler flavor='RoKAI'): (sum w x / sum|w|) * sqrt(n) / sd     (not mean-centred)
  ksea   (Casado 2013 / KSEAapp / decoupler flavor='KSEA'):  (sum w x / sum|w| - mu) * sqrt(n) / sd
  wmean  (decoupler waggr fun='wmean', times=0):            sum w x / sum|w|
  norm_wmean (decoupler waggr times=1000): z of wmean vs. the permutation null of x over the N sites; computed
         analytically (infinite-permutation limit): sum w (x-mu) / (sd * sqrt(sum w^2 - (sum w)^2/N))
  ulm    (decoupler ulm): t-value of the slope of x ~ 1 + w over the N sites
  mlm    (decoupler mlm): t-values of x ~ 1 + W (all kinases with n >= min_substrates in that condition)
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.linalg as sla

from ._common import DEFAULT_MIN_SUB, condition_stats, frame, ks_matrix


def _sums(X, lib, binary=False):
    K = ks_matrix(lib, X.index, binary=binary)
    N, mu, sd, ss = condition_stats(X)
    V = X.to_numpy()[K.site_rows].astype(np.float64)
    M = ~np.isnan(V)
    Vc = np.where(M, V - mu[None, :], 0.0)
    V0 = np.where(M, V, 0.0)
    Mf = M.astype(np.float64)
    W = K.W
    Wb = W.copy(); Wb.data[:] = 1.0
    W2 = W.multiply(W).tocsr()
    return dict(K=K, N=N, mu=mu, sd=sd, ss=ss, n=Wb @ Mf, Sw=abs(W) @ Mf, Swr=W @ Mf, Sww=W2 @ Mf,
                Swx=W @ V0, Swxc=W @ Vc)


def zscore(X, lib, min_substrates=DEFAULT_MIN_SUB, return_n=False, flavor="RoKAI"):
    s = _sums(X, lib)
    with np.errstate(all="ignore"):
        mean = s["Swx"] / s["Sw"]
        if flavor == "KSEA":
            mean = mean - s["mu"][None, :]
        z = mean * np.sqrt(s["n"]) / s["sd"][None, :]
    return frame(z, X, s["K"].kinases, s["n"], min_substrates, return_n)


def ksea(X, lib, min_substrates=DEFAULT_MIN_SUB, return_n=False):
    return zscore(X, lib, min_substrates, return_n, flavor="KSEA")


def wmean(X, lib, min_substrates=DEFAULT_MIN_SUB, return_n=False):
    s = _sums(X, lib)
    with np.errstate(all="ignore"):
        z = s["Swx"] / s["Sw"]
    return frame(z, X, s["K"].kinases, s["n"], min_substrates, return_n)


def norm_wmean(X, lib, min_substrates=DEFAULT_MIN_SUB, return_n=False):
    s = _sums(X, lib)
    N = s["N"][None, :]
    with np.errstate(all="ignore"):
        var_w = s["Sww"] - s["Swr"] ** 2 / N
        z = s["Swxc"] / (s["sd"][None, :] * np.sqrt(var_w))
    return frame(z, X, s["K"].kinases, s["n"], min_substrates, return_n)


def ulm(X, lib, min_substrates=DEFAULT_MIN_SUB, return_n=False):
    s = _sums(X, lib)
    N = s["N"][None, :]
    with np.errstate(all="ignore"):
        var_w = s["Sww"] - s["Swr"] ** 2 / N           # (N-1) * var(w)
        r = s["Swxc"] / np.sqrt(var_w * s["ss"][None, :])
        df = N - 2
        t = r * np.sqrt(df / ((1.0 - r + 2.2e-16) * (1.0 + r + 2.2e-16)))
    return frame(t, X, s["K"].kinases, s["n"], min_substrates, return_n)


def mlm(X, lib, min_substrates=DEFAULT_MIN_SUB, return_n=False, rcond=1e-10):
    """Per condition OLS via the (centred) normal equations restricted to kinases with n >= min_substrates.

    Rank-deficient designs (e.g. two kinases with identical measured substrate sets) make decoupler's mlm raise;
    here the pseudo-inverse (eigenvalues < rcond*max dropped) is used instead, so such kinases share the
    minimum-norm coefficient (documented deviation; full-rank conditions are identical to decoupler)."""
    K = ks_matrix(lib, X.index)
    N, mu, sd, ss = condition_stats(X)
    V = X.to_numpy()[K.site_rows].astype(np.float64)
    W = K.W.tocsc()
    nk, nc = W.shape[0], V.shape[1]
    out = np.full((nk, nc), np.nan)
    Wb = K.W.copy(); Wb.data[:] = 1.0
    nmat = Wb @ (~np.isnan(V)).astype(np.float64)
    for j in range(nc):
        m = ~np.isnan(V[:, j])
        ks = np.flatnonzero(nmat[:, j] >= min_substrates)
        if len(ks) == 0:
            continue
        Wc = W[:, np.flatnonzero(m)].tocsr()[ks]       # k x n_lib_measured
        xc = V[m, j] - mu[j]
        Sw = np.asarray(Wc.sum(1)).ravel()
        G = (Wc @ Wc.T).toarray() - np.outer(Sw, Sw) / N[j]
        b = Wc @ xc
        evals, evecs = sla.eigh(G)
        keep = evals > rcond * max(evals.max(), 1e-300)
        Ginv = (evecs[:, keep] / evals[keep]) @ evecs[:, keep].T
        beta = Ginv @ b
        df = N[j] - len(ks) - 1
        if df <= 0:
            continue
        sse = max(ss[j] - beta @ b, 0.0) / df
        with np.errstate(all="ignore"):
            out[ks, j] = beta / np.sqrt(sse * np.diag(Ginv))
    return frame(out, X, K.kinases, nmat, min_substrates, return_n)
