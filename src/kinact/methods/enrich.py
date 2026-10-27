"""Rank/permutation-based baselines: VIPER (decoupler 2.3.1 used directly) and GSEA (numba reimplementation of
decoupler's gsea, checked against dc.mt.gsea in scripts/check_baseline_refs.py).

Both are run per condition on the condition's measured sites only (decoupler refuses NaN input)."""
from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor

import numba as nb
import numpy as np
import pandas as pd

from ._common import DEFAULT_MIN_SUB, frame, ks_matrix


# ----------------------------------------------------------------------------------------------- VIPER
_G = {}


def _viper_one(j):
    from decoupler.mt._viper import _func_viper
    V, W, min_sub, kw = _G["V"], _G["W"], _G["min_sub"], _G["kw"]
    x_all = _G["X"][:, j]
    m_all = ~np.isnan(x_all)
    x = x_all[m_all].astype(np.float64)
    # library rows that are measured, mapped to positions inside x
    pos_all = np.cumsum(m_all) - 1
    v = V[:, j]
    mr = ~np.isnan(v)
    Wm = W[:, np.flatnonzero(mr)]
    n = np.asarray((Wm != 0).sum(1)).ravel()
    ks = np.flatnonzero(n >= min_sub)
    res = np.full(W.shape[0], np.nan)
    if len(ks) == 0:
        return j, res
    adj = np.zeros((len(x), len(ks)))
    rows = pos_all[_G["site_rows"][np.flatnonzero(mr)]]
    Wk = Wm[ks].toarray().T                          # measured lib sites x ks
    adj[rows, :] = Wk
    es, _ = _func_viper(x.reshape(1, -1), adj, **kw)
    res[ks] = es[0]
    return j, res


def viper(X, lib, min_substrates=DEFAULT_MIN_SUB, return_n=False, n_jobs=None, pleiotropy=True):
    """decoupler.mt viper (aREA + pleiotropy correction, decoupler defaults) per condition."""
    K = ks_matrix(lib, X.index)
    Xa = X.to_numpy()
    V = Xa[K.site_rows].astype(np.float64)
    Wb = K.W.copy(); Wb.data[:] = 1.0
    nmat = Wb @ (~np.isnan(V)).astype(np.float64)
    _G.update(X=Xa, V=V, W=K.W.tocsr(), site_rows=K.site_rows, min_sub=min_substrates,
              kw={"pleiotropy": pleiotropy})
    out = np.full((len(K.kinases), V.shape[1]), np.nan)
    n_jobs = n_jobs or min(32, os.cpu_count())
    if n_jobs == 1 or V.shape[1] < 4:
        for j in range(V.shape[1]):
            out[:, j] = _viper_one(j)[1]
    else:
        import multiprocessing as mp
        with ProcessPoolExecutor(n_jobs, mp_context=mp.get_context("fork")) as ex:
            for j, r in ex.map(_viper_one, range(V.shape[1]), chunksize=max(1, V.shape[1] // (4 * n_jobs))):
                out[:, j] = r
    _G.clear()
    return frame(out, X, K.kinases, nmat, min_substrates, return_n)


# ----------------------------------------------------------------------------------------------- GSEA
@nb.njit(cache=True)
def _es_from_hits(pos, w, nvar):
    """decoupler _esrank on a set given its sorted hit positions pos (0-based ranks) and |value| weights w."""
    m = pos.size
    tot = w.sum()
    if tot == 0.0:
        return 0.0
    dec = 1.0 / (nvar - m)
    mx_pos = 0.0
    mx_neg = 0.0
    h = 0.0
    for j in range(m):
        before = h - (pos[j] - j) * dec           # running sum just before hit j (after the miss at pos-1)
        if pos[j] - j > 0 and before < mx_neg:
            mx_neg = before
        h += w[j] / tot
        after = h - (pos[j] - j) * dec
        if after > mx_pos:
            mx_pos = after
        if after < mx_neg:
            mx_neg = after
    if mx_pos > -mx_neg:
        return mx_pos
    return mx_neg


@nb.njit(cache=True)
def _gsea_condition(row_sorted_abs, set_ptr, set_pos, times, seed):
    nvar = row_sorted_abs.size
    nset = set_ptr.size - 1
    out = np.zeros(nset)
    np.random.seed(seed)
    flag = np.zeros(nvar, dtype=np.bool_)
    for k in range(nset):
        pos = np.sort(set_pos[set_ptr[k]:set_ptr[k + 1]])
        m = pos.size
        if m == 0 or m >= nvar:
            out[k] = np.nan
            continue
        es = _es_from_hits(pos, row_sorted_abs[pos], nvar)
        sp_ = 0.0; npos = 0; sn = 0.0; nneg = 0
        rp = np.empty(m, dtype=np.int64)
        for t in range(times):
            # Floyd's algorithm: uniform m-subset of 0..nvar-1
            c = 0
            for jj in range(nvar - m, nvar):
                r = np.random.randint(0, jj + 1)
                if flag[r]:
                    r = jj
                flag[r] = True
                rp[c] = r
                c += 1
            for q in range(m):
                flag[rp[q]] = False
            ps = np.sort(rp)
            e = _es_from_hits(ps, row_sorted_abs[ps], nvar)
            if e >= 0:
                sp_ += e; npos += 1
            else:
                sn += e; nneg += 1
        if es >= 0 and npos > 0:
            out[k] = es / (sp_ / npos)
        elif es < 0 and nneg > 0:
            out[k] = -es / (sn / nneg)
        else:
            out[k] = 0.0
    return out


def _gsea_cols(args):
    cols, = args
    Xa, site_rows, W, min_sub, times, seed = (_G[k] for k in ["X", "site_rows", "W", "min_sub", "times", "seed"])
    res = []
    for j in cols:
        x_all = Xa[:, j]
        m_all = ~np.isnan(x_all)
        x = x_all[m_all].astype(np.float64)
        order = np.argsort(-x, kind="stable")
        rank = np.empty(len(x), np.int64); rank[order] = np.arange(len(x))
        pos_all = np.cumsum(m_all) - 1
        v = x_all[site_rows]
        mr = ~np.isnan(v)
        Wm = W[:, np.flatnonzero(mr)].tocsr()
        rows = rank[pos_all[site_rows[np.flatnonzero(mr)]]]
        n = np.diff(Wm.indptr)
        o = _gsea_condition(np.abs(x[order]), Wm.indptr.astype(np.int64), rows[Wm.indices].astype(np.int64),
                            times, seed)
        o[n < min_sub] = np.nan
        res.append((j, o, n))
    return res


def gsea(X, lib, min_substrates=DEFAULT_MIN_SUB, return_n=False, times=1000, seed=42, n_jobs=None):
    """GSEA NES as in decoupler.mt.gsea (unweighted set membership, |x|-weighted running sum, NES = ES / mean of
    same-sign permutation ES; set-membership permutations, `times` per kinase). Permutation draws differ from
    decoupler's (per-kinase Floyd sampling, seed per condition) -> agreement checked by Spearman."""
    K = ks_matrix(lib, X.index, binary=True)
    Xa = X.to_numpy()
    _G.update(X=Xa, site_rows=K.site_rows, W=K.W.tocsr(), min_sub=min_substrates, times=times, seed=seed)
    nc = Xa.shape[1]
    out = np.full((len(K.kinases), nc), np.nan)
    nmat = np.zeros((len(K.kinases), nc))
    n_jobs = n_jobs or min(32, os.cpu_count())
    chunks = [list(c) for c in np.array_split(np.arange(nc), max(1, min(nc, 4 * n_jobs)))]
    if n_jobs == 1 or nc < 4:
        results = [r for c in chunks for r in _gsea_cols((c,))]
    else:
        import multiprocessing as mp
        _gsea_cols(([0],))  # compile before fork
        with ProcessPoolExecutor(n_jobs, mp_context=mp.get_context("fork")) as ex:
            results = [r for rr in ex.map(_gsea_cols, [(c,) for c in chunks]) for r in rr]
    for j, o, n in results:
        out[:, j] = o; nmat[:, j] = n
    _G.clear()
    return frame(out, X, K.kinases, nmat, min_substrates, return_n)
