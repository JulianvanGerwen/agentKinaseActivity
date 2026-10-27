"""RoKAI (Yilmaz et al. 2021, Nat Commun 12:1177) network propagation + z-score inference.

Re-implementation of src/rokai/{rokai_networks,rokai_core,rokai_circuit,rokai_inference}.m
(github.com/serhan-yilmaz/RoKAI, master, fetched 2026-10-27) with the published functional network
(data/rokai_network_data_uniprotkb.mat -> external/rokai/rokai_network.npz via scripts/prep_rokai_network.py):
kinase-substrate (RoKAI's PSP), kinase-kinase STRING PPI (x1e-3), site-site PTMcode co-evolution + structure
distance. Defaults as in rokai.m: IncludePPI/Coevolution/StructureDistance/MissingSites = true, damping 0.5
(r = 1), tau = 1e8. Circuit per condition: R = diag(Rd + rowsum(A)) - A, Rd = 1 (measured) or 1/tau (missing),
v = R^-1 c (c = log2FC of measured sites, 0 elsewhere); refined values are kept for measured sites only.
Inference: RoKAI z-score on the refined values of measured network sites, sum(v_S) / (sd(v) * sqrt(n)).
(rokai.m defaults to 'mean'; the paper / web app report z-score, which is used here.)

Library: lib='rokai_ks' uses the network's own KS edges (the published default); any other library replaces
the kinase-substrate layer (sites mapped by UniProt_position as in rokai_input_mapping.m; kinases mapped to
the STRING kinase layer by UniProt accession), keeping the PPI and site-site layers.
"""
from __future__ import annotations

import os
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from ._common import DEFAULT_MIN_SUB, frame, load_lib_df

ROOT = Path(__file__).resolve().parents[3]
NET = ROOT / "external/rokai"
TAU = 1e8


@lru_cache(maxsize=1)
def load_network():
    z = np.load(NET / "rokai_network.npz")
    mats = {}
    for k in ["Wkin2site", "Wkin2kin", "Wsite2site_coev", "Wsite2site_sd"]:
        mats[k] = sp.csr_matrix((z[f"{k}_data"], z[f"{k}_indices"], z[f"{k}_indptr"]), shape=tuple(z[f"{k}_shape"]))
    sites = pd.read_csv(NET / "rokai_sites.tsv", sep="\t")
    kins = pd.read_csv(NET / "rokai_kinases.tsv", sep="\t")
    return mats, sites, kins


def _site_id(keys):
    """'P28482_T185' -> 'P28482_185' (RoKAI protein-position identifier, upper-case)."""
    s = pd.Series(np.asarray(keys, dtype=object)).astype(str)
    acc = s.str.rsplit("_", n=1).str[0]
    pos = s.str.rsplit("_", n=1).str[1].str.replace(r"[^0-9]", "", regex=True)
    return (acc + "_" + pos).str.upper().to_numpy()


def _first_indexer(index, keys):
    """Position of the first occurrence of each key in `index` (-1 if absent), like MATLAB ismember."""
    ser = pd.Series(np.arange(len(index)), index=np.asarray(index))
    ser = ser[~ser.index.duplicated()]
    return ser.reindex(np.asarray(keys)).fillna(-1).astype(int).to_numpy()


def _universe_acc():
    U = pd.read_csv(ROOT / "results/ground_truth/kinase_universe.tsv", sep="\t", keep_default_na=False)
    return dict(zip(U["uniprot"], U["kinase"]))


@lru_cache(maxsize=8)
def _build(lib_key: str):
    """Network for a library: node order = [sites..., kinases...]; returns A (sym, no diag), site ids, kinase names,
    Wk2s (kinase x site, binary)."""
    mats, rs, rk = load_network()
    rsite = rs["Identifier"].str.upper().to_numpy()
    acc2k = _universe_acc()
    rk_gene = [acc2k.get(a, None) for a in rk["KinaseID"]]
    if lib_key == "rokai_ks":
        # kinases named by U symbol when mappable, else RoKAI's gene label
        kn = [g if g is not None else f"ROKAI:{lab}" for g, lab in zip(rk_gene, rk["Gene"])]
        sites = rsite
        Wk2s = mats["Wkin2site"].tocsr()
        Wk2k = mats["Wkin2kin"]
        kinases = pd.Index(kn)
        Ws2s = (mats["Wsite2site_coev"] + mats["Wsite2site_sd"]).tocsr()
        # duplicate kinase labels (unlikely) are kept as separate nodes; output keeps first
    else:
        lib = load_lib_df(lib_key)
        lsite = _site_id(lib["site"])
        extra = pd.Index(pd.unique(lsite)).difference(pd.Index(rsite))
        sites = np.concatenate([rsite, extra.to_numpy()])
        kinases = pd.Index(sorted(lib["kinase"].unique()))
        e = pd.DataFrame({"k": kinases.get_indexer(lib["kinase"]), "s": _first_indexer(sites, lsite)}).drop_duplicates()
        Wk2s = sp.csr_matrix((np.ones(len(e)), (e["k"], e["s"])), shape=(len(kinases), len(sites)))
        # PPI layer mapped by UniProt -> U symbol
        rk_pos = {g: i for i, g in enumerate(rk_gene) if g is not None}
        kmap = np.array([rk_pos.get(k, -1) for k in kinases])
        P = sp.csr_matrix((np.ones((kmap >= 0).sum()), (np.flatnonzero(kmap >= 0), kmap[kmap >= 0])),
                          shape=(len(kinases), len(rk_gene)))
        Wk2k = (P @ mats["Wkin2kin"] @ P.T).tocsr()
        ns = len(sites)
        Ws = (mats["Wsite2site_coev"] + mats["Wsite2site_sd"]).tocoo()
        Ws2s = sp.csr_matrix((Ws.data, (Ws.row, Ws.col)), shape=(ns, ns))
    Wk2k = Wk2k * 1e-3
    A = sp.bmat([[Ws2s, Wk2s.T], [Wk2s, Wk2k]], format="csr")
    A = A - sp.diags(A.diagonal())
    A.eliminate_zeros()
    return A.tocsr(), pd.Index(sites), kinases, Wk2s


_G = {}


def _solve_cols(cols):
    A, deg, C, ns = _G["A"], _G["deg"], _G["C"], _G["ns"]
    out = []
    for j in cols:
        c = C[:, j]
        meas = ~np.isnan(c)
        Rd = np.where(meas, 1.0, 1.0 / TAU)
        rhs = np.concatenate([np.where(meas, c, 0.0), np.zeros(A.shape[0] - ns)])
        Rd = np.concatenate([Rd, np.full(A.shape[0] - ns, 1.0 / TAU)])
        R = (sp.diags(Rd + deg) - A).tocsc()
        v = spla.spsolve(R, rhs, permc_spec="MMD_AT_PLUS_A", use_umfpack=False)
        out.append((j, np.where(meas, v[:ns], np.nan)))
    return out


def propagate(X: pd.DataFrame, lib_key="rokai_ks", n_jobs=None):
    """Refined site values (network sites x conditions; NaN where not measured)."""
    A, sites, kinases, Wk2s = _build(lib_key)
    xid = _site_id(X.index)
    # input mapping: duplicate ids (same protein-position, different residue label) -> first occurrence
    xi = pd.Index(xid)
    first = ~xi.duplicated()
    pos = _first_indexer(sites, xi[first])
    rows = np.flatnonzero(first)[pos >= 0]
    C = np.full((len(sites), X.shape[1]), np.nan)
    C[pos[pos >= 0]] = X.to_numpy()[rows].astype(np.float64)
    # drop isolated measured/unmeasured nodes from the solve? not needed: isolated measured -> v = x exactly.
    _G.update(A=A, deg=np.asarray(A.sum(1)).ravel(), C=C, ns=len(sites))
    nc = X.shape[1]
    n_jobs = n_jobs or min(32, os.cpu_count())
    chunks = [list(c) for c in np.array_split(np.arange(nc), max(1, min(nc, 2 * n_jobs)))]
    Vs = np.full_like(C, np.nan)
    if n_jobs == 1 or nc < 4:
        res = [r for c in chunks for r in _solve_cols(c)]
    else:
        import multiprocessing as mp
        with ProcessPoolExecutor(n_jobs, mp_context=mp.get_context("fork")) as ex:
            res = [r for rr in ex.map(_solve_cols, chunks) for r in rr]
    for j, v in res:
        Vs[:, j] = v
    _G.clear()
    return Vs, sites, kinases, Wk2s


def rokai(X, lib="rokai_ks", min_substrates=DEFAULT_MIN_SUB, return_n=False, n_jobs=None):
    if not isinstance(lib, str):
        raise TypeError("rokai needs a library name (network is built and cached per library)")
    Vs, sites, kinases, Wk2s = propagate(X, lib, n_jobs=n_jobs)
    M = ~np.isnan(Vs)
    sd = np.nanstd(Vs, axis=0, ddof=1)
    n = Wk2s @ M.astype(np.float64)
    s = Wk2s @ np.where(M, Vs, 0.0)
    with np.errstate(all="ignore"):
        z = s / (sd[None, :] * np.sqrt(n))
    out = frame(z, X, kinases, n, min_substrates, return_n)
    if lib == "rokai_ks":  # keep only U-mapped kinase columns (first occurrence)
        keep = ~pd.Index(kinases).duplicated() & ~pd.Index(kinases).str.startswith("ROKAI:")
        out = (out[0].loc[:, keep], out[1].loc[:, keep]) if return_n else out.loc[:, keep]
    return out
