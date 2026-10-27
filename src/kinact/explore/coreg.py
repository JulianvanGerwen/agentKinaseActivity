"""Branch B4 — data-driven substrate/reporter expansion from co-regulation across the compendium (DESIGN §7.4).

Everything here is label-free (no GT information) and fold-safe: every statistic that is later used to score a
set of held-out conditions H is computed only on ALLOWED conditions = DEV conditions with qc_status == 'ok' that are
NOT in H and do NOT share a control key (PMID, Control, Sample) with any condition of H (DESIGN §4,
docs/reviews/r1_data.md item 11). The DEV cache already excludes TEST and DEV-confirm groups.

Components
  allowed_mask(meta, conds, heldout)       boolean mask over `conds` of conditions allowed for fitting.
  standardise(V)                           per-condition robust z (median / 1.4826 MAD), clipped at +-clip.
  seed_profiles(Y, Wseed)                  kinase x condition 'activity profiles' = weighted-z of seed substrates
                                           (NaN where no seed is measured).
  profile_site_corr(P, Y, nmin)            pairwise-complete Pearson r between every kinase profile and every site
                                           over the allowed conditions (GPU matmuls), + overlap n.
  site_site_topk(Y, rows, nmin, k)         for each (library) site its k most co-regulated sites (signed r).
  expansion_edges(r, n, ...)               top-m co-regulated sites per kinase (shared rule, no per-kinase params):
                                           rank by shrunk |r| (r * n/(n+n0)), keep |r| >= rmin, drop seed edges;
                                           weight = signed shrunk r (or |r| if sign='pos' keeps only r > 0).
  impute_from_neighbors(...)               idea 3: per condition, value of an unmeasured (library) site predicted
                                           from its measured co-regulated neighbours (signed weighted mean).
Ideas 1/2 (kinase-anchored factor = seed profile; expansion = sites loading on it) use profile_site_corr;
idea 3 (RoKAI-like propagation on a data-driven site graph) uses site_site_topk + impute_from_neighbors.
Shared hyperparameters only (m, rmin, n0, nmin, beta, sign); no per-kinase or per-edge free weights.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp

MAD_K = 1.4826


# ----------------------------------------------------------------------------- allowed conditions
def control_key(meta: pd.DataFrame) -> pd.Series:
    return meta["PMID"].astype(str) + "|" + meta["Control"].astype(str) + "|" + meta["Sample"].astype(str)


def allowed_mask(meta: pd.DataFrame, conds, heldout, require_qc_ok: bool = True) -> np.ndarray:
    """Mask over `conds`: allowed for unsupervised fitting when scoring `heldout` (DESIGN §4).
    Excludes held-out conditions, any condition sharing (PMID, Control, Sample) with a held-out condition, and
    (require_qc_ok) conditions with qc_status != 'ok'."""
    m = meta.reindex(list(conds))
    assert m.index.notna().all() and not m["split"].isna().any(), "unknown condition in conds"
    if (m["split"] == "test").any() or (m["split"] == "devconfirm").any():
        raise ValueError("TEST / DEV-confirm conditions may never enter a co-regulation fit")
    held = set(heldout)
    key = control_key(m)
    held_keys = set(control_key(meta.reindex([c for c in held if c in meta.index])))
    ok = ~m.index.isin(list(held)) & ~key.isin(held_keys).to_numpy()
    if require_qc_ok:
        ok &= (m["qc_status"] == "ok").to_numpy()
    return np.asarray(ok, bool)


def fold_heldout(meta: pd.DataFrame, fold: str):
    """All DEV conditions of a CV fold (GT or not) are held out together."""
    return list(meta.index[meta["split"] == fold])


# ----------------------------------------------------------------------------- standardisation
def standardise(V: np.ndarray, clip: float = 6.0) -> np.ndarray:
    """Per column (condition) robust z; NaN kept."""
    med = np.nanmedian(V, axis=0)
    sc = MAD_K * np.nanmedian(np.abs(V - med[None, :]), axis=0)
    sc = np.where(sc > 0, sc, 1.0)
    return np.clip((V - med[None, :]) / sc[None, :], -clip, clip)


def seed_profiles(Y: np.ndarray, Wseed: sp.csr_matrix, min_meas: int = 1) -> np.ndarray:
    """K x C: sum_s w_ks y_sc / sqrt(sum_s w_ks^2 m_sc) over measured seeds; NaN if < min_meas seeds measured."""
    M = ~np.isnan(Y)
    Y0 = np.where(M, Y, 0.0)
    Mf = M.astype(np.float64)
    num = Wseed @ Y0
    den = Wseed.multiply(Wseed).tocsr() @ Mf
    nm = (Wseed != 0).astype(np.float64) @ Mf
    with np.errstate(all="ignore"):
        P = num / np.sqrt(den)
    return np.where((nm >= min_meas) & (den > 0), P, np.nan)


# ----------------------------------------------------------------------------- correlations (GPU)
def _pairwise_corr(A, MA, B, MB, torch):
    """Pairwise-complete Pearson between rows of A (a x C) and rows of B (b x C); masks MA, MB (float 0/1),
    values already zeroed where missing. Returns r, n (a x b). Inputs are torch tensors."""
    n = MA @ MB.T
    sa = A @ MB.T
    sb = MA @ B.T
    saa = (A * A) @ MB.T
    sbb = MA @ (B * B).T
    sab = A @ B.T
    with torch.no_grad():
        nn_ = n.clamp(min=1)
        cov = sab - sa * sb / nn_
        va = (saa - sa * sa / nn_).clamp(min=0)
        vb = (sbb - sb * sb / nn_).clamp(min=0)
        r = cov / torch.sqrt(va * vb).clamp(min=1e-12)
        r = torch.where((n >= 3) & (va > 1e-9) & (vb > 1e-9), r, torch.zeros_like(r))
    return r.clamp(-1, 1), n


def profile_site_corr(P: np.ndarray, Y: np.ndarray, device="cuda", dtype="float64", chunk=40000):
    """r, n (K x S, float32 / int32) between profiles P (K x C) and sites Y (S x C) over the given columns
    (caller passes only allowed conditions)."""
    import torch
    dt = getattr(torch, dtype)
    MP = ~np.isnan(P)
    Pt = torch.as_tensor(np.where(MP, P, 0.0), dtype=dt, device=device)
    MPt = torch.as_tensor(MP, dtype=dt, device=device)
    R = np.zeros((P.shape[0], Y.shape[0]), np.float32)
    N = np.zeros((P.shape[0], Y.shape[0]), np.int32)
    for i in range(0, Y.shape[0], chunk):
        Yi = Y[i:i + chunk]
        M = ~np.isnan(Yi)
        Yt = torch.as_tensor(np.where(M, Yi, 0.0), dtype=dt, device=device)
        Mt = torch.as_tensor(M, dtype=dt, device=device)
        r, n = _pairwise_corr(Pt, MPt, Yt, Mt, torch)
        R[:, i:i + chunk] = r.float().cpu().numpy()
        N[:, i:i + chunk] = n.round().int().cpu().numpy()
    return R, N


def site_site_topk(Y: np.ndarray, rows: np.ndarray, k: int = 20, nmin: int = 30, n0: float = 30.0,
                   device="cuda", chunk=2000, cand: np.ndarray | None = None):
    """For each site in `rows`: its k most co-regulated other sites (by |r * n/(n+n0)|, n >= nmin) among `cand`
    (default all sites). Returns (nbr idx [len(rows) x k] into Y rows, -1 = none; r; n)."""
    import torch
    cand = np.arange(Y.shape[0]) if cand is None else np.asarray(cand)
    Mc = ~np.isnan(Y[cand])
    Yc = torch.as_tensor(np.where(Mc, Y[cand], 0.0), dtype=torch.float32, device=device)
    Mct = torch.as_tensor(Mc, dtype=torch.float32, device=device)
    # centre candidate sites by their mean over measured conditions (reduces float32 cancellation)
    mu = (Yc.sum(1) / Mct.sum(1).clamp(min=1)).unsqueeze(1)
    Yc = (Yc - mu) * Mct
    cand_t = torch.as_tensor(cand, device=device)
    out_i = np.full((len(rows), k), -1, np.int64)
    out_r = np.zeros((len(rows), k), np.float32)
    out_n = np.zeros((len(rows), k), np.int32)
    for i in range(0, len(rows), chunk):
        rr = rows[i:i + chunk]
        Ms = ~np.isnan(Y[rr])
        Ys = torch.as_tensor(np.where(Ms, Y[rr], 0.0), dtype=torch.float32, device=device)
        Mst = torch.as_tensor(Ms, dtype=torch.float32, device=device)
        mus = (Ys.sum(1) / Mst.sum(1).clamp(min=1)).unsqueeze(1)
        Ys = (Ys - mus) * Mst
        r, n = _pairwise_corr(Ys, Mst, Yc, Mct, torch)
        score = (r * n / (n + n0)).abs()
        score = torch.where(n >= nmin, score, torch.zeros_like(score))
        self_ = cand_t.unsqueeze(0) == torch.as_tensor(rr, device=device).unsqueeze(1)
        score[self_] = 0
        v, j = torch.topk(score, k, dim=1)
        rk = torch.gather(r, 1, j); nk = torch.gather(n, 1, j)
        ok = (v > 0).cpu().numpy()
        out_i[i:i + chunk] = np.where(ok, cand[j.cpu().numpy()], -1)
        out_r[i:i + chunk] = np.where(ok, rk.cpu().numpy(), 0)
        out_n[i:i + chunk] = np.where(ok, nk.round().int().cpu().numpy(), 0)
    return out_i, out_r, out_n


# ----------------------------------------------------------------------------- expansion
def topm_per_kinase(R: np.ndarray, N: np.ndarray, mmax: int = 200, nmin: int = 30, n0: float = 30.0,
                    exclude: sp.csr_matrix | None = None):
    """Store, per kinase, the mmax sites with largest |shrunk r| (n >= nmin), excluding seed edges.
    Returns idx (K x mmax, -1 = none), r, n, rs (shrunk r)."""
    import torch
    Rt = torch.as_tensor(R, device="cuda"); Nt = torch.as_tensor(N, device="cuda").float()
    rs = Rt * Nt / (Nt + n0)
    sc = torch.where(Nt >= nmin, rs.abs(), torch.zeros_like(rs))
    if exclude is not None:
        ex = exclude.tocoo()
        sc[torch.as_tensor(ex.row, device="cuda").long(), torch.as_tensor(ex.col, device="cuda").long()] = 0
    v, j = torch.topk(sc, mmax, dim=1)
    ok = (v > 0).cpu().numpy()
    j = j.cpu().numpy()
    out = dict(idx=np.where(ok, j, -1), r=np.take_along_axis(R, j, 1) * ok,
               n=np.take_along_axis(N, j, 1) * ok, rs=np.take_along_axis(rs.cpu().numpy(), j, 1) * ok)
    return out


def expansion_edges(top: dict, kinases, sites, m: int, rmin: float = 0.0, sign: str = "signed",
                    weight: str = "rs") -> pd.DataFrame:
    """Kinase-site edges from stored top lists: first m entries with |r| >= rmin (sign='pos': r > 0 only, applied
    before truncation). weight: 'rs' signed shrunk r, 'r' signed r, 'bin' sign only."""
    idx, r, rs = top["idx"], top["r"], top["rs"]
    keep = idx >= 0
    keep &= np.abs(r) >= rmin
    if sign == "pos":
        keep &= r > 0
    rank = np.cumsum(keep, axis=1)
    keep &= rank <= m
    ki, col = np.nonzero(keep)
    w = {"rs": rs, "r": r, "bin": np.sign(r)}[weight][ki, col]
    return pd.DataFrame({"kinase": np.asarray(kinases)[ki], "site": np.asarray(sites)[idx[ki, col]],
                         "weight": w.astype(float)})


# ----------------------------------------------------------------------------- imputation (idea 3)
def impute_from_neighbors(Y: np.ndarray, rows: np.ndarray, nbr: np.ndarray, rn: np.ndarray, min_nbr: int = 2,
                          rmin: float = 0.3):
    """For sites `rows` (len L) and conditions (columns of Y): yhat = sum_j r_j y_j / sum_j |r_j| over measured
    neighbours j with |r_j| >= rmin; returns yhat (L x C, NaN if < min_nbr measured neighbours) and the effective
    reliability q = (sum |r_j| y-weighted)^2-free proxy: mean |r_j| over used neighbours (L x C)."""
    L, C = len(rows), Y.shape[1]
    use = (nbr >= 0) & (np.abs(rn) >= rmin)
    j = np.where(use, nbr, 0)
    Yn = Y[j.ravel()].reshape(L, nbr.shape[1], C)                  # L x k x C
    Mn = ~np.isnan(Yn) & use[:, :, None]
    w = np.where(Mn, rn[:, :, None], 0.0)
    num = (w * np.where(Mn, Yn, 0.0)).sum(1)
    den = np.abs(w).sum(1)
    cnt = Mn.sum(1)
    with np.errstate(all="ignore"):
        yhat = num / den
        q = den / np.maximum(cnt, 1)
    ok = cnt >= min_nbr
    return np.where(ok, yhat, np.nan), np.where(ok, q, 0.0)


# ----------------------------------------------------------------------------- global-factor residualisation
def global_factors(P: np.ndarray, q: int) -> np.ndarray:
    """Top-q condition factors (q x C, orthonormal rows) of the kinase-profile matrix P (K x C; NaN -> 0 after
    per-kinase centring). Label-free; computed on allowed conditions only (caller passes allowed columns)."""
    M = ~np.isnan(P)
    mu = np.where(M, P, 0).sum(1) / np.maximum(M.sum(1), 1)
    A = np.where(M, P - mu[:, None], 0.0)
    _, _, Vt = np.linalg.svd(A, full_matrices=False)
    return Vt[:q]


def residualize(Y: np.ndarray, F: np.ndarray, device="cuda", chunk=50000, ridge=1e-3) -> np.ndarray:
    """Per row y with mask m: remove the least-squares fit on factors F (q x C) using measured entries only.
    NaN pattern kept."""
    import torch
    if F is None or len(F) == 0:
        return Y
    Ft = torch.as_tensor(F, dtype=torch.float64, device=device)                # q x C
    q = Ft.shape[0]
    FF = (Ft.unsqueeze(0) * Ft.unsqueeze(1)).reshape(q * q, -1)                 # q^2 x C
    out = np.empty_like(Y)
    for i in range(0, Y.shape[0], chunk):
        Yi = Y[i:i + chunk]
        M = ~np.isnan(Yi)
        Yt = torch.as_tensor(np.where(M, Yi, 0.0), dtype=torch.float64, device=device)
        Mt = torch.as_tensor(M, dtype=torch.float64, device=device)
        G = (Mt @ FF.T).reshape(-1, q, q) + ridge * torch.eye(q, dtype=torch.float64, device=device)
        b = Yt @ Ft.T                                                            # n x q
        beta = torch.linalg.solve(G, b.unsqueeze(-1)).squeeze(-1)
        R = (Yt - beta @ Ft) * Mt
        out[i:i + chunk] = np.where(M, R.cpu().numpy(), np.nan)
    return out


# ----------------------------------------------------------------------------- design for the B2 core
def make_design(site_index, kinases, base: pd.DataFrame, exp: pd.DataFrame | None):
    """bayes.Design with Wc = base library (weights), Wp = expansion edges (signed weights; edges already in base
    dropped). W(beta) = Wc + beta * Wp."""
    from .bayes import Design
    site_index = pd.Index(site_index); kidx = pd.Index(list(kinases))

    def edges(df, agg):
        ki = kidx.get_indexer(df["kinase"]); si = site_index.get_indexer(df["site"])
        keep = (ki >= 0) & (si >= 0)
        e = pd.DataFrame({"s": si[keep], "k": ki[keep], "w": df["weight"].to_numpy(float)[keep]})
        return e.groupby(["s", "k"], as_index=False)["w"].agg(agg)
    eb = edges(base, "max")
    ep = edges(exp, "mean") if exp is not None and len(exp) else eb.iloc[:0]
    ep = ep.merge(eb[["s", "k"]], how="left", indicator=True)
    ep = ep[ep["_merge"] == "left_only"].drop(columns="_merge")
    rows = np.union1d(eb["s"].to_numpy(), ep["s"].to_numpy()).astype(int)

    def mat(e):
        return sp.csr_matrix((e["w"].to_numpy(float), (np.searchsorted(rows, e["s"].to_numpy()), e["k"].to_numpy())),
                             shape=(len(rows), len(kidx)))
    return Design(kidx, rows, mat(eb), mat(ep))
