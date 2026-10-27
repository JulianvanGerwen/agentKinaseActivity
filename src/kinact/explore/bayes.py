"""Branch 2 — Bayesian linear-Gaussian core (DESIGN §7.2; docs/reviews/r1_algo.md item 2).

Model for one condition c (single-sample):
    y_s = (x_s - med_c) / sigma_c = sum_k W_sk a_k + b + eps_s,   eps_s ~ N(0, v_s)
    a ~ N(0, tau^2 K),  b ~ N(0, 1e6) (intercept='joint') or b := 0 (intercept='median': x centred by its median)
  v_s     per-site relative noise variance, estimated label-free from allowed compendium conditions
          (robust: median of squared robust-standardised values / 0.4549, shrunk towards 1 with n0 pseudo-obs,
          normalised to median 1). v_s = 1 for sites without enough data, and in the 'no noise model' ablation.
  sigma_c robust per-condition scale: 1.4826 * MAD over all measured sites of (x - med_c)/sqrt(v_s).
  W       W_cur + alpha * W_pred: curated library edges (weight 1) plus predicted edges not in the curated
          library (binary, scaled by alpha). Only measured sites enter (site-specific measured mask per condition).
  K       = (I + rho_f F + rho_g G) / (1 + rho_f + rho_g) on the diagonal-normalised scale, F = same Manning family,
          G = same Manning group (heterogeneous groups 'Other'/'Atypical' get no group coupling).
Posterior (closed form): Lambda = W_M' D W_M + (tau^2 K)^-1, mean = Lambda^-1 W_M' D y, sd = sqrt(diag Lambda^-1),
z = mean / sd. Hyperparameters: tau, rho_f, rho_g, alpha, intercept (5).
Reporting: kinase k is scored in c iff it has >= min_sub measured library substrates; with borrow=True also
kinases without measured substrates that have a K-relative (K_kl > 0) with >= min_sub measured substrates.
Nothing here uses GT labels; no per-kinase / per-edge free parameters.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
import scipy.linalg as sla
import scipy.sparse as sp

NO_GROUP_COUPLING = ("Other", "Atypical")
MAD_K = 1.4826
CHI2_1_MEDIAN = 0.454936


# ----------------------------------------------------------------------------- kinase kernel
def kinase_kernel(universe: pd.DataFrame, kinases, rho_f: float, rho_g: float,
                  no_group=NO_GROUP_COUPLING) -> np.ndarray:
    U = universe.set_index("kinase").reindex(list(kinases))
    fam = U["manning_family"].fillna("").astype(str).to_numpy()
    grp = U["manning_group"].fillna("").astype(str).to_numpy()
    F = (fam[:, None] == fam[None, :]) & (fam[:, None] != "")
    G = (grp[:, None] == grp[None, :]) & (grp[:, None] != "") & ~np.isin(grp, list(no_group))[:, None]
    K = np.eye(len(U)) + rho_f * F + rho_g * G
    d = np.sqrt(np.diag(K))
    return K / d[:, None] / d[None, :]


# ----------------------------------------------------------------------------- noise model
def robust_center_scale(V: np.ndarray, w: np.ndarray | None = None):
    """Per column median and 1.4826*MAD of V (or of (V-med)*w row-weights). NaN = missing."""
    med = np.nanmedian(V, axis=0)
    R = V - med[None, :]
    if w is not None:
        R = R * w[:, None]
    sc = MAD_K * np.nanmedian(np.abs(R), axis=0)
    return med, sc


def site_noise_variance(V: np.ndarray, n0: float = 10.0, min_obs: int = 3, chunk: int = 20000) -> np.ndarray:
    """Label-free relative noise variance per site (rows of V = allowed conditions only).

    r = (x - med_c)/scale_c per condition; v_raw = median(r^2)/0.4549 over the site's measured conditions;
    v = (n v_raw + n0 * 1) / (n + n0) (shrink to 1); sites with < min_obs measurements -> 1; normalised so the
    median over sites with data is 1."""
    med, sc = robust_center_scale(V)
    v = np.ones(V.shape[0])
    n = np.zeros(V.shape[0])
    for i in range(0, V.shape[0], chunk):
        R = (V[i:i + chunk] - med[None, :]) / sc[None, :]
        n[i:i + chunk] = np.sum(~np.isnan(R), axis=1)
        with np.errstate(all="ignore"):
            v[i:i + chunk] = np.nanmedian(R * R, axis=1) / CHI2_1_MEDIAN
    ok = n >= min_obs
    v = np.where(ok, (n * np.where(ok, v, 1.0) + n0) / (n + n0), 1.0)
    v = v / np.median(v[ok]) if ok.any() else v
    return np.clip(v, 0.05, 50.0)


# ----------------------------------------------------------------------------- library design
@dataclass
class Design:
    kinases: pd.Index
    lib_rows: np.ndarray        # row positions (into the site index) of library sites
    Wc: sp.csr_matrix           # lib sites x kinases, curated part (weights)
    Wp: sp.csr_matrix           # lib sites x kinases, predicted part (binary; edges not in curated)

    def W(self, alpha: float) -> sp.csr_matrix:
        return (self.Wc + alpha * self.Wp).tocsr() if alpha else self.Wc


def build_design(site_index, kinases, lib_cur: pd.DataFrame, lib_pred: pd.DataFrame | None = None) -> Design:
    site_index = pd.Index(site_index)
    kidx = pd.Index(list(kinases))

    def edges(df, binary):
        ki = kidx.get_indexer(df["kinase"]); si = site_index.get_indexer(df["site"])
        keep = (ki >= 0) & (si >= 0)
        w = np.ones(keep.sum()) if binary else df["weight"].to_numpy(float)[keep]
        return pd.DataFrame({"s": si[keep], "k": ki[keep], "w": w}).groupby(["s", "k"], as_index=False)["w"].max()

    ec = edges(lib_cur, False)
    ep = edges(lib_pred, True) if lib_pred is not None else ec.iloc[:0]
    ep = ep.merge(ec[["s", "k"]], how="left", indicator=True)
    ep = ep[ep["_merge"] == "left_only"].drop(columns="_merge")
    rows = np.union1d(ec["s"].to_numpy(), ep["s"].to_numpy()).astype(int)

    def mat(e):
        return sp.csr_matrix((e["w"].to_numpy(float), (np.searchsorted(rows, e["s"].to_numpy()), e["k"].to_numpy())),
                             shape=(len(rows), len(kidx)))
    return Design(kidx, rows, mat(ec), mat(ep))


# ----------------------------------------------------------------------------- sufficient statistics
def standardise(V: np.ndarray, v: np.ndarray):
    """y = (x - med)/sigma with sigma = robust scale of (x-med)/sqrt(v); returns y, d = 1/v (0 if missing)."""
    w = 1.0 / np.sqrt(v)
    med, sc = robust_center_scale(V, w)
    Y = (V - med[None, :]) / sc[None, :]
    M = ~np.isnan(V)
    D = np.where(M, (1.0 / v)[:, None], 0.0)
    return np.where(M, Y, 0.0), D, med, sc


def suff_stats(design: Design, alpha: float, Y: np.ndarray, D: np.ndarray, device="cuda", chunk=24):
    """Per condition (columns of Y, D: all sites x B, standardised; D = 0 where missing):
    G (B, p+1, p+1) with index 0 = intercept, h (B, p+1), n (B, p) measured-substrate counts. torch tensors."""
    import torch
    W = design.W(alpha)
    Wd = torch.as_tensor(W.toarray(), dtype=torch.float32, device=device)              # S x p
    Wb = (Wd != 0).float()
    rows = design.lib_rows
    B = Y.shape[1]
    p = Wd.shape[1]
    Gs, hs, ns = [], [], []
    for i in range(0, B, chunk):
        Dl = torch.as_tensor(D[rows, i:i + chunk], dtype=torch.float32, device=device)  # S x b
        Yl = torch.as_tensor(Y[rows, i:i + chunk], dtype=torch.float32, device=device)
        b = Dl.shape[1]
        G = torch.empty((b, p + 1, p + 1), dtype=torch.float32, device=device)
        A = Wd.T.unsqueeze(0) * Dl.T.unsqueeze(1)                                        # b x p x S
        G[:, 1:, 1:] = torch.bmm(A, Wd.unsqueeze(0).expand(b, -1, -1))
        g0 = (Wd.T @ Dl).T                                                               # b x p
        G[:, 0, 1:] = g0
        G[:, 1:, 0] = g0
        G[:, 0, 0] = torch.as_tensor(D[:, i:i + chunk].sum(0), dtype=torch.float32, device=device)
        h = torch.empty((b, p + 1), dtype=torch.float32, device=device)
        h[:, 1:] = (Wd.T @ (Dl * Yl)).T
        h[:, 0] = torch.as_tensor((D[:, i:i + chunk] * Y[:, i:i + chunk]).sum(0), dtype=torch.float32, device=device)
        Gs.append(G); hs.append(h)
        ns.append((Wb.T @ (Dl > 0).float()).T)
    return torch.cat(Gs), torch.cat(hs), torch.cat(ns)


def posterior(G, h, K: np.ndarray, tau: float, intercept: str = "median", chunk=256):
    """Batched closed-form posterior. Returns mean, sd (B x p, float64 numpy)."""
    import torch
    dev = G.device
    p = K.shape[0]
    P0 = torch.linalg.inv(torch.as_tensor(K, dtype=torch.float64, device=dev)) / tau ** 2   # (numpy inv: BLAS oversubscription)
    means, sds = [], []
    for i in range(0, G.shape[0], chunk):
        g = G[i:i + chunk].double(); hh = h[i:i + chunk].double()
        if intercept == "joint":
            L = g.clone()
            L[:, 1:, 1:] += P0
            L[:, 0, 0] += 1e-6
            rhs = hh
        elif intercept in ("median", "none"):
            L = g[:, 1:, 1:] + P0
            rhs = hh[:, 1:]
        else:
            raise ValueError(intercept)
        C = torch.linalg.cholesky(L)
        mu = torch.cholesky_solve(rhs.unsqueeze(-1), C).squeeze(-1)
        var = torch.diagonal(torch.cholesky_inverse(C), dim1=1, dim2=2)
        if intercept == "joint":
            mu, var = mu[:, 1:], var[:, 1:]
        means.append(mu.cpu().numpy()); sds.append(var.sqrt().cpu().numpy())
    return np.concatenate(means), np.concatenate(sds)


def report_mask(n: np.ndarray, K: np.ndarray, min_sub: int = 1, borrow: bool = False) -> np.ndarray:
    """B x p boolean: kinase scored in condition."""
    has = n >= min_sub
    if not borrow:
        return has
    rel = (K > 1e-12) & ~np.eye(K.shape[0], dtype=bool)
    return has | ((has.astype(float) @ rel.astype(float)) > 0)


# ----------------------------------------------------------------------------- one-call scorer
@dataclass
class BayesConfig:
    tau: float = 1.0
    rho_f: float = 1.0
    rho_g: float = 0.0
    alpha: float = 0.0
    intercept: str = "median"
    min_sub: int = 1
    borrow: bool = False
    noise_model: bool = True
    output: str = "z"           # 'z' or 'mean'


def score(X: pd.DataFrame, design: Design, universe: pd.DataFrame, cfg: BayesConfig, v: np.ndarray | None = None,
          device="cuda", return_all=False):
    """Score conditions (columns of X; index must equal the site index used for `design`). v: per-site noise
    variances aligned to X.index (None or cfg.noise_model False -> 1). Returns conditions x kinases (NaN = unscored)."""
    V = X.to_numpy(np.float64)
    vv = np.ones(V.shape[0]) if (v is None or not cfg.noise_model) else v
    Y, D, _, _ = standardise(V, vv)
    G, h, n = suff_stats(design, cfg.alpha, Y, D, device=device)
    K = kinase_kernel(universe, design.kinases, cfg.rho_f, cfg.rho_g)
    mu, sd = posterior(G, h, K, cfg.tau, cfg.intercept)
    n = n.cpu().numpy()
    rep = report_mask(n, K, cfg.min_sub, cfg.borrow)
    z = mu / sd
    out = np.where(rep, z if cfg.output == "z" else mu, np.nan)
    S = pd.DataFrame(out, index=X.columns, columns=design.kinases)
    if return_all:
        return S, dict(mean=pd.DataFrame(mu, index=X.columns, columns=design.kinases),
                       sd=pd.DataFrame(sd, index=X.columns, columns=design.kinases),
                       z=pd.DataFrame(z, index=X.columns, columns=design.kinases),
                       n=pd.DataFrame(n, index=X.columns, columns=design.kinases), report=rep)
    return S


# ----------------------------------------------------------------------------- CPU single-condition path
class SingleSample:
    """CPU scorer for one condition at a time (runtime benchmark / CLI path). Precomputes P0 and the CSR design."""

    def __init__(self, design: Design, universe: pd.DataFrame, cfg: BayesConfig, v_lib: np.ndarray | None,
                 v_median_all: float = 1.0):
        self.cfg = cfg
        self.W = design.W(cfg.alpha).tocsr()
        self.K = kinase_kernel(universe, design.kinases, cfg.rho_f, cfg.rho_g)
        self.P0 = np.linalg.inv(self.K) / cfg.tau ** 2
        self.v_lib = np.ones(self.W.shape[0]) if (v_lib is None or not cfg.noise_model) else v_lib
        self.kinases = design.kinases

    def score_one(self, x_lib: np.ndarray, med: float, sigma: float):
        """x_lib: values at library sites (NaN = missing); med/sigma: condition centre/scale over all sites."""
        m = np.flatnonzero(~np.isnan(x_lib))
        Wm = self.W[m]
        d = 1.0 / self.v_lib[m]
        y = (x_lib[m] - med) / sigma
        WD = Wm.T.multiply(d[None, :]).tocsr()
        G = (WD @ Wm).toarray()
        h = WD @ y
        L = G + self.P0
        C = sla.cholesky(L, lower=True, check_finite=False)
        mu = sla.cho_solve((C, True), h, check_finite=False)
        Ci = sla.solve_triangular(C, np.eye(len(L)), lower=True, check_finite=False)
        sd = np.sqrt((Ci * Ci).sum(0))
        n = np.asarray((Wm != 0).sum(0)).ravel()
        rep = report_mask(n[None, :], self.K, self.cfg.min_sub, self.cfg.borrow)[0]
        return np.where(rep, mu / sd, np.nan), mu, sd, n
