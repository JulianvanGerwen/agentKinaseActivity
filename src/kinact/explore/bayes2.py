"""Branch B2.1 + B3 — refinements and calibration of the Bayesian linear-Gaussian core (bayes.py).

Same model as bayes.py (y = (x - med_c)/sigma_c = W a + eps, a ~ N(0, tau^2 K), median-centred, no site-noise model;
borrow on), with these additions:

Multi-source library (item 3)
  W = W_cur + sum_j alpha_j W_j, W_j = binary global top-k_j predicted edges of source j (r2np, r3, cantley) that
  are not curated edges. An edge predicted by several sources gets the SUM of its source weights. Because
  G = W'MW and h = W'My are bilinear in the alphas, per-source blocks (G_ij, h_i, g0_i) are computed once and
  combined for any alpha vector (`combine`).

Output statistics (item 1), all from the same posterior (Lambda = G + (tau^2 K)^-1, A = Lambda^-1 W'M):
  mean      mu = A y
  z_post    mu / sqrt(diag Lambda^-1)                     (B2 output; posterior z)
  z_freq    mu / sqrt(q),  q = diag(A A') = diag(Lambda^-1 G Lambda^-1)
            frequentist SD of mu under the site-exchangeable null y ~ N(0, I): exact under site permutation.
  hybrid_l  mu / q^(l/2), l in (0,1): l=0 ranks like `mean`, l=1 like z_freq.

Calibration (item 2 / branch B3); m_eff,k = (sum_s A_ks)^2 / sum_s A_ks^2 (effective number of substrates,
  computed analytically: sum_s A_ks = (Lambda^-1 g0)_k with g0 = W'M1).
  vif      compound-symmetry variance inflation: Var(mu_k) = q_k (1 + rho_e (m_eff,k - 1)); one global rho_e fitted
           on DEV role=null conditions only (`fit_rho_null`).
  empnull  per-condition empirical null: z_freq / s_c(m), s_c(m)^2 = s0_c^2 (1 + rho_c (m_eff - 1)), with
           (s0_c, rho_c) fitted label-free on the condition's own scored kinases by robust (MAD) variances in
           m_eff quartiles (Efron-style; most kinases are unperturbed in any one condition).
  Note: a per-condition scalar rescaling never changes M1 (within-condition ranking); only the m-dependent part does.

Kinase kernel (item 4): K = (I + rho_f F + rho_g G + rho_s S) / diag-normalised, S = label-free substrate-profile
  similarity of predicted libraries (cosine of top-k membership vectors, Hadamard power p; PSD by Schur).

Robustness (item 5): per-condition scale 'mad' (default) or 'sd'; optional winsorising |y| <= c.
Nothing here uses GT labels; no per-kinase / per-edge free parameters.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .bayes import MAD_K, kinase_kernel, report_mask

SOURCES = ("r2np", "r3", "cantley")


# ----------------------------------------------------------------------------- library design
@dataclass
class MultiDesign:
    kinases: pd.Index
    lib_rows: np.ndarray                 # row positions into the site index
    mats: dict = field(default_factory=dict)   # name -> csr (len(lib_rows) x p); 'cur' + predicted sources

    @property
    def names(self):
        return list(self.mats)

    def W(self, alphas: dict) -> sp.csr_matrix:
        W = self.mats["cur"].copy()
        for n, a in alphas.items():
            if a:
                W = W + a * self.mats[n]
        return W.tocsr()


def build_multi_design(site_index, kinases, lib_cur: pd.DataFrame, preds: dict) -> MultiDesign:
    """preds: name -> predicted library DataFrame (already restricted to top-k). Predicted edges that are curated
    edges are dropped from the predicted matrices (binary weight 1)."""
    site_index = pd.Index(site_index)
    kidx = pd.Index(list(kinases))

    def edges(df, binary):
        ki = kidx.get_indexer(df["kinase"]); si = site_index.get_indexer(df["site"])
        keep = (ki >= 0) & (si >= 0)
        w = np.ones(keep.sum()) if binary else df["weight"].to_numpy(float)[keep]
        return pd.DataFrame({"s": si[keep], "k": ki[keep], "w": w}).groupby(["s", "k"], as_index=False)["w"].max()

    E = {"cur": edges(lib_cur, False)}
    cur_codes = E["cur"]["s"].to_numpy(np.int64) * len(kidx) + E["cur"]["k"].to_numpy()
    for n, df in preds.items():
        e = edges(df, True)
        code = e["s"].to_numpy(np.int64) * len(kidx) + e["k"].to_numpy()
        E[n] = e[~np.isin(code, cur_codes)]
    rows = np.unique(np.concatenate([e["s"].to_numpy() for e in E.values()])).astype(int)
    mats = {n: sp.csr_matrix((e["w"].to_numpy(float), (np.searchsorted(rows, e["s"].to_numpy()), e["k"].to_numpy())),
                             shape=(len(rows), len(kidx))) for n, e in E.items()}
    return MultiDesign(kidx, rows, mats)


# ----------------------------------------------------------------------------- standardisation
def standardise(V: np.ndarray, scale: str = "mad", winsor: float | None = None):
    """y = (x - med_c)/s_c per column, then optional clip |y| <= winsor (in units of s_c).
    scale: 'mad' (1.4826 MAD; B2), 'sd' (root mean square of x - med), or 'wsd' (MAD-standardise, clip, then divide
    by the RMS of the clipped values so that the per-site variance is 1 — makes z_freq exact under site permutation).
    Returns Y (0 where missing), mask M (float) and the scale."""
    med = np.nanmedian(V, axis=0)
    R = V - med[None, :]
    if scale in ("mad", "wsd"):
        s = MAD_K * np.nanmedian(np.abs(R), axis=0)
    elif scale == "sd":
        s = np.sqrt(np.nanmean(R * R, axis=0))
    else:
        raise ValueError(scale)
    Y = R / s[None, :]
    if winsor is not None:
        Y = np.clip(Y, -winsor, winsor)
    if scale == "wsd":
        r = np.sqrt(np.nanmean(Y * Y, axis=0))
        Y = Y / r[None, :]
        s = s * r
    M = ~np.isnan(V)
    return np.where(M, Y, 0.0), M.astype(np.float64), s


# ----------------------------------------------------------------------------- sufficient statistics (blocks)
def block_stats(md: MultiDesign, Y: np.ndarray, M: np.ndarray, device="cuda", chunk=16, names=None):
    """Per condition (columns): G[(i,j)] = W_i' M W_j (B,p,p), h[i] = W_i' M y, g0[i] = W_i' M 1, n[i] = # measured
    edges per kinase. torch tensors (float32) on `device`."""
    import torch
    names = names or md.names
    rows = md.lib_rows
    Wd = {n: torch.as_tensor(md.mats[n].toarray(), dtype=torch.float32, device=device) for n in names}
    B = Y.shape[1]
    out_G = {(a, b): [] for i, a in enumerate(names) for b in names[i:]}
    out_h = {n: [] for n in names}; out_g = {n: [] for n in names}; out_n = {n: [] for n in names}
    for c0 in range(0, B, chunk):
        Ml = torch.as_tensor(M[rows, c0:c0 + chunk], dtype=torch.float32, device=device)       # S x b
        Yl = torch.as_tensor(Y[rows, c0:c0 + chunk], dtype=torch.float32, device=device)
        b = Ml.shape[1]
        for i, a in enumerate(names):
            A = Wd[a].T.unsqueeze(0) * Ml.T.unsqueeze(1)                                         # b x p x S
            for bn in names[i:]:
                out_G[(a, bn)].append(torch.bmm(A, Wd[bn].unsqueeze(0).expand(b, -1, -1)))
            out_h[a].append((Wd[a].T @ (Ml * Yl)).T)
            out_g[a].append((Wd[a].T @ Ml).T)
            out_n[a].append(((Wd[a] != 0).float().T @ Ml).T)
    cat = lambda d: {k: torch.cat(v) for k, v in d.items()}
    return dict(G=cat(out_G), h=cat(out_h), g0=cat(out_g), n=cat(out_n), names=list(names))


def combine(blocks: dict, alphas: dict, with_G=True):
    """alphas: name -> weight (cur fixed at 1). Returns G, h, g0 (torch) and n (numpy; edges with nonzero weight)."""
    w = {n: (1.0 if n == "cur" else float(alphas.get(n, 0.0))) for n in blocks["names"]}
    act = [n for n in blocks["names"] if w[n] != 0]
    G = None
    if with_G:
        for (a, b), M in blocks["G"].items():
            if a in act and b in act:
                t = w[a] * w[b] * M
                if a != b:
                    t = t + t.transpose(1, 2)
                G = t if G is None else G + t
    h = sum(w[n] * blocks["h"][n] for n in act)
    g0 = sum(w[n] * blocks["g0"][n] for n in act)
    n = sum(blocks["n"][n] for n in act).cpu().numpy()     # upper bound if sources overlap (only used for >=1 / bins)
    return G, h, g0, n


# ----------------------------------------------------------------------------- posterior with frequentist SD
def posterior2(G, h, g0, K: np.ndarray, tau: float, chunk=128):
    """Returns numpy (B x p): mu, sd_post = sqrt(diag Lambda^-1), q = diag(Lambda^-1 G Lambda^-1), suma = Lambda^-1 g0."""
    import torch
    dev = G.device
    P0 = torch.linalg.inv(torch.as_tensor(K, dtype=torch.float64, device=dev)) / tau ** 2
    out = {k: [] for k in ("mu", "sd_post", "q", "suma")}
    for i in range(0, G.shape[0], chunk):
        g = G[i:i + chunk].double()
        L = g + P0
        C = torch.linalg.cholesky(L)
        Li = torch.cholesky_inverse(C)
        mu = (Li @ h[i:i + chunk].double().unsqueeze(-1)).squeeze(-1)
        out["mu"].append(mu.cpu().numpy())
        out["sd_post"].append(torch.diagonal(Li, dim1=1, dim2=2).sqrt().cpu().numpy())
        out["q"].append(((Li @ g) * Li).sum(-1).clamp_min(0).cpu().numpy())
        out["suma"].append((Li @ g0[i:i + chunk].double().unsqueeze(-1)).squeeze(-1).cpu().numpy())
    return {k: np.concatenate(v) for k, v in out.items()}


def m_eff(q, suma):
    with np.errstate(all="ignore"):
        return np.where(q > 0, suma ** 2 / q, 0.0)


def output_stat(post: dict, kind: str, lam: float | None = None, rho_e: float = 0.0):
    """kind: mean | z_post | z_freq | hybrid (lam) | z_vif (rho_e)."""
    mu, q = post["mu"], post["q"]
    with np.errstate(all="ignore"):
        if kind == "mean":
            return mu
        if kind == "z_post":
            return mu / post["sd_post"]
        if kind == "z_freq":
            return mu / np.sqrt(q)
        if kind == "hybrid":
            return mu / q ** (lam / 2)
        if kind == "z_vif":
            me = m_eff(q, post["suma"])
            return mu / np.sqrt(q * (1 + rho_e * np.maximum(me - 1, 0)))
    raise ValueError(kind)


# ----------------------------------------------------------------------------- calibration
def empirical_null(Z: np.ndarray, ME: np.ndarray, scored: np.ndarray, n_q: int = 4, rho_max: float = 1.0,
                   fit_rho: bool = True):
    """Per-condition (row) empirical null. In each row, scored kinases are split into n_q quantile groups of
    u = m_eff - 1; robust variance v_g = (1.4826 MAD(z))^2 per group (centre 0) is regressed on mean u_g:
    v_g = s0^2 (1 + rho u_g), rho in [0, rho_max] (least squares on v, non-negative). Returns z / s(u), s0, rho."""
    out = np.full_like(Z, np.nan, dtype=float)
    s0s = np.full(Z.shape[0], np.nan); rhos = np.zeros(Z.shape[0])
    for c in range(Z.shape[0]):
        ok = scored[c] & np.isfinite(Z[c])
        if ok.sum() < 4 * n_q:
            continue
        z = Z[c, ok]; u = np.maximum(ME[c, ok] - 1, 0)
        edges = np.quantile(u, np.linspace(0, 1, n_q + 1))
        gi = np.clip(np.searchsorted(edges, u, side="right") - 1, 0, n_q - 1)
        ug, vg = [], []
        for g in range(n_q):
            m = gi == g
            if m.sum() >= 3:
                ug.append(u[m].mean()); vg.append((MAD_K * np.median(np.abs(z[m]))) ** 2)
        ug, vg = np.array(ug), np.array(vg)
        if fit_rho and len(ug) >= 2 and np.ptp(ug) > 0:
            # v = a + b u with a = s0^2, b = s0^2 rho ; LS then clamp
            Xd = np.c_[np.ones_like(ug), ug]
            a, b = np.linalg.lstsq(Xd, vg, rcond=None)[0]
            if a <= 0:
                a = max(np.min(vg), 1e-6)
            rho = float(np.clip(b / a, 0, rho_max))
            if rho == 0 or rho == rho_max:
                a = float(np.median(vg / (1 + rho * ug)))
        else:
            rho, a = 0.0, float(np.median(vg)) if len(vg) else 1.0
        s0s[c], rhos[c] = np.sqrt(a), rho
        uu = np.maximum(ME[c] - 1, 0)
        out[c] = Z[c] / np.sqrt(a * (1 + rho * uu))
    return out, s0s, rhos


def fit_rho_null(Zf: np.ndarray, ME: np.ndarray, scored: np.ndarray, grid=None):
    """Global compound-symmetry rho_e from null conditions (rows of Zf = z_freq, scored mask): choose rho_e making
    the pooled robust variance of z_freq/sqrt(VIF) equal across m_eff quartiles (min squared log-ratio to the
    lowest quartile's). Returns rho_e and the table."""
    grid = grid if grid is not None else [0, 0.002, 0.005, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3]
    ok = scored & np.isfinite(Zf)
    z = Zf[ok]; u = np.maximum(ME[ok] - 1, 0)
    qs = np.quantile(u, [0.25, 0.5, 0.75])
    gi = np.searchsorted(qs, u, side="right")
    groups = [g for g in range(4) if (gi == g).sum() >= 10]
    tab = []
    for r in grid:
        zz = z / np.sqrt(1 + r * u)
        v = np.array([(MAD_K * np.median(np.abs(zz[gi == g]))) ** 2 for g in groups])
        tab.append((r, float(np.sum(np.log(v / v.mean()) ** 2)), v.tolist()))
    best = min(tab, key=lambda t: t[1])
    return best[0], tab


# ----------------------------------------------------------------------------- label-free kinase similarity
def topk_similarity(lib: pd.DataFrame, kinases, k: int = 500, power: int = 2) -> np.ndarray:
    """Cosine similarity of binary top-k membership vectors (kinases x sites), Hadamard power `power`, unit diagonal;
    kinases without predictions get a zero row/col (diag 1 kept by the caller's I term)."""
    kidx = pd.Index(list(kinases))
    l = lib[lib["rank"] <= k] if "rank" in lib.columns else lib
    ki = kidx.get_indexer(l["kinase"]); keep = ki >= 0
    si, _ = pd.factorize(l["site"][keep])
    A = sp.csr_matrix((np.ones(keep.sum()), (ki[keep], si)), shape=(len(kidx), si.max() + 1))
    Gm = (A @ A.T).toarray()
    d = np.sqrt(np.diag(Gm))
    with np.errstate(all="ignore"):
        S = np.where(np.outer(d, d) > 0, Gm / np.outer(d, d), 0.0)
    np.fill_diagonal(S, 0.0)
    return S ** power


def kernel2(universe, kinases, rho_f, rho_g, rho_s=0.0, S=None):
    """(I + rho_f F + rho_g G + rho_s S) diag-normalised. S: off-diagonal PSD-completable similarity (0 diag)."""
    K0 = kinase_kernel(universe, kinases, rho_f, rho_g)
    if not rho_s or S is None:
        return K0
    # undo normalisation of K0, add S (+ its diagonal so S+I is a Gram matrix), renormalise
    U = universe.set_index("kinase").reindex(list(kinases))
    fam = U["manning_family"].fillna("").astype(str).to_numpy()
    grp = U["manning_group"].fillna("").astype(str).to_numpy()
    nf = (fam != "").astype(float); ng = ((grp != "") & ~np.isin(grp, ["Other", "Atypical"])).astype(float)
    raw_d = 1 + rho_f * nf + rho_g * ng
    Kraw = K0 * np.sqrt(np.outer(raw_d, raw_d))
    has = (S.sum(1) > 0).astype(float)
    Kraw = Kraw + rho_s * (S + np.diag(has))
    d = np.sqrt(np.diag(Kraw))
    return Kraw / np.outer(d, d)


# ----------------------------------------------------------------------------- one-call scorer
@dataclass
class Bayes2Config:
    tau: float = 0.125
    rho_f: float = 3.0
    rho_g: float = 0.25
    rho_s: float = 0.0
    alphas: dict = field(default_factory=lambda: {"r2np": 1.0})
    output: str = "z_freq"         # mean | z_post | z_freq | hybrid | z_vif | empnull
    lam: float = 0.5
    rho_e: float = 0.0
    scale: str = "mad"
    winsor: float | None = None
    borrow: bool = True


def score(X: pd.DataFrame, md: MultiDesign, universe: pd.DataFrame, cfg: Bayes2Config, S=None, device="cuda",
          return_all=False):
    """Score conditions (columns of X; index = site index of `md`). Returns conditions x kinases (NaN = unscored)."""
    Y, M, _ = standardise(X.to_numpy(np.float64), cfg.scale, cfg.winsor)
    names = ["cur"] + [n for n, a in cfg.alphas.items() if a]
    bl = block_stats(md, Y, M, device=device, names=names)
    G, h, g0, n = combine(bl, cfg.alphas)
    K = kernel2(universe, md.kinases, cfg.rho_f, cfg.rho_g, cfg.rho_s, S)
    post = posterior2(G, h, g0, K, cfg.tau)
    rep = report_mask(n, K, 1, cfg.borrow)
    if cfg.output == "empnull":
        zf = output_stat(post, "z_freq")
        out, _, _ = empirical_null(zf, m_eff(post["q"], post["suma"]), rep)
    else:
        out = output_stat(post, cfg.output, cfg.lam, cfg.rho_e)
    Sx = pd.DataFrame(np.where(rep, out, np.nan), index=X.columns, columns=md.kinases)
    if return_all:
        return Sx, dict(post=post, n=n, report=rep, K=K)
    return Sx


# ----------------------------------------------------------------------------- CPU single-condition path
class SingleSample2:
    """CPU scorer for one condition at a time (runtime / CLI path). Returns score (cfg.output), calibrated z
    (empirical null on z_freq), mu, q, n."""

    def __init__(self, md: MultiDesign, universe: pd.DataFrame, cfg: Bayes2Config, S=None):
        import scipy.linalg as sla
        self.sla = sla
        self.cfg = cfg
        self.W = md.W(cfg.alphas).tocsr()
        self.Wb = self.W.copy(); self.Wb.data[:] = 1.0
        self.K = kernel2(universe, md.kinases, cfg.rho_f, cfg.rho_g, cfg.rho_s, S)
        self.P0 = np.linalg.inv(self.K) / cfg.tau ** 2
        self.rows = md.lib_rows
        self.rel = ((self.K > 1e-12) & ~np.eye(len(self.K), dtype=bool)).astype(float)

    def score_one(self, x: np.ndarray):
        """x: values at ALL sites of the site index (NaN = missing)."""
        sla = self.sla
        Y, M, _ = standardise(x[:, None].astype(np.float64), self.cfg.scale, self.cfg.winsor)
        xl = Y[self.rows, 0]; ml = M[self.rows, 0] > 0
        m = np.flatnonzero(ml)
        Wm = self.W[m]
        G = (Wm.T @ Wm).toarray()
        h = Wm.T @ xl[m]
        g0 = np.asarray(Wm.sum(0)).ravel()
        C = sla.cholesky(G + self.P0, lower=True, check_finite=False)
        Ci = sla.solve_triangular(C, np.eye(len(G)), lower=True, check_finite=False)
        Li = Ci.T @ Ci
        mu = Li @ h
        q = np.maximum(((Li @ G) * Li).sum(1), 0)
        post = dict(mu=mu[None], q=q[None], suma=(Li @ g0)[None], sd_post=np.sqrt(np.diag(Li))[None])
        n = np.asarray(self.Wb[m].sum(0)).ravel()
        has = n >= 1
        rep = has | ((has.astype(float) @ self.rel) > 0) if self.cfg.borrow else has
        zf = output_stat(post, "z_freq")
        zc, _, _ = empirical_null(zf, m_eff(post["q"], post["suma"]), rep[None])
        sc = zc if self.cfg.output == "empnull" else output_stat(post, self.cfg.output, self.cfg.lam, self.cfg.rho_e)
        return np.where(rep, sc[0], np.nan), np.where(rep, zc[0], np.nan), mu, q, n
