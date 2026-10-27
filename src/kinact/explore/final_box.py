"""Branch B8 — final development box (pre-declared in results/b8/b8_report.md).

E1  perm_scales / fuse_calibrated : per-signature-kinase permutation-null scale of the fused numerator
    F_k = z_core + a z_sig (v_k = mean F_k^2 over site-permuted TRAINING conditions); z = F_k / sqrt(v_k).
E2  fit_factors / project          : top-r site-space singular vectors of the allowed (label-free) DEV matrix;
    each input condition is projected off them (least squares on its measured factor sites), then rescaled to RMS 1.
E3  bag_z                           : average of the top-3 inner-CV core configs' z, scaled by the pooled uncentred
    correlation (permuted training conditions) so the null scale stays 1.
Everything here is label-free; labels enter only through the caller's nested selection.
"""
from __future__ import annotations

import numpy as np

CLIP = 4.0
MIN_PERM_CELLS = 30
MIN_FACTOR_SITES = 50
SITE_FRAC = 0.20


# ----------------------------------------------------------------------------- E1
def perm_scales(zc_perm: np.ndarray, zs_perm: np.ndarray, alpha: float, sig_cols: np.ndarray):
    """zc_perm, zs_perm: (n_perm_conds x p) core z and (rho_sig-scaled) z_sig on permuted training conditions.
    sig_cols: bool (p,) kinases with a signature. Returns (v_both, v_sigonly) arrays (p,), NaN for non-signature."""
    both = np.isfinite(zc_perm) & np.isfinite(zs_perm)
    F = np.where(both, zc_perm + alpha * np.nan_to_num(zs_perm), 0.0)
    n = both.sum(0)
    with np.errstate(all="ignore"):
        v = (F ** 2).sum(0) / n
    pooled = (F[:, sig_cols] ** 2).sum() / max(n[sig_cols].sum(), 1)
    v = np.where(n >= MIN_PERM_CELLS, v, pooled)
    so = np.isfinite(zs_perm)
    S2 = np.where(so, np.nan_to_num(zs_perm) ** 2, 0.0)
    ns = so.sum(0)
    with np.errstate(all="ignore"):
        vs = S2.sum(0) / ns
    pooled_s = S2[:, sig_cols].sum() / max(ns[sig_cols].sum(), 1)
    vs = np.where(ns >= MIN_PERM_CELLS, vs, pooled_s)
    return np.where(sig_cols, v, np.nan), np.where(sig_cols, vs, np.nan)


def fuse_calibrated(z_core, z_sig, alpha, v_both, v_sig):
    zc, zs = np.asarray(z_core, np.float64), np.asarray(z_sig, np.float64)
    both = np.isfinite(zc) & np.isfinite(zs)
    out = np.where(np.isfinite(zc), zc, zs / np.sqrt(v_sig)[None, :])
    with np.errstate(all="ignore"):
        fused = (zc + alpha * zs) / np.sqrt(v_both)[None, :]
    return np.where(both, fused, out)


# ----------------------------------------------------------------------------- E2
def fit_factors(Y: np.ndarray, M: np.ndarray, cols: np.ndarray, r_max: int = 3, site_frac: float = SITE_FRAC,
                device="cuda"):
    """Y, M: sites x conditions (wsd standardised, 0 where missing; mask). cols: allowed condition columns.
    Returns (site_rows, V) with V (len(site_rows) x r_max) orthonormal right singular vectors (site space)."""
    import torch
    Ma = M[:, cols]
    rows = np.flatnonzero(Ma.mean(1) >= site_frac)
    A = np.clip(Y[np.ix_(rows, cols)], -CLIP, CLIP) * Ma[rows]          # sites x n, missing -> 0
    At = torch.as_tensor(A.T, dtype=torch.float64, device=device)       # n x sites
    _, s, Vh = torch.linalg.svd(At, full_matrices=False)
    return rows, Vh[:r_max].T.cpu().numpy(), s[:r_max + 2].cpu().numpy()


def project(Y: np.ndarray, M: np.ndarray, rows: np.ndarray, V: np.ndarray, r: int):
    """Project columns of Y (sites x B, 0 where missing) off the first r factors. Returns Y' (same shape, 0 where
    missing), rescaled per column to RMS 1 over measured sites; columns with < MIN_FACTOR_SITES measured factor sites
    are returned unchanged."""
    Yp = Y.astype(np.float64, copy=True)
    Vr = V[:, :r]
    for b in range(Y.shape[1]):
        m = M[rows, b] > 0
        if m.sum() < MIN_FACTOR_SITES:
            continue
        Vm = Vr[m]
        yv = Y[rows[m], b].astype(np.float64)
        beta, *_ = np.linalg.lstsq(Vm, np.clip(yv, -CLIP, CLIP), rcond=None)
        Yp[rows[m], b] = yv - Vm @ beta
        mm = M[:, b] > 0
        rms = np.sqrt(np.mean(Yp[mm, b] ** 2))
        Yp[:, b] = np.where(mm, Yp[:, b] / rms, 0.0)
    return Yp


# ----------------------------------------------------------------------------- E3
def bag_z(Z_list, Zperm_list):
    """Z_list: list of (B x p) z of the bagged configs; Zperm_list: same configs on permuted training conditions.
    z_bag = sum_i z_i / sqrt(sum_ij c_ij), c_ij = uncentred correlation over permuted cells defined for all."""
    P = np.stack(Zperm_list)                                            # n x C x p
    ok = np.isfinite(P).all(0)
    v = [P[i][ok] for i in range(len(P))]
    c = np.array([[np.mean(a * b) / np.sqrt(np.mean(a * a) * np.mean(b * b)) for b in v] for a in v])
    Z = np.stack(Z_list)
    S = Z.sum(0)                                                        # NaN if any config NaN
    # a kinase reported by some configs only: fall back to the mean of the defined z (rare; report sets differ by K)
    n_def = np.isfinite(Z).sum(0)
    S = np.where(n_def == len(Z), S, np.where(n_def > 0, np.nansum(Z, 0) / np.maximum(n_def, 1) * len(Z), np.nan))
    return S / np.sqrt(c.sum()), c
