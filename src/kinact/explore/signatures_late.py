"""Branch B7b — late fusion of perturbation-response signatures (DESIGN §12 C10).

Pre-declared design: results/b7b/b7b_report.md, "PRE-DECLARED DESIGN". Signatures are the B7 fits
(kinact.explore.signatures). They are used OUTSIDE the Bayesian kernel:

  z_sig(k,c) = sum_s w_ks clip(y_cs) / sqrt(sum_s w_ks^2)   over k's top-n_top signature sites measured in c (>= 5)
  z_sig     /= rho_sig  (one scalar per signature fit = RMS of z_sig over the DEV null conditions)
  z          = (z_core + a z_sig) / sqrt(1 + a^2)   if both are defined; else whichever is defined.

Odd in y and with no intercept, so there is no kinase-level offset (C10 ii). Kinases without a signature keep z_core.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp

CLIP = 4.0
MIN_SITES = 5


def sig_project(Ws: sp.csr_matrix, Y: np.ndarray, M: np.ndarray, min_sites: int = MIN_SITES) -> np.ndarray:
    """Ws: sites x p signed signature weights (full site index). Y: sites x B standardised data (0 where missing),
    M: sites x B mask. Returns B x p raw z_sig (NaN for kinases without a signature / < min_sites measured)."""
    R = np.unique(Ws.nonzero()[0])
    p, B = Ws.shape[1], Y.shape[1]
    if len(R) == 0:
        return np.full((B, p), np.nan)
    W = Ws[R].tocsc().astype(np.float64)
    Yr = np.clip(Y[R].astype(np.float64), -CLIP, CLIP) * M[R]
    Mr = M[R].astype(np.float64)
    num = np.asarray(W.T @ Yr)                       # p x B
    den = np.sqrt(np.asarray(W.multiply(W).T @ Mr))
    cnt = np.asarray((W != 0).astype(np.float64).T @ Mr)
    with np.errstate(all="ignore"):
        z = np.where((cnt >= min_sites) & (den > 0), num / den, np.nan)
    return z.T


def rho_sig(z_null: np.ndarray) -> float:
    """Scalar standardiser: RMS of the raw z_sig over the null conditions (all defined cells)."""
    v = z_null[np.isfinite(z_null)]
    return float(np.sqrt(np.mean(v * v))) if len(v) else 1.0


def fuse(z_core: np.ndarray, z_sig: np.ndarray, alpha: float) -> np.ndarray:
    zc, zs = np.asarray(z_core, np.float64), np.asarray(z_sig, np.float64)
    both = np.isfinite(zc) & np.isfinite(zs)
    out = np.where(np.isfinite(zc), zc, zs)
    return np.where(both, (zc + alpha * zs) / np.sqrt(1.0 + alpha ** 2), out)
