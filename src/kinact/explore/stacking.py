"""Branch B6: stacking / ensembling of kinase-activity score matrices (DEV only).

All component matrices are cond x kinase (NaN = unscored). Combination rules are odd functions of the input
(sign-symmetric: flipping the input vector flips every component and therefore the stacked score), carry no
kinase-specific parameters and no additive offsets (an additive missingness indicator would not be sign-symmetric and
would act as a label-sign prior).

Pipeline
  1. z-normalise each component j by a single label-free scale sigma_j = SD of the component on site-permuted DEV GT
     conditions (so every component is ~N(0,1) under the permutation null); NaN -> 0 ("no evidence" on z-scale).
     Alternative: per-condition rank-normal transform (Phi^-1 of mid-rank among the condition's scored kinases).
  2. stacked s = sum_j w_j z_j (unscored by every component -> NaN).
  3. recalibration: s / sqrt(w_o' R_o w_o), R = correlation of normalised components on site-permuted conditions
     (pairwise complete), o = components observed in that cell. Label-free; 0 extra fitted parameters.
  4. Learned w (linear stacker): smooth-AUROC pairwise loss sigmoid(-(t_k - t_j)/(T sd_c)) within training M1-sel
     conditions over scored evaluable negatives, pair weights reproducing the harness macro
     (pairs -> (kinase, group) -> kinase), ridge penalty lam * ||w - w0||^2 shrinking toward w0 = core-only.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.special import ndtri


def null_scale(P: np.ndarray) -> float:
    """SD (about 0) of a component on null/permuted data."""
    v = P[np.isfinite(P)]
    return float(np.sqrt(np.mean(v * v)))


def rank_normal(S: np.ndarray) -> np.ndarray:
    """Per-row Phi^-1((rank - 0.5)/n) over finite entries; NaN kept."""
    out = np.full(S.shape, np.nan)
    for i in range(S.shape[0]):
        ok = np.isfinite(S[i])
        n = ok.sum()
        if n < 2:
            continue
        r = pd.Series(S[i, ok]).rank(method="average").to_numpy()
        out[i, ok] = ndtri((r - 0.5) / n)
    return out


def null_corr(Ps: list[np.ndarray]) -> np.ndarray:
    """Pairwise-complete correlation (about 0, i.e. E[z_i z_j]/sqrt(E z_i^2 E z_j^2)) of normalised null components."""
    J = len(Ps)
    R = np.eye(J)
    for a in range(J):
        for b in range(a + 1, J):
            ok = np.isfinite(Ps[a]) & np.isfinite(Ps[b])
            if ok.sum() < 50:
                continue
            x, y = Ps[a][ok], Ps[b][ok]
            R[a, b] = R[b, a] = float((x * y).mean() / np.sqrt((x * x).mean() * (y * y).mean()))
    return R


def combine(Zs: list[np.ndarray], w: np.ndarray, R: np.ndarray | None = None) -> np.ndarray:
    """Zs: normalised components (NaN = unscored). Returns stacked score, recalibrated by the per-cell null SD
    sqrt(w_o' R_o w_o) when R is given. Cells with no observed component (or zero null SD) -> NaN."""
    Z = np.stack(Zs, -1)                       # B x p x J
    O = np.isfinite(Z)
    s = np.where(O, Z, 0.0) @ w
    if R is not None:
        Wo = O * w                              # B x p x J
        var = np.einsum("bpj,jk,bpk->bp", Wo, R, Wo)
        with np.errstate(all="ignore"):
            s = s / np.sqrt(var)
        s[~(var > 1e-12)] = np.nan
    s[~O.any(-1)] = np.nan
    return s


# ------------------------------------------------------------------ learned linear stacker
def pair_data(ev, conds: list, train_conds, Zs: list[np.ndarray]):
    """Training targets (M1-sel, scored by >=1 component) for conditions `train_conds`.
    Returns dict with per-target tensors (numpy): X (n_t x p x J, NaN->0), O (any component observed, n_t x p),
    ki, sign, w (macro weights), neg (n_t x p, evaluable negatives & observed), nN (all evaluable negatives)."""
    pos = {c: i for i, c in enumerate(conds)}
    Z = np.stack(Zs, -1)
    obs = np.isfinite(Z).any(-1)
    tg = ev.targets
    tg = tg[tg["conditionID"].isin(set(train_conds)) & tg["in_M1sel"]].copy()
    tg["col"] = tg["conditionID"].map(pos).astype(int)
    tg = tg[obs[tg["col"].to_numpy(), tg["ki"].to_numpy()]]
    npg = tg.groupby(["kinase", "group"])["ki"].transform("count")
    ngr = tg.groupby("kinase")["group"].transform("nunique")
    tg["w"] = 1.0 / npg / ngr / tg["kinase"].nunique()
    col = tg["col"].to_numpy()
    neg = ev.neg_mask[tg["ci"].to_numpy()] & obs[col]
    nN = ev.neg_mask[tg["ci"].to_numpy()].sum(1)
    return dict(X=np.nan_to_num(Z[col]), O=obs[col], Ob=np.isfinite(Z[col]), ki=tg["ki"].to_numpy(), sign=tg["sign"].to_numpy(float),
                w=tg["w"].to_numpy(), neg=neg, nN=nN)


def fit_weights(pd_: dict, w0: np.ndarray, lam: float, R: np.ndarray | None = None, T: float = 0.5,
                max_iter: int = 200, device="cpu"):
    """Minimise sum_t w_t mean_j sigmoid(-(t_k - t_j)/(T sd_c)) + lam ||w - w0||^2 (L-BFGS, torch).
    sd_c = SD of the stacked score over the condition's observed kinases (makes the loss scale-free; the penalty
    fixes the scale). With R, the stacked score is recalibrated per cell exactly as in `combine` (trained on the
    output that is evaluated)."""
    import torch
    t = lambda a, dt=torch.float64: torch.as_tensor(a, dtype=dt, device=device)
    X, O, neg = t(pd_["X"]), t(pd_["O"]), t(pd_["neg"])
    ki = torch.as_tensor(pd_["ki"], device=device)
    sg, wt, nN = t(pd_["sign"]), t(pd_["w"]), t(pd_["nN"])
    W0 = t(w0)
    Ob = t(pd_["Ob"]) if R is not None else None
    Rt = t(R) if R is not None else None
    w = W0.clone().requires_grad_(True)
    opt = torch.optim.LBFGS([w], lr=1, max_iter=max_iter, history_size=20, line_search_fn="strong_wolfe")

    def loss_fn():
        s = X @ w                                            # n_t x p
        if Rt is not None:
            Wo = Ob * w
            var = torch.einsum("tpj,jk,tpk->tp", Wo, Rt, Wo)
            s = torch.where(O > 0, s / torch.sqrt(var.clamp_min(1e-8)), torch.zeros_like(s))
        cnt = O.sum(1).clamp_min(2)
        mu = (s * O).sum(1) / cnt
        sd = torch.sqrt((((s - mu[:, None]) ** 2) * O).sum(1) / (cnt - 1)).clamp_min(1e-8)
        Zs = s * sg[:, None]
        tk = Zs.gather(1, ki[:, None])
        d = (tk - Zs) / (T * sd[:, None])
        L = torch.sigmoid(-d) * neg
        return (wt * L.sum(1) / nN).sum() + lam * ((w - W0) ** 2).sum()

    def closure():
        opt.zero_grad()
        l = loss_fn()
        l.backward()
        return l

    opt.step(closure)
    return w.detach().cpu().numpy()
