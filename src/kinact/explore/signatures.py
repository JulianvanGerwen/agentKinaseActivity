"""Branch B7 — perturbation-derived kinase response signatures as extra edges in the Bayesian core (DESIGN §12 C10).

PRE-DECLARED DESIGN (written before any B7 result was computed)
----------------------------------------------------------------
Signature fit (one call per training set; `fit_signatures`):
  * Training conditions T: role == 'gt', qc_status == 'ok', split in the TRAINING folds only. The caller passes the
    evaluated conditions; the fit asserts that no evaluated condition, and no condition sharing a drug group with an
    evaluated condition, is in T (C10 i / "drug-held-out").
  * Input y_c = per-condition wsd standardisation (bayes2.standardise, no winsor) of the log2FC, clipped at |y| <= 4
    (fit only; robustness).
  * Positives: GT direct pairs, confidence high/medium, sign s in {-1,+1}. Variant 'clean' uses only conditions with
    n_targets <= 10 (PRIMARY, pre-declared: same label quality as M1-sel, fewer co-targets sharing a response);
    variant 'all' uses every training GT condition.
  * Sign alignment and group averaging: v_cs = s_ck * y_cs; averaged within drug group first (dose series of one drug
    are not independent), then across groups: m_ks = mean_g v_gks, n_ks = # groups with site s measured.
  * Generic-response contrast (PRIMARY: on): subtract the generic perturbation response
    r_s = mean over ALL training groups of the group mean of s_c * y_cs (s_c = majority sign of c's direct positives).
    It removes responses shared by any perturbation (stress, global pY loss, ...), i.e. popularity-like information.
  * Site noise sigma_s = SD over training groups of (group-level aligned vector - r_s); sites measured in < 5 training
    groups are not eligible; sigma floored at its 10% quantile.
  * Kinase eligible iff >= 2 training drug groups as a direct positive (fixed n_min = 2); site eligible for k iff
    n_ks >= 2.
  * Shrinkage (HP kappa, pseudo-groups): w_ks = (m_ks / sigma_s) * n_ks / (n_ks + kappa). Edges = top n_top sites by
    |w_ks| (HP n_top). Weights are SIGNED (a site may fall when k is activated).
Scoring: W = W_cur + 1*W_r2np@100 + alpha_sig * W_sig (HP alpha_sig); same posterior / z_vif as core-v1 (bayes2);
  rho_e refitted on the DEV role=null conditions for each design. Core HPs (tau, rho_f, rho_g) per outer fold = the
  core-v1 nested picks (so alpha_sig = 0 reproduces core-v1 exactly).
C10 ii: signatures carry no kinase offset — they enter only as design-matrix columns multiplying the evaluated
  condition's data; score(-x) = -score(x) (median-centring and the linear posterior are odd in x).
HPs (<= 3, nested CV): alpha_sig in {0, .25, .5, 1, 2}, n_top in {20, 50, 100, 200}, kappa in {0, 1, 4}.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import scipy.sparse as sp

CLIP = 4.0
MIN_SITE_GROUPS = 5
N_MIN_GROUPS = 2


class SignatureLeak(AssertionError):
    pass


# ----------------------------------------------------------------------------- training sets
def training_conditions(C: pd.DataFrame, meta: pd.DataFrame, train_folds, eval_conds) -> list:
    """GT conditions usable for signature fitting: role gt, qc ok, split in train_folds, in the DEV cache, and not
    sharing a drug group with any evaluated condition. Raises SignatureLeak if a train fold overlaps eval conds."""
    train_folds = set(train_folds)
    Cg = C[(C["role"] == "gt") & C["split"].isin(train_folds)]
    Cg = Cg[Cg["conditionID"].isin(meta.index)]
    Cg = Cg[meta.loc[Cg["conditionID"], "qc_status"].to_numpy() == "ok"]
    grp_eval = set(meta.reindex(list(eval_conds))["group_id"].dropna())
    tr = [c for c, g in zip(Cg["conditionID"], Cg["group_id"]) if g not in grp_eval]
    if set(tr) & set(eval_conds):
        raise SignatureLeak("evaluated condition in training set")
    return tr


def positives(P: pd.DataFrame, C: pd.DataFrame, train_conds, nt_max: int | None = None) -> pd.DataFrame:
    """Direct high/medium positives in train_conds (optionally n_targets <= nt_max): conditionID, kinase, sign, group."""
    Ci = C.set_index("conditionID")
    p = P[(P["cls"] == "direct") & P["confidence"].isin(["high", "medium"]) & P["sign"].isin([-1, 1])
          & P["conditionID"].isin(set(train_conds))].copy()
    p["group"] = Ci.loc[p["conditionID"], "group_id"].to_numpy()
    p["n_targets"] = Ci.loc[p["conditionID"], "n_targets"].to_numpy()
    if nt_max is not None:
        p = p[p["n_targets"] <= nt_max]
    return p[["conditionID", "kinase", "sign", "group"]].reset_index(drop=True)


def _group_mean(Y: np.ndarray, cols: np.ndarray, signs: np.ndarray, groups: np.ndarray):
    """Group-level nanmean of sign-aligned columns. Returns (n_groups x S) array, group labels."""
    ug = pd.unique(groups)
    out = np.full((len(ug), Y.shape[0]), np.nan, dtype=np.float32)
    for i, g in enumerate(ug):
        m = groups == g
        V = Y[:, cols[m]] * signs[m][None, :].astype(np.float32)
        with np.errstate(all="ignore"), warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            out[i] = np.nanmean(V, axis=1)
    return out, ug


@dataclass
class Signatures:
    """kinase -> (site row indices into the site index, signed weights), sorted by |w| descending (up to n_max)."""
    sites: dict = field(default_factory=dict)
    weights: dict = field(default_factory=dict)
    n_groups: dict = field(default_factory=dict)
    train_conds: list = field(default_factory=list)
    eligible_sites: np.ndarray | None = None
    generic: np.ndarray | None = None

    def matrix(self, n_sites: int, kinases, n_top: int, alpha: float = 1.0) -> sp.csr_matrix:
        kidx = pd.Index(list(kinases))
        r, c, v = [], [], []
        for k, s in self.sites.items():
            j = kidx.get_loc(k)
            r.append(s[:n_top]); c.append(np.full(min(n_top, len(s)), j)); v.append(self.weights[k][:n_top])
        if not r:
            return sp.csr_matrix((n_sites, len(kidx)))
        return sp.csr_matrix((alpha * np.concatenate(v), (np.concatenate(r), np.concatenate(c))),
                             shape=(n_sites, len(kidx)))

    def shuffled(self, seed: int = 0) -> "Signatures":
        """Control: within each signature, permute site labels — weights kept, sites replaced by a random draw (without
        replacement) from the eligible-site pool. Destroys site identity, keeps weight distribution and size."""
        rng = np.random.default_rng(seed)
        out = Signatures(n_groups=dict(self.n_groups), train_conds=self.train_conds,
                         eligible_sites=self.eligible_sites, generic=self.generic)
        for k, s in self.sites.items():
            out.sites[k] = rng.choice(self.eligible_sites, size=len(s), replace=False)
            out.weights[k] = self.weights[k].copy()
        return out

    def generic_control(self) -> "Signatures":
        """Control: every kinase that has a signature gets the SAME generic-response signature (top sites of the
        standardised generic response r_s / sigma_s), with its own weight norm. Measures popularity carried by the
        generic perturbation response."""
        out = Signatures(n_groups=dict(self.n_groups), train_conds=self.train_conds,
                         eligible_sites=self.eligible_sites, generic=self.generic)
        g = self.generic
        order = np.argsort(-np.abs(g))
        for k, s in self.sites.items():
            n = len(s)
            sel = order[:n]
            w = g[sel]
            nrm = np.sqrt(np.sum(self.weights[k] ** 2)) / max(np.sqrt(np.sum(w ** 2)), 1e-12)
            out.sites[k] = sel
            out.weights[k] = w * nrm
        return out


def fit_signatures(Y: np.ndarray, cond_index, train_conds, pos: pd.DataFrame, cond_group: pd.Series,
                   eval_conds=(), kappa=1.0, n_max: int = 200, contrast: bool = True,
                   site_idx: np.ndarray | None = None):
    """Y: standardised sites x conditions array (NaN = missing), columns = cond_index. Only columns in train_conds are
    read. pos: positives() table restricted to train_conds. cond_group: conditionID -> group_id (for leakage checks).
    Returns Signatures with up to n_max sites per kinase (indices into Y's rows, or into site_idx if given); if kappa
    is a list/tuple, a dict kappa -> Signatures (one pass over the data)."""
    kappas = list(kappa) if isinstance(kappa, (list, tuple)) else [kappa]
    cond_index = pd.Index(cond_index)
    train_conds = list(train_conds)
    eval_conds = list(eval_conds)
    if set(train_conds) & set(eval_conds):
        raise SignatureLeak("evaluated condition in training set")
    if set(cond_group.reindex(train_conds)) & set(cond_group.reindex(eval_conds).dropna()):
        raise SignatureLeak("training condition shares a drug group with an evaluated condition")
    if not set(pos["conditionID"]) <= set(train_conds):
        raise SignatureLeak("positive table contains non-training conditions")
    ci = cond_index.get_indexer(train_conds)
    assert (ci >= 0).all()
    Yt = np.clip(Y[:, ci].astype(np.float32), -CLIP, CLIP)
    tpos = pd.Index(train_conds)
    # condition-level majority sign (for the generic response)
    sc = pos.groupby("conditionID")["sign"].sum()
    sc = np.sign(sc).replace(0, 1)
    gen_conds = [c for c in train_conds if c in sc.index]
    gcols = tpos.get_indexer(gen_conds)
    GM, _ = _group_mean(Yt, gcols, sc.loc[gen_conds].to_numpy(), cond_group.loc[gen_conds].to_numpy())
    with np.errstate(all="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        n_site = np.sum(~np.isnan(GM), axis=0)
        generic = np.nanmean(GM, axis=0)
        sig = np.nanstd(GM - generic[None, :], axis=0, ddof=1)
    elig = (n_site >= MIN_SITE_GROUPS) & np.isfinite(sig) & (sig > 0)
    floor = np.quantile(sig[elig], 0.10)
    sig = np.where(elig, np.maximum(sig, floor), np.nan)
    generic = np.where(elig, generic, 0.0)
    outs = {kp: Signatures(train_conds=train_conds, eligible_sites=np.flatnonzero(elig),
                           generic=np.where(elig, generic / sig, 0.0)) for kp in kappas}
    for k, pk in pos.groupby("kinase"):
        ng = pk["group"].nunique()
        if ng < N_MIN_GROUPS:
            continue
        D, _ = _group_mean(Yt, tpos.get_indexer(pk["conditionID"]), pk["sign"].to_numpy(), pk["group"].to_numpy())
        if contrast:
            D = D - generic[None, :]
        with np.errstate(all="ignore"), warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            n = np.sum(~np.isnan(D), axis=0)
            m = np.nanmean(D, axis=0)
        for kp, out in outs.items():
            w = (m / sig) * n / (n + kp)
            ok = elig & (n >= 2) & np.isfinite(w)
            idx = np.flatnonzero(ok)
            order = idx[np.argsort(-np.abs(w[idx]), kind="stable")][:n_max]
            out.sites[k] = order if site_idx is None else site_idx[order]
            out.weights[k] = w[order].astype(np.float64)
            out.n_groups[k] = int(ng)
    return outs if isinstance(kappa, (list, tuple)) else outs[kappas[0]]


# ----------------------------------------------------------------------------- scoring blocks
def sig_blocks(Wb: sp.csr_matrix, Ws: sp.csr_matrix, Y: np.ndarray, M: np.ndarray, device="cuda", chunk=16):
    """Cross/self blocks of the signature columns for the conditions in Y/M (sites x B, full site index).
    Wb: base design (sites x p, full site index); Ws: signature design (sites x p). Returns torch dict
    G_bs (B,p,p) = Wb' M Ws, G_ss = Ws' M Ws, h_s = Ws' M y, g0_s = Ws' M 1, n_s (numpy) = # measured sig edges."""
    import torch
    R = np.unique(Ws.nonzero()[0])
    p = Ws.shape[1]
    B = Y.shape[1]
    if len(R) == 0:
        z = torch.zeros((B, p, p), device=device)
        zv = torch.zeros((B, p), device=device)
        return dict(G_bs=z, G_ss=z.clone(), h_s=zv, g0_s=zv.clone(), n_s=np.zeros((B, p)))
    Wbd = torch.as_tensor(Wb[R].toarray(), dtype=torch.float32, device=device)
    Wsd = torch.as_tensor(Ws[R].toarray(), dtype=torch.float32, device=device)
    Wsb = (Wsd != 0).float()
    out = {k: [] for k in ("G_bs", "G_ss", "h_s", "g0_s", "n_s")}
    for c0 in range(0, B, chunk):
        Ml = torch.as_tensor(M[R, c0:c0 + chunk], dtype=torch.float32, device=device)
        Yl = torch.as_tensor(Y[R, c0:c0 + chunk], dtype=torch.float32, device=device)
        b = Ml.shape[1]
        A = Wsd.T.unsqueeze(0) * Ml.T.unsqueeze(1)                      # b x p x R
        out["G_ss"].append(torch.bmm(A, Wsd.unsqueeze(0).expand(b, -1, -1)))
        Ab = Wbd.T.unsqueeze(0) * Ml.T.unsqueeze(1)
        out["G_bs"].append(torch.bmm(Ab, Wsd.unsqueeze(0).expand(b, -1, -1)))
        out["h_s"].append((Wsd.T @ (Ml * Yl)).T)
        out["g0_s"].append((Wsd.T @ Ml).T)
        out["n_s"].append((Wsb.T @ Ml).T.cpu().numpy())
    res = {k: torch.cat(v) for k, v in out.items() if k != "n_s"}
    res["n_s"] = np.concatenate(out["n_s"])
    return res


def combine_sig(base: dict, sig: dict, alpha: float):
    """base: dict(G, h, g0, n) torch/numpy for the same conditions. G(alpha) = G_bb + a(G_bs + G_bs') + a^2 G_ss."""
    if alpha == 0:
        return base["G"], base["h"], base["g0"], base["n"]
    Gbs = sig["G_bs"]
    G = base["G"] + alpha * (Gbs + Gbs.transpose(1, 2)) + alpha ** 2 * sig["G_ss"]
    h = base["h"] + alpha * sig["h_s"]
    g0 = base["g0"] + alpha * sig["g0_s"]
    n = base["n"] + sig["n_s"]
    return G, h, g0, n
