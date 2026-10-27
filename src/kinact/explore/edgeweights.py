"""Branch B4b — edge weights as a SHARED function of label-free edge/site features (DESIGN §7.4).

Model: the B2.1 Bayesian core (bayes2; y = W a + eps, a ~ N(0, tau^2 K), median-centred, z_post output) with
    W_sk = base_sk * exp(5 tanh(theta . phi(k, s) / 5))   (soft clamp, ~exp(theta . phi) for |.| < 2),
base = 1 for curated_union and r2np@100 edges (core-v1, alpha = 1). theta = 0 reproduces core-v1 exactly.
phi (centred and scaled over edges; no GT labels, no per-kinase / per-edge free parameters; <= 12 params):
    is_r2np        predicted-only edge (r2np@100, not curated)
    psp_vivo       curated edge with >= 1 PSP in-vivo row
    psp_vitro      curated PSP edge with in-vitro rows only
    n_src          log(1 + number of curated source tags)                  (0 for predicted-only edges)
    r2np_pct       r2np within-kinase percentile (top-500 list; 0 if absent)
    r3_pct         r3 within-kinase percentile (top-500; 0 if absent)
    cantley_pct    Cantley within-kinase percentile (top-500; 0 if absent)
    log_prom       log # kinases whose library (curated + r2np@100) contains the site
    res_mismatch   (kinase in Manning TK group) XOR (site residue is Y)
    log_sd         log SD of the standardised site value over ALLOWED conditions         (fold-dependent)
    meas_freq      fraction of ALLOWED conditions in which the site is measured          (fold-dependent)
ALLOWED = DEV conditions allowed for the held-out fold (cache/b4/allowed.json: excludes the held-out fold and
conditions sharing a control with it); inside inner CV additionally excluding the inner validation fold.

Training objective (aligned to M1-sel): within every training GT condition (M1-sel only), pairwise smooth-AUROC loss
sigmoid(-(t_k - t_j) / (T * sd_c)) (or logistic softplus(-.)), t = sign * z (z_vif by default), over the condition's evaluable negatives j that are scored;
pair weights reproduce the kinase-macro averaging (pairs -> (kinase, group) -> kinase -> macro) of kinact.eval;
sd_c = SD of z over scored kinases in c (scale-free); + lam * ||theta||^2. Full-batch L-BFGS on GPU, gradients
through Cholesky (torch autograd).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

FEATURES = ["is_r2np", "psp_vivo", "psp_vitro", "n_src", "r2np_pct", "r3_pct", "cantley_pct", "log_prom",
            "res_mismatch", "log_sd", "meas_freq"]
STATIC = FEATURES[:9]
DYNAMIC = FEATURES[9:]
LOGW_MAX = 5.0      # log-weights soft-clamped to (-5, 5) via 5 tanh(x / 5) (numerical safety in line searches)


# ----------------------------------------------------------------------------- edges + static features
def edge_table(md, site_index, universe: pd.DataFrame, lib_cur: pd.DataFrame, psp_all: pd.DataFrame,
               psp_vivo: pd.DataFrame, pred500: dict) -> pd.DataFrame:
    """One row per nonzero edge of md ('cur' + every other matrix, which bayes2 keeps disjoint from 'cur').
    Columns: row (position in md.lib_rows), k (kinase position), base, site, kinase + STATIC features (raw)."""
    site_index = pd.Index(site_index)
    parts = []
    for name, Mx in md.mats.items():
        C = Mx.tocoo()
        parts.append(pd.DataFrame({"row": C.row, "k": C.col, "base": C.data, "is_r2np": float(name != "cur")}))
    E = pd.concat(parts, ignore_index=True)
    E = E[E["base"] != 0].groupby(["row", "k"], as_index=False).agg(base=("base", "sum"), is_r2np=("is_r2np", "min"))
    E["site"] = np.asarray(site_index[md.lib_rows[E["row"].to_numpy()]])
    E["kinase"] = np.asarray(md.kinases[E["k"].to_numpy()])
    key = E["kinase"] + "|" + E["site"]
    kk = lambda df: set(df["kinase"] + "|" + df["site"])
    cur = E["is_r2np"] == 0
    vivo, pspa = kk(psp_vivo), kk(psp_all)
    E["psp_vivo"] = (cur & key.isin(vivo)).astype(float)
    E["psp_vitro"] = (cur & key.isin(pspa) & ~key.isin(vivo)).astype(float)
    nsrc = (lib_cur.assign(key=lib_cur["kinase"] + "|" + lib_cur["site"])
            .set_index("key")["source"].astype(str).str.split(";").map(len))
    nsrc = nsrc.groupby(level=0).max()
    E["n_src"] = np.where(cur, np.log1p(key.map(nsrc).fillna(0).to_numpy()), 0.0)
    for n, lib in pred500.items():
        w = lib.assign(key=lib["kinase"] + "|" + lib["site"]).groupby("key")["weight"].max()
        E[f"{n}_pct"] = key.map(w).fillna(0.0).to_numpy()
    E["log_prom"] = np.log(E.groupby("site")["k"].transform("count").to_numpy())
    grp = universe.set_index("kinase").reindex(list(md.kinases))["manning_group"].fillna("").to_numpy()
    is_tk = grp[E["k"].to_numpy()] == "TK"
    is_y = E["site"].str.extract(r"_([STY])\d+$")[0].to_numpy() == "Y"
    E["res_mismatch"] = (is_tk ^ is_y).astype(float)
    return E.reset_index(drop=True)


def site_stats(Y: np.ndarray, M: np.ndarray, cols: np.ndarray):
    """Per lib row: log SD of standardised values and measurement frequency over the condition columns `cols`
    (boolean or index array into the columns of Y/M). Y: standardised values (0 where missing), M: 0/1 mask."""
    Ys, Ms = Y[:, cols], M[:, cols]
    n = Ms.sum(1)
    with np.errstate(all="ignore"):
        mean = (Ys * Ms).sum(1) / n
        var = (Ys * Ys * Ms).sum(1) / n - mean ** 2
    sd = np.sqrt(np.clip(var, 1e-4, None))
    med_sd = np.nanmedian(sd[n >= 3])
    log_sd = np.where(n >= 3, np.log(sd), np.log(med_sd))
    return log_sd, n / Ms.shape[1]


def feature_matrix(E: pd.DataFrame, log_sd: np.ndarray | None = None, meas_freq: np.ndarray | None = None,
                   features=FEATURES) -> np.ndarray:
    """Edges x features, each column centred and scaled (SD 1) over edges; constant columns -> 0."""
    cols = []
    for f in features:
        if f == "log_sd":
            v = log_sd[E["row"].to_numpy()]
        elif f == "meas_freq":
            v = meas_freq[E["row"].to_numpy()]
        else:
            v = E[f].to_numpy(float)
        sd = v.std()
        cols.append((v - v.mean()) / sd if np.ptp(v) > 0 else np.zeros_like(v))
    return np.column_stack(cols).astype(np.float32)


# ----------------------------------------------------------------------------- differentiable core
class DiffCore:
    """z_post of the bayes2 core with edge weights base * exp(Phi theta), differentiable in theta.

    Y, M: standardised values / 0-1 mask on md.lib_rows for the conditions to score (S x B).
    K: kinase kernel (p x p); tau: prior SD."""

    def __init__(self, E: pd.DataFrame, n_rows: int, p: int, K: np.ndarray, tau: float, device="cuda"):
        import torch
        self.t = torch
        self.dev = device
        self.r = torch.as_tensor(E["row"].to_numpy().copy(), device=device)
        self.k = torch.as_tensor(E["k"].to_numpy(), device=device)
        self.base = torch.as_tensor(E["base"].to_numpy(np.float32), device=device)
        self.S, self.p = n_rows, p
        self.P0 = (torch.linalg.inv(torch.as_tensor(K, dtype=torch.float64, device=device)) / tau ** 2)

    def W(self, theta, Phi):
        t = self.t
        w = self.base if theta is None else self.base * t.exp(LOGW_MAX * t.tanh((Phi @ theta) / LOGW_MAX))
        W = t.zeros(self.S, self.p, device=self.dev, dtype=t.float32)
        return W.index_put((self.r, self.k), w)

    def z(self, W, Y, M, out="z_vif", rho_e=0.0):
        """Y, M: torch (S x b). out: z_post | mean | z_freq | z_vif (rho_e) | all (dict mu, q, me, z_post).
        Returns b x p float64 (z_vif = mu / sqrt(q (1 + rho_e (m_eff - 1))), as bayes2.output_stat)."""
        t = self.t
        Gs, hs, gs = [], [], []
        for j in range(Y.shape[1]):
            idx = t.nonzero(M[:, j] > 0).squeeze(1)
            Ws = W.index_select(0, idx).double()
            Gs.append(Ws.T @ Ws); hs.append(Ws.T @ Y.index_select(0, idx)[:, j].double()); gs.append(Ws.sum(0))
        G = t.stack(Gs).double(); h = t.stack(hs).double(); g0 = t.stack(gs).double()
        L = t.linalg.cholesky(G + self.P0)
        Li = t.cholesky_inverse(L)
        mu = (Li @ h.unsqueeze(-1)).squeeze(-1)
        if out == "mean":
            return mu
        if out == "z_post":
            return mu / t.sqrt(t.diagonal(Li, dim1=1, dim2=2))
        q = ((Li @ G) * Li).sum(-1).clamp_min(1e-12)
        suma = (Li @ g0.unsqueeze(-1)).squeeze(-1)
        me = suma ** 2 / q
        if out == "all":
            return dict(mu=mu, q=q, me=me, z_post=mu / t.sqrt(t.diagonal(Li, dim1=1, dim2=2)))
        if out == "z_freq":
            return mu / t.sqrt(q)
        return mu / t.sqrt(q * (1 + rho_e * (me - 1).clamp_min(0)))


def pair_tables(ev, conds_train, rep: np.ndarray, cond_pos: dict, m1sel_only=True):
    """Training targets from the frozen evaluator (ev.targets, ev.neg_mask) restricted to `conds_train`.
    Returns DataFrame (col = column of the score batch, ki, sign, w) and negative masks (n_t x p, already & rep).
    Pair weights reproduce the macro: 1/(#pairs of (k,group)) / (#groups of k) / (#kinases)."""
    tg = ev.targets
    tg = tg[tg["conditionID"].isin(set(conds_train))]
    if m1sel_only:
        tg = tg[tg["in_M1sel"]]
    tg = tg.copy()
    tg["col"] = tg["conditionID"].map(cond_pos).astype(int)
    tg = tg[rep[tg["col"].to_numpy(), tg["ki"].to_numpy()]]
    npg = tg.groupby(["kinase", "group"])["ki"].transform("count")
    ngr = tg.groupby("kinase")["group"].transform("nunique")
    tg["w"] = 1.0 / npg / ngr / tg["kinase"].nunique()
    neg = ev.neg_mask[tg["ci"].to_numpy()] & rep[tg["col"].to_numpy()]
    nN = ev.neg_mask[tg["ci"].to_numpy()].sum(1)
    return tg.reset_index(drop=True), neg, nN


def pairwise_loss(z, cols_local, ki, sign, w, neg, nN, T=0.5, kind="logistic"):
    """z: b x p (torch). Targets index into rows of z by cols_local. Returns weighted sum of mean pair losses.
    kind: 'logistic' softplus(-d) (convex surrogate) | 'sigmoid' sigmoid(-d) (smooth 1 - AUROC)."""
    import torch
    rep = torch.isfinite(z)
    zc = torch.where(rep, z, torch.zeros_like(z))
    cnt = rep.sum(1).clamp_min(2)
    mu = zc.sum(1) / cnt
    sd = torch.sqrt(((zc - mu[:, None]) ** 2 * rep).sum(1) / (cnt - 1)).clamp_min(1e-8)
    Z = zc[cols_local] * sign[:, None]
    tk = Z.gather(1, ki[:, None])
    d = (tk - Z) / (T * sd[cols_local][:, None])
    L = (torch.nn.functional.softplus(-d) if kind == "logistic" else torch.sigmoid(-d)) * neg
    return (w * L.sum(1) / nN).sum()


def fit_theta(core: DiffCore, Phi, Y, M, rep, ev, conds, train_conds, lam: float, T=0.5, chunk=48, max_iter=30,
              theta0=None, verbose=False, kind="sigmoid", out="z_vif", rho_e=0.0):
    """Full-batch L-BFGS on pairwise loss over training conditions (M1-sel targets). Y/M: torch S x B for `conds`
    (all conditions; only train columns are used). Returns theta (numpy) and history."""
    import torch
    pos = {c: i for i, c in enumerate(conds)}
    tg, neg, nN = pair_tables(ev, train_conds, rep, pos)
    cols = np.unique(tg["col"].to_numpy())
    d = Phi.shape[1]
    theta = torch.zeros(d, device=core.dev, dtype=torch.float32, requires_grad=True)
    if theta0 is not None:
        with torch.no_grad():
            theta.copy_(torch.as_tensor(theta0))
    negT = torch.as_tensor(neg, device=core.dev, dtype=torch.float64)
    nNT = torch.as_tensor(nN, device=core.dev, dtype=torch.float64)
    wT = torch.as_tensor(tg["w"].to_numpy(), device=core.dev, dtype=torch.float64)
    kiT = torch.as_tensor(tg["ki"].to_numpy(), device=core.dev)
    sgT = torch.as_tensor(tg["sign"].to_numpy(), device=core.dev, dtype=torch.float64)
    colT = tg["col"].to_numpy()
    repT = torch.as_tensor(rep, device=core.dev)
    chunks = [cols[i:i + chunk] for i in range(0, len(cols), chunk)]
    tsel = [np.flatnonzero(np.isin(colT, ch)) for ch in chunks]
    opt = torch.optim.LBFGS([theta], lr=1, max_iter=max_iter, history_size=10, line_search_fn="strong_wolfe")
    hist = []

    def closure():
        opt.zero_grad()
        tot = 0.0
        for ch, ts in zip(chunks, tsel):
            W = core.W(theta, Phi)
            chT = torch.as_tensor(ch, device=core.dev)
            z = core.z(W, Y[:, chT], M[:, chT], out, rho_e)
            z = torch.where(repT[chT], z, torch.full_like(z, float("nan")))
            loc = torch.as_tensor(np.searchsorted(ch, colT[ts]), device=core.dev)
            l = pairwise_loss(z, loc, kiT[ts], sgT[ts], wT[ts], negT[ts], nNT[ts], T, kind)
            l.backward()
            tot += float(l)
        reg = lam * (theta ** 2).sum()
        reg.backward()
        tot += float(reg)
        hist.append(tot)
        if verbose:
            print(f"  loss {tot:.5f} theta {theta.detach().cpu().numpy().round(3)}", flush=True)
        return torch.tensor(tot)

    opt.step(closure)
    return theta.detach().cpu().numpy(), hist


def score_theta(core: DiffCore, theta, Phi, Y, M, rep, chunk=64, out="z_vif", rho_e=0.0):
    """Score all columns of Y/M with fixed theta (None = core-v1). Returns numpy B x p with NaN outside rep
    (out='all': dict of numpy arrays mu, q, me, z_post)."""
    import torch
    with torch.no_grad():
        th = None if theta is None else torch.as_tensor(np.asarray(theta, np.float32), device=core.dev)
        W = core.W(th, Phi)
        parts = [core.z(W, Y[:, i:i + chunk], M[:, i:i + chunk], out, rho_e) for i in range(0, Y.shape[1], chunk)]
    if out == "all":
        return {k: np.concatenate([p[k].cpu().numpy() for p in parts]) for k in parts[0]}
    Z = np.concatenate([p.cpu().numpy() for p in parts])
    return np.where(rep, Z, np.nan)


def fit_rho_e(core: DiffCore, theta, Phi, Yn, Mn, rep_n):
    """Global compound-symmetry rho_e for the weighted design, fitted on null conditions only (bayes2.fit_rho_null)."""
    from . import bayes2 as B2
    a = score_theta(core, theta, Phi, Yn, Mn, rep_n, out="all")
    zf = a["mu"] / np.sqrt(a["q"])
    rho, _ = B2.fit_rho_null(zf, a["me"], rep_n)
    return rho


def promiscuity_weights(E: pd.DataFrame, beta: float) -> np.ndarray:
    """Hand-specified 1-parameter variant: per-edge log-weight = -beta * log(promiscuity)."""
    return -beta * E["log_prom"].to_numpy(np.float32)
