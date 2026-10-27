"""Branch B1 — kinase–substrate LIBRARY variants scored with simple fixed statistics (DESIGN §7.1).

All scoring is single-sample (each DEV condition on its own measured sites); nothing here uses GT labels.

Library specs (a config's `lib` dict):
  {"type": "static", "name": "psp_all"|...}                      curated library, weight 1
  {"type": "pred", "name": "r2np"|"r3"|"cantley", "k": 50, "mode": "measured"|"global", "gamma": 1}
       measured: per condition, the k best-scoring sites among the sites MEASURED in that condition
       global:   the k best sites among all candidates (stored rank <= k), then restricted to measured
       weight = pct**gamma with pct = 1 - (r-1)/k, r = rank inside the kinase's selected list (r=1..k).
       (The stored global percentile is >= 0.997 for k <= 500, so pct**gamma with it would be degenerate.)
  {"type": "mix", "cur": "curated_union", "pred": {...measured pred spec...}, "alpha": 0.3, "m": 3|10|None}
       union of curated (weight 1) and predicted (weight alpha*pct**gamma) edges; an edge in both keeps
       max weight (=1, curated). Predicted edges are added only for kinases with < m curated substrates
       measured in that condition (m=None: always, plain union).

Statistics (from per-(kinase, condition) sums over measured substrates):
  wz    = sum w x / sqrt(sum w^2) / sd_c        (= RoKAI/benchmarKIN z for binary w; sd_c over all measured sites)
  mean  = sum w x / sum w                       (per-condition standardisation (x-mu_c)/sd_c is a monotone
                                                 map inside a condition, so M1 is identical; not recomputed)
  ridge = ridge-regression t: x_c (library sites, centred) ~ W^T a, lambda = 1 (squared-weight units)
A kinase is scored iff it has >= min_sub measured substrates (NaN otherwise).
"""
from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from ..data import load_dev
from ..eval import DEFAULT_ROOT, Evaluator, load_gt, read_gt_tsv
from ..library import EXT, PREDICTED_FILES, SiteMapper, load_library, verify_predicted_file

ROOT = DEFAULT_ROOT
FOLDS = [f"fold{i}" for i in range(5)]
B1_CACHE = ROOT / "cache/b1"
KMAX = 200


def gt_freeze_hash() -> str:
    return hashlib.sha256((ROOT / "results/ground_truth/GT_FREEZE.json").read_bytes()).hexdigest()


# ------------------------------------------------------------------------------------------- data
class DevData:
    def __init__(self):
        P, C, U, _ = load_gt()
        bins = read_gt_tsv(ROOT / "results/ground_truth/harmony_bins_group.tsv")
        self.P, self.C, self.U = P, C, U
        self.kin = pd.Index(list(U["kinase"]))
        self.ev = Evaluator(P, C, U, split_filter=FOLDS, bins=bins, zero_policy="ignore")
        self.conds = list(self.ev.cond_index)
        X = load_dev(conditions=self.conds)
        X = X[X.notna().any(axis=1)]
        self.X = X
        self.sites = X.index
        V = X.to_numpy(np.float64)
        self.M = ~np.isnan(V)
        self.V0 = np.where(self.M, V, 0.0)
        N = self.M.sum(0)
        mu = self.V0.sum(0) / N
        self.sd = np.sqrt((np.where(self.M, V - mu, 0.0) ** 2).sum(0) / (N - 1))
        self.fold = self.C.set_index("conditionID")["split"].reindex(self.conds).to_numpy()


def static_matrix(dd: DevData, lib_df: pd.DataFrame) -> sp.csr_matrix:
    """kinase(U) x site(dd.sites) CSR of weights (max over duplicate edges)."""
    ki = dd.kin.get_indexer(lib_df["kinase"])
    si = dd.sites.get_indexer(lib_df["site"])
    keep = (ki >= 0) & (si >= 0)
    e = pd.DataFrame({"k": ki[keep], "s": si[keep], "w": lib_df["weight"].to_numpy(np.float64)[keep]})
    e = e.groupby(["k", "s"], as_index=False)["w"].max()
    return sp.csr_matrix((e["w"], (e["k"], e["s"])), shape=(len(dd.kin), len(dd.sites)))


def pred_dense(dd: DevData, name: str) -> np.ndarray:
    """Dense kinase(U) x DEV-site predicted score matrix (NaN = no score). Cached in cache/b1/."""
    B1_CACHE.mkdir(parents=True, exist_ok=True)
    key = hashlib.sha256("\n".join(dd.sites).encode()).hexdigest()[:12]
    f = B1_CACHE / f"pred_{name}_{key}.npy"
    if f.exists():
        return np.load(f)
    verify_predicted_file(name)
    sm = SiteMapper()
    ks = pd.read_parquet(EXT / "ksm_index/sites.parquet", columns=["site_id"])
    ksm_pos = dd.sites.get_indexer(sm.from_keys(ks["site_id"]))  # ksm site idx -> DEV site position
    kin = pd.read_parquet(EXT / "ksm_index/kinases.parquet")
    acc2k = dict(zip(dd.U["uniprot"], dd.U["kinase"]))
    kpos = np.array([dd.kin.get_loc(acc2k[a]) if a in acc2k else -1 for a in kin["uniprot_acc"]])
    D = np.full((len(dd.kin), len(dd.sites)), -np.inf, np.float32)
    spec = PREDICTED_FILES[name]
    if name == "cantley":
        A = np.load(spec["path"], mmap_mode="r")
        cols = np.flatnonzero(ksm_pos >= 0)
        for ki in np.flatnonzero(kpos >= 0):
            row = np.asarray(A[ki, cols], np.float32)
            ok = np.isfinite(row)
            np.maximum.at(D[kpos[ki]], ksm_pos[cols[ok]], row[ok])
    else:
        import pyarrow.parquet as pq
        pf = pq.ParquetFile(spec["path"])
        for b in pf.iter_batches(batch_size=8_000_000, columns=["kinase_idx", "site_idx", "score"]):
            k = b.column(0).to_numpy(); s = b.column(1).to_numpy(); v = b.column(2).to_numpy().astype(np.float32)
            r, c = kpos[k], ksm_pos[s]
            keep = (r >= 0) & (c >= 0) & np.isfinite(v)
            np.maximum.at(D, (r[keep], c[keep]), v[keep])
    D[~np.isfinite(D)] = np.nan
    np.save(f, D)
    return D


def measured_topk(dd: DevData, D: np.ndarray, kmax=KMAX):
    """Per condition and kinase: DEV-site indices of the kmax best predicted sites among measured sites,
    sorted by score (C x K x kmax int32; -1 = fewer valid sites)."""
    C, K = len(dd.conds), D.shape[0]
    out = np.full((C, K, kmax), -1, np.int32)
    Dn = np.where(np.isnan(D), -np.inf, D)
    for c in range(C):
        idx = np.flatnonzero(dd.M[:, c])
        sub = Dn[:, idx]
        kk = min(kmax, sub.shape[1])
        part = np.argpartition(-sub, kk - 1, axis=1)[:, :kk]
        ps = np.take_along_axis(sub, part, 1)
        o = np.argsort(-ps, axis=1, kind="stable")
        part = np.take_along_axis(part, o, 1)
        ps = np.take_along_axis(ps, o, 1)
        sel = idx[part].astype(np.int32)
        sel[~np.isfinite(ps)] = -1
        out[c, :, :kk] = sel
    return out


# ------------------------------------------------------------------------------------------- sums
def sums_static(dd: DevData, W: sp.csr_matrix):
    """K x C arrays: num = sum w x, sw = sum w, sw2 = sum w^2, n = # measured substrates."""
    Mf = dd.M.astype(np.float64)
    Wb = W.copy(); Wb.data[:] = 1.0
    return dict(num=W @ dd.V0, sw=W @ Mf, sw2=W.multiply(W).tocsr() @ Mf, n=Wb @ Mf)


def pred_weights(k, gamma):
    pct = 1.0 - np.arange(k) / k
    return pct ** gamma


def sums_measured(dd: DevData, sel: np.ndarray, k: int, gamma: float, alpha=1.0, exclude=None):
    """Sums for per-condition measured top-k selection. exclude: bool C x K x k mask of edges to drop."""
    s = sel[:, :, :k]
    valid = s >= 0
    if exclude is not None:
        valid &= ~exclude[:, :, :k]
    w = alpha * pred_weights(k, gamma)[None, None, :] * valid
    C = s.shape[0]
    x = dd.V0[np.where(valid, s, 0), np.arange(C)[:, None, None]]
    num = (w * x).sum(2).T
    return dict(num=num, sw=w.sum(2).T, sw2=(w ** 2).sum(2).T, n=valid.sum(2).T.astype(np.float64))


def stat_from_sums(dd: DevData, S: dict, stat: str, min_sub: int) -> pd.DataFrame:
    with np.errstate(all="ignore"):
        if stat == "wz":
            z = S["num"] / np.sqrt(S["sw2"]) / dd.sd[None, :]
        elif stat == "mean":
            z = S["num"] / S["sw"]
        else:
            raise ValueError(stat)
    z = np.where((S["n"] >= min_sub) & np.isfinite(z), z, np.nan)
    return pd.DataFrame(z.T, index=dd.conds, columns=dd.kin)


# ------------------------------------------------------------------------------------------- ridge
def _ridge_one(x, rows, cols, w, K, lam):
    """x: values at the condition's library sites (len S); edges (rows=kinase, cols=site pos)."""
    W = sp.csr_matrix((w, (rows, cols)), shape=(K, len(x)))
    ks = np.flatnonzero(np.diff(W.indptr) > 0)
    out = np.full(K, np.nan)
    if len(ks) == 0 or len(x) < 3:
        return out
    Wk = W[ks]
    xc = x - x.mean()
    G = (Wk @ Wk.T).toarray()
    evals, Q = np.linalg.eigh(G)
    inv = 1.0 / (evals + lam)
    b = Wk @ xc
    a = Q @ (inv * (Q.T @ b))
    res = xc - Wk.T @ a
    df = np.sum(evals * inv)
    s2 = res @ res / max(len(x) - df - 1, 1.0)
    var = s2 * ((Q ** 2) @ (evals * inv ** 2))
    out[ks] = a / np.sqrt(var)
    return out


def ridge_scores(dd: DevData, edges_fn, min_sub: int, lam=1.0, n_jobs=32) -> pd.DataFrame:
    """edges_fn(c) -> (kinase idx, DEV site idx, weight) arrays for condition c (measured sites only)."""
    from joblib import Parallel, delayed

    def job(c):
        k, s, w = edges_fn(c)
        m = dd.M[s, c]
        k, s, w = k[m], s[m], w[m]
        cnt = np.bincount(k, minlength=len(dd.kin))
        keep = cnt[k] >= min_sub
        k, s, w = k[keep], s[keep], w[keep]
        us, pos = np.unique(s, return_inverse=True)
        return _ridge_one(dd.V0[us, c], k, pos, w, len(dd.kin), lam)

    from threadpoolctl import threadpool_limits
    with threadpool_limits(1):
        res = Parallel(n_jobs=n_jobs, backend="threading")(delayed(job)(c) for c in range(len(dd.conds)))
    return pd.DataFrame(np.vstack(res), index=dd.conds, columns=dd.kin)


# ------------------------------------------------------------------------------------------- engine
class Engine:
    """Builds and caches library components; scores a config dict -> conditions x kinases frame."""

    def __init__(self, dd: DevData, perm_seed: int | None = None):
        self.dd = dd
        self.perm_seed = perm_seed  # shuffled-library control: site-label bijection applied to every library
        self._W, self._D, self._sel, self._cur_codes = {}, {}, {}, {}

    # shuffled control: relabel sites. Curated: kinact.controls.shuffled_library (bijection on the library's own
    # sites). Predicted: bijection on DEV sites (columns of the score matrix), preserving per-kinase score lists.
    def W(self, name):
        if name not in self._W:
            df = load_library(name)
            if self.perm_seed is not None:
                from ..controls import shuffled_library
                df = shuffled_library(df, seed=self.perm_seed)
            self._W[name] = static_matrix(self.dd, df)
        return self._W[name]

    def D(self, name):
        if name not in self._D:
            D = pred_dense(self.dd, name)
            if self.perm_seed is not None:
                rng = np.random.default_rng(self.perm_seed + 1)
                D = D[:, rng.permutation(D.shape[1])]
            self._D[name] = D
        return self._D[name]

    def sel(self, name):
        if name not in self._sel:
            self._sel[name] = measured_topk(self.dd, self.D(name))
        return self._sel[name]

    def global_W(self, name, k, gamma):
        D = self.D(name)
        Dn = np.where(np.isnan(D), -np.inf, D)
        # global rank among DEV sites (approximates the stored rank among all candidates; restricted to measured
        # at scoring time)
        part = np.argsort(-Dn, axis=1)[:, :k]
        ok = np.isfinite(np.take_along_axis(Dn, part, 1))
        w = np.broadcast_to(pred_weights(k, gamma), part.shape)
        r = np.repeat(np.arange(D.shape[0]), k).reshape(part.shape)
        return sp.csr_matrix((w[ok], (r[ok], part[ok])), shape=D.shape)

    def overlap(self, cur, pred):
        """C x K x KMAX bool: predicted measured-top edge is also a curated edge."""
        W = self.W(cur)
        codes = set((W.tocoo().row.astype(np.int64) * len(self.dd.sites) + W.tocoo().col).tolist())
        sel = self.sel(pred)
        kk = np.arange(sel.shape[1], dtype=np.int64)[None, :, None]
        code = kk * len(self.dd.sites) + sel
        return np.isin(code, np.fromiter(codes, np.int64)) & (sel >= 0)

    def sums(self, lib: dict):
        t = lib["type"]
        if t == "static":
            return sums_static(self.dd, self.W(lib["name"]))
        if t == "pred" and lib["mode"] == "global":
            return sums_static(self.dd, self.global_W(lib["name"], lib["k"], lib["gamma"]))
        if t == "pred":
            return sums_measured(self.dd, self.sel(lib["name"]), lib["k"], lib["gamma"])
        if t == "mix":
            p = lib["pred"]
            Sc = sums_static(self.dd, self.W(lib["cur"]))
            ov = self.overlap(lib["cur"], p["name"])
            Sp = sums_measured(self.dd, self.sel(p["name"]), p["k"], p["gamma"], alpha=lib["alpha"], exclude=ov)
            gate = np.ones_like(Sc["n"]) if lib["m"] is None else (Sc["n"] < lib["m"]).astype(np.float64)
            return {key: Sc[key] + gate * Sp[key] for key in Sc}
        raise ValueError(t)

    def edges_fn(self, lib: dict):
        """Per-condition edge lists for ridge."""
        dd = self.dd
        t = lib["type"]
        if t == "static" or (t == "pred" and lib["mode"] == "global"):
            W = self.W(lib["name"]) if t == "static" else self.global_W(lib["name"], lib["k"], lib["gamma"])
            Wc = W.tocoo()
            arr = (Wc.row, Wc.col, Wc.data)
            return lambda c: arr

        def pred_edges(c, p, alpha, excl=None):
            s = self.sel(p["name"])[c, :, :p["k"]]
            w = np.broadcast_to(alpha * pred_weights(p["k"], p["gamma"]), s.shape)
            ok = s >= 0
            if excl is not None:
                ok &= ~excl[c, :, :p["k"]]
            r = np.repeat(np.arange(s.shape[0]), s.shape[1]).reshape(s.shape)
            return r[ok], s[ok], w[ok]

        if t == "pred":
            return lambda c: pred_edges(c, lib, 1.0)
        p = lib["pred"]
        Wc = self.W(lib["cur"]).tocoo()
        ov = self.overlap(lib["cur"], p["name"])
        ncur = sums_static(dd, self.W(lib["cur"]))["n"]

        def mix(c):
            r, s, w = pred_edges(c, p, lib["alpha"], ov)
            if lib["m"] is not None:
                g = ncur[r, c] < lib["m"]
                r, s, w = r[g], s[g], w[g]
            return (np.r_[Wc.row, r], np.r_[Wc.col, s], np.r_[Wc.data, w])
        return mix

    def score(self, cfg: dict) -> pd.DataFrame:
        if cfg["stat"] == "ridge":
            return ridge_scores(self.dd, self.edges_fn(cfg["lib"]), cfg["min_sub"], lam=cfg.get("lam", 1.0))
        return stat_from_sums(self.dd, self.sums(cfg["lib"]), cfg["stat"], cfg["min_sub"])


def cfg_name(cfg: dict) -> str:
    L = cfg["lib"]
    if L["type"] == "static":
        s = L["name"]
    elif L["type"] == "pred":
        s = f"{L['name']}@{L['k']}/{L['mode'][:4]}/g{L['gamma']:g}"
    else:
        p = L["pred"]
        s = f"{L['cur']}+{p['name']}@{p['k']}/g{p['gamma']:g}*a{L['alpha']:g}/m{L['m'] if L['m'] else 'inf'}"
    return f"{cfg['stat']}|{s}|min{cfg['min_sub']}"
