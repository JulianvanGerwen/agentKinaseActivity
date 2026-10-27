"""core-v1 — user-facing kinase activity scorer (RESEARCH_TREE D6/D7; DESIGN §1 output contract).

Model (identical to `kinact.explore.bayes2` with the core-v1 config; reimplemented for speed, no GT inputs):
  y   = standardised input: per condition x - median, divided by the RMS of the MAD-standardised residuals ('wsd')
  y_s = sum_k W_sk a_k + eps_s over the measured library sites s;   a ~ N(0, tau^2 K)
  W   = curated edges (weight) + sum_j alpha_j * (predicted edges of source j that are not curated edges)
  K   = (I + rho_f F + rho_g G) diagonal-normalised (F: same Manning family; G: same Manning group except
        Other/Atypical)
Posterior: Lambda = W_m'W_m + P0 (P0 = (tau^2 K)^-1), Li = Lambda^-1, mu = Li W_m'y, sd = sqrt(diag Li),
  q = diag(Li G Li) = diag(Li) - diag(Li P0 Li), suma = Li W_m'1, m_eff = suma^2 / q,
  z = z_vif = mu / sqrt(q (1 + rho_e max(m_eff - 1, 0))).
A kinase is scored if it has >= min_substrates measured library substrates, or (borrow) if a direct K-relative
has them. Unscored kinases are NaN in z / mean / sd; n_substrates is always given.

Speed (DESIGN §1 desideratum 4):
  * conditions with identical measured-library-site patterns share one factorisation (only W'y differs);
  * a pattern within r_max sites of an already factorised pattern (same study, a few missing sites) gets an exact
    rank-r Woodbury update of that pattern's inverse (~4 p^2 r flops instead of ~p^3);
  * G = W_m'W_m is assembled from precomputed per-site kinase pairs (index_add), W'y / W'1 / counts for all
    conditions in one sparse product;
  * diag(Li P0 Li) uses P0 = (D - Vt C Vt')/tau^2 (K is identity + family/group indicators, rank ~140), O(p^2 r_K);
  * CPU: torch float64, patterns spread over a thread pool (torch releases the GIL; intra-op threads 1).
    The thread budget is the live BLAS threadpool limit (threadpoolctl), or `n_jobs`;
  * GPU: batched float64 Cholesky over patterns; G for all patterns as one sparse (p^2 x sites) @ mask product.
Label-free fitted quantities: rho_e (config/core_v1.json; refitted and checked by scripts/build_core_artifacts.py on
DEV null conditions only). core-v1 has no site-noise model.
"""
from __future__ import annotations

import json
import os
import re
import warnings
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = ROOT / "config" / "core_v1.json"
MAD_K = 1.4826
NO_GROUP_COUPLING = ("Other", "Atypical")
SITE_KEY_RE = re.compile(r"^(?P<acc>[A-Z0-9]+(?:-\d+)?)_(?P<res>[STY])(?P<pos>\d+)$")
LONG_COLS = ("uniprot", "residue", "position", "value")


# ============================================================================ building blocks
def kinase_indicators(fam, grp, rho_f: float, rho_g: float, no_group=NO_GROUP_COUPLING):
    """K_raw = I + V diag(gamma) V' with V = [family indicators, group indicators] (empty family/group and
    Other/Atypical groups give no coupling). Returns V (p x r), gamma (r)."""
    fam = np.asarray(fam, dtype=object).astype(str); grp = np.asarray(grp, dtype=object).astype(str)
    cols, gam = [], []
    if rho_f:
        for f in sorted(set(fam) - {""}):
            cols.append(fam == f); gam.append(rho_f)
    if rho_g:
        for g in sorted(set(grp) - {""} - set(no_group)):
            cols.append(grp == g); gam.append(rho_g)
    V = np.array(cols, dtype=np.float64).T if cols else np.zeros((len(fam), 0))
    return V, np.asarray(gam, dtype=np.float64)


def kinase_kernel(fam, grp, rho_f: float, rho_g: float, no_group=NO_GROUP_COUPLING) -> np.ndarray:
    """(I + rho_f F + rho_g G) diagonal-normalised; same definition as kinact.explore.bayes.kinase_kernel."""
    fam = np.asarray(fam, dtype=object).astype(str); grp = np.asarray(grp, dtype=object).astype(str)
    F = (fam[:, None] == fam[None, :]) & (fam[:, None] != "")
    G = (grp[:, None] == grp[None, :]) & (grp[:, None] != "") & ~np.isin(grp, list(no_group))[:, None]
    K = np.eye(len(fam)) + rho_f * F + rho_g * G
    d = np.sqrt(np.diag(K))
    return K / d[:, None] / d[None, :]


def column_scales(V: np.ndarray, scale: str = "wsd", winsor: float | None = None, chunk: int = 32):
    """Per column (condition) over its finite values: centre med and scale s (MAD / SD), and for 'wsd' the RMS r of
    the (winsorised) MAD-standardised values. Standardised y = ((x - med)/s [clipped]) / r, identical arithmetic to
    kinact.explore.bayes2.standardise. Returns med, s, r (r = 1 unless 'wsd')."""
    V = np.asarray(V)
    B = V.shape[1]
    med = np.full(B, np.nan); s = np.full(B, np.nan); rr = np.ones(B)
    for j0 in range(0, B, chunk):
        C = np.ascontiguousarray(V[:, j0:j0 + chunk].T, dtype=np.float64)
        for jj in range(C.shape[0]):
            v = C[jj][np.isfinite(C[jj])]
            if not len(v):
                continue
            j = j0 + jj
            med[j] = np.median(v)
            r = v - med[j]
            if scale in ("mad", "wsd"):
                s[j] = MAD_K * np.median(np.abs(r))
            elif scale == "sd":
                s[j] = np.sqrt(np.mean(r * r))
            else:
                raise ValueError(scale)
            if scale == "wsd":
                with np.errstate(all="ignore"):
                    y = r / s[j]
                    if winsor is not None:
                        y = np.clip(y, -winsor, winsor)
                    rr[j] = np.sqrt(np.mean(y * y))
    return med, s, rr


def apply_scales(V: np.ndarray, med, s, rr, winsor: float | None = None) -> np.ndarray:
    with np.errstate(all="ignore"):
        Y = (np.asarray(V, dtype=np.float64) - med[None, :]) / s[None, :]
        if winsor is not None:
            Y = np.clip(Y, -winsor, winsor)
        return Y / rr[None, :]


def fit_rho_e(Zf: np.ndarray, ME: np.ndarray, scored: np.ndarray, grid=None):
    """Global compound-symmetry rho_e from null conditions (same procedure as bayes2.fit_rho_null): choose rho_e
    making the pooled robust variance of z_freq/sqrt(1 + rho_e (m_eff - 1)) equal across m_eff quartiles."""
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


# ============================================================================ input parsing
def to_wide(data, condition: str | None = None) -> pd.DataFrame:
    """User input -> sites x conditions float DataFrame indexed by site key `{UniProt}_{res}{pos}`.

    Accepted: Series (index = site keys; one condition, named `condition` or the Series name); wide DataFrame (index
    = site keys, or a 'site'/'site_key' text column; other columns = conditions); long DataFrame with columns
    uniprot, residue, position, value (+ optional 'condition'). Non-numeric values -> NaN (missing). Duplicate site
    keys are averaged (NaN-skipping), as in the compendium cache."""
    if isinstance(data, pd.Series):
        name = condition or (data.name if data.name is not None else "condition")
        df = data.to_frame(str(name))
    elif isinstance(data, pd.DataFrame):
        df = data
        low = {str(c).lower(): c for c in df.columns}
        if all(c in low for c in LONG_COLS):
            pos = pd.to_numeric(df[low["position"]], errors="coerce")
            keys = (df[low["uniprot"]].astype(str).str.strip().str.upper() + "_"
                    + df[low["residue"]].astype(str).str.strip().str.upper() + pos.astype("Int64").astype(str))
            cond = (df[low["condition"]].astype(str) if "condition" in low
                    else pd.Series(condition or "condition", index=df.index))
            t = pd.DataFrame({"site": keys.to_numpy(), "condition": cond.to_numpy(),
                              "value": pd.to_numeric(df[low["value"]], errors="coerce").to_numpy()})[pos.notna().to_numpy()]
            df = t.pivot_table(index="site", columns="condition", values="value", aggfunc="mean", sort=False,
                               dropna=False)
            df.columns.name = None
        else:
            for c in ("site", "site_key"):
                if c in low and not pd.api.types.is_numeric_dtype(df[low[c]]):
                    df = df.set_index(low[c])
                    break
    else:
        raise TypeError("input must be a pandas Series or DataFrame")
    if not all(pd.api.types.is_float_dtype(t) for t in df.dtypes):
        df = df.apply(pd.to_numeric, errors="coerce").astype(np.float64)
    if any(not isinstance(c, str) for c in df.columns):
        df = df.set_axis([str(c) for c in df.columns], axis=1)
    if df.index.has_duplicates:
        df = df.groupby(level=0, sort=False).mean()
    return df


def _normalise_keys(keys: pd.Index) -> pd.Index:
    s = pd.Index(keys).astype(str).str.strip().str.upper()
    return pd.Index(s.str.replace(r"\s+", "", regex=True))


# ============================================================================ result
@dataclass
class ScoreResult:
    """conditions x kinases frames: z = z_vif (calibrated z, positive = activated; NaN = unscored); mean / sd =
    posterior mean / SD of the activity change (units of the standardised input; NaN = unscored);
    n_substrates = measured library substrates (int, also for unscored kinases)."""
    z: pd.DataFrame
    mean: pd.DataFrame
    sd: pd.DataFrame
    n_substrates: pd.DataFrame
    info: dict = field(default_factory=dict)

    def long(self, dropna: bool = False) -> pd.DataFrame:
        """One row per (condition, kinase): condition, kinase, z, mean, sd, n_substrates."""
        cond = np.repeat(self.z.index.to_numpy(), self.z.shape[1])
        kin = np.tile(self.z.columns.to_numpy(), self.z.shape[0])
        out = pd.DataFrame({"condition": cond, "kinase": kin, "z": self.z.to_numpy().ravel(),
                            "mean": self.mean.to_numpy().ravel(), "sd": self.sd.to_numpy().ravel(),
                            "n_substrates": self.n_substrates.to_numpy().ravel()})
        return out.dropna(subset=["z"]).reset_index(drop=True) if dropna else out

    def wide(self) -> pd.DataFrame:
        """kinase rows; one condition -> columns z, mean, sd, n_substrates; several -> `<stat>:<condition>`."""
        stats = {"z": self.z, "mean": self.mean, "sd": self.sd, "n_substrates": self.n_substrates}
        if self.z.shape[0] == 1:
            out = pd.DataFrame({k: v.iloc[0] for k, v in stats.items()})
            out["n_substrates"] = out["n_substrates"].astype(int)
        else:
            out = pd.concat({k: v.T for k, v in stats.items()}, axis=1)
            out.columns = [f"{s}:{c}" for s, c in out.columns]
        out.index.name = "kinase"
        return out


# ============================================================================ model
class KinaseActivityModel:
    """core-v1 scorer: `KinaseActivityModel.from_config()` then `.score(df_or_series)` -> ScoreResult.

    Benchmark protocol (kinact.runtime.bench_runtime): `prepare(X, lib=None)` = load/mapping + standardisation,
    `score(state)` with a prepared state = posterior."""

    def __init__(self, cfg: dict, edges: pd.DataFrame, kinases: pd.DataFrame):
        from threadpoolctl import threadpool_limits
        with threadpool_limits(limits=1, user_api="blas"):   # p~500 linear algebra: BLAS threads oversubscribe
            self._init(cfg, edges, kinases)

    def _init(self, cfg: dict, edges: pd.DataFrame, kinases: pd.DataFrame):
        self.cfg = dict(cfg)
        self.kinases = pd.Index(kinases["kinase"].astype(str))
        p = len(self.kinases)
        fam = kinases["manning_family"].fillna("").astype(str).to_numpy()
        grp = kinases["manning_group"].fillna("").astype(str).to_numpy()
        tau2 = float(cfg["tau"]) ** 2
        self.K = kinase_kernel(fam, grp, cfg["rho_f"], cfg["rho_g"])
        self._P0 = np.linalg.inv(self.K) / tau2
        self._P0 = (self._P0 + self._P0.T) / 2
        rel = (self.K > 1e-12) & ~np.eye(p, dtype=bool)
        self._rel = rel.astype(np.float64)
        # low-rank form of P0 for diag(Li P0 Li): K_raw = I + V g V', K = Dn^-1/2 K_raw Dn^-1/2
        #   => P0 = (Dn - Vt C Vt') / tau^2, Vt = Dn^1/2 V, C = (diag(1/g) + V'V)^-1
        V, gam = kinase_indicators(fam, grp, cfg["rho_f"], cfg["rho_g"])
        Dn = 1.0 + (V * V) @ gam
        C = np.linalg.inv(np.diag(1.0 / gam) + V.T @ V) if len(gam) else np.zeros((0, 0))
        Vt = np.sqrt(Dn)[:, None] * V
        assert np.allclose((np.diag(Dn) - Vt @ C @ Vt.T) / tau2, self._P0, rtol=1e-9, atol=1e-9 * np.abs(self._P0).max())
        self._lr = (Dn / tau2, Vt, C / tau2)

        # ---- design: library sites x kinases
        kpos = self.kinases.get_indexer(edges["kinase"].astype(str))
        e = edges[kpos >= 0]
        kpos = kpos[kpos >= 0]
        alpha = {"cur": 1.0, **{n: float(v.get("alpha", 1.0)) for n, v in cfg["library"].get("predicted", {}).items()}}
        w = e["weight"].to_numpy(float) * e["source"].map(alpha).fillna(0.0).to_numpy(float)
        nz = w != 0
        self.sites = pd.Index(np.unique(e["site"].astype(str).to_numpy()[nz]))
        s = self.sites.get_indexer(e["site"].astype(str).to_numpy()[nz])
        W = sp.coo_matrix((w[nz], (s, kpos[nz])), shape=(len(self.sites), p)).tocsr()
        W.sum_duplicates()
        self._W = W
        self._WT = W.T.tocsr()
        Wb = W.copy(); Wb.data[:] = 1.0
        self._WbT = Wb.T.tocsr()
        # per-site kinase pairs (sorted by site): G_flat[k1 p + k2] = sum_{measured s} W_sk1 W_sk2
        cnt = np.diff(W.indptr)
        pc = cnt ** 2
        ptr = np.concatenate([[0], np.cumsum(pc)])
        row_of = np.repeat(np.arange(W.shape[0]), pc)
        loc = np.arange(ptr[-1]) - ptr[row_of]
        i1 = W.indptr[row_of] + loc // cnt[row_of]; i2 = W.indptr[row_of] + loc % cnt[row_of]
        self._pairs = dict(ptr=ptr, cnt=pc, k=W.indices[i1].astype(np.int64) * p + W.indices[i2],
                           w=W.data[i1] * W.data[i2], site=row_of)
        self._cache = {}
        self.rho_e = float(cfg["rho_e"])
        self.min_sub = int(cfg.get("min_substrates", 1))

    # ------------------------------------------------------------------ construction
    @classmethod
    def from_config(cls, path: str | os.PathLike | None = None) -> "KinaseActivityModel":
        """Load config (default config/core_v1.json) and the artefacts it points to (`artifacts`, relative to the
        config file)."""
        path = Path(path) if path is not None else DEFAULT_CONFIG
        cfg = json.loads(path.read_text())
        art = Path(cfg.get("artifacts", "../models/core_v1"))
        art = art if art.is_absolute() else (path.resolve().parent / art).resolve()
        if not (art / "edges.tsv.gz").exists():
            raise FileNotFoundError(f"model artefacts missing in {art}; run scripts/build_core_artifacts.py")
        edges = pd.read_csv(art / "edges.tsv.gz", sep="\t")
        kin = pd.read_csv(art / "kinases.tsv", sep="\t", keep_default_na=False, na_values=[""])
        man = json.loads((art / "manifest.json").read_text()) if (art / "manifest.json").exists() else {}
        fit = man.get("rho_e_refit")
        if fit is not None and abs(fit - float(cfg["rho_e"])) > 1e-12:
            warnings.warn(f"config rho_e={cfg['rho_e']} differs from the artefact refit {fit}")
        m = cls(cfg, edges, kin)
        m.artifacts_dir, m.manifest, m.config_path = art, man, path
        return m

    # ------------------------------------------------------------------ input -> standardised library matrix
    def prepare(self, data, lib=None, condition: str | None = None) -> dict:
        """Load/mapping + standardisation (`lib` is ignored; bench_runtime signature). Every finite input value
        enters its condition's centre/scale; only sites in the library enter the likelihood."""
        X = to_wide(data, condition)
        li = self.sites.get_indexer(X.index)
        miss = li < 0
        if miss.any():                                   # retry unmatched keys after normalisation (case/whitespace)
            sample = X.index[miss][:2000]
            if not _normalise_keys(sample).equals(pd.Index(sample).astype(str)):
                li[miss] = self.sites.get_indexer(_normalise_keys(X.index[miss]))
        ok = li >= 0
        V = X.to_numpy()
        wz, sc = self.cfg.get("winsor"), self.cfg.get("scale", "wsd")
        nj = max(1, min(_current_threads(), V.shape[1] // 16))
        if nj > 1:                                       # numpy partition/reductions release the GIL
            with ThreadPoolExecutor(nj) as ex:
                parts = list(ex.map(lambda j: column_scales(V[:, j:j + 16], sc, wz), range(0, V.shape[1], 16)))
            med, s, rr = (np.concatenate(z) for z in zip(*parts))
        else:
            med, s, rr = column_scales(V, sc, wz)
        Yk = apply_scales(V[ok], med, s, rr, wz)
        fin = np.isfinite(Yk)
        Ylib = np.zeros((len(self.sites), X.shape[1])); Mlib = np.zeros((len(self.sites), X.shape[1]), bool)
        Ylib[li[ok]] = np.where(fin, Yk, 0.0)
        Mlib[li[ok]] = fin
        info = {"n_input_sites": int(len(X)), "n_library_sites_matched": int(ok.sum()),
                "n_unknown_sites": int((~ok).sum()), "n_conditions": int(X.shape[1])}
        if len(X) and not ok.any():
            warnings.warn("no input site matches the kinase-substrate library (keys must be `{UniProt}_{res}{pos}`, "
                          "e.g. P31749_S473); all kinases are NaN")
        return {"conditions": list(X.columns), "Y": Ylib, "M": Mlib, "info": info}

    # ------------------------------------------------------------------ scoring
    def score(self, data, device: str = "cpu", n_jobs: int | None = None, condition: str | None = None) -> ScoreResult:
        """Score sites x conditions input (see `to_wide` for accepted formats; a prepared state is accepted too).
        device: 'cpu' or 'cuda[:i]'. n_jobs: CPU threads (default: live BLAS threadpool limit)."""
        st = data if isinstance(data, dict) and "M" in data else self.prepare(data, condition=condition)
        po = self.posterior(st, device, n_jobs)
        mu, q, suma = po["mu"], po["q"], po["suma"]
        with np.errstate(all="ignore"):
            me = np.where(q > 0, suma ** 2 / q, 0.0)
            z = mu / np.sqrt(q * (1 + self.rho_e * np.maximum(me - 1, 0)))
        rep = po["report"] & np.isfinite(z)
        mk = lambda A: pd.DataFrame(np.where(rep, A, np.nan), index=st["conditions"], columns=self.kinases)
        return ScoreResult(z=mk(z), mean=mk(mu), sd=mk(po["sd"]),
                           n_substrates=pd.DataFrame(po["n"], index=st["conditions"], columns=self.kinases),
                           info={**st["info"], "n_patterns": po["n_patterns"], "n_woodbury": po["n_woodbury"],
                                 "device": str(device), "model": self.cfg.get("name"), "rho_e": self.rho_e})

    def posterior(self, st: dict, device: str = "cpu", n_jobs: int | None = None) -> dict:
        """Raw posterior quantities (conditions x kinases): mu, sd, q, suma, n (int), report (bool)."""
        Y, M = st["Y"], st["M"]
        B, p = Y.shape[1], len(self.kinases)
        Mf = M.astype(np.float64)
        H = np.asarray(self._WT @ (Y * Mf))                 # p x B   W'y
        G0 = np.asarray(self._WT @ Mf)                      # p x B   W'1
        N = np.asarray(self._WbT @ Mf).T                    # B x p   measured substrates
        if B:
            _, first, pat = np.unique(np.packbits(M, axis=0).T, axis=0, return_index=True, return_inverse=True)
            pat = pat.ravel()
        else:
            first, pat = np.zeros(0, int), np.zeros(0, int)
        out = {k: np.zeros((B, p)) for k in ("mu", "sd", "q", "suma")}
        if str(device).startswith("cuda"):
            nw = self._run_gpu(M, H, G0, first, pat, out, device)
        else:
            nw = self._run_cpu(M, H, G0, first, pat, out, n_jobs)
        has = N >= self.min_sub
        if self.cfg.get("borrow", True):
            from threadpoolctl import threadpool_limits
            with threadpool_limits(limits=min(_current_threads(), 8), user_api="blas"):
                rep = has | ((has.astype(np.float64) @ self._rel) > 0)
        else:
            rep = has
        return {**out, "n": np.rint(N).astype(np.int64), "report": rep, "n_patterns": int(len(first)),
                "n_woodbury": int(nw)}

    # ------------------------------------------------------------------ CPU path
    def _tconst(self):
        import torch
        if "cpu" not in self._cache:
            Dn, Vt, C = self._lr
            pr = self._pairs
            self._cache["cpu"] = dict(
                P0=torch.from_numpy(self._P0), Dn=torch.from_numpy(Dn), C=torch.from_numpy(C),
                VtT=_csr(torch.from_numpy(np.ascontiguousarray(Vt.T))),
                ptr=torch.from_numpy(pr["ptr"][:-1]), cnt=torch.from_numpy(pr["cnt"]),
                pk=torch.from_numpy(pr["k"]), pw=torch.from_numpy(pr["w"]))
        return self._cache["cpu"]

    def _G_cpu(self, mrows):
        import torch
        c = self._tconst()
        p = len(self.kinases)
        sel = torch.from_numpy(mrows)
        cnt = c["cnt"][sel]
        start = c["ptr"][sel]
        tot = int(cnt.sum())
        idx = torch.arange(tot) + torch.repeat_interleave(start - (torch.cumsum(cnt, 0) - cnt), cnt)
        return torch.zeros(p * p, dtype=torch.float64).index_add_(0, c["pk"][idx], c["pw"][idx]).view(p, p)

    def _qd(self, Li):
        """diag(Li P0 Li) = rowsum((Li^2) Dn) - rowsum((LV C) * LV), LV = Li Vt (Li symmetric)."""
        c = self._tconst()
        LV = (c["VtT"] @ Li).T
        return (Li * Li) @ c["Dn"] - ((LV @ c["C"]) * LV).sum(1)

    def _run_cpu_numpy(self, M, H, G0, first, pat, out):
        """Fallback without torch: scipy Cholesky per pattern, sequential, no Woodbury."""
        import scipy.linalg as sla
        Dn, Vt, C = self._lr
        for i in range(len(first)):
            cols = np.flatnonzero(pat == i)
            Wm = self._W[np.flatnonzero(M[:, first[i]])]
            cf = sla.cho_factor((Wm.T @ Wm).toarray() + self._P0, lower=True, check_finite=False)
            Li = sla.cho_solve(cf, np.eye(len(Dn)), check_finite=False)
            dLi = np.diag(Li).copy()
            LV = Li @ Vt
            qd = (Li * Li) @ Dn - ((LV @ C) * LV).sum(1)
            R = Li @ np.column_stack([H[:, cols], G0[:, first[i]]])
            out["mu"][cols] = R[:, :-1].T; out["suma"][cols] = R[:, -1]
            out["sd"][cols] = np.sqrt(dLi); out["q"][cols] = np.maximum(dLi - qd, 0)
        return 0

    def _run_cpu(self, M, H, G0, first, pat, out, n_jobs, window=32):
        try:
            import torch
        except ImportError:                                                   # pragma: no cover
            return self._run_cpu_numpy(M, H, G0, first, pat, out)
        from threadpoolctl import threadpool_limits
        P, p = len(first), len(self.kinases)
        if P == 0:
            return 0
        order = np.argsort(pat, kind="stable")
        bounds = np.searchsorted(pat[order], np.arange(P + 1))
        groups = [order[bounds[i]:bounds[i + 1]] for i in range(P)]
        Mp = M[:, first]
        c = self._tconst()
        P0 = c["P0"]
        Ht = torch.from_numpy(np.ascontiguousarray(H)); G0t = torch.from_numpy(np.ascontiguousarray(G0))
        nj = n_jobs if n_jobs is not None else _current_threads()
        nj = max(1, min(int(nj), P))
        # ---- phase 1 (cheap, sequential): Woodbury bases. Pattern i reuses the closest of the last `window` bases if
        #      it differs in <= r_max library sites; update cost ~4 p^2 r vs ~1.3-2 p^3 for a factorisation.
        r_max = int(self.cfg.get("woodbury_r_max", p // 4 if nj == 1 else 0))   # threads: bandwidth-bound -> off
        base_of = np.full(P, -1)
        if r_max > 0 and P > 1:
            packed = np.packbits(Mp, axis=0).T.copy()
            bases = []
            for i in range(P):
                if bases:
                    cand = np.asarray(bases[-window:])
                    d = np.bitwise_count(packed[cand] ^ packed[i]).sum(1)
                    j = int(d.argmin())
                    if d[j] <= r_max:
                        base_of[i] = cand[j]
                        continue
                bases.append(i)
        has_dep = np.zeros(P, bool); has_dep[base_of[base_of >= 0]] = True
        store = {}
        inv = _spd_inverse if nj == 1 else (lambda A: torch.cholesky_inverse(torch.linalg.cholesky(A)))

        def write(i, R, dLi, qd):
            cols = groups[i]
            R = R.numpy()
            out["mu"][cols] = R[:, :-1].T
            out["suma"][cols] = R[:, -1]
            out["sd"][cols] = dLi.sqrt().numpy()
            out["q"][cols] = (dLi - qd).clamp_min(0).numpy()

        def rhs(i):
            cols = torch.from_numpy(groups[i])
            return torch.cat([Ht.index_select(1, cols), G0t[:, int(groups[i][0])][:, None]], 1)

        def job_base(i):
            Li = inv(self._G_cpu(np.flatnonzero(Mp[:, i])).add_(P0))
            dLi = torch.diagonal(Li).clone()
            qd = self._qd(Li)
            write(i, Li @ rhs(i), dLi, qd)
            if has_dep[i]:
                store[i] = (Li, dLi, qd)

        def job_update(i):
            b = base_of[i]
            Li, dLi_b, qd_b = store[b]
            mi, mb = Mp[:, i], Mp[:, b]
            dif = mi != mb
            rows = np.concatenate([np.flatnonzero(dif & mi), np.flatnonzero(dif & mb)])
            n_add = int((dif & mi).sum())
            U = torch.from_numpy(self._W[rows].toarray())                                # r x p
            sgn = torch.ones(len(rows), dtype=torch.float64); sgn[n_add:] = -1.0
            F = Li @ U.T                                                                 # p x r
            E = torch.linalg.solve(torch.diag(sgn) + U @ F, F.T).T                       # F (S + U'Li U)^-1
            PF = P0 @ F
            dLi = dLi_b - (E * F).sum(1)
            qd = qd_b - 2 * ((Li @ PF) * E).sum(1) + ((E @ (F.T @ PF)) * E).sum(1)
            Rh = rhs(i)
            write(i, Li @ Rh - E @ (F.T @ Rh), dLi, qd)

        b_ids = np.flatnonzero(base_of < 0); u_ids = np.flatnonzero(base_of >= 0)
        old = torch.get_num_threads()
        torch.set_num_threads(1)
        try:
            with threadpool_limits(limits=1):
                if nj == 1:
                    for i in b_ids:
                        job_base(i)
                    for i in u_ids:
                        job_update(i)
                else:
                    # torch intra-op (OpenMP) thread counts are per OS thread: pin every worker to 1
                    with ThreadPoolExecutor(nj, initializer=torch.set_num_threads, initargs=(1,)) as ex:
                        list(ex.map(job_base, b_ids))
                        list(ex.map(job_update, u_ids))
        finally:
            torch.set_num_threads(old)
        return len(u_ids)

    # ------------------------------------------------------------------ GPU path
    def _gconst(self, dev):
        import torch
        key = str(dev)
        if key not in self._cache:
            p = len(self.kinases)
            pr = self._pairs
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                Pm = torch.sparse_coo_tensor(torch.from_numpy(np.vstack([pr["k"], pr["site"]])), torch.from_numpy(pr["w"]),
                                             size=(p * p, self._W.shape[0])).coalesce().to_sparse_csr().to(dev)
            Dn, Vt, C = self._lr
            self._cache[key] = dict(Pm=Pm, P0=torch.as_tensor(self._P0, device=dev),
                                    Dn=torch.as_tensor(Dn, device=dev), Vt=torch.as_tensor(Vt, device=dev),
                                    C=torch.as_tensor(C, device=dev))
        return self._cache[key]

    def _run_gpu(self, M, H, G0, first, pat, out, device, chunk=192):
        import torch
        p = len(self.kinases)
        dev = torch.device(device)
        c = self._gconst(dev)
        Ht = torch.as_tensor(H, device=dev)
        G0t = torch.as_tensor(G0, device=dev)
        pat_t = torch.as_tensor(pat, device=dev)
        for c0 in range(0, len(first), chunk):
            ids = np.arange(c0, min(c0 + chunk, len(first)))
            Mp = torch.as_tensor(M[:, first[ids]], dtype=torch.float64, device=dev)       # sites x b
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                G = (c["Pm"] @ Mp).T.reshape(len(ids), p, p)
            Li = torch.cholesky_inverse(torch.linalg.cholesky(G.add_(c["P0"])))
            del G
            dLi = torch.diagonal(Li, dim1=1, dim2=2)
            LV = Li @ c["Vt"]
            qq = (dLi - ((Li * Li) @ c["Dn"] - ((LV @ c["C"]) * LV).sum(-1))).clamp_min(0)
            del LV
            cols = np.flatnonzero((pat >= ids[0]) & (pat <= ids[-1]))
            for k0 in range(0, len(cols), chunk):
                cc = cols[k0:k0 + chunk]
                loc = pat_t[cc] - int(ids[0])
                R = torch.stack([Ht[:, cc].T, G0t[:, cc].T], -1)                          # b x p x 2
                R = torch.bmm(Li[loc], R)
                out["mu"][cc] = R[..., 0].cpu().numpy(); out["suma"][cc] = R[..., 1].cpu().numpy()
                out["sd"][cc] = dLi[loc].sqrt().cpu().numpy(); out["q"][cc] = qq[loc].cpu().numpy()
            del Li
        return 0


def _csr(t):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return t.to_sparse_csr()


def _spd_inverse(A, base: int = 64):
    """Inverse of an SPD matrix (torch float64): recursive 2x2 Schur-complement blocking (gemm-dominated, ~4/3 n^3
    flops; faster than potrf+potri at n~500 single-threaded); Cholesky + cholesky_inverse at the leaves."""
    import torch
    n = A.shape[0]
    if n <= base:
        return torch.cholesky_inverse(torch.linalg.cholesky(A))
    k = n // 2
    A12 = A[:k, k:]
    I11 = _spd_inverse(A[:k, :k], base)
    X = I11 @ A12
    I22 = _spd_inverse(A[k:, k:] - A12.T @ X, base)
    Y = X @ I22
    out = torch.empty_like(A)
    out[:k, :k] = I11 + Y @ X.T
    out[:k, k:] = -Y
    out[k:, :k] = -Y.T
    out[k:, k:] = I22
    return out


def _current_threads() -> int:
    """Thread budget = the live BLAS threadpool limit (so threadpool_limits(n) / bench_runtime settings apply)."""
    try:
        from threadpoolctl import threadpool_info
        n = [d["num_threads"] for d in threadpool_info() if d.get("user_api") == "blas"]
        if n:
            return int(min(n))
    except Exception:                                           # pragma: no cover
        pass
    return os.cpu_count() or 1


def score(data, config: str | os.PathLike | None = None, device: str = "cpu") -> ScoreResult:
    """One-call API: KinaseActivityModel.from_config(config).score(data, device=device)."""
    return KinaseActivityModel.from_config(config).score(data, device=device)
