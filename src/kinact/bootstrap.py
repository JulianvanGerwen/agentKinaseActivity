"""Stratified cluster bootstrap over GT groups (DESIGN §5, §13 B1/B6; r1_stats #5, r2_stats N7).

Resample unit = GT group; strata = group's majority GT source (by M1-sel pairs); B=2000,
seed 0; ONE resample-count matrix W (B × n_groups) shared by all methods and endpoints.
Because whole groups are resampled, (kinase, group) means are fixed and the full metric
(kinase set, bin membership, macro) is recomputed per resample as weighted sums:
    kinase_b = Σ_g W[b,g]·v_kg / Σ_g W[b,g]   (kinase present iff denominator > 0)
    metric_b = mean over present kinases (harmony: per bin, then mean over non-empty bins).
Inference (DESIGN §12 C8, bias-corrected because the kinase composition of the macro average
changes across resamples):
    CI  = basic (reverse-percentile) interval [2·θ̂ − q97.5, 2·θ̂ − q2.5]
    p   = one-sided shifted-null p for Δ > 0: (1 + #{Δ*_b − mean(Δ*) ≥ Δ̂}) / (B + 1)
    raw percentile interval [q2.5, q97.5], mean(θ*) and bias = mean(θ*) − θ̂ are kept as diagnostics.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import sparse

from .eval import Evaluator


@dataclass
class BootstrapPlan:
    groups: pd.Index
    strata: pd.Series
    W: np.ndarray  # (B, n_groups) float32 counts
    seed: int

    @classmethod
    def from_evaluator(cls, ev: Evaluator, B: int = 2000, seed: int = 0) -> "BootstrapPlan":
        return cls.build(ev.group_index, ev.group_stratum, B, seed)

    @classmethod
    def build(cls, groups, strata, B: int = 2000, seed: int = 0) -> "BootstrapPlan":
        groups = pd.Index(groups)
        strata = pd.Series(strata).reindex(groups)
        rng = np.random.default_rng(seed)
        W = np.zeros((B, len(groups)), dtype=np.float32)
        for s in sorted(strata.unique()):
            idx = np.flatnonzero((strata == s).values)
            draws = rng.integers(0, idx.size, size=(B, idx.size))
            for b in range(B):
                W[b] += np.bincount(idx[draws[b]], minlength=len(groups)).astype(np.float32)
        return cls(groups, strata, W, seed)

    @property
    def B(self) -> int:
        return self.W.shape[0]


def _kg_from(ev: Evaluator, pairs: pd.DataFrame, endpoint: str):
    rows, value, weight, top = ev.endpoint_rows(pairs, endpoint)
    if len(rows) == 0:
        return None, top
    return ev.kg_table(rows, value, weight, top), top


def metric_replicates(ev: Evaluator, pairs: pd.DataFrame, endpoint: str, plan: BootstrapPlan,
                      chunk: int = 500) -> np.ndarray:
    """Vector (B,) of the endpoint recomputed on every resample."""
    kg, top = _kg_from(ev, pairs, endpoint)
    if kg is None:
        return np.full(plan.B, np.nan)
    gpos = pd.Series(np.arange(len(plan.groups)), index=plan.groups)
    g_e = gpos.reindex(kg["group"].astype(str)).to_numpy()
    if np.isnan(g_e.astype(float)).any():
        raise ValueError("pairs contain groups not in the bootstrap plan")
    g_e = g_e.astype(int)
    unit_key = (kg[top].astype(str) + "\x00" + kg["kinase"].astype(str)) if top else kg["kinase"]
    u_codes, u_uni = pd.factorize(unit_key)
    E = len(kg)
    U = sparse.csr_matrix((np.ones(E), (np.arange(E), u_codes)), shape=(E, len(u_uni)))
    if top:
        t_codes, _ = pd.factorize(pd.Series(u_uni).str.split("\x00").str[0])
        T = sparse.csr_matrix((np.ones(len(u_uni)), (np.arange(len(u_uni)), t_codes)),
                              shape=(len(u_uni), t_codes.max() + 1))
    v = kg["v"].to_numpy(float)
    out = np.empty(plan.B)
    for a in range(0, plan.B, chunk):
        Wg = plan.W[a:a + chunk][:, g_e].astype(float)       # (b, E)
        num = (U.T @ (Wg * v).T).T                          # (b, units)
        den = (U.T @ Wg.T).T
        present = den > 0
        um = np.where(present, num / np.where(present, den, 1), 0.0)
        if top:
            tsum = (T.T @ um.T).T
            tcnt = (T.T @ present.T.astype(float)).T
            tp = tcnt > 0
            tm = np.where(tp, tsum / np.where(tp, tcnt, 1), 0.0)
            out[a:a + chunk] = tm.sum(1) / np.maximum(tp.sum(1), 1)
        else:
            out[a:a + chunk] = um.sum(1) / np.maximum(present.sum(1), 1)
    return out


def _summ(point: float, reps: np.ndarray, alpha: float = 0.05) -> dict:
    """C8: basic (reverse-percentile) CI; raw percentile CI and bootstrap bias as diagnostics."""
    lo, hi = np.nanpercentile(reps, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    m = float(np.nanmean(reps))
    return {"estimate": float(point), "se": float(np.nanstd(reps, ddof=1)),
            "ci_low": float(2 * point - hi), "ci_high": float(2 * point - lo),
            "ci_low_pct": float(lo), "ci_high_pct": float(hi),
            "mean_replicates": m, "bias": m - float(point)}


def shifted_null_p(point: float, reps: np.ndarray) -> float:
    """C8 one-sided p for Δ > 0: (1 + #{Δ*_b − mean(Δ*) ≥ Δ̂}) / (B + 1); NaN replicates count
    against (as exceedances)."""
    reps = np.asarray(reps, float)
    centred = reps - np.nanmean(reps)
    exceed = ~(centred < point)          # ≥ point, or NaN
    return float((1 + exceed.sum()) / (reps.size + 1))


def percentile_p(reps: np.ndarray) -> float:
    """Pre-C8 diagnostic: (1 + #{Δ*_b ≤ 0}) / (B + 1) (biased when replicates are shifted)."""
    reps = np.asarray(reps, float)
    return float((1 + np.sum(~(reps > 0))) / (reps.size + 1))


def bootstrap_ci(ev: Evaluator, scores: pd.DataFrame, endpoint: str = "M1-sel",
                 plan: BootstrapPlan | None = None, pairs: pd.DataFrame | None = None) -> dict:
    plan = plan or BootstrapPlan.from_evaluator(ev)
    P = ev.pairs(scores) if pairs is None else pairs
    reps = metric_replicates(ev, P, endpoint, plan)
    return _summ(ev.metric(P, endpoint)["value"], reps) | {"replicates": reps}


def paired_delta(ev: Evaluator, scores_a: pd.DataFrame, scores_b: pd.DataFrame,
                 endpoint: str = "M1-sel", plan: BootstrapPlan | None = None) -> dict:
    """Δ = metric(A) − metric(B) on shared resamples; one-sided shifted-null p for Δ > 0 (C8),
    basic CI in ci_low/ci_high, raw percentile CI / p in *_pct (diagnostic).

    endpoint 'M2' (or 'M2:<filter>') evaluates both on the matched universe S(A)∩S(B).
    """
    plan = plan or BootstrapPlan.from_evaluator(ev)
    if endpoint.startswith("M2"):
        mask = ev.matched_mask(scores_a, scores_b)
        pa, pb = ev.pairs(scores_a, mask), ev.pairs(scores_b, mask)
        ep = endpoint.replace("M2", "M1-sel-scored", 1)
    else:
        pa, pb = ev.pairs(scores_a), ev.pairs(scores_b)
        ep = endpoint
    ra, rb = metric_replicates(ev, pa, ep, plan), metric_replicates(ev, pb, ep, plan)
    d = ra - rb
    point = ev.metric(pa, ep)["value"] - ev.metric(pb, ep)["value"]
    return _summ(point, d) | {"p_one_sided": shifted_null_p(point, d),
                              "p_one_sided_pct": percentile_p(d), "B": plan.B, "replicates": d}


def holm(pvals) -> np.ndarray:
    """Holm step-down adjusted p-values (monotone, capped at 1)."""
    p = np.asarray(pvals, float)
    m = p.size
    o = np.argsort(p)
    adj = np.maximum.accumulate((m - np.arange(m)) * p[o])
    out = np.empty(m)
    out[o] = np.minimum(adj, 1.0)
    return out


def beats(delta_m1: dict, holm_p: float, delta_m2_point: float, delta_unseen_point: float,
          delta_min: float = 0.01, alpha: float = 0.05) -> bool:
    """TEST 'beats' rule (§5 + B6): Holm p<α AND ΔM1≥δ_min AND ΔM1(T-unseen)>0 AND ΔM2≥0."""
    return bool(holm_p < alpha and delta_m1["estimate"] >= delta_min
                and delta_unseen_point > 0 and delta_m2_point >= 0)
