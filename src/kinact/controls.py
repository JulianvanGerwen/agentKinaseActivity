"""Negative / positive controls (DESIGN §5, §13 B3/B7, r3 S10).

All control score frames are conditions × kinases. Controls that legitimately emit exact
zeros (oracle, prior-only) must be evaluated with Evaluator(..., zero_policy="ignore").
"""
from __future__ import annotations

from typing import Callable, Iterable, Sequence

import numpy as np
import pandas as pd

from .eval import EVAL_CONFIDENCE


def _direct(gt_pairs: pd.DataFrame, confidence=EVAL_CONFIDENCE) -> pd.DataFrame:
    d = gt_pairs[(gt_pairs["cls"].astype(str) == "direct")
                 & gt_pairs["confidence"].astype(str).str.lower().isin(confidence)]
    return d.drop_duplicates(["conditionID", "kinase"])


def random_scores(conditions: Sequence[str], kinases: Sequence[str], seed: int = 0,
                  coverage: float = 1.0) -> pd.DataFrame:
    """N(0,1) scores; a fraction (1-coverage) set to NaN (unscored). Expected M1 = 0.5."""
    rng = np.random.default_rng(seed)
    z = rng.standard_normal((len(conditions), len(kinases)))
    if coverage < 1:
        z[rng.random(z.shape) >= coverage] = np.nan
    return pd.DataFrame(z, index=list(conditions), columns=list(kinases))


def oracle(gt_pairs: pd.DataFrame, conditions: Sequence[str], kinases: Sequence[str],
           anti: bool = False, confidence=EVAL_CONFIDENCE) -> pd.DataFrame:
    """score = s for direct positives, 0 elsewhere (M1 = 1; anti-oracle: M1 = 0)."""
    S = pd.DataFrame(index=list(conditions), columns=list(kinases))
    d = _direct(gt_pairs, confidence)
    d = d[d["conditionID"].isin(S.index) & d["kinase"].isin(S.columns)]
    ri = S.index.get_indexer(d["conditionID"])
    ci = S.columns.get_indexer(d["kinase"])
    v = np.zeros(S.shape)
    v[ri, ci] = pd.to_numeric(d["sign"]).to_numpy(float) * (-1 if anti else 1)
    return pd.DataFrame(v, index=S.index, columns=S.columns)


def training_frequency(gt_pairs: pd.DataFrame, gt_conditions: pd.DataFrame,
                       train_conditions: Iterable[str], kinases: Sequence[str],
                       confidence=EVAL_CONFIDENCE) -> pd.Series:
    """freq_k = fraction of training GT conditions (role=='gt', n_targets>0) in which k is a
    direct positive. Uses only the given training conditions (never TEST / held-out fold)."""
    tc = set(train_conditions)
    gc = gt_conditions[gt_conditions["conditionID"].isin(tc)
                       & (gt_conditions["role"].astype(str) == "gt")
                       & (pd.to_numeric(gt_conditions["n_targets"]) > 0)]
    if "split" in gc.columns and (gc["split"].astype(str) == "test").any():
        raise ValueError("training conditions include split=='test'")
    n = gc["conditionID"].nunique()
    d = _direct(gt_pairs, confidence)
    d = d[d["conditionID"].isin(set(gc["conditionID"]))]
    cnt = d.groupby("kinase")["conditionID"].nunique()
    return (cnt / max(n, 1)).reindex(list(kinases)).fillna(0.0)


def condition_sign(gt_pairs: pd.DataFrame, conditions: Sequence[str],
                   confidence=EVAL_CONFIDENCE) -> pd.Series:
    """s_c = majority sign of c's direct positives (tie -> +1; none -> 0)."""
    d = _direct(gt_pairs, confidence)
    m = pd.to_numeric(d["sign"]).groupby(d["conditionID"]).sum()
    s = np.sign(m).replace(0, 1)
    return s.reindex(list(conditions)).fillna(0.0)


def prior_only(gt_pairs, gt_conditions, train_conditions, eval_conditions, kinases,
               confidence=EVAL_CONFIDENCE) -> pd.DataFrame:
    """B3 prior-only control: score_ck = s_c · freq_k (oracle condition sign, no phospho)."""
    f = training_frequency(gt_pairs, gt_conditions, train_conditions, kinases, confidence)
    s = condition_sign(gt_pairs, eval_conditions, confidence)
    return pd.DataFrame(np.outer(s.values, f.values), index=list(eval_conditions),
                        columns=list(kinases))


def shift_prior(X: pd.DataFrame, freq: pd.Series) -> pd.DataFrame:
    """B3 shift-prior control (label-free): score_ck = freq_k · mean_sites(x_c)."""
    mu = X.mean(axis=0, skipna=True)
    return pd.DataFrame(np.outer(mu.values, freq.values), index=X.columns, columns=freq.index)


def site_permutation(X: pd.DataFrame, seed: int = 0, mode: str = "within_condition") -> pd.DataFrame:
    """Site-permuted input (sites × conditions).
    within_condition: shuffle observed values among observed sites, independently per column
    (missingness pattern kept). rows: permute site labels jointly (keeps co-variation)."""
    rng = np.random.default_rng(seed)
    if mode == "rows":
        return pd.DataFrame(X.to_numpy()[rng.permutation(X.shape[0])], index=X.index,
                            columns=X.columns)
    v = X.to_numpy(float).copy()
    for j in range(v.shape[1]):
        ok = np.flatnonzero(~np.isnan(v[:, j]))
        v[ok, j] = v[rng.permutation(ok), j]
    return pd.DataFrame(v, index=X.index, columns=X.columns)


def shuffled_library(lib: pd.DataFrame, seed: int = 0, mode: str = "relabel") -> pd.DataFrame:
    """Degree-preserving shuffled-library control (r3 S10).
    relabel: random bijection on the library's site labels -> per-kinase substrate counts,
    weights and kinase–kinase substrate overlap preserved exactly, no duplicate edges.
    edges: permute the site column across edges (counts preserved; may create duplicates,
    which are dropped, so counts can shrink slightly)."""
    rng = np.random.default_rng(seed)
    out = lib.copy()
    if mode == "relabel":
        sites = pd.unique(lib["site"])
        mp = dict(zip(sites, sites[rng.permutation(len(sites))]))
        out["site"] = lib["site"].map(mp)
        return out
    out["site"] = lib["site"].to_numpy()[rng.permutation(len(lib))]
    return out.drop_duplicates(["kinase", "site"])


def sign_flip_check(method_fn: Callable, X: pd.DataFrame, lib=None, atol: float = 1e-8,
                    rtol: float = 1e-6) -> dict:
    """Negating the input must negate the scores (NaN pattern identical)."""
    a = method_fn(X, lib)
    b = method_fn(-X, lib)
    b = b.reindex(index=a.index, columns=a.columns)
    va, vb = a.to_numpy(float), b.to_numpy(float)
    nan_ok = bool(np.array_equal(np.isnan(va), np.isnan(vb)))
    ok = ~np.isnan(va) & ~np.isnan(vb)
    dev = float(np.max(np.abs(va[ok] + vb[ok]))) if ok.any() else 0.0
    scale = float(np.max(np.abs(va[ok]))) if ok.any() else 0.0
    return {"passed": nan_ok and dev <= atol + rtol * scale, "max_abs_dev": dev,
            "nan_pattern_equal": nan_ok}
