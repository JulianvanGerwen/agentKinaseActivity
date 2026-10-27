"""Evaluation on the external (TEST-only) benchmarKIN sets prepared by scripts/prep_external.py.

Sets: perturbation `tyrosine` (18 exps), `hernandez_new` (15 exps, 3 studies not in our compendium) and the CPTAC
tumour benchmark (7 cohorts). Score convention everywhere: DataFrame index = condition (experiment id / tumour sample),
columns = HGNC kinase symbols, unscored = NaN (never 0).

Metrics
  (i)  M1 ("ours"): perturbation = kinact.eval.Evaluator on synthetic GT tables (negatives = all other U kinases, unscored
       count 1/2, unscored target 1/2; pairs -> (kinase, group) -> kinase -> macro). Tumour = kinase-centric: per
       (cohort, kinase) AUROC of score over GS+ vs GS- samples (1/2 if the kinase is unscored), mean over cohorts, macro
       over kinases.
  (ii) benchmarKIN's own: perturbation = per-experiment RMS scaling, sign, pooled positives vs all other scored cells,
       decoupler mcauroc (1000 draws, n_pos negatives with replacement; experiments without a scored target dropped);
       tumour = per-cohort kinase z-score over samples, GS+/GS- pooled across cohorts, 1000 x 80% subsamples without
       replacement (kinases with >= 30 scored samples).
Real methods must only be run via `run_external(..., confirm_test=True)` at TEST time.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd

from .eval import Evaluator

ROOT = Path(os.environ.get("KINACT_ROOT", "/workspace"))
EXT = ROOT / "cache" / "external"
GT = ROOT / "results" / "ground_truth"
PERT_SETS = ("tyrosine", "hernandez_new")
COHORTS = ("BRCA", "CCRCC", "GBM", "HNSCC", "LSCC", "LUAD", "UCEC")


# ------------------------------------------------------------------ loading
def load_universe() -> list:
    return list(pd.read_csv(GT / "kinase_universe.tsv", sep="\t", keep_default_na=False)["kinase"])


def load_pert(name: str, ext: Path = EXT):
    """-> (site x experiment matrix, experiments table, targets table restricted to U)."""
    d = pd.read_parquet(ext / f"pert_{name}.parquet")
    e = pd.read_csv(ext / "pert_experiments.tsv", sep="\t", keep_default_na=False)
    t = pd.read_csv(ext / "pert_targets.tsv", sep="\t", keep_default_na=False)
    e, t = e[e.set == name].reset_index(drop=True), t[(t.set == name) & (t.in_U.astype(str) == "True")]
    return d, e, t.drop_duplicates(["experiment", "kinase"]).reset_index(drop=True)


def load_tumour(ext: Path = EXT):
    data = {c: pd.read_parquet(ext / f"tumour_{c}.parquet") for c in COHORTS}
    gs = pd.read_csv(ext / "tumour_gs_long.tsv", sep="\t", keep_default_na=False)
    gs = gs[gs.in_U.astype(str) == "True"].reset_index(drop=True)
    return data, gs


# ------------------------------------------------------------------ (i) M1 on perturbation sets
def pert_gt_tables(exps: pd.DataFrame, targets: pd.DataFrame):
    """Synthetic gt_pairs / gt_conditions for kinact.eval.Evaluator. Group = study (PMID) x drug."""
    gp = pd.DataFrame({"conditionID": targets.experiment, "kinase": targets.kinase, "cls": "direct",
                       "sign": targets.sign.astype(int), "confidence": "high", "rule_id": "ext", "note": ""})
    n = targets.groupby("experiment").size()
    gc = pd.DataFrame({"conditionID": exps.experiment, "PMID": exps.pmid, "source": "X",
                       "group_id": exps.pmid.astype(str) + "_" + exps.drug.astype(str),
                       "n_targets": exps.experiment.map(n).fillna(0).astype(int), "role": "gt", "split": "external"})
    gc["in_M1sel"] = (gc.n_targets >= 1) & (gc.n_targets <= 10)
    return gp, gc


def m1_perturbation(scores: pd.DataFrame, name: str, ext: Path = EXT, universe=None) -> dict:
    _, exps, tg = load_pert(name, ext)
    gp, gc = pert_gt_tables(exps, tg)
    ev = Evaluator(gp, gc, universe or load_universe(), zero_policy="ignore")
    P = ev.pairs(scores)
    seen = exps.set_index("experiment").drug_in_dev.astype(str).eq("True")
    P["drug_seen"] = P.conditionID.map(seen)
    out = {"M1-sel": ev.metric(P, "M1-sel"), "M1-all": ev.metric(P, "M1-all"),
           "coverage": float(P.target_scored.mean()), "n_pairs": int(len(P)),
           "n_groups": int(P.group.nunique())}
    for lab, m in (("seen", P.drug_seen), ("unseen", ~P.drug_seen)):
        sub = P[m]
        if len(sub):
            kg = ev.kg_table(sub, "auroc", None, None)
            out[f"M1-all_drug_{lab}"] = {"value": float(kg.groupby("kinase").v.mean().mean()),
                                         "n_pairs": int(len(sub))}
    out["pairs"] = P
    return out


# ------------------------------------------------------------------ (ii) benchmarKIN metrics
def auroc(pos, neg) -> float:
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    ps = np.sort(pos)
    lo, hi = np.searchsorted(ps, neg, "left"), np.searchsorted(ps, neg, "right")
    return float(((len(ps) - hi) + 0.5 * (hi - lo)).sum() / (len(ps) * len(neg)))


def mc_auroc(pos, neg, n_iter=1000, seed=42, n_pos=None, n_neg=None, replace=True) -> np.ndarray:
    """Default = decoupler mcauroc (all positives, len(pos) negatives drawn with replacement). With n_pos/n_neg and
    replace=False = benchmarKIN tumour ROC sampler."""
    rng = np.random.default_rng(seed)
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    if len(pos) == 0 or len(neg) == 0:
        return np.array([np.nan])
    out = np.empty(n_iter)
    for i in range(n_iter):
        p = pos if n_pos is None else rng.choice(pos, n_pos, replace=replace)
        q = rng.choice(neg, len(pos) if n_neg is None else n_neg, replace=replace)
        out[i] = auroc(p, q)
    return out


def rms_scale(A: pd.DataFrame) -> pd.DataFrame:
    """R scale(x, center=FALSE) per column (divide by sqrt(sum x^2/(n-1)) over scored kinases)."""
    n = A.notna().sum(axis=0)
    rms = np.sqrt((A ** 2).sum(axis=0) / (n - 1).clip(lower=1))
    return A / rms.replace(0, np.nan)


def benchmarkin_perturbation(scores: pd.DataFrame, name: str | None = None, ext: Path = EXT, exps=None, targets=None,
                             rm_bg: bool = True, n_iter: int = 1000, seed: int = 42) -> dict:
    """scores: experiment x kinase. Pooled AUROC (benchmarKIN run_perturbBench)."""
    if exps is None:
        _, exps, targets = load_pert(name, ext)
    ids = [e for e in exps.experiment if e in scores.index]
    A = scores.loc[ids].T                                   # kinase x experiment
    sign = exps.set_index("experiment").sign.reindex(ids).astype(float).to_numpy()
    S = rms_scale(A) * sign[None, :]
    S = S.mask(S == 0)
    G = pd.DataFrame(False, index=S.index, columns=S.columns)
    for e, k in zip(targets.experiment, targets.kinase):
        if e in G.columns and k in G.index:
            G.loc[k, e] = True
    if rm_bg:
        has = (G & S.notna()).any(axis=0)
        S, G = S.loc[:, has], G.loc[:, has]
    v, g = S.to_numpy(float).ravel(), G.to_numpy(bool).ravel()
    ok = np.isfinite(v)
    pos, neg = v[ok & g], v[ok & ~g]
    dist = mc_auroc(pos, neg, n_iter=n_iter, seed=seed)
    return {"auroc_mean": float(np.nanmean(dist)), "auroc_full": auroc(pos, neg), "n_pos": int(len(pos)),
            "n_neg": int(len(neg)), "n_experiments": int(S.shape[1])}


# ------------------------------------------------------------------ tumour
def _align(scores_by_cohort: dict, gs: pd.DataFrame, min_num: int = 30, zscore: bool = True) -> pd.DataFrame:
    rows = []
    for c, S in scores_by_cohort.items():
        A = S.T                                              # kinase x sample
        if min_num:
            A = A[A.notna().sum(axis=1) >= min_num]
        if zscore:
            A = A.sub(A.mean(axis=1), axis=0).div(A.std(axis=1, ddof=1), axis=0)
        g = gs[gs.cohort.eq(c) & gs.kinase.isin(A.index) & gs["sample"].isin(A.columns)]
        if g.empty:
            continue
        z = A.to_numpy(float)[A.index.get_indexer(g.kinase), A.columns.get_indexer(g["sample"])]
        rows.append(g.assign(z=z))
    if not rows:
        return pd.DataFrame(columns=list(gs.columns) + ["z"])
    return pd.concat(rows, ignore_index=True).dropna(subset=["z"])


def benchmarkin_tumour(scores_by_cohort: dict, gs: pd.DataFrame, min_num: int = 30, sub_frac: float = 0.8,
                       n_iter: int = 1000, seed: int = 1) -> dict:
    d = _align(scores_by_cohort, gs, min_num)
    both = d.groupby(["cohort", "kinase"]).label.agg(lambda s: s.nunique() == 2)
    keep = set(both[both].index)
    d = d[[(c, k) in keep for c, k in zip(d.cohort, d.kinase)]]
    pos, neg = d.z[d.label == 1].to_numpy(), d.z[d.label == 0].to_numpy()
    if not len(pos) or not len(neg):
        return {"auroc_mean": np.nan, "n_pos": 0, "n_neg": 0}
    dist = mc_auroc(pos, neg, n_iter, seed, n_pos=int(np.floor(sub_frac * len(pos))),
                    n_neg=int(np.floor(sub_frac * len(neg))), replace=False)
    return {"auroc_mean": float(np.mean(dist)), "auroc_full": auroc(pos, neg), "n_pos": int(len(pos)),
            "n_neg": int(len(neg)), "n_kinases": int(d.kinase.nunique())}


def m1_tumour(scores_by_cohort: dict, gs: pd.DataFrame) -> dict:
    """Kinase-centric macro AUROC: (cohort, kinase) with GS+ and GS- samples; unscored kinase/sample -> 1/2 for the
    kinase-cohort; mean over cohorts, then macro over kinases."""
    rec = []
    for (c, k), g in gs.groupby(["cohort", "kinase"]):
        if g.label.nunique() < 2 or c not in scores_by_cohort:
            continue
        S = scores_by_cohort[c]
        v = S[k].reindex(g["sample"].values).to_numpy(float) if k in S.columns else np.full(len(g), np.nan)
        lab = g.label.to_numpy()
        ok = np.isfinite(v)
        a = auroc(v[ok & (lab == 1)], v[ok & (lab == 0)]) if (ok & (lab == 1)).any() and (ok & (lab == 0)).any() else np.nan
        rec.append((c, k, 0.5 if np.isnan(a) else a, not np.isnan(a)))
    r = pd.DataFrame(rec, columns=["cohort", "kinase", "auroc", "scored"])
    if r.empty:
        return {"M1-kinase": np.nan, "n_kinases": 0}
    km = r.groupby("kinase").auroc.mean()
    return {"M1-kinase": float(km.mean()), "n_kinases": int(len(km)), "coverage": float(r.scored.mean()),
            "per_kinase": km}


# ------------------------------------------------------------------ synthetic controls (validation only)
def random_scores(index, universe=None, seed=0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    cols = universe or load_universe()
    return pd.DataFrame(rng.standard_normal((len(index), len(cols))), index=list(index), columns=cols)


def oracle_pert_scores(name: str, ext: Path = EXT, seed=0, universe=None) -> pd.DataFrame:
    """Targets get sign*(5+noise)... only a harness check; labels are the answer."""
    _, exps, tg = load_pert(name, ext)
    S = random_scores(exps.experiment, universe, seed)
    for e, k, s in zip(tg.experiment, tg.kinase, tg.sign):
        S.loc[e, k] = s * 10.0
    return S


def oracle_tumour_scores(gs: pd.DataFrame, data: dict, seed=0, universe=None) -> dict:
    out = {}
    for c, d in data.items():
        S = random_scores(d.columns, universe, seed)
        g = gs[gs.cohort == c]
        for k, s, lab in zip(g.kinase, g["sample"], g.label):
            if s in S.index:
                S.loc[s, k] = 10.0 if lab == 1 else -10.0
        out[c] = S
    return out


# ------------------------------------------------------------------ TEST-time driver
def run_external(score_fn, confirm_test: bool = False, ext: Path = EXT, n_iter: int = 1000) -> dict:
    """score_fn(site x condition DataFrame) -> condition x kinase DataFrame (NaN = unscored).
    TEST-time only: refuses to run without confirm_test=True (these sets are part of the final TEST)."""
    if not confirm_test:
        raise RuntimeError("external sets are TEST data; pass confirm_test=True from the PI-launched TEST script")
    res = {}
    for name in PERT_SETS:
        d, exps, tg = load_pert(name, ext)
        S = score_fn(d)
        res[name] = {"M1": {k: v for k, v in m1_perturbation(S, name, ext).items() if k != "pairs"},
                     "benchmarkin": benchmarkin_perturbation(S, exps=exps, targets=tg, n_iter=n_iter)}
    data, gs = load_tumour(ext)
    sc = {c: score_fn(d) for c, d in data.items()}
    res["tumour"] = {"M1": {k: v for k, v in m1_tumour(sc, gs).items() if k != "per_kinase"},
                     "benchmarkin": benchmarkin_tumour(sc, gs, n_iter=n_iter)}
    return res
