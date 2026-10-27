"""Frozen evaluation harness (DESIGN.md v2 §5 + §13; docs/harness.md).

Core quantities
---------------
t = s * score (s = GT sign of the target; larger t = more correct).
Pair AUROC for target k in condition c with evaluable negatives N_c (U minus every kinase
listed for c in gt_pairs, whatever its class/confidence) and scored set S_c:
    k in S_c : [ sum_{j in N_c∩S_c} (1[t_j<t_k] + ½·1[t_j=t_k]) + ½·|N_c \\ S_c| ] / |N_c|
    k not in S_c : ½
Averaging: pairs -> (weighted) mean within (kinase, GT group) -> mean over groups per kinase
-> macro over kinases.
"""
from __future__ import annotations

import json
import os
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
from scipy import sparse, stats

DEFAULT_ROOT = Path(os.environ.get("KINACT_ROOT", "/workspace"))
EVAL_CONFIDENCE = ("high", "medium")
BIN_EDGES = (1, 5, 20)  # bins: 0 | 1-4 | 5-19 | >=20
BIN_LABELS = ("0", "1-4", "5-19", ">=20")


class TestSplitLocked(RuntimeError):
    """Raised when TEST conditions are touched before FREEZE.json exists."""
    __test__ = False


def freeze_path(path: str | os.PathLike | None = None) -> Path:
    if path is not None:
        return Path(path)
    return Path(os.environ.get("KINACT_FREEZE_PATH", DEFAULT_ROOT / "FREEZE.json"))


def check_split_access(splits: Iterable[str], freeze: str | os.PathLike | None = None) -> None:
    """Raise TestSplitLocked if 'test' is among `splits` and FREEZE.json is absent."""
    if "test" in {str(s) for s in splits} and not freeze_path(freeze).exists():
        raise TestSplitLocked(
            f"evaluation touches split=='test' but {freeze_path(freeze)} does not exist")


def psp_bin(n) -> np.ndarray:
    """Map (possibly fractional, e.g. median) PSP-substrate counts to bin labels."""
    n = np.asarray(n, dtype=float)
    idx = np.searchsorted(np.asarray(BIN_EDGES, float), n, side="right")
    out = np.asarray(BIN_LABELS, dtype=object)[idx]
    out[np.isnan(n)] = "NA"
    return out


def read_gt_tsv(path) -> pd.DataFrame:
    """Read a GT/split TSV without pandas' default NA strings ("null", "NA", "None", "nan" stay
    literal strings, e.g. role == "null"); only truly empty fields become NaN."""
    return pd.read_csv(path, sep="\t", keep_default_na=False, na_values=[""])


def load_gt(gt_dir: str | os.PathLike = DEFAULT_ROOT / "results/ground_truth"):
    """Agent-visible (DEV + devconfirm) GT: (gt_pairs, gt_conditions, universe, bins|None).
    Accepts gt_pairs.tsv or gt_pairs_dev.tsv naming; all readers keep literal "null"/"NA" strings. TEST rows: see load_sealed_test()."""
    gt_dir = Path(gt_dir)

    def rd(stem, required=True):
        for n in (f"{stem}.tsv", f"{stem}_dev.tsv"):
            if (gt_dir / n).exists():
                return read_gt_tsv(gt_dir / n)
        if required:
            raise FileNotFoundError(gt_dir / f"{stem}.tsv")
        return None
    return rd("gt_pairs"), rd("gt_conditions"), rd("kinase_universe"), rd("harmony_bins", False)


def _as_bool(s: pd.Series) -> pd.Series:
    if s.dtype == bool:
        return s
    return s.astype(str).str.strip().str.lower().isin(["true", "1", "yes", "t"])


# --------------------------------------------------------------------------------------
# Evaluator
# --------------------------------------------------------------------------------------
@dataclass
class Evaluator:
    """Pre-indexed GT for a set of splits; evaluate any number of score matrices on it.

    gt_pairs: conditionID, kinase, cls ∈ {direct, downstream, masked}, sign, confidence
    gt_conditions: conditionID, source, group_id, n_targets, in_M1sel, role, split
    universe: iterable of kinase symbols (U) or the kinase_universe DataFrame
    split_filter: splits to include (None = all present -> includes 'test' if present)
    bins: optional DataFrame conditionID, kinase, n_psp (quantified PSP substrates of
          kinase in condition). Harmony bin per (kinase, group) = bin(median n_psp).
    train_pos_kinases: kinases that are direct positives in the training data (T-unseen).
    """

    gt_pairs: pd.DataFrame
    gt_conditions: pd.DataFrame
    universe: Sequence[str] | pd.DataFrame
    split_filter: Sequence[str] | None = None
    bins: pd.DataFrame | None = None
    train_pos_kinases: Iterable[str] | None = None
    confidence: Sequence[str] = EVAL_CONFIDENCE
    freeze: str | os.PathLike | None = None
    conditions: Sequence[str] | None = None  # optional further restriction (e.g. one fold)
    positive_cls: str = "direct"  # "downstream" -> secondary downstream-positive AUROC
    zero_policy: str = "raise"  # NaN->0 bug check: "raise" | "warn" | "ignore"

    def __post_init__(self):
        if isinstance(self.universe, pd.DataFrame):
            self.universe = list(self.universe["kinase"])
        self.kinases = pd.Index(pd.unique(np.asarray(self.universe, dtype=object)))
        kpos = {k: i for i, k in enumerate(self.kinases)}

        gc = self.gt_conditions.copy()
        if self.split_filter is not None:
            sf = [self.split_filter] if isinstance(self.split_filter, str) else list(self.split_filter)
            gc = gc[gc["split"].astype(str).isin(sf)]
        if self.conditions is not None:
            gc = gc[gc["conditionID"].isin(set(self.conditions))]
        check_split_access(gc["split"].astype(str).unique(), self.freeze)
        gc = gc[gc["role"].astype(str) == "gt"]
        gc = gc[pd.to_numeric(gc["n_targets"]) > 0]
        gc = gc.drop_duplicates("conditionID").set_index("conditionID")
        gc["in_M1sel"] = _as_bool(gc["in_M1sel"])

        gp = self.gt_pairs[self.gt_pairs["conditionID"].isin(gc.index)].copy()
        gp["cls"] = gp["cls"].astype(str)
        gp["confidence"] = gp["confidence"].astype(str).str.lower()
        direct = gp[(gp["cls"] == self.positive_cls) & gp["confidence"].isin(self.confidence)]
        bad = direct[~direct["kinase"].isin(kpos)]
        if len(bad):
            raise ValueError(f"{len(bad)} direct GT pairs with kinase outside U: "
                             f"{sorted(bad['kinase'].unique())[:10]}")
        if self.positive_cls != "direct":
            direct = direct[pd.to_numeric(direct["sign"]) != 0]
        sgn = pd.to_numeric(direct["sign"])
        if not set(np.unique(sgn)) <= {-1, 1}:
            raise ValueError("direct positives must have sign in {-1,+1} (B5)")
        direct = direct.assign(sign=sgn.astype(int)).drop_duplicates(["conditionID", "kinase"])
        conflict = direct.duplicated(["conditionID", "kinase"], keep=False)
        assert not conflict.any()

        # conditions that end up with no evaluable direct positives are dropped
        gc = gc[gc.index.isin(direct["conditionID"].unique())]
        self.cond_table = gc
        self.cond_index = pd.Index(gc.index)
        cpos = {c: i for i, c in enumerate(self.cond_index)}

        # excluded-from-negatives: any kinase listed for the condition (all cls/confidence)
        excl = gp[gp["kinase"].isin(kpos)]
        n_c, n_k = len(self.cond_index), len(self.kinases)
        neg = np.ones((n_c, n_k), dtype=bool)
        ci = excl["conditionID"].map(cpos)
        ok = ci.notna()
        neg[ci[ok].astype(int).values, excl.loc[ok, "kinase"].map(kpos).astype(int).values] = False
        self.neg_mask = neg

        direct = direct[direct["conditionID"].isin(cpos)]
        groups = gc["group_id"].astype(str)
        self.group_index = pd.Index(sorted(groups.unique()))
        gpos = {g: i for i, g in enumerate(self.group_index)}
        tgt = pd.DataFrame({
            "conditionID": direct["conditionID"].values,
            "kinase": direct["kinase"].values,
            "sign": direct["sign"].values,
        })
        tgt["ci"] = tgt["conditionID"].map(cpos).astype(int)
        tgt["ki"] = tgt["kinase"].map(kpos).astype(int)
        tgt["group"] = tgt["conditionID"].map(groups).values
        tgt["gi"] = tgt["group"].map(gpos).astype(int)
        tgt["source"] = tgt["conditionID"].map(gc["source"].astype(str)).values
        tgt["in_M1sel"] = tgt["conditionID"].map(gc["in_M1sel"]).values.astype(bool)
        tgt["n_targets"] = tgt["conditionID"].map(pd.to_numeric(gc["n_targets"])).values
        tgt["w_all"] = 1.0 / tgt["n_targets"]
        tgt["n_neg"] = neg[tgt["ci"].values].sum(1)
        tgt = tgt[tgt["n_neg"] > 0].sort_values(["ci", "ki"]).reset_index(drop=True)
        self.targets = tgt

        # harmony bins per (kinase, group)
        if self.bins is not None:
            lab = harmony_bin_labels(self.bins, tgt)
            key = pd.MultiIndex.from_arrays([tgt["kinase"], tgt["group"]])
            tgt["bin"] = lab.reindex(key).fillna("NA").values
            if (tgt["bin"] == "NA").any():
                warnings.warn(f"{int((tgt['bin'] == 'NA').sum())} GT pairs without a PSP bin")
        else:
            tgt["bin"] = "NA"

        # bootstrap stratum per group = majority source by M1-sel pairs (fallback all pairs)
        def _major(df):
            sel = df[df["in_M1sel"]]
            use = sel if len(sel) else df
            return use["source"].value_counts().sort_index().idxmax()
        self.group_stratum = (tgt.groupby("group").apply(_major, include_groups=False)
                              .reindex(self.group_index))
        self.train_pos = None if self.train_pos_kinases is None else set(self.train_pos_kinases)

    # ------------------------------------------------------------------ pair AUROCs
    def _matrix(self, scores: pd.DataFrame) -> np.ndarray:
        S = scores.reindex(index=self.cond_index, columns=self.kinases)
        return S.to_numpy(dtype=float)

    def pairs(self, scores: pd.DataFrame, mask: np.ndarray | None = None) -> pd.DataFrame:
        """Per-pair table with auroc (½-rule) and auroc_scored (scored negatives only)."""
        if mask is None:
            check_unscored_encoding(scores.reindex(index=self.cond_index.intersection(scores.index)),
                                    self.zero_policy)
        S = self._matrix(scores)
        if mask is not None:
            S = np.where(mask, S, np.nan)
        tgt = self.targets
        n = len(tgt)
        auroc = np.full(n, 0.5)
        aur_s = np.full(n, np.nan)
        scored = np.zeros(n, bool)
        n_neg_sc = np.zeros(n, int)
        tie = np.zeros(n, bool)
        ci_arr, ki_arr, sg_arr = tgt["ci"].values, tgt["ki"].values, tgt["sign"].values
        bounds = np.flatnonzero(np.r_[True, ci_arr[1:] != ci_arr[:-1], True])
        for a, b in zip(bounds[:-1], bounds[1:]):
            c = ci_arr[a]
            v = S[c, self.neg_mask[c]]
            N = v.size
            sv = np.sort(v[~np.isnan(v)])
            m = sv.size
            nu = N - m
            sk = S[c, ki_arr[a:b]]
            ok = ~np.isnan(sk)
            lo = np.searchsorted(sv, sk, "left")
            hi = np.searchsorted(sv, sk, "right")
            ties = hi - lo
            wins = np.where(sg_arr[a:b] > 0, lo, m - hi) + 0.5 * ties
            auroc[a:b] = np.where(ok, (wins + 0.5 * nu) / N, 0.5)
            if m > 0:
                aur_s[a:b] = np.where(ok, wins / max(m, 1), np.nan)
            scored[a:b] = ok
            n_neg_sc[a:b] = m
            tie[a:b] = ok & (ties > 0)
        out = tgt[["conditionID", "kinase", "group", "source", "sign", "in_M1sel",
                   "n_targets", "w_all", "bin", "n_neg"]].copy()
        out["n_neg_scored"] = n_neg_sc
        out["target_scored"] = scored
        out["auroc"] = auroc
        out["auroc_scored"] = aur_s
        out["tie"] = tie
        return out

    # ------------------------------------------------------------------ endpoints
    def endpoint_rows(self, pairs: pd.DataFrame, endpoint: str):
        """Return (rows, value column, weight column or None, top-level column or None).

        endpoint grammar: base[:filter]   base ∈ {M1-sel, M1-all, M1-sel-scored,
        M1-all-scored, harmony}; filter ∈ {source=H, sign=-1, unseen, bin=1-4}.
        """
        base, _, filt = endpoint.partition(":")
        p = pairs
        weight, top = None, None
        if base.startswith("M1-sel") or base == "harmony" or base.startswith("M2"):
            p = p[p["in_M1sel"]]
        elif base.startswith("M1-all"):
            weight = "w_all"
        else:
            raise ValueError(f"unknown endpoint {endpoint}")
        value = "auroc"
        if base.endswith("scored") or base.startswith("M2"):
            value = "auroc_scored"
            p = p[p["target_scored"] & (p["n_neg_scored"] > 0)]
        if base == "harmony":
            top = "bin"
            p = p[p["bin"] != "NA"]
        for f in filter(None, filt.split(",")):
            if f == "unseen":
                if self.train_pos is None:
                    raise ValueError("T-unseen needs train_pos_kinases")
                p = p[~p["kinase"].isin(self.train_pos)]
            else:
                k, _, v = f.partition("=")
                col = p[k]
                p = p[col.astype(str) == (str(int(v)) if k == "sign" else v)]
        return p, value, weight, top

    @staticmethod
    def kg_table(rows: pd.DataFrame, value: str, weight: str | None, top: str | None):
        """(kinase, group[, top]) weighted means of the pair values."""
        w = rows[weight] if weight else pd.Series(1.0, index=rows.index)
        keys = ["kinase", "group"] + ([top] if top else [])
        d = pd.DataFrame({k: rows[k] for k in keys})
        d["wv"] = w * rows[value]
        d["w"] = w
        g = d.groupby(keys, sort=True)[["wv", "w"]].sum()
        g["v"] = g["wv"] / g["w"]
        return g.reset_index()

    def metric(self, pairs: pd.DataFrame, endpoint: str) -> dict:
        rows, value, weight, top = self.endpoint_rows(pairs, endpoint)
        if len(rows) == 0:
            return {"value": np.nan, "n_kinases": 0, "n_pairs": 0}
        kg = self.kg_table(rows, value, weight, top)
        if top:
            km = kg.groupby([top, "kinase"])["v"].mean()
            per_bin = km.groupby(level=0).mean()
            res = {"value": float(per_bin.mean()), "n_kinases": int(kg["kinase"].nunique()),
                   "n_pairs": int(len(rows)),
                   "per_bin": {b: float(v) for b, v in per_bin.items()},
                   "n_kinases_per_bin": {b: int(v) for b, v in km.groupby(level=0).size().items()}}
            return res
        km = kg.groupby("kinase")["v"].mean()
        return {"value": float(km.mean()), "n_kinases": int(len(km)), "n_pairs": int(len(rows))}

    # ------------------------------------------------------------------ M2
    def matched_mask(self, scores_a: pd.DataFrame, scores_b: pd.DataFrame) -> np.ndarray:
        return ~np.isnan(self._matrix(scores_a)) & ~np.isnan(self._matrix(scores_b))

    def m2(self, scores_a, scores_b, base: str = "M1-sel") -> dict:
        """M2: universe in c = S_c(A) ∩ S_c(B) for targets and negatives."""
        mask = self.matched_mask(scores_a, scores_b)
        pa, pb = self.pairs(scores_a, mask), self.pairs(scores_b, mask)
        ep = base + "-scored"
        ra, rb = self.metric(pa, ep), self.metric(pb, ep)
        return {"A": ra["value"], "B": rb["value"], "delta": ra["value"] - rb["value"],
                "n_pairs": ra["n_pairs"], "n_kinases": ra["n_kinases"],
                "pair_fraction": ra["n_pairs"] / max(1, int(self.targets["in_M1sel"].sum())
                                                     if base == "M1-sel" else len(self.targets))}

    # ------------------------------------------------------------------ coverage-significant
    def coverage_significant(self, pairs: pd.DataFrame, endpoint: str = "M1-sel",
                             n_perm: int = 2000, q: float = 0.10, seed: int = 0) -> pd.DataFrame:
        """Per-kinase permutation p (target rank uniform among scored negatives ∪ {k},
        unscored negatives fixed at ½, group structure kept) and BH q-values."""
        rows, value, weight, _ = self.endpoint_rows(pairs, endpoint)
        if len(rows) == 0:
            return pd.DataFrame(columns=["kinase", "auroc", "p_perm", "q_bh", "significant"])
        rng = np.random.default_rng(seed)
        m = rows["n_neg_scored"].values
        N = rows["n_neg"].values
        sc = rows["target_scored"].values
        if value == "auroc_scored":
            null = rng.integers(0, m[:, None] + 1, size=(len(rows), n_perm)) / np.maximum(m, 1)[:, None]
        else:
            r = rng.integers(0, m[:, None] + 1, size=(len(rows), n_perm))
            null = np.where(sc[:, None], (r + 0.5 * (N - m)[:, None]) / N[:, None], 0.5)
        w = rows[weight].values if weight else np.ones(len(rows))
        kg_key = pd.MultiIndex.from_arrays([rows["kinase"].values, rows["group"].values])
        kg_codes, kg_uni = pd.factorize(kg_key)
        A = sparse.csr_matrix((w, (kg_codes, np.arange(len(rows)))), shape=(len(kg_uni), len(rows)))
        wsum = np.asarray(A.sum(1)).ravel()
        kg_obs = (A @ rows[value].values) / wsum
        kg_null = (A @ null) / wsum[:, None]
        k_codes, k_uni = pd.factorize(kg_uni.get_level_values(0))
        K = sparse.csr_matrix((np.ones(len(k_codes)), (k_codes, np.arange(len(k_codes)))),
                              shape=(len(k_uni), len(k_codes)))
        cnt = np.asarray(K.sum(1)).ravel()
        k_obs = (K @ kg_obs) / cnt
        k_null = (K @ kg_null) / cnt[:, None]
        p = (1 + (k_null >= k_obs[:, None] - 1e-12).sum(1)) / (n_perm + 1)
        any_scored = pd.Series(sc).groupby(rows["kinase"].values).any().reindex(k_uni).values
        p = np.where(any_scored, p, 1.0)
        qv = bh(p)
        return pd.DataFrame({"kinase": k_uni, "auroc": k_obs, "p_perm": p, "q_bh": qv,
                             "significant": (qv <= q) & (k_obs > 0.5)})

    # ------------------------------------------------------------------ full evaluation
    def evaluate(self, scores: pd.DataFrame, n_perm: int = 2000, q: float = 0.10) -> dict:
        if len(self.targets) == 0:
            raise ValueError("no evaluable GT pairs for this split/confidence selection")
        P = self.pairs(scores)
        out: dict = {}
        for ep in ["M1-sel", "M1-all", "M1-sel-scored", "M1-all-scored"]:
            out[ep] = self.metric(P, ep)
        sel = P[P["in_M1sel"]]
        out["coverage_M1-sel"] = float(sel["target_scored"].mean()) if len(sel) else np.nan
        out["coverage_M1-all"] = float(P["target_scored"].mean()) if len(P) else np.nan
        out["tie_share_M1-sel"] = float(sel.loc[sel["target_scored"], "tie"].mean()) \
            if sel["target_scored"].any() else np.nan
        out["by_source"] = {s: self.metric(P, f"M1-sel:source={s}") for s in sorted(P["source"].unique())}
        out["by_source_all"] = {s: self.metric(P, f"M1-all:source={s}") for s in sorted(P["source"].unique())}
        out["by_sign"] = {int(s): self.metric(P, f"M1-sel:sign={int(s)}") for s in sorted(P["sign"].unique())}
        if (P["bin"] != "NA").any():
            out["harmony"] = self.metric(P, "harmony")
        if self.train_pos is not None:
            out["T-unseen"] = self.metric(P, "M1-sel:unseen")
            out["T-unseen_all"] = self.metric(P, "M1-all:unseen")
        cs = self.coverage_significant(P, "M1-sel", n_perm=n_perm, q=q)
        nsig = int(cs["significant"].sum())
        out["coverage_significant"] = {"n_significant": nsig, "n_kinases": int(len(cs)),
                                       "expected_false_discoveries": q * nsig}
        if (P["bin"] != "NA").any() and len(cs):
            kb = (P[P["in_M1sel"]].groupby("kinase")["bin"].agg(lambda s: s.mode().iloc[0]))
            cs["bin"] = cs["kinase"].map(kb)
            out["coverage_significant"]["per_bin"] = {
                b: int(v) for b, v in cs.groupby("bin")["significant"].sum().items()}
        out["per_kinase"] = self.per_kinase(P, cs)
        out["pairs"] = P
        return out

    def per_kinase(self, P: pd.DataFrame, cs: pd.DataFrame | None = None) -> pd.DataFrame:
        cols = {}
        for ep in ["M1-sel", "M1-all", "M1-sel-scored"]:
            rows, value, weight, _ = self.endpoint_rows(P, ep)
            if len(rows):
                kg = self.kg_table(rows, value, weight, None)
                cols[ep] = kg.groupby("kinase")["v"].mean()
        g = P.groupby("kinase")
        tab = pd.DataFrame({
            "n_pairs": g.size(), "n_conditions": g["conditionID"].nunique(),
            "n_groups": g["group"].nunique(), "n_pairs_M1sel": g["in_M1sel"].sum(),
            "pair_coverage": g["target_scored"].mean(),
            "sources": g["source"].agg(lambda s: ",".join(sorted(set(s)))),
            "bin": g["bin"].agg(lambda s: s.mode().iloc[0]),
        })
        for k, v in cols.items():
            tab[k] = v
        if cs is not None and len(cs):
            tab = tab.join(cs.set_index("kinase")[["p_perm", "q_bh", "significant"]])
        if self.train_pos is not None:
            tab["seen_in_train"] = tab.index.isin(self.train_pos)
        return tab.reset_index().rename(columns={"index": "kinase"})


def evaluate(scores, gt_pairs, gt_conditions, split_filter=None, universe=None, bins=None,
             train_pos_kinases=None, freeze=None, n_perm=2000, **kw) -> dict:
    """Contract entry point (docs/interfaces.md)."""
    if universe is None:
        universe = sorted(set(gt_pairs["kinase"]) | set(scores.columns))
    ev = Evaluator(gt_pairs, gt_conditions, universe, split_filter, bins=bins,
                   train_pos_kinases=train_pos_kinases, freeze=freeze, **kw)
    return ev.evaluate(scores, n_perm=n_perm)


def harmony_bin_labels(bins: pd.DataFrame, tgt: pd.DataFrame) -> pd.Series:
    """Bin label per (kinase, group): accepts results/ground_truth/harmony_bins.tsv keyed by
    (kinase, group_id) with a frozen `bin` column, or per (conditionID, kinase) counts
    (n_psp_quantified | n_psp) -> median over the group's positive conditions -> bin."""
    b = bins.copy()
    ncol = next((c for c in ("n_psp_quantified", "n_psp") if c in b.columns), None)
    if "conditionID" in b.columns:
        b = b.merge(tgt[["conditionID", "kinase", "group"]], on=["conditionID", "kinase"])
        if ncol is None:
            raise ValueError("per-condition bin table needs n_psp_quantified")
        med = b.groupby(["kinase", "group"])[ncol].median()
        return pd.Series(psp_bin(med.values), index=med.index)
    b["group"] = b["group_id"].astype(str)
    b = b.drop_duplicates(["kinase", "group"]).set_index(["kinase", "group"])
    if "bin" in b.columns:
        return b["bin"].astype(str)
    return pd.Series(psp_bin(b[ncol].values), index=b.index)


def check_unscored_encoding(scores: pd.DataFrame, policy: str = "raise",
                            frac: float = 0.5) -> list:
    """Unscored must be NaN (r3 S2). Flags conditions with no NaN and >frac exact 0.0
    (likely NaN->0 bug) and kinases with a constant non-NaN score across ≥2 conditions."""
    if policy == "ignore" or scores.shape[0] == 0:
        return []
    v = scores.to_numpy(float)
    nan_c = np.isnan(v).any(1)
    zero_c = (v == 0).mean(1) > frac
    bad = list(scores.index[~nan_c & zero_c])
    const_k = []
    if v.shape[0] >= 2:
        with np.errstate(all="ignore"):
            n_ok = (~np.isnan(v)).sum(0)
            const = (np.nanmax(v, 0) == np.nanmin(v, 0)) & (n_ok >= 2)
        const_k = list(scores.columns[const])
    if bad or const_k:
        msg = (f"possible NaN->0 encoding: {len(bad)} conditions with no NaN and >{frac:.0%} "
               f"exact zeros {bad[:5]}; {len(const_k)} kinases constant across conditions "
               f"{const_k[:5]}. Unscored must be NaN (use zero_policy='ignore' for controls).")
        if policy == "raise" and bad:
            raise ValueError(msg)
        warnings.warn(msg)
    return bad


# --------------------------------------------------------------------------------------
# sealed TEST labels (B10, r3 S5)
# --------------------------------------------------------------------------------------
SEALED_DIR = DEFAULT_ROOT / "cache/sealed"


def sealed_gt_hash(sealed_dir: str | os.PathLike = SEALED_DIR) -> str:
    """sha256 over sorted '<filename>\t<sha256(file)>' lines of all *.tsv in sealed_dir."""
    import hashlib
    lines = []
    for f in sorted(Path(sealed_dir).glob("*.tsv")):
        lines.append(f"{f.name}\t{hashlib.sha256(f.read_bytes()).hexdigest()}")
    return hashlib.sha256("\n".join(lines).encode()).hexdigest()


def load_sealed_test(sealed_dir: str | os.PathLike = SEALED_DIR,
                     freeze: str | os.PathLike | None = None):
    """Only path to TEST labels. Requires FREEZE.json with gt_hash == sealed_gt_hash()."""
    fp = freeze_path(freeze)
    if not fp.exists():
        raise TestSplitLocked(f"{fp} missing; TEST labels are sealed")
    meta = json.loads(fp.read_text())
    want = meta.get("gt_hash") or meta.get("gt_sha256")
    got = sealed_gt_hash(sealed_dir)
    if want != got:
        raise TestSplitLocked(f"sealed GT hash mismatch: FREEZE {want} vs files {got}")
    sd = Path(sealed_dir)
    rd = lambda n: read_gt_tsv(sd / n) if (sd / n).exists() else None
    return {"gt_pairs": rd("gt_pairs.tsv"), "gt_conditions": rd("gt_conditions.tsv"),
            "harmony_bins": rd("harmony_bins.tsv"),
            "harmony_bins_group": rd("harmony_bins_group.tsv"),
            "condition_groups": rd("condition_groups.tsv")}


# --------------------------------------------------------------------------------------
# multiple testing / calibration
# --------------------------------------------------------------------------------------
def bh(p) -> np.ndarray:
    p = np.asarray(p, float)
    n = p.size
    if n == 0:
        return p
    o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1.0)
    return out


def calibration(scores: pd.DataFrame, null_conditions: Sequence[str],
                n_substrates: pd.DataFrame | None = None, thr: float = 1.96) -> dict:
    """Fraction |z|>thr on null conditions: overall, per kinase (summary), per substrate bin.

    n_substrates: conditions×kinases substrate counts (method's own) used for binning.
    Also returns KS distance of pooled two-sided p-values from Uniform(0,1).
    """
    Z = scores.reindex(index=[c for c in null_conditions if c in scores.index])
    z = Z.to_numpy(float)
    ok = ~np.isnan(z)
    exc = np.abs(z) > thr
    res = {"n_conditions": int(Z.shape[0]), "n_scored": int(ok.sum()),
           "frac_exceed": float(exc[ok].mean()) if ok.any() else np.nan}
    per_k = pd.Series(np.where(ok, exc, np.nan).T.tolist(), index=Z.columns).map(
        lambda r: np.nanmean(r) if np.any(~np.isnan(r)) else np.nan)
    res["per_kinase"] = per_k.dropna()
    res["per_kinase_median"] = float(per_k.median()) if per_k.notna().any() else np.nan
    if ok.any():
        pv = 2 * stats.norm.sf(np.abs(z[ok]))
        res["ks_uniform"] = float(stats.kstest(pv, "uniform").statistic)
    if n_substrates is not None:
        ns = n_substrates.reindex(index=Z.index, columns=Z.columns).to_numpy(float)
        lab = psp_bin(ns.ravel()).reshape(ns.shape)
        res["per_bin"] = {b: float(exc[ok & (lab == b)].mean())
                          for b in BIN_LABELS + ("NA",) if (ok & (lab == b)).any()}
    return res


def to_jsonable(metrics: dict) -> dict:
    """Strip DataFrames for ledger/JSON output."""
    out = {}
    for k, v in metrics.items():
        if isinstance(v, (pd.DataFrame, pd.Series)):
            continue
        if isinstance(v, dict):
            out[k] = to_jsonable(v)
        elif isinstance(v, (np.floating, np.integer)):
            out[k] = v.item()
        else:
            out[k] = v
    return out


def dumps(metrics: dict) -> str:
    return json.dumps(to_jsonable(metrics), sort_keys=True, default=str)
