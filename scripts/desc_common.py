"""Shared setup for the DESCRIPTIVE analyses (DESIGN §5 'Descriptive'; never used for model selection).

DEV only: cache/dev compendium + results/ground_truth DEV rows (folds 0-4). No cache/sealed, no external/.

Score matrices (conditions x kinases, NaN = unscored) for ALL DEV conditions, cached in cache/descriptive/:
  core-v1    : kinact.core package (config/core_v1.json), z_vif
  core-v2    : core-v1 + B7b late-fusion signatures, alpha = 1 (D11). Condition in fold f uses the B7 PRIMARY
               ('clean+contrast') signature fit on the 4 folds != f (cache/b7/sig_clean+contrast_fold{f}.pkl) with
               B7's nested (n_top, kappa) pick for fold f; rho_sig = RMS of z_sig over the 7 DEV null conditions.
  psp-z-min3 : RoKAI/benchmarKIN z-score, psp_all, >= 3 measured substrates
  cur-z-min1 : z-score, curated_union, >= 1 measured substrate
  prior-only : s_c * freq_k, freq from GT conditions of the other 4 folds (GT conditions only; oracle sign)
`masked=True`: every self-site (selfsites.self_site_keys, psp+aloop + all PSP 'activity' regulatory sites on U
kinases) is set to NaN in the input before scoring (GT-A, DESIGN §3).
"""
import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path("/workspace")
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))
from kinact import controls as ctl  # noqa: E402
from kinact.core import KinaseActivityModel  # noqa: E402
from kinact.data import load_dev, load_dev_meta  # noqa: E402
from kinact.eval import Evaluator, load_gt, read_gt_tsv  # noqa: E402
from kinact.explore import bayes2 as B2  # noqa: E402
from kinact.explore import selfsites as SS  # noqa: E402
from kinact.explore import signatures_late as SL  # noqa: E402
from kinact.methods.linear import zscore  # noqa: E402

FOLDS = [f"fold{i}" for i in range(5)]
CACHE = ROOT / "cache/descriptive"; CACHE.mkdir(parents=True, exist_ok=True)
OUT = ROOT / "results/descriptive"; OUT.mkdir(parents=True, exist_ok=True)
METHODS = ["core-v1", "core-v2", "psp-z-min3", "cur-z-min1", "prior-only"]
ALPHA_V2 = 1.0
GT = ROOT / "results/ground_truth"


class Data:
    def __init__(self):
        self.P, self.C, self.U, _ = load_gt()
        self.bins = read_gt_tsv(GT / "harmony_bins_group.tsv")
        self.Udf = read_gt_tsv(GT / "kinase_universe.tsv")
        self.kin = list(self.Udf["kinase"])
        self.meta = load_dev_meta()
        self.X = load_dev()
        assert set(self.meta["split"]) <= set(FOLDS), "DEV cache must contain folds 0-4 only"
        self.ev = Evaluator(self.P, self.C, self.U, split_filter=FOLDS, bins=self.bins, zero_policy="ignore")
        self.gt_conds = list(self.ev.cond_index)
        self.fold_of = self.meta["split"]
        self.nulls = [c for c in self.C.loc[(self.C["role"] == "null") & self.C["split"].isin(FOLDS), "conditionID"]
                      if c in self.X.columns]


def b7_picks():
    d = json.loads((ROOT / "results/b7/b7_eval.json").read_text())
    ch = d["results"]["B7 clean+contrast (nested)"]["chosen"]
    return {f: (int(c["n_top"]), float(c["kappa"])) for f, c in ch.items()}


def zsig_oof(d: Data, X: pd.DataFrame) -> pd.DataFrame:
    """Fold-OOF standardised z_sig for all columns of X (sites index must equal d.X.index)."""
    assert X.index.equals(d.X.index)
    Y, M, _ = B2.standardise(X.to_numpy(np.float64), "wsd", None)
    cols = pd.Index(X.columns)
    nul = cols.get_indexer(d.nulls)
    picks = b7_picks()
    Z = np.full((X.shape[1], len(d.kin)), np.nan)
    fo = d.fold_of.reindex(cols).to_numpy()
    for f in FOLDS:
        n_top, kp = picks[f]
        sig = pickle.load(open(ROOT / f"cache/b7/sig_clean+contrast_{f}.pkl", "rb"))[kp]
        Ws = sig.matrix(X.shape[0], d.kin, n_top)
        rho = SL.rho_sig(SL.sig_project(Ws, Y[:, nul], M[:, nul]))
        idx = np.flatnonzero(fo == f)
        Z[idx] = SL.sig_project(Ws, Y[:, idx], M[:, idx]) / rho
    return pd.DataFrame(Z, index=cols, columns=d.kin)


def prior_frame(d: Data) -> pd.DataFrame:
    parts = []
    gc = pd.Series(d.gt_conds)
    for f in FOLDS:
        ev = list(gc[gc.map(d.fold_of) == f]); tr = list(gc[gc.map(d.fold_of) != f])
        parts.append(ctl.prior_only(d.P, d.C, tr, ev, d.kin))
    return pd.concat(parts).reindex(index=d.gt_conds, columns=d.kin)


def score_all(d: Data, masked: bool = False, methods=METHODS, model=None) -> dict:
    """{method: conditions x kinases DataFrame} for all DEV conditions (prior-only: GT conditions only). Cached."""
    tag = "masked" if masked else "full"
    out, X = {}, None
    for m in methods:
        fp = CACHE / f"scores_{tag}_{m}.parquet"
        if fp.exists():
            out[m] = pd.read_parquet(fp); continue
        if X is None:
            X = d.X.astype(np.float64)
            if masked:
                X = SS.mask_self_sites(X, SS.self_site_keys(d.Udf))
        if m == "core-v1":
            mdl = model or KinaseActivityModel.from_config()
            S = mdl.score(X).z.reindex(columns=d.kin)
        elif m == "core-v2":
            zc = out.get("core-v1")
            if zc is None:
                zc = score_all(d, masked, ["core-v1"])["core-v1"]
            zs = zsig_oof(d, X).reindex(index=zc.index)
            S = pd.DataFrame(SL.fuse(zc.to_numpy(), zs.to_numpy(), ALPHA_V2), index=zc.index, columns=d.kin)
        elif m == "psp-z-min3":
            S = zscore(X, "psp_all", min_substrates=3).reindex(columns=d.kin)
        elif m == "cur-z-min1":
            S = zscore(X, "curated_union", min_substrates=1).reindex(columns=d.kin)
        elif m == "prior-only":
            S = prior_frame(d)
        S = S.astype(np.float64)
        S.columns = S.columns.astype(str)
        S.to_parquet(fp)
        out[m] = S
    return out


def inhibitor_agents() -> set:
    t = pd.read_csv(GT / "drug_targets_manual.tsv", sep="\t", keep_default_na=False, comment="#")
    return set(t["agent"]) | {"MAPKi_unspecified", "genistein", "staurosporine"}


def non_inhibitor_conditions(d: Data) -> list:
    """DEV conditions (qc ok) with no kinase-inhibitor agent (drug_targets_manual agents, incl. lipid-kinase-only)
    and not Hijazi (source H)."""
    g = read_gt_tsv(GT / "condition_groups.tsv").set_index("conditionID")
    g = g.reindex(d.X.columns)
    inh = inhibitor_agents()
    has_inh = g["drugs"].fillna("").map(lambda s: bool(set(filter(None, s.split(";"))) & inh))
    src = d.C.set_index("conditionID")["source"].reindex(d.X.columns).fillna("")
    ok = d.meta["qc_status"].reindex(d.X.columns) == "ok"
    return list(d.X.columns[(~has_inh.to_numpy()) & (src != "H").to_numpy() & ok.to_numpy()])


def qc_ok_conditions(d: Data) -> list:
    return list(d.meta.index[d.meta["qc_status"] == "ok"])
