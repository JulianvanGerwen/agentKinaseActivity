"""Data loaders: DEV compendium cache and kinase-substrate libraries (docs/interfaces.md).

DEV cache (cache/dev/): quant_dev.npy float32 [peptides x conditions], rows = sites_dev.txt
(peptide_UID), columns = conditions_dev.txt. Peptides map to canonical site keys
(`{UniProt}_{residue}{pos}`, cache/sites.parquet); several peptides (multiplicities) may share a site.
"""
from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .eval import DEFAULT_ROOT, read_gt_tsv

DEV_DIR = DEFAULT_ROOT / "cache/dev"
LIB_DIR = DEFAULT_ROOT / "cache/libraries"
SITES_PQ = DEFAULT_ROOT / "cache/sites.parquet"


def load_dev_meta(dev_dir=DEV_DIR) -> pd.DataFrame:
    """conditions_dev.tsv indexed by conditionID (columns incl. group_id, split, qc_status)."""
    return read_gt_tsv(Path(dev_dir) / "conditions_dev.tsv").set_index("conditionID", drop=False)


@lru_cache(maxsize=4)
def _load_dev_raw(dev_dir: str):
    d = Path(dev_dir)
    q = np.load(d / "quant_dev.npy", mmap_mode="r")
    peps = (d / "sites_dev.txt").read_text().split()
    conds = (d / "conditions_dev.txt").read_text().split()
    assert q.shape == (len(peps), len(conds)), (q.shape, len(peps), len(conds))
    pep2site = pd.read_parquet(SITES_PQ, columns=["peptide_UID", "site_key"]).set_index("peptide_UID")["site_key"]
    sites = pep2site.reindex(peps)
    assert sites.notna().all()
    return q, peps, conds, sites.to_numpy()


def load_dev(dev_dir=DEV_DIR, conditions=None, collapse=True, dtype=np.float32) -> pd.DataFrame:
    """Sites x conditions DataFrame (log2FC, NaN = not quantified).

    collapse=True: peptides sharing a canonical site key are averaged (nanmean) -> index = site_key
    (unique, first-occurrence order). collapse=False: index = peptide_UID. `conditions` restricts
    columns (loaded lazily from the memmap, so subsets are cheap).
    """
    q, peps, conds, sites = _load_dev_raw(str(dev_dir))
    if conditions is not None:
        cidx = pd.Index(conds).get_indexer(list(conditions))
        if (cidx < 0).any():
            raise KeyError("conditions not in DEV cache: " + str(list(np.asarray(conditions)[cidx < 0])[:5]))
        cols = list(conditions)
        M = np.asarray(q[:, cidx], dtype=np.float32)
    else:
        cols, M = conds, np.asarray(q, dtype=np.float32)
    if not collapse:
        return pd.DataFrame(M.astype(dtype, copy=False), index=peps, columns=cols)
    codes, uniq = pd.factorize(sites)
    if len(uniq) == len(sites):
        return pd.DataFrame(M.astype(dtype, copy=False), index=uniq, columns=cols)
    ok = ~np.isnan(M)
    S = sp.csr_matrix((np.ones(len(codes), np.float32), (codes, np.arange(len(codes)))),
                      shape=(len(uniq), len(codes)))
    num = S @ np.where(ok, M, 0.0)
    cnt = S @ ok.astype(np.float32)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = np.where(cnt > 0, num / cnt, np.nan)
    return pd.DataFrame(out.astype(dtype, copy=False), index=uniq, columns=cols)


def load_library_df(name: str, lib_dir=LIB_DIR) -> pd.DataFrame:
    return pd.read_parquet(Path(lib_dir) / f"{name}.parquet")


def library_matrix(lib: str | pd.DataFrame, site_index, kinases=None, binary=False):
    """Sparse kinase x site matrix (CSR, float32) aligned to `site_index` (columns).

    lib: library name or DataFrame(kinase, site, weight). Edges whose site is not in site_index are
    dropped; duplicate (kinase, site) edges keep the max weight. kinases: row order (default sorted).
    Returns (W, kinase_index, site_index).
    """
    df = load_library_df(lib) if isinstance(lib, str) else lib
    site_index = pd.Index(site_index)
    if kinases is None:
        kinases = sorted(df["kinase"].unique())
    kidx = pd.Index(kinases)
    ki = kidx.get_indexer(df["kinase"])
    si = site_index.get_indexer(df["site"])
    keep = (ki >= 0) & (si >= 0)
    w = np.ones(keep.sum(), np.float32) if binary else df["weight"].to_numpy(np.float32)[keep]
    e = pd.DataFrame({"k": ki[keep], "s": si[keep], "w": w}).groupby(["k", "s"], as_index=False)["w"].max()
    W = sp.csr_matrix((e["w"].to_numpy(np.float32), (e["k"], e["s"])), shape=(len(kidx), len(site_index)))
    return W, kidx, site_index
