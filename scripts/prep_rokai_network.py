#!/usr/bin/env python
"""Convert the RoKAI functional network (github.com/serhan-yilmaz/RoKAI, data/rokai_network_data_uniprotkb.mat,
downloaded 2026-10-27 into external/rokai/) into numpy/parquet files usable without MATLAB.

The .mat stores MATLAB `table` objects (Site, Kinase) that scipy cannot read; they are parsed with the
`mat-io` package (pip install mat-io; set MATIO_PATH to its install dir if not on sys.path).
Output: external/rokai/rokai_network.npz (Wkin2site, Wkin2kin, Wsite2site_coev, Wsite2site_sd as CSR),
external/rokai/rokai_sites.tsv (Identifier, Protein, Position), external/rokai/rokai_kinases.tsv (KinaseID, Gene).
"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / "external/rokai"
if os.environ.get("MATIO_PATH"):
    sys.path.insert(0, os.environ["MATIO_PATH"])
import warnings  # noqa: E402

warnings.filterwarnings("ignore")
import matio  # noqa: E402


def _tab(o):
    p = o.properties
    names = [str(np.asarray(v).ravel()[0]) for v in p["varnames"].ravel()]
    cols = {}
    for n, col in zip(names, p["data"].ravel()):
        a = np.asarray(col).ravel()
        if a.dtype == object:
            a = np.array([str(np.asarray(x).ravel()[0]) if np.asarray(x).size else "" for x in a])
        cols[n] = a
    return pd.DataFrame(cols)


def main():
    d = matio.load_from_mat(str(D / "data/rokai_network_data_uniprotkb.mat"))
    nd = d["NetworkData"]
    site = _tab(nd["Site"][0, 0])
    kin = _tab(nd["Kinase"][0, 0])
    mats = {k: sp.csr_matrix(nd[k][0, 0]) for k in ["Wkin2site", "Wkin2kin", "Wsite2site_coev", "Wsite2site_sd"]}
    assert mats["Wkin2site"].shape == (len(kin), len(site)), (mats["Wkin2site"].shape, len(kin), len(site))
    np.savez_compressed(D / "rokai_network.npz", **{k: None for k in []},
                        **{f"{k}_{a}": getattr(m, a) for k, m in mats.items() for a in ["data", "indices", "indptr"]},
                        **{f"{k}_shape": np.array(m.shape) for k, m in mats.items()})
    site[[c for c in ["Identifier", "Protein", "Position"] if c in site]].to_csv(D / "rokai_sites.tsv", sep="\t", index=False)
    kin[[c for c in ["KinaseID", "Gene"] if c in kin]].to_csv(D / "rokai_kinases.tsv", sep="\t", index=False)
    print(len(site), "sites;", len(kin), "kinases;", {k: m.nnz for k, m in mats.items()})


if __name__ == "__main__":
    main()
