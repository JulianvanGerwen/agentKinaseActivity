#!/usr/bin/env python
"""Branch B4 precompute: fold-safe co-regulation statistics (label-free).

For each DEV fold f, ALLOWED = DEV conditions with qc ok, not in fold f, not sharing (PMID, Control, Sample) with a
fold-f condition (kinact.explore.coreg.allowed_mask). On the allowed columns only:
  (1) kinase seed profiles (weighted z of seed substrates; seeds = curated_union, or curated_union + r2np@100)
      and the pairwise-complete Pearson r between each profile and every site with >= NMIN allowed observations;
      the top-MMAX sites per kinase (excluding seed edges) are stored -> cache/b4/top_<seed>_<fold>.npz
  (2) site-site top-K neighbours of every base-library site (curated_union + r2np@100) -> cache/b4/nbr_<fold>.npz
Also writes cache/b4/allowed.json (allowed condition lists per fold) and timing.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from kinact.data import load_dev, load_dev_meta  # noqa: E402
from kinact.eval import read_gt_tsv  # noqa: E402
from kinact.explore import coreg as CR  # noqa: E402
from kinact.library import load_library  # noqa: E402

FOLDS = [f"fold{i}" for i in range(5)]
CACHE = ROOT / "cache/b4"
NMIN, N0, MMAX, KNB = 30, 30.0, 200, 20
QS = [0, 1, 3]


def kin_site_matrix(lib, kinases, sites):
    ki = pd.Index(kinases).get_indexer(lib["kinase"]); si = pd.Index(sites).get_indexer(lib["site"])
    keep = (ki >= 0) & (si >= 0)
    e = pd.DataFrame({"k": ki[keep], "s": si[keep], "w": lib["weight"].to_numpy(float)[keep]})
    e = e.groupby(["k", "s"], as_index=False)["w"].max()
    return sp.csr_matrix((e["w"], (e["k"], e["s"])), shape=(len(kinases), len(sites)))


def seed_libs():
    cur = load_library("curated_union")[["kinase", "site", "weight"]]
    r2 = load_library("r2np", k=100)[["kinase", "site"]].assign(weight=1.0)
    return {"cur": cur, "cur+r2np": pd.concat([cur, r2]).groupby(["kinase", "site"], as_index=False)["weight"].max()}


def main():
    t0 = time.time()
    CACHE.mkdir(parents=True, exist_ok=True)
    meta = load_dev_meta()
    X = load_dev()
    X = X[X.notna().any(axis=1)]
    sites = X.index
    Udf = read_gt_tsv(ROOT / "results/ground_truth/kinase_universe.tsv")
    kin = list(Udf["kinase"])
    Y = CR.standardise(X.to_numpy(np.float64))
    libs = seed_libs()
    Wseed = {k: kin_site_matrix(v, kin, sites) for k, v in libs.items()}
    base_rows = np.unique(Wseed["cur+r2np"].tocoo().col)
    print(f"load {time.time()-t0:.0f}s; sites {len(sites)}; base lib rows {len(base_rows)}", flush=True)
    allowed_info, timing = {}, {}
    for f in FOLDS:
        tf = time.time()
        held = CR.fold_heldout(meta, f)
        al = CR.allowed_mask(meta, X.columns, held)
        allowed_info[f] = list(X.columns[al])
        Ya = Y[:, al]
        nobs = (~np.isnan(Ya)).sum(1)
        cand = np.flatnonzero(nobs >= NMIN)
        Pall = {sname: CR.seed_profiles(Ya, W) for sname, W in Wseed.items()}
        t1 = 0.0
        for q in QS:
            F = CR.global_factors(Pall["cur+r2np"], q) if q else None
            Yc = CR.residualize(Ya[cand], F) if q else Ya[cand]
            for sname, W in Wseed.items():
                P = CR.residualize(Pall[sname], F) if q else Pall[sname]
                R, N = CR.profile_site_corr(P, Yc)
                top = CR.topm_per_kinase(R, N, MMAX, NMIN, N0, exclude=W[:, cand])
                idx = np.where(top["idx"] >= 0, cand[np.maximum(top["idx"], 0)], -1)
                np.savez_compressed(CACHE / f"top_{sname}_q{q}_{f}.npz", idx=idx, r=top["r"], n=top["n"],
                                    rs=top["rs"], kinases=np.array(kin), prof_cov=(~np.isnan(P)).sum(1))
                del R, N
            ts = time.time()
            Yfull = Ya.copy(); Yfull[cand] = Yc
            nbr, rn, nn = CR.site_site_topk(Yfull, base_rows, k=KNB, nmin=NMIN, n0=N0, cand=cand)
            np.savez_compressed(CACHE / f"nbr_q{q}_{f}.npz", rows=base_rows, nbr=nbr, r=rn, n=nn)
            t1 += time.time() - ts
        timing[f] = {"total_s": time.time() - tf, "site_site_s": t1, "n_allowed": int(al.sum()),
                     "n_heldout": len(held), "n_cand_sites": int(len(cand))}
        print(f, timing[f], flush=True)
    np.save(CACHE / "sites.npy", np.array(sites))
    json.dump({"allowed": allowed_info, "timing": timing, "NMIN": NMIN, "N0": N0, "MMAX": MMAX, "KNB": KNB},
              open(CACHE / "allowed.json", "w"))
    print(f"done {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
