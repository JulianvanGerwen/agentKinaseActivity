#!/usr/bin/env python
"""Equivalence + runtime of the user-facing core-v1 scorer (kinact.core) — DEV data only.

Equivalence: kinact.core z vs kinact.explore.bayes2 core-v1 z_vif (same config) on the 200 fixed bench conditions
(seed 0, as scripts/run_baselines.py) — against (a) bayes2.SingleSample2 (float64 CPU path) and (b) bayes2.score
(batched GPU path; float32 sufficient statistics) — and on all 1,358 DEV conditions against (b).
Runtime: kinact.runtime.bench_runtime on the same 200 conditions (load/mapping = prepare, scoring = posterior;
warm, median of 5) at 1 and 16 threads and on the GPU; plus full DEV compendium wall time (prepare + score, median
of 3 warm runs). Writes results/core_v1/{equivalence,runtime}.json.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from kinact.core import KinaseActivityModel  # noqa: E402
from kinact.data import load_dev  # noqa: E402
from kinact.eval import read_gt_tsv  # noqa: E402
from kinact.explore import bayes2 as B2  # noqa: E402
from kinact.library import load_library  # noqa: E402
from kinact.runtime import bench_runtime  # noqa: E402

OUT = ROOT / "results/core_v1"


def compare(A: pd.DataFrame, B: pd.DataFrame) -> dict:
    a, b = A.to_numpy(), B.to_numpy()
    nan_eq = bool(np.array_equal(np.isnan(a), np.isnan(b)))
    ok = ~np.isnan(a) & ~np.isnan(b)
    sp_min, sp_min_r = 1.0, 1.0
    for x, y in zip(a, b):
        m = ~np.isnan(x) & ~np.isnan(y)
        if m.sum() > 2:
            sp_min = min(sp_min, spearmanr(x[m], y[m])[0])
            # exact ties (e.g. symmetric borrowed kinases) are broken by 1e-12 rounding: compare at 1e-8 resolution
            sp_min_r = min(sp_min_r, spearmanr(np.round(x[m], 8), np.round(y[m], 8))[0])
    return {"n_conditions": int(a.shape[0]), "n_scored_cells": int(ok.sum()), "nan_pattern_equal": nan_eq,
            "max_abs_diff": float(np.abs(a[ok] - b[ok]).max()), "pooled_spearman": float(spearmanr(a[ok], b[ok])[0]),
            "min_per_condition_spearman": float(sp_min), "min_per_condition_spearman_round1e-8": float(sp_min_r)}


class GPU:
    def __init__(self, m):
        self.m = m

    def prepare(self, X, lib=None):
        return self.m.prepare(X)

    def score(self, st):
        return self.m.score(st, device="cuda")


def wall(fn, n=3):
    fn()
    ts = []
    for _ in range(n):
        t0 = time.perf_counter(); fn(); ts.append(time.perf_counter() - t0)
    return float(np.median(ts))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    X = load_dev().astype(np.float64)
    rng = np.random.default_rng(0)
    conds = sorted(rng.choice(sorted(X.columns), 200, replace=False))
    Xb = X[conds]
    m = KinaseActivityModel.from_config()

    # ------------------------------------------------------------------ equivalence
    U = read_gt_tsv(ROOT / "results/ground_truth/kinase_universe.tsv")
    md = B2.build_multi_design(X.index, list(U["kinase"]), load_library("curated_union"),
                               {"r2np": load_library("r2np", k=100)})
    cfg = B2.Bayes2Config(tau=m.cfg["tau"], rho_f=m.cfg["rho_f"], rho_g=m.cfg["rho_g"], alphas={"r2np": 1.0},
                          output="z_vif", rho_e=m.cfg["rho_e"], scale="wsd", winsor=None, borrow=True)
    R = m.score(Xb)
    Rg = m.score(Xb, device="cuda")
    ss = B2.SingleSample2(md, U, cfg)
    Zs = pd.DataFrame(np.vstack([ss.score_one(Xb[c].to_numpy())[0] for c in conds]), index=conds,
                      columns=list(U["kinase"])).reindex(columns=m.kinases)
    Sb = B2.score(Xb, md, U, cfg).reindex(columns=m.kinases)
    Ra = m.score(X)
    Sa = B2.score(X, md, U, cfg).reindex(columns=m.kinases)
    eq = {"bench200_cpu_vs_bayes2_SingleSample2_float64": compare(R.z, Zs),
          "bench200_gpu_vs_bayes2_SingleSample2_float64": compare(Rg.z, Zs),
          "bench200_cpu_vs_bayes2_score_batched": compare(R.z, Sb),
          "bench200_bayes2_SingleSample2_vs_bayes2_score (bayes2 internal float32 error)": compare(Zs, Sb),
          "dev1358_cpu_vs_bayes2_score_batched": compare(Ra.z, Sa),
          "bench200_mean_vs_bayes2_posterior_mean_maxabs": None}
    _, aux = B2.score(Xb, md, U, cfg, return_all=True)
    mu_ref = np.where(aux["report"], aux["post"]["mu"], np.nan)
    sd_ref = np.where(aux["report"], aux["post"]["sd_post"], np.nan)
    eq["bench200_mean_vs_bayes2_posterior_mean_maxabs"] = float(np.nanmax(np.abs(R.mean.to_numpy() - mu_ref)))
    eq["bench200_sd_vs_bayes2_posterior_sd_maxabs"] = float(np.nanmax(np.abs(R.sd.to_numpy() - sd_ref)))
    eq["bench200_n_substrates_equal"] = bool(np.array_equal(R.n_substrates.to_numpy(), aux["n"].astype(int)))
    (OUT / "equivalence.json").write_text(json.dumps(eq, indent=1))
    print(json.dumps(eq, indent=1), flush=True)

    # ------------------------------------------------------------------ runtime
    rt = {"bench": "kinact.runtime.bench_runtime, 200 fixed DEV conditions (seed 0), warm, median of 5",
          "cpu": bench_runtime(m, Xb, threads=(1, 16)),
          "gpu": bench_runtime(GPU(m), Xb, threads=(16,))}
    full = {}
    for nt in (1, 16):
        with threadpool_limits(nt):
            full[f"cpu_{nt}t_prepare_s"] = wall(lambda: m.prepare(X))
            st = m.prepare(X)
            full[f"cpu_{nt}t_score_s"] = wall(lambda: m.score(st))
            full[f"cpu_{nt}t_total_s"] = full[f"cpu_{nt}t_prepare_s"] + full[f"cpu_{nt}t_score_s"]
    with threadpool_limits(16):
        full["gpu_score_s"] = wall(lambda: m.score(st, device="cuda"))
        full["gpu_total_s"] = full["cpu_16t_prepare_s"] + full["gpu_score_s"]
    full["n_conditions"] = int(X.shape[1])
    full["n_patterns"] = int(m.score(st).info["n_patterns"])
    with threadpool_limits(1):
        full["n_woodbury_1t"] = int(m.score(st).info["n_woodbury"])
    rt["dev_compendium_1358"] = full
    import torch
    rt["hardware"] = {"cpu_count": __import__("os").cpu_count(), "gpu": torch.cuda.get_device_name(0),
                      "loadavg_at_end": list(__import__("os").getloadavg())}
    (OUT / "runtime.json").write_text(json.dumps(rt, indent=1))
    print(json.dumps(rt, indent=1))


if __name__ == "__main__":
    main()
