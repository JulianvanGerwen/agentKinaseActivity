"""Runtime benchmark (DESIGN §1 desideratum 4; reported, not a gate).

A method is either a callable score(X, lib) -> conditions×kinases, or an object with
prepare(X, lib) -> state (load/mapping: site matching, matrix building) and
score(state) -> frame (scoring). Optionally null(state) is timed separately.
Each phase: one warm-up run, then median of n_rep runs, under threadpool limits 1 and 16.
"""
from __future__ import annotations

import os
import time
from typing import Callable

import numpy as np

N_COMPENDIUM = 1937


def _time(fn, n_rep):
    fn()  # warm
    ts = []
    for _ in range(n_rep):
        t0 = time.perf_counter()
        fn()
        ts.append(time.perf_counter() - t0)
    return float(np.median(ts)), ts


def bench_runtime(method, X_subset, lib=None, n_rep: int = 5, threads=(1, 16)) -> dict:
    from threadpoolctl import threadpool_limits

    n_cond = X_subset.shape[1]
    has_phases = hasattr(method, "prepare") and hasattr(method, "score")
    out = {"n_conditions": int(n_cond), "n_rep": n_rep, "phases": has_phases}
    for nt in threads:
        old = {k: os.environ.get(k) for k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS",
                                              "OPENBLAS_NUM_THREADS")}
        for k in old:
            os.environ[k] = str(nt)  # affects only libraries initialised later; limits below are live
        try:
            with threadpool_limits(limits=nt):
                r = {}
                if has_phases:
                    r["load_mapping_s"], _ = _time(lambda: method.prepare(X_subset, lib), n_rep)
                    state = method.prepare(X_subset, lib)
                    r["scoring_s"], _ = _time(lambda: method.score(state), n_rep)
                    if hasattr(method, "null"):
                        r["null_s"], _ = _time(lambda: method.null(state), n_rep)
                    total = r["load_mapping_s"] + r["scoring_s"] + r.get("null_s", 0.0)
                else:
                    total, _ = _time(lambda: method(X_subset, lib), n_rep)
                r["total_s"] = total
                r["ms_per_condition"] = 1e3 * total / max(n_cond, 1)
                r["compendium_extrapolated_s"] = total / max(n_cond, 1) * N_COMPENDIUM
                out[f"threads_{nt}"] = r
        finally:
            for k, v in old.items():
                if v is None:
                    os.environ.pop(k, None)
                else:
                    os.environ[k] = v
    return out
