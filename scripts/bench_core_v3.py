#!/usr/bin/env python
"""Equivalence summary + runtime of the user-facing core-v3 scorer (config/core_v3.json) — DEV data only.

Equivalence (numbers for the docs; the asserting test is tests/test_core_v3.py): deployed core-v3 z vs the research
path (bayes2 batched z_vif + signatures_late.sig_project / rho_sig + final_box.fuse_calibrated with the artefact
signatures and E1 scales) on the 200 fixed bench conditions (seed 0) and all 1,358 DEV conditions; plus the
fused-vs-core-v1 increment.
Runtime: kinact.runtime.bench_runtime on the same 200 conditions (load/mapping = prepare, scoring = posterior + fusion;
warm, median of 5) at 1 and 16 threads and on the GPU, for core-v1 and core-v3; full DEV compendium (1,358 conditions)
wall time (prepare + score, median of 3 warm runs). Writes results/core_v3/{equivalence,runtime}.json.
"""
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))
from bench_core import GPU, compare, wall  # noqa: E402
from kinact.core import KinaseActivityModel  # noqa: E402
from kinact.data import load_dev  # noqa: E402
from kinact.eval import read_gt_tsv  # noqa: E402
from kinact.explore import bayes2 as B2, final_box as FB, signatures_late as SL  # noqa: E402
from kinact.library import load_library  # noqa: E402
from kinact.runtime import bench_runtime  # noqa: E402

OUT = ROOT / "results/core_v3"
CFG3 = ROOT / "config/core_v3.json"


def research_z(m, X32, md, U, bcfg, kin):
    """Research path on artefact signatures: bayes2 z_vif fused with SL.sig_project (float32 input as in B7b)."""
    art = m.artifacts_dir
    sig = pd.read_csv(art / "signatures.tsv.gz", sep="\t")
    sc = pd.read_csv(art / "e1_scales.tsv", sep="\t").set_index("kinase").reindex(kin)
    import scipy.sparse as sp
    si = pd.Index(X32.index).get_indexer(sig.site); ki = pd.Index(kin).get_indexer(sig.kinase)
    Ws = sp.csr_matrix((sig.weight.to_numpy(), (si, ki)), shape=(len(X32.index), len(kin)))
    Y, M, _ = B2.standardise(X32.to_numpy(np.float64), "wsd", None)
    zs = SL.sig_project(Ws, Y.astype(np.float32), M.astype(np.float32)) / m.sig["rho"]
    zc = B2.score(X32, md, U, bcfg).reindex(columns=kin).to_numpy()
    return FB.fuse_calibrated(zc, zs, m.sig["alpha"], sc.v_both.to_numpy(), sc.v_sig.to_numpy())


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    X32 = load_dev()
    X = X32.astype(np.float64)
    rng = np.random.default_rng(0)
    conds = sorted(rng.choice(sorted(X.columns), 200, replace=False))
    Xb = X[conds]
    m3 = KinaseActivityModel.from_config(CFG3)
    m1 = KinaseActivityModel.from_config()
    U = read_gt_tsv(ROOT / "results/ground_truth/kinase_universe.tsv"); kin = list(U["kinase"])
    md = B2.build_multi_design(X.index, kin, load_library("curated_union"), {"r2np": load_library("r2np", k=100)})
    bcfg = B2.Bayes2Config(tau=m3.cfg["tau"], rho_f=m3.cfg["rho_f"], rho_g=m3.cfg["rho_g"], alphas={"r2np": 1.0},
                           output="z_vif", rho_e=m3.cfg["rho_e"], scale="wsd", winsor=None, borrow=True)
    # ------------------------------------------------------------------ equivalence
    eq = {}
    for lab, cs in (("bench200", conds), ("dev1358", list(X.columns))):
        R = m3.score(X[cs])
        ref = pd.DataFrame(research_z(m3, X32[cs], md, U, bcfg, kin), index=cs, columns=kin)
        eq[f"{lab}_cpu_vs_research_path"] = compare(R.z, ref)
        Rg = m3.score(X[cs], device="cuda")
        eq[f"{lab}_gpu_vs_cpu"] = compare(Rg.z, R.z)
        has = m3.sig["has"]
        d = (R.z.to_numpy() - R.z_core.to_numpy())[:, has]
        eq[f"{lab}_signature_kinases"] = {"n": int(has.sum()), "cells_changed_vs_core_v1": int((np.abs(d) > 1e-12).sum()),
                                          "z_sig_defined_frac": float(np.isfinite(R.z_sig.to_numpy()[:, has]).mean())}
        eq[f"{lab}_non_signature_kinases_equal_core_v1"] = bool(np.array_equal(
            np.nan_to_num(R.z.to_numpy()[:, ~has], nan=9e9), np.nan_to_num(m1.score(X[cs]).z.to_numpy()[:, ~has], nan=9e9)))
    (OUT / "equivalence.json").write_text(json.dumps(eq, indent=1))
    print(json.dumps(eq, indent=1), flush=True)

    # ------------------------------------------------------------------ runtime
    rt = {"bench": "kinact.runtime.bench_runtime, 200 fixed DEV conditions (seed 0), warm, median of 5",
          "core_v3": {"cpu": bench_runtime(m3, Xb, threads=(1, 16)), "gpu": bench_runtime(GPU(m3), Xb, threads=(16,))},
          "core_v1": {"cpu": bench_runtime(m1, Xb, threads=(1, 16)), "gpu": bench_runtime(GPU(m1), Xb, threads=(16,))}}
    st = m3.prepare(Xb)
    with threadpool_limits(1):
        zc = m3.score(st).z_core.to_numpy()
        rt["core_v3"]["fusion_only_200_1t_s"] = wall(lambda: m3.fuse(zc, m3.sig_z(st)), n=5)
    full = {}
    for nt in (1, 16):
        with threadpool_limits(nt):
            full[f"cpu_{nt}t_prepare_s"] = wall(lambda: m3.prepare(X))
            st = m3.prepare(X)
            full[f"cpu_{nt}t_score_s"] = wall(lambda: m3.score(st))
            full[f"cpu_{nt}t_total_s"] = full[f"cpu_{nt}t_prepare_s"] + full[f"cpu_{nt}t_score_s"]
    with threadpool_limits(16):
        full["gpu_score_s"] = wall(lambda: m3.score(st, device="cuda"))
        full["gpu_total_s"] = full["cpu_16t_prepare_s"] + full["gpu_score_s"]
    full["n_conditions"] = int(X.shape[1])
    full["ms_per_condition_1t"] = 1e3 * full["cpu_1t_total_s"] / X.shape[1]
    full["compendium_1937_extrapolated_1t_s"] = full["cpu_1t_total_s"] / X.shape[1] * 1937
    rt["core_v3"]["dev_compendium_1358"] = full
    import torch
    rt["hardware"] = {"cpu_count": os.cpu_count(), "gpu": torch.cuda.get_device_name(0),
                      "loadavg_at_end": list(os.getloadavg())}
    (OUT / "runtime.json").write_text(json.dumps(rt, indent=1))
    print(json.dumps(rt, indent=1))


if __name__ == "__main__":
    main()
