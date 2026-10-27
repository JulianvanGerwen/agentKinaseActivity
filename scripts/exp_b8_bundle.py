#!/usr/bin/env python
"""B8 final bundle table + runtime. Bundle = core-v2 + adopted experiments (results/b8/b8_eval.json `adopt` flags;
only E1 adopted -> bundle == E1 OOF frame, no new evaluation/choice). Writes results/b8/b8_bundle.json.
Runtime: CPU single-thread, 200 fixed DEV conditions (first 200 GT conds), core-v1 (KinaseActivityModel.score)
and the late-fusion step (z_sig projection + calibrated fusion) timed separately, median of 5 after one warm-up."""
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

from exp_b8_eval import B8, OUT, CACHE, FOLDS
from kinact.eval import to_jsonable
from kinact.explore import signatures_late as SL, final_box as FB
from kinact import bootstrap as bs


def med_time(fn, n=5):
    fn()
    ts = []
    for _ in range(n):
        t = time.perf_counter(); fn(); ts.append(time.perf_counter() - t)
    return float(np.median(ts))


def main():
    ev = json.loads((OUT / "b8_eval.json").read_text())["results"]
    adopted = [e for e in ("E1", "E2", "E3") if ev[e]["adopt"]]
    assert adopted == ["E1"], adopted            # combination rule only needed for >1 adopted
    b = B8()
    S_v2 = pd.read_pickle(CACHE / "S_core-v2.pkl"); S_b = pd.read_pickle(CACHE / "S_E1.pkl")
    c = b.c
    rep = {"adopted": adopted, "core-v2": b.full(S_v2), "bundle": b.full(S_b)}
    rep["bundle"]["delta_vs_core_v2"] = ev["E1"]["delta_vs_core_v2"]
    rep["calibration"] = ev["E1"]["calibration"]
    # signature-kinase pairs
    has = lambda P: np.array([bool(b.sigmask[c.fold_of[cid]][c.kin.index(k)]) for k, cid in zip(P["kinase"], P["conditionID"])])
    for lab, S in (("core-v2", S_v2), ("bundle", S_b)):
        P = c.ev.pairs(S); m = has(P)
        rep[lab]["has-signature"] = c.ev.metric(P[m], "M1-sel")["value"]
        rep[lab]["no-signature"] = c.ev.metric(P[~m], "M1-sel")["value"]
    # ---------------- runtime (CPU, 1 thread)
    from kinact.core import KinaseActivityModel
    model = KinaseActivityModel.from_config()
    conds = c.conds[:200]
    X = c.X[conds]
    f = "fold0"; a, nt = b.alpha[f], b.picks[f][0]
    Ws = b.fits[(f, f)].matrix(b.nsites, c.kin, nt)
    Y, M = b.ctx._YM("gt", np.arange(200))
    zc = np.nan_to_num(np.zeros((200, len(c.kin))))
    v = np.ones(len(c.kin)); vs = np.ones(len(c.kin))
    rt = {}
    with threadpool_limits(1):
        rt["core_v1_score_s"] = med_time(lambda: model.score(X, device="cpu", n_jobs=1))
        rt["sig_project_s"] = med_time(lambda: SL.sig_project(Ws, Y, M))
        rt["fuse_calibrated_s"] = med_time(lambda: FB.fuse_calibrated(zc, SL.sig_project(Ws, Y, M), a, v, vs))
    rt = {k: v_ for k, v_ in rt.items()}
    rt["ms_per_cond"] = {"core_v1": 1e3 * rt["core_v1_score_s"] / 200, "late_fusion_E1": 1e3 * rt["fuse_calibrated_s"] / 200}
    rt["ms_per_cond"]["bundle_total"] = rt["ms_per_cond"]["core_v1"] + rt["ms_per_cond"]["late_fusion_E1"]
    rt["compendium_1937_s_extrapolated"] = rt["ms_per_cond"]["bundle_total"] * 1937 / 1e3
    rep["runtime_cpu_1thread"] = rt
    print(json.dumps(to_jsonable(rt), indent=1), flush=True)
    (OUT / "b8_bundle.json").write_text(json.dumps(to_jsonable(rep), indent=1, default=str))


if __name__ == "__main__":
    main()
