#!/usr/bin/env python
"""B1 library-variant grid on pooled DEV folds 0-4 (no devconfirm/test).

Grid fixed a priori (before any result): see build_grid(). Every config is appended to results/dev_ledger.tsv
(branch "B1"). Outputs: results/b1/grid.tsv (one row per config), cache/b1/grid_pairs.npz (per-pair AUROC
vectors for the nested-CV selection in exp_b1_report.py).
Usage: python scripts/exp_b1_grid.py [--no-ledger]
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from kinact import ledger  # noqa: E402
from kinact.eval import to_jsonable  # noqa: E402
from kinact.explore.lib_variants import B1_CACHE, FOLDS, DevData, Engine, cfg_name, gt_freeze_hash  # noqa: E402

OUT = ROOT / "results/b1"
PREDS = ["r2np", "r3", "cantley"]


def build_grid():
    g = []
    for name in ["psp_all", "psp_invivo", "curated_union", "invitro_sugiyama"]:
        for ms in [1, 3]:
            for st in ["wz", "mean", "ridge"]:
                g.append({"stat": st, "lib": {"type": "static", "name": name}, "min_sub": ms})
    for p in PREDS:
        for k in [25, 50, 100, 200]:
            for mode in ["measured", "global"]:
                for gam in [0, 1, 4]:
                    for st in ["wz", "mean"]:
                        g.append({"stat": st, "lib": {"type": "pred", "name": p, "k": k, "mode": mode, "gamma": gam},
                                  "min_sub": 1})
                    if mode == "measured" and k in (50, 100) and gam in (0, 1):
                        g.append({"stat": "ridge", "lib": {"type": "pred", "name": p, "k": k, "mode": mode,
                                                           "gamma": gam}, "min_sub": 1})
    for cur in ["psp_all", "curated_union"]:
        for p in PREDS:
            for k in [50, 100]:
                for a in [0.1, 0.3, 1.0]:
                    for m in [3, 10, None]:
                        lib = {"type": "mix", "cur": cur, "pred": {"name": p, "k": k, "gamma": 1, "mode": "measured"},
                               "alpha": a, "m": m}
                        for st in ["wz", "mean"]:
                            g.append({"stat": st, "lib": lib, "min_sub": 1})
                        if k == 50 and a == 0.3 and m in (3, None):
                            g.append({"stat": "ridge", "lib": lib, "min_sub": 1})
    return g


def summarise(r):
    v = lambda x: x["value"] if isinstance(x, dict) else x
    hb = r["harmony"]["per_bin"]
    d = {"M1-sel": v(r["M1-sel"]), "M1-all": v(r["M1-all"]), "M1-sel-scored": v(r["M1-sel-scored"]),
         "M1-all-scored": v(r["M1-all-scored"]), "cov_sel": r["coverage_M1-sel"], "cov_all": r["coverage_M1-all"],
         "harmony": v(r["harmony"])}
    for b in ["0", "1-4", "5-19", ">=20"]:
        d[f"h_{b}"] = hb.get(b, np.nan)
    for s in ["H", "D", "O"]:
        d[f"src_{s}"] = v(r["by_source"].get(s, np.nan))
    for s in [-1, 1]:
        d[f"sign_{s}"] = v(r["by_sign"].get(s, np.nan))
    d["n_sig"] = r["coverage_significant"]["n_significant"]
    return d


def fold_metrics(dd, P):
    fold = P["conditionID"].map(dict(zip(dd.conds, dd.fold)))
    out = {}
    for f in FOLDS:
        out[f"tr_{f}"] = dd.ev.metric(P[fold != f], "M1-sel")["value"]
        out[f"te_{f}"] = dd.ev.metric(P[fold == f], "M1-sel")["value"]
    return out


def main():
    use_ledger = "--no-ledger" not in sys.argv
    t0 = time.time()
    gth = gt_freeze_hash()
    assert gth.startswith("1e6f791c"), gth
    dd = DevData()
    eng = Engine(dd)
    grid = build_grid()
    print(f"{len(grid)} configs; GT_FREEZE {gth[:16]}", flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    rows, A, AS, TS, NN = [], [], [], [], []
    for i, cfg in enumerate(grid):
        t1 = time.time()
        S = eng.score(cfg)
        r = dd.ev.evaluate(S)
        P = r["pairs"]
        d = {"i": i, "name": cfg_name(cfg), "stat": cfg["stat"], "lib_type": cfg["lib"]["type"],
             **summarise(r), **fold_metrics(dd, P), "sec": time.time() - t1, "cfg": json.dumps(cfg)}
        rows.append(d)
        A.append(P["auroc"].to_numpy(np.float32)); AS.append(P["auroc_scored"].to_numpy(np.float32))
        TS.append(P["target_scored"].to_numpy()); NN.append(P["n_neg_scored"].to_numpy(np.int32))
        if use_ledger:
            lm = {k: d[k] for k in ["M1-sel", "M1-all", "M1-sel-scored", "M1-all-scored", "cov_sel", "harmony",
                                    "h_0", "h_1-4", "h_5-19", "h_>=20", "src_H", "src_D", "src_O", "sign_-1",
                                    "sign_1", "n_sig"]}
            lm["b1_config_counter"] = i + 1
            lm["gt_freeze"] = gth[:16]
            ledger.append(f"b1_{cfg['stat']}", cfg, "dev", to_jsonable(lm), branch="B1", folds=range(5))
        print(f"[{i+1}/{len(grid)}] {d['name']}: M1-sel {d['M1-sel']:.3f} cov {d['cov_sel']:.2f} "
              f"scored {d['M1-sel-scored']:.3f} harm {d['harmony']:.3f} ({d['sec']:.1f}s)", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "grid.tsv", sep="\t", index=False, float_format="%.4f")
    np.savez_compressed(B1_CACHE / "grid_pairs.npz", auroc=np.vstack(A), auroc_scored=np.vstack(AS),
                        target_scored=np.vstack(TS), n_neg_scored=np.vstack(NN))
    (OUT / "grid_meta.json").write_text(json.dumps({"gt_freeze_sha256": gth, "n_configs": len(grid),
                                                    "runtime_s": time.time() - t0, "folds": FOLDS}, indent=1))
    print(f"done {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
