#!/usr/bin/env python
"""Build the core-v1 model artefacts (models/core_v1/) for kinact.core.KinaseActivityModel.

Inputs (label-free / allowed DEV data only):
  * kinase universe + Manning group/family: results/ground_truth/kinase_universe.tsv (no labels);
  * libraries: cache/libraries/<curated>.parquet and the predicted sources at their top-k (config/core_v1.json);
  * rho_e refit: the DEV `role == null` conditions of folds 0-4 (decryptM 0 vs 0 nM), read from cache/dev only.
    The refit uses the same procedure as bayes2.fit_rho_null (cross-checked) and must reproduce the config value.
Nothing from cache/sealed/, no DEV-confirm or TEST values, no GT labels.

Outputs: models/core_v1/edges.tsv.gz (site, kinase, weight, source; predicted edges that are curated are dropped,
as in bayes2.build_multi_design), kinases.tsv, manifest.json (hashes, counts, rho_e refit table).
"""
import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from kinact import core  # noqa: E402
from kinact.data import load_dev  # noqa: E402
from kinact.eval import load_gt, read_gt_tsv  # noqa: E402
from kinact.library import LIBDIR, load_library, sha256_file  # noqa: E402

FOLDS = [f"fold{i}" for i in range(5)]


def build_edges(cfg, kinases):
    kset = set(kinases)
    cur = load_library(cfg["library"]["curated"])
    cur = cur[cur["kinase"].isin(kset)]
    ec = cur.groupby(["site", "kinase"], as_index=False)["weight"].max().assign(source="cur")
    cur_keys = set(zip(ec["site"], ec["kinase"]))
    parts = [ec]
    for name, spec in cfg["library"].get("predicted", {}).items():
        lib = load_library(name, k=int(spec["k"]))
        lib = lib[lib["kinase"].isin(kset)][["site", "kinase"]].drop_duplicates()
        keep = np.array([(s, k) not in cur_keys for s, k in zip(lib["site"], lib["kinase"])], bool)
        parts.append(lib[keep].assign(weight=1.0, source=name))
    e = pd.concat(parts, ignore_index=True)[["site", "kinase", "weight", "source"]]
    return e.sort_values(["source", "kinase", "site"]).reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=str(ROOT / "config/core_v1.json"))
    ap.add_argument("--check-bayes2", action="store_true", help="cross-check rho_e fit with explore.bayes2")
    a = ap.parse_args()
    cfg_path = Path(a.config)
    cfg = json.loads(cfg_path.read_text())
    out = (cfg_path.resolve().parent / cfg["artifacts"]).resolve()
    out.mkdir(parents=True, exist_ok=True)

    U = read_gt_tsv(ROOT / "results/ground_truth/kinase_universe.tsv")
    kin = U[["kinase", "uniprot", "manning_group", "manning_family"]]
    kin.to_csv(out / "kinases.tsv", sep="\t", index=False)
    edges = build_edges(cfg, list(kin["kinase"]))
    edges.to_csv(out / "edges.tsv.gz", sep="\t", index=False)
    print(f"edges {len(edges)} ({edges['source'].value_counts().to_dict()}), sites {edges['site'].nunique()}")

    # ---- rho_e refit on DEV nulls (label-free calibration input)
    _, C, _, _ = load_gt()
    nulls = list(C.loc[(C["role"] == "null") & C["split"].isin(FOLDS), "conditionID"])
    Xn = load_dev(conditions=nulls)
    model = core.KinaseActivityModel(cfg, edges, kin)
    po = model.posterior(model.prepare(Xn.astype(np.float64)))
    with np.errstate(all="ignore"):
        zf = po["mu"] / np.sqrt(po["q"])
        me = np.where(po["q"] > 0, po["suma"] ** 2 / po["q"], 0.0)
    rho, tab = core.fit_rho_e(zf, me, po["report"])
    print(f"rho_e refit on {len(nulls)} DEV nulls: {rho} (config {cfg['rho_e']})")
    check = None
    if a.check_bayes2:
        from kinact.explore import bayes2 as B2
        rb, _ = B2.fit_rho_null(zf, me, po["report"])
        check = float(rb)
        assert rb == rho
    if abs(rho - cfg["rho_e"]) > 1e-12:
        raise SystemExit(f"rho_e refit {rho} != config {cfg['rho_e']}")

    man = {
        "name": cfg["name"], "built": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "config_sha256": hashlib.sha256(cfg_path.read_bytes()).hexdigest(),
        "libraries_sha256": {n: sha256_file(LIBDIR / f"{n}.parquet")
                             for n in [cfg["library"]["curated"], *cfg["library"].get("predicted", {})]},
        "universe_sha256": sha256_file(ROOT / "results/ground_truth/kinase_universe.tsv"),
        "n_edges": int(len(edges)), "n_edges_by_source": edges["source"].value_counts().to_dict(),
        "n_sites": int(edges["site"].nunique()), "n_kinases": int(len(kin)),
        "site_noise_model": None,
        "rho_e_refit": float(rho), "rho_e_fit_conditions": nulls, "rho_e_fit_table": tab,
        "rho_e_bayes2_crosscheck": check,
        "data": "kinase_universe.tsv (labels-free), cache/libraries, cache/dev (DEV role=null conditions only)",
    }
    (out / "manifest.json").write_text(json.dumps(man, indent=1, default=float))
    print("wrote", out)


if __name__ == "__main__":
    main()
