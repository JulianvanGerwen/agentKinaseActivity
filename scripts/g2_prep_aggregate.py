#!/usr/bin/env python
"""Aggregate results/g2_prep/resplit_seed{0..4}.json -> results/g2_prep/resplit_summary.tsv/json (per seed, mean±sd)."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
OUT = Path("/workspace/results/g2_prep")
M = ["core-v3", "core-v2", "core-v1", "B7 joint (nested)", "PSP z-score min3", "curated z-score min1", "prior-only", "z-score nested lib/min"]
rows, drows, picks = [], [], {}
for s in range(5):
    d = json.loads((OUT / f"resplit_seed{s}.json").read_text())
    for m in M:
        r = d["results"][m]
        rows.append({"seed": s, "method": m, **{k: r[k] for k in ["M1-sel", "M1-all", "harmony", "T-unseen", "M1-sel-scored", "coverage_M1-sel"]}})
    for comp, dd in d["delta_core_v3"].items():
        for ep, x in dd.items():
            drows.append({"seed": s, "comparison": comp, "endpoint": ep, **{k: x[k] for k in ["estimate", "ci_low", "ci_high", "p_one_sided"]}})
    picks[s] = {k: d[k] for k in ["core_picks", "b7_chosen", "ntop_kappa", "alpha", "zscore_nested_chosen", "fold_m1sel_pairs", "n_sig_kinases_outer"]}
R = pd.DataFrame(rows); D = pd.DataFrame(drows)
R.to_csv(OUT / "resplit_per_seed.tsv", sep="\t", index=False); D.to_csv(OUT / "resplit_deltas.tsv", sep="\t", index=False)
agg = {}
for m in M:
    x = R[R.method == m]
    agg[m] = {ep: {"seed0": float(x[x.seed == 0][ep].iloc[0]), "mean_1_4": float(x[x.seed > 0][ep].mean()),
                   "mean_0_4": float(x[ep].mean()), "sd_0_4": float(x[ep].std(ddof=1))} for ep in ["M1-sel", "M1-all", "harmony", "T-unseen"]}
dagg = {}
for (comp, ep), x in D.groupby(["comparison", "endpoint"]):
    dagg[f"{comp}|{ep}"] = {"seed0": float(x[x.seed == 0].estimate.iloc[0]), "mean_1_4": float(x[x.seed > 0].estimate.mean()),
                            "mean_0_4": float(x.estimate.mean()), "sd_0_4": float(x.estimate.std(ddof=1)),
                            "min": float(x.estimate.min()), "max_p": float(x.p_one_sided.max())}
(OUT / "resplit_summary.json").write_text(json.dumps({"methods": agg, "deltas": dagg, "picks": picks}, indent=1, default=str))
pd.set_option("display.width", 250)
print(R.pivot(index="method", columns="seed", values="M1-sel").round(4))
for ep in ["M1-all", "harmony", "T-unseen"]:
    print(ep); print(R.pivot(index="method", columns="seed", values=ep).round(4))
print(D[D.endpoint == "M1-sel"].pivot(index="comparison", columns="seed", values="estimate").round(4))
print(D[D.endpoint == "M1-sel"].pivot(index="comparison", columns="seed", values="p_one_sided").round(4))
print(D[D.endpoint == "T-unseen"].pivot(index="comparison", columns="seed", values="estimate").round(4))
for s, p in picks.items(): print(s, p["ntop_kappa"], p["alpha"], p["zscore_nested_chosen"], p["fold_m1sel_pairs"])
