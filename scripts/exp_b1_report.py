#!/usr/bin/env python
"""B1 report: reads results/b1/grid.tsv (from exp_b1_grid.py), computes
  - nested-CV estimate of the selected config (choose on 4 folds by M1-sel, score the 5th; pooled OOF scores),
  - shuffled-library control (3 seeds) for the overall best config,
  - M2 vs PSP z-score, paired group bootstrap (B=2000, seed 0) vs PSP z-score and vs prior-only (OOF freq),
and writes results/b1/b1_report.md + results/b1/b1_report.json.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
from kinact import bootstrap as bs, controls as ctl, ledger  # noqa: E402
from kinact.eval import to_jsonable  # noqa: E402
from kinact.explore.lib_variants import FOLDS, DevData, Engine, cfg_name, gt_freeze_hash  # noqa: E402
from exp_b1_grid import summarise  # noqa: E402

OUT = ROOT / "results/b1"
REF = {"stat": "wz", "lib": {"type": "static", "name": "psp_all"}, "min_sub": 3}
f3 = lambda x: "NA" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{x:.3f}"


def family(r):
    c = json.loads(r["cfg"])["lib"]
    if c["type"] == "static":
        return c["name"]
    if c["type"] == "pred":
        return f"{c['name']}/{c['mode']}"
    return f"mix:{c['cur']}+{c['pred']['name']}"


def main():
    t0 = time.time()
    gth = gt_freeze_hash()
    G = pd.read_csv(OUT / "grid.tsv", sep="\t")
    G["family"] = G.apply(family, axis=1)
    dd = DevData()
    eng = Engine(dd)
    ev = dd.ev
    cfgs = [json.loads(c) for c in G["cfg"]]

    # ---- controls / reference
    S_ref = eng.score(REF)
    po = []
    for f in FOLDS:
        ec = [c for c, fo in zip(dd.conds, dd.fold) if fo == f]
        tc = list(dd.C.loc[dd.C["split"].isin([x for x in FOLDS if x != f]), "conditionID"])
        po.append(ctl.prior_only(dd.P, dd.C, tc, ec, list(dd.kin)))
    S_prior = pd.concat(po).loc[dd.conds]

    # ---- best apparent + nested CV (pooled OOF)
    ibest = int(G["M1-sel"].idxmax())
    chosen, parts = {}, []
    for f in FOLDS:
        j = int(G[f"tr_{f}"].idxmax())
        chosen[f] = {"config": G.loc[j, "name"], "train_M1sel": G.loc[j, f"tr_{f}"], "test_M1sel": G.loc[j, f"te_{f}"]}
        S = eng.score(cfgs[j])
        parts.append(S.loc[[c for c, fo in zip(dd.conds, dd.fold) if fo == f]])
    S_nest = pd.concat(parts).loc[dd.conds]
    S_best = eng.score(cfgs[ibest])

    res = {"ref_psp_z": ev.evaluate(S_ref), "prior_only": ev.evaluate(S_prior), "best": ev.evaluate(S_best),
           "nested_OOF": ev.evaluate(S_nest)}
    # nested within families that are pre-specified alternatives (curated only; predicted-only)
    fam_nest = {}
    for lab, mask in {"curated_static_only": G["lib_type"] == "static",
                      "predicted_only": G["lib_type"] == "pred", "mixtures_only": G["lib_type"] == "mix",
                      "wz_only": G["stat"] == "wz"}.items():
        sub = G[mask]
        pp = []
        for f in FOLDS:
            j = int(sub[f"tr_{f}"].idxmax())
            pp.append(eng.score(cfgs[j]).loc[[c for c, fo in zip(dd.conds, dd.fold) if fo == f]])
        fam_nest[lab] = ev.evaluate(pd.concat(pp).loc[dd.conds])["M1-sel"]["value"]

    # ---- shuffled-library control for the best config
    shuf = []
    for seed in [0, 1, 2]:
        e2 = Engine(dd, perm_seed=seed)
        r = ev.evaluate(e2.score(cfgs[ibest]))
        shuf.append(summarise(r))
        res[f"shuffled_s{seed}"] = r

    # ---- M2 + bootstrap
    m2 = {"best": ev.m2(S_best, S_ref), "nested_OOF": ev.m2(S_nest, S_ref)}
    plan = bs.BootstrapPlan.from_evaluator(ev, B=2000, seed=0)
    deltas = {}
    for a_lab, A in [("best", S_best), ("nested_OOF", S_nest)]:
        for b_lab, B in [("ref_psp_z", S_ref), ("prior_only", S_prior)]:
            d = bs.paired_delta(ev, A, B, "M1-sel", plan)
            deltas[f"{a_lab}_vs_{b_lab}"] = {k: v for k, v in d.items() if k != "replicates"}

    summ = {k: summarise(v) for k, v in res.items()}
    for lab in ["nested_OOF", *[f"shuffled_s{s}" for s in range(3)]]:
        lm = {k: summ[lab][k] for k in ["M1-sel", "M1-all", "M1-sel-scored", "cov_sel", "harmony"]}
        lm["gt_freeze"] = gth[:16]
        ledger.append(f"b1_{lab.split('_s')[0]}", {"b1_control": lab, "base": G.loc[ibest, "name"]}, "dev",
                      to_jsonable(lm), branch="B1", folds=range(5))

    # ---- write report
    cols = ["M1-sel", "M1-all", "M1-sel-scored", "cov_sel", "harmony", "h_0", "h_1-4", "h_5-19", "h_>=20",
            "src_H", "src_D", "src_O", "sign_-1", "sign_1", "n_sig"]
    hdr = "| config | M1-sel | M1-all | scored | cov | harm | b0 | b1-4 | b5-19 | b≥20 | H | D | O | s− | s+ | nsig |"
    sep = "|" + "---|" * 16

    def line(name, d):
        return f"| {name} | " + " | ".join(str(int(d[c])) if c == "n_sig" else f3(d[c]) for c in cols) + " |"

    L = ["# B1 — kinase–substrate library variants (DEV folds 0–4 pooled)", "",
         f"GT_FREEZE sha256 `{gth[:16]}`; {len(dd.conds)} GT conditions, {len(G)} configs (all in "
         "`results/dev_ledger.tsv`, branch B1). Scripts: `scripts/exp_b1_grid.py` → `scripts/exp_b1_report.py`; "
         "code `src/kinact/explore/lib_variants.py` (definitions in its docstring). Single-sample scoring, no GT "
         "information in any method. `wz` = Σwx/√Σw²/σ_c, `mean` = Σwx/Σw, `ridge` = ridge t (λ=1). "
         "Predicted libs: `meas` = top-k among sites measured in the condition, `glob` = global top-k then "
         "restricted to measured; weight (1−(r−1)/k)^γ. Mixtures: curated (w=1) + predicted top-k measured "
         "(γ=1, w·α) for kinases with < m measured curated substrates.", "",
         "## 1. Reference, controls, best per library family (M1-sel; harmony bins b*; source; sign)", "",
         hdr, sep,
         line("PSP z-score (psp_all, min3)", summ["ref_psp_z"]),
         line("prior-only (OOF freq, oracle sign)", summ["prior_only"])]
    fb = G.loc[G.groupby("family")["M1-sel"].idxmax()].sort_values("M1-sel", ascending=False)
    for _, r in fb.iterrows():
        L.append(line(r["name"], r))
    L += ["", "## 2. Best per statistic and top-10 overall", "", hdr, sep]
    for st in ["wz", "mean", "ridge"]:
        r = G.loc[G[G["stat"] == st]["M1-sel"].idxmax()]
        L.append(line("best " + st + ": " + r["name"], r))
    for _, r in G.sort_values("M1-sel", ascending=False).head(10).iterrows():
        L.append(line(r["name"], r))
    # library effects on harmony bins: best per family for each bin
    L += ["", "## 3. Factor effects (mean M1-sel over the grid cells sharing the factor; wz+mean, predicted only)", ""]
    pr = G[G["lib_type"] == "pred"].copy()
    pr["pred"] = pr["cfg"].map(lambda c: json.loads(c)["lib"]["name"])
    pr["k"] = pr["cfg"].map(lambda c: json.loads(c)["lib"]["k"])
    pr["mode"] = pr["cfg"].map(lambda c: json.loads(c)["lib"]["mode"])
    pr["gamma"] = pr["cfg"].map(lambda c: json.loads(c)["lib"]["gamma"])
    pr = pr[pr["stat"] != "ridge"]
    for fac in ["pred", "mode", "k", "gamma", "stat"]:
        t = pr.groupby(fac)[["M1-sel", "M1-sel-scored", "cov_sel", "h_0", "h_1-4", "h_5-19", "h_>=20"]].mean()
        L.append(f"- **{fac}**: " + "; ".join(f"{i}: {f3(v['M1-sel'])} (b0 {f3(v['h_0'])}, b1-4 {f3(v['h_1-4'])}, "
                                               f"b5-19 {f3(v['h_5-19'])}, b≥20 {f3(v['h_>=20'])})"
                                               for i, v in t.iterrows()))
    mx = G[G["lib_type"] == "mix"].copy()
    mx["alpha"] = mx["cfg"].map(lambda c: json.loads(c)["lib"]["alpha"])
    mx["m"] = mx["cfg"].map(lambda c: str(json.loads(c)["lib"]["m"]))
    mx = mx[mx["stat"] != "ridge"]
    for fac in ["alpha", "m"]:
        t = mx.groupby(fac)[["M1-sel", "cov_sel"]].mean()
        L.append(f"- **mixture {fac}**: " + "; ".join(f"{i}: {f3(v['M1-sel'])}" for i, v in t.iterrows()))
    L += ["", "## 4. Selection optimism: nested CV (choose config by M1-sel on 4 folds, score the 5th; pooled OOF)", "",
          f"- Apparent best (all 5 folds): **{G.loc[ibest, 'name']}** M1-sel {f3(G.loc[ibest, 'M1-sel'])}.",
          f"- **Nested OOF M1-sel {f3(summ['nested_OOF']['M1-sel'])}** (M1-all {f3(summ['nested_OOF']['M1-all'])}, "
          f"scored {f3(summ['nested_OOF']['M1-sel-scored'])}, cov {f3(summ['nested_OOF']['cov_sel'])}, harmony "
          f"{f3(summ['nested_OOF']['harmony'])}); optimism = {G.loc[ibest, 'M1-sel'] - summ['nested_OOF']['M1-sel']:+.3f}.",
          "- Chosen per held-out fold: " + "; ".join(f"{f}: {v['config']} (train {f3(v['train_M1sel'])}, "
                                                    f"held-out {f3(v['test_M1sel'])})" for f, v in chosen.items()),
          "- Nested OOF restricted to one family: " + "; ".join(f"{k} {f3(v)}" for k, v in fam_nest.items()), "",
          "## 5. Shuffled-library control (best config; site-label bijection, per-kinase counts/weights kept)", "",
          hdr, sep, line("best: " + G.loc[ibest, "name"], summ["best"])]
    for s in range(3):
        L.append(line(f"shuffled seed {s}", summ[f"shuffled_s{s}"]))
    L += ["", "## 6. M2 (matched universe) vs PSP z-score and paired group bootstrap (B=2000, seed 0, C8 basic CI)", ""]
    for k, v in m2.items():
        L.append(f"- M2 {k} vs PSP z: {f3(v['A'])} vs {f3(v['B'])} (Δ {v['delta']:+.3f}; {v['n_pairs']} pairs, "
                 f"{v['n_kinases']} kinases, pair fraction {f3(v['pair_fraction'])})")
    for k, d in deltas.items():
        L.append(f"- ΔM1-sel {k}: {d['estimate']:+.3f} [95% CI {d['ci_low']:+.3f}, {d['ci_high']:+.3f}], "
                 f"one-sided p={d['p_one_sided']:.4f}")
    (OUT / "b1_report_auto.md").write_text("\n".join(L) + "\n")
    js = {"gt_freeze_sha256": gth, "best": G.loc[ibest, "name"], "summ": summ, "chosen": chosen,
          "family_nested": fam_nest, "m2": m2, "deltas": deltas, "runtime_s": time.time() - t0}
    (OUT / "b1_report.json").write_text(json.dumps(to_jsonable(js), indent=1, default=str))
    print("\n".join(L))


if __name__ == "__main__":
    main()
