#!/usr/bin/env python
"""B5 part 1: self-site coverage, standalone predictive value of s_k, paradoxical increases by drug/target class.
Writes results/b5/signal.json and results/b5/signal_tables.md."""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import exp_b5_common as CM  # noqa: E402
from kinact import ledger  # noqa: E402
from kinact.eval import to_jsonable  # noqa: E402


def main():
    ctx = CM.Ctx()
    ev, kin, conds = ctx.ev, ctx.kin, ctx.conds
    out = {"gt_freeze_sha256": CM.E2.gt_freeze_hash()}
    L = []
    allX = ~np.isnan(ctx.X.to_numpy())                       # all DEV conditions
    tgt = ev.targets
    for w, tb in ctx.tables.items():
        s, n = ctx.signal(w)
        inc = tb[tb["site"].isin(ctx.site_index)]
        meas_any = inc[inc["site"].isin(ctx.site_index[allX.any(1)])]
        nk_cond = (n > 0).sum(1)
        tp = tgt.assign(n_self=n.to_numpy()[tgt["ci"], tgt["ki"]] if False else
                        [n.at[c, k] for c, k in zip(tgt["conditionID"], tgt["kinase"])])
        sel = tp[tp["in_M1sel"]]
        cov = {"n_sites": len(tb), "n_kinases_table": int(tb["kinase"].nunique()),
               "n_inhib_sites": int((tb["sign"] < 0).sum()),
               "n_kinases_site_in_compendium_DEV": int(inc["kinase"].nunique()),
               "n_kinases_measured_any_DEV_cond": int(meas_any["kinase"].nunique()),
               "n_kinases_measured_in_GT_conds": int(((n > 0).sum(0) > 0).sum()),
               "median_kinases_per_GT_cond": float(nk_cond.median()),
               "M1sel_pairs": int(len(sel)), "M1sel_pairs_measured": int((sel["n_self"] > 0).sum()),
               "M1sel_kinases": int(sel["kinase"].nunique()),
               "M1sel_kinases_measured_any": int(sel.loc[sel["n_self"] > 0, "kinase"].nunique()),
               "all_pairs": int(len(tp)), "all_pairs_measured": int((tp["n_self"] > 0).sum()),
               "by_sign_measured_frac": {str(g): float((d["n_self"] > 0).mean()) for g, d in sel.groupby("sign")},
               "by_source_measured_frac": {g: float((d["n_self"] > 0).mean()) for g, d in sel.groupby("source")}}
        # standalone: harness, scored-only
        r = ev.evaluate(s, n_perm=200)
        Pp = ev.pairs(s)
        stand = {"M1-sel": r["M1-sel"]["value"], "M1-sel-scored": r["M1-sel-scored"]["value"],
                 "coverage": r["coverage_M1-sel"]}
        for ep in ["M1-sel-scored:sign=-1", "M1-sel-scored:sign=1", "M1-sel-scored:source=H",
                   "M1-sel-scored:source=D", "M1-sel-scored:source=O", "M1-all-scored"]:
            m = ev.metric(Pp, ep)
            stand[ep] = {"value": m["value"], "n_pairs": m["n_pairs"], "n_kinases": m["n_kinases"]}
        # direction: sign(s_k) vs GT sign on measured target pairs (all pairs, unweighted)
        tp["s"] = [s.at[c, k] for c, k in zip(tp["conditionID"], tp["kinase"])]
        m = tp[tp["n_self"] > 0]
        stand["direction_agree"] = {str(g): {"frac_agree": float((np.sign(d["s"]) == d["sign"]).mean()),
                                             "mean_s": float(d["s"].mean()), "n": int(len(d))}
                                    for g, d in m.groupby("sign")}
        out[w] = {"coverage": cov, "standalone": stand}
        ledger.append("b5_selfsite_alone", {"set": w, "signal": "mean(sign*y)/rms"}, "dev",
                      {"M1-sel": stand["M1-sel"], "M1-sel-scored": stand["M1-sel-scored"], "coverage": stand["coverage"]},
                      branch="B5", folds=range(5))
        if w == "psp":
            # paradox by target family / drug (inhibitor pairs, s=-1)
            Cd = ctx.meta["Condition"]
            fam = ctx.Udf.set_index("kinase")["manning_family"]
            grp = ctx.Udf.set_index("kinase")["manning_group"]
            inh = m[m["sign"] < 0].assign(fam=lambda d: d["kinase"].map(fam), grp=lambda d: d["kinase"].map(grp),
                                         drug=lambda d: d["conditionID"].map(Cd))
            byfam = inh.groupby("fam").agg(n=("s", "size"), n_kin=("kinase", "nunique"), mean_s=("s", "mean"),
                                           frac_up=("s", lambda x: float((x > 0).mean())),
                                           kinases=("kinase", lambda x: ",".join(sorted(set(x))[:6])))
            out["paradox_by_family"] = byfam.sort_values("n", ascending=False).reset_index().to_dict("records")
            out["paradox_by_group"] = inh.groupby("grp").agg(n=("s", "size"), mean_s=("s", "mean"),
                                                             frac_up=("s", lambda x: float((x > 0).mean()))).reset_index().to_dict("records")
            out["paradox_by_source"] = inh.groupby("source").agg(n=("s", "size"), mean_s=("s", "mean"),
                                                                 frac_up=("s", lambda x: float((x > 0).mean()))).reset_index().to_dict("records")
            # within-target-kinase: per target kinase most-up examples
            byk = inh.groupby("kinase").agg(n=("s", "size"), mean_s=("s", "mean"), frac_up=("s", lambda x: float((x > 0).mean())),
                                           drugs=("drug", lambda x: ";".join(sorted(set(map(str, x)))[:4])))
            out["paradox_by_kinase"] = byk[byk["n"] >= 3].sort_values("mean_s", ascending=False).reset_index().to_dict("records")
            stim = m[m["sign"] > 0].assign(fam=lambda d: d["kinase"].map(fam))
            out["stim_by_family"] = stim.groupby("fam").agg(n=("s", "size"), mean_s=("s", "mean"),
                                                            frac_up=("s", lambda x: float((x > 0).mean()))).sort_values("n", ascending=False).reset_index().to_dict("records")
            inh.to_csv(CM.OUT / "inhibitor_target_selfsignal.tsv", sep="\t", index=False)
    (CM.OUT / "signal.json").write_text(json.dumps(to_jsonable(out), indent=1, default=str))
    print(json.dumps(to_jsonable(out), indent=1, default=str)[:12000])


if __name__ == "__main__":
    main()
