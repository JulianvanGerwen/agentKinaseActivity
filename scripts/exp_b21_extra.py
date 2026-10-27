#!/usr/bin/env python
"""B2.1 follow-up: nested variants for the recommended family (out grid, winsor 3 only): z_vif / empnull / mean;
paired bootstrap empnull vs mean and vs 'B2 config, new code'; fixed recommended config (non-nested, optimistic)."""
import numpy as np, pandas as pd
import exp_b21_eval as EV
import exp_b21_grid as G
from exp_b21_common import Ctx, nested, full_eval, bs, ledger_row, b2_reference, dump, calib_bins, FOLDS

ctx = Ctx(ks={"r2np": [100]})
g = EV.load_grids(); g = g[g["exp"] == "out"]
sc = EV.Scorer(ctx, g, None)
w3 = g["winsor"] == 3.0; a1 = g["alphas"] == "{'r2np@100': 1.0}"
V = {"w3: mean": w3 & (g.output == "mean"), "w3: empnull": w3 & (g.output == "empnull"),
     "w3: z_vif": w3 & (g.output == "z_vif"), "w3: z_freq": w3 & (g.output == "z_freq"),
     "w0: B2 config (z_post)": g["winsor"].isna() & (g.output == "z_post") & a1}
S, ch = {}, {}
for n, m in V.items():
    S[n], ch[n] = nested(g, sc, ctx, m, hp=["output", "winsor", "tau", "rho_f", "rho_g", "alphas"])
ref, _ = b2_reference(ctx)
plan = bs.BootstrapPlan.from_evaluator(ctx.ev, B=2000, seed=0)
res = {}
for n in S:
    r = full_eval(ctx, S[n]); r["calibration"] = EV.calib_variant(ctx, g, ch[n], None)
    r["d_vs_B2ref"] = bs.paired_delta(ctx.ev, S[n], ref, "M1-sel", plan)
    r["d_vs_B2cfg_newcode"] = bs.paired_delta(ctx.ev, S[n], S["w0: B2 config (z_post)"], "M1-sel", plan)
    r["chosen"] = ch[n]; res[n] = r
    ledger_row("bayes2_nested_extra", {"variant": n, "chosen": {f: c["key"] for f, c in ch[n].items()}}, r)
    d, d2 = r["d_vs_B2ref"], r["d_vs_B2cfg_newcode"]
    print(f"{n}: M1-sel {r['M1-sel']:.4f} M1-all {r['M1-all']:.4f} harm {r['harmony_per_bin']} T-uns {r['T-unseen']:.3f} "
          f"src {r['by_source']} | vsB2 {d['estimate']:+.3f} [{d['ci_low']:+.3f},{d['ci_high']:+.3f}] p={d['p_one_sided']:.3f} "
          f"| vsB2cfg {d2['estimate']:+.3f} [{d2['ci_low']:+.3f},{d2['ci_high']:+.3f}] p={d2['p_one_sided']:.3f} | calib {r['calibration']}")
    print("   chosen", {f: (c['output'], c['tau'], c['rho_f'], c['rho_g']) for f, c in ch[n].items()})
d = bs.paired_delta(ctx.ev, S["w3: empnull"], S["w3: mean"], "M1-sel", plan)
print("empnull - mean", {k: round(d[k], 4) for k in ["estimate", "ci_low", "ci_high", "p_one_sided"]})
res["empnull_minus_mean"] = d
dump("b21_extra.json", res)
