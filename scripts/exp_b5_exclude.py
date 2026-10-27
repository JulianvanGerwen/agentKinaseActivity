#!/usr/bin/env python
"""B5 part 3: GT-A exclusion helper check — B2 reference core with all kinase self-/regulatory-activity sites
removed from the input (selfsites.mask_self_sites) vs the full input. Writes results/b5/exclude.json."""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import exp_b5_common as CM  # noqa: E402
from kinact import bootstrap as bs, ledger  # noqa: E402
from kinact.explore import selfsites as SS  # noqa: E402

ctx = CM.Ctx()
ref = ctx.core()
out = {"gt_freeze_sha256": CM.E2.gt_freeze_hash()}
plan = bs.BootstrapPlan.from_evaluator(ctx.ev, B=2000, seed=0)
for name, kw in {"psp_only": dict(which="psp", include_all_regulatory=False),
                 "conservative": dict(which="psp+aloop", include_all_regulatory=True)}.items():
    keys = SS.self_site_keys(ctx.Udf, **kw)
    Xm = SS.mask_self_sites(ctx.X[ctx.conds], keys)
    ctx.Y, ctx.D = CM.E2.standardise_by_fold(Xm, ctx.conds, ctx.meta, {}, noise=False)
    S = ctx.core()["z"]
    d = bs.paired_delta(ctx.ev, S, ref["z"], "M1-sel", plan)
    r = ctx.ev.evaluate(S, n_perm=200)
    out[name] = {"n_keys": len(keys), "n_keys_in_DEV": int(ctx.X.index.isin(list(keys)).sum()),
                 "M1-sel": r["M1-sel"]["value"], "by_source": {k: v["value"] for k, v in r["by_source"].items()},
                 "by_sign": {str(k): v["value"] for k, v in r["by_sign"].items()},
                 "delta_vs_full": {x: float(d[x]) for x in ["estimate", "ci_low", "ci_high"]}}
    ledger.append("b5_core_selfsites_masked", {"mask": name}, "dev", {"M1-sel": out[name]["M1-sel"]}, branch="B5", folds=range(5))
(CM.OUT / "exclude.json").write_text(json.dumps(out, indent=1))
print(json.dumps(out, indent=1))
