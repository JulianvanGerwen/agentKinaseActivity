#!/usr/bin/env python
"""B5 part 2: combine self-site signal with the B2 reference core (fixed B2 config; z or posterior mean).
 (a) late additive: core/rms_c + w * T(s), T by mode (add | agree | mag | up), set (psp | aloop | psp+aloop)
 (b) Bayesian extra observation: self-sites enter the likelihood as y_i = lam * sign_i * a_k + eps_i
Hyperparameters per family (<=3): late (set, mode, w); bayes (set, lam). Nested CV over DEV folds 0-4 (per outer
fold: argmax inner M1-sel on the other folds, w=0 / lam=0 = reference included). Every config -> ledger (B5).
Writes results/b5/grid.tsv, results/b5/combine.json, results/b5/combine_tables.md."""
import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import exp_b5_common as CM  # noqa: E402
import exp_b2_eval as E2E  # noqa: E402
from kinact import bootstrap as bs, ledger  # noqa: E402
from kinact.eval import to_jsonable  # noqa: E402
from kinact.explore import selfsites as SS  # noqa: E402

MODES = ["add", "agree", "mag", "up"]
WS = [0.1, 0.25, 0.5, 1.0, 2.0]
LAMS = [0.25, 0.5, 1.0, 2.0, 4.0]


def main():
    t0 = time.time()
    ctx = CM.Ctx()
    ev, kin, conds, fold_of = ctx.ev, ctx.kin, ctx.conds, ctx.fold_of
    ref = ctx.core()
    sig = {w: ctx.signal(w)[0].to_numpy() for w in CM.SETS}
    rows, store = [], {}

    def add(cfg, S):
        Pp = ev.pairs(S)
        m = CM.E2.inner_metrics(ev, Pp, fold_of)
        key = ledger.config_hash(cfg)
        store[key] = S
        rows.append({**cfg, "config_hash": key, **m})
        ledger.append("b5_" + cfg["family"], cfg, "dev", {k: m[k] for k in ["M1-sel", "M1-all", "M1-sel-scored", "coverage"]},
                      branch="B5", folds=range(5))

    for out in ["z", "mean"]:
        add({"family": "ref", "core": out}, ref[out])
        c = ref[out].to_numpy()
        for st, mode, w in itertools.product(CM.SETS, MODES, WS):
            S = pd.DataFrame(SS.late_combine(c, sig[st], w, mode), index=conds, columns=kin)
            add({"family": "late", "core": out, "set": st, "mode": mode, "w": w}, S)
    for st, lam in itertools.product(CM.SETS, LAMS):
        d = SS.design_with_self(ctx.design, ctx.tables[st], ctx.site_index, lam)
        r = ctx.core(d)
        for out in ["z", "mean"]:
            add({"family": "bayes", "core": out, "set": st, "lam": lam}, r[out])
    g = pd.DataFrame(rows)
    g.to_csv(CM.OUT / "grid.tsv", sep="\t", index=False)
    print(f"grid {len(g)} configs {time.time()-t0:.0f}s", flush=True)

    def nested(mask):
        sub = g[mask]
        parts, chosen = [], {}
        cf = pd.Series(conds).map(fold_of).to_numpy()
        for f in CM.FOLDS:
            b = sub.loc[sub[f"inner_{f}"].idxmax()]
            chosen[f] = {k: b[k] for k in ["family", "core", "set", "mode", "w", "lam"] if k in b and pd.notna(b[k])}
            parts.append(store[b["config_hash"]].loc[np.asarray(conds)[cf == f]])
        return pd.concat(parts).loc[conds], chosen

    S, chosen = {"ref_z": ref["z"], "ref_mean": ref["mean"]}, {}
    isref = g["family"] == "ref"
    for out in ["z", "mean"]:
        co = g["core"] == out
        S[f"late_{out}"], chosen[f"late_{out}"] = nested(co & (isref | (g["family"] == "late")))
        S[f"bayes_{out}"], chosen[f"bayes_{out}"] = nested(co & (isref | (g["family"] == "bayes")))
        for mode in MODES:     # forced mode (no w=0 option) — what each rule does when switched on
            S[f"late_{out}_{mode}"], chosen[f"late_{out}_{mode}"] = nested(co & (g["family"] == "late") & (g["mode"] == mode))
        S[f"bayes_{out}_forced"], chosen[f"bayes_{out}_forced"] = nested(co & (g["family"] == "bayes"))
    S["any"], chosen["any"] = nested(g["family"].notna())
    res = {}
    for k, v in S.items():
        r = E2E.summarise(ev.evaluate(v))
        r["T-unseen"] = E2E.t_unseen(ev, v, ctx.P, ctx.C, fold_of)
        res[k] = r
        if k not in ("ref_z", "ref_mean"):
            ledger.append("b5_nested", {"variant": k, "chosen": chosen[k]}, "dev",
                          {x: r[x] for x in ["M1-sel", "M1-all", "M1-sel-scored", "coverage"]}, branch="B5", folds=range(5))
    plan = bs.BootstrapPlan.from_evaluator(ev, B=2000, seed=0)
    boot = {}
    for k in S:
        if k.startswith("ref"):
            continue
        b = "ref_mean" if ("mean" in k) else "ref_z"
        for ep in ["M1-sel", "M1-all", "harmony"]:
            d = bs.paired_delta(ev, S[k], S[b], ep, plan)
            boot[f"{k}_vs_{b}_{ep}"] = {x: float(d[x]) for x in ["estimate", "se", "ci_low", "ci_high", "p_one_sided"]}
        boot[f"{k}_vs_{b}_M2"] = ev.m2(S[k], S[b])
    d = bs.paired_delta(ev, S["ref_mean"], S["ref_z"], "M1-sel", plan)
    boot["ref_mean_vs_ref_z_M1-sel"] = {x: float(d[x]) for x in ["estimate", "se", "ci_low", "ci_high", "p_one_sided"]}
    # pooled (non-nested, optimistic) best per family
    pooled = {f"{fam}_{out}": g[(g["family"] == fam) & (g["core"] == out)].sort_values("M1-sel").iloc[-1][
        ["set", "mode", "w", "lam", "M1-sel", "M1-all"]].to_dict() for fam in ["late", "bayes"] for out in ["z", "mean"]}
    out = {"gt_freeze_sha256": CM.E2.gt_freeze_hash(), "n_configs": int(len(g)), "results": res, "chosen": chosen,
           "bootstrap": boot, "pooled_best": pooled, "wall_s": time.time() - t0}
    (CM.OUT / "combine.json").write_text(json.dumps(to_jsonable(out), indent=1, default=str))
    L = ["| variant | M1-sel | M1-all | cov | harm | b0 | b1-4 | b5-19 | b>=20 | H | D | O | s- | s+ | T-uns | nsig |",
         "|" + "---|" * 16]
    for k, r in res.items():
        hb = r["harmony_per_bin"]
        L.append(f"| {k} | {r['M1-sel']:.3f} | {r['M1-all']:.3f} | {r['coverage']:.3f} | {r['harmony']:.3f} | "
                 + " | ".join(f"{hb.get(b, np.nan):.3f}" for b in ["0", "1-4", "5-19", ">=20"]) + " | "
                 + " | ".join(f"{r['by_source'].get(s, np.nan):.3f}" for s in "HDO") + " | "
                 + f"{r['by_sign'].get('-1', np.nan):.3f} | {r['by_sign'].get('1', np.nan):.3f} | {r['T-unseen']['value']:.3f} | {r['n_sig']} |")
    L += ["", "| comparison | Δ | 95% CI | p |", "|---|---|---|---|"]
    for k, d in boot.items():
        if k.endswith("_M2"):
            L.append(f"| {k} | {d} | | |")
        else:
            L.append(f"| {k} | {d['estimate']:+.3f} | [{d['ci_low']:+.3f}, {d['ci_high']:+.3f}] | {d['p_one_sided']:.4f} |")
    L += ["", "Chosen per fold:"] + [f"- {k}: " + "; ".join(f"{f}: {json.dumps(c, default=str)}" for f, c in v.items())
                                     for k, v in chosen.items()]
    L += ["", "Pooled best (non-nested): " + json.dumps(pooled, default=str)]
    (CM.OUT / "combine_tables.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))
    print(f"done {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
