#!/usr/bin/env python
"""Branch B4b evaluation: nested OOF variants from results/b4b/fit_<fold>.npz vs core-v1 (z_vif, D7).
Writes results/b4b/b4b_eval.json and results/b4b/b4b_tables.md; ledger branch B4b."""
import json

import numpy as np
import pandas as pd

from exp_b21_common import calib_bins, full_eval, t_unseen  # noqa: F401
from exp_b4b_common import FOLDS, OUT, Ctx, bs, gt_freeze_hash, ledger_row
from kinact.bootstrap import metric_replicates, shifted_null_p
from kinact.eval import to_jsonable

LAMS = ["0.1", "0.01", "0.001"]
BETAS = ["0.0", "0.25", "0.5", "0.75", "1.0"]


def unseen_mask(ctx, Pp):
    P, C = ctx.P, ctx.C
    gc = C[C["split"].isin(FOLDS) & (C["role"] == "gt")].set_index("conditionID")["split"]
    pos = P[(P["cls"] == "direct") & P["confidence"].isin(["high", "medium"]) & P["conditionID"].isin(gc.index)]
    pos = pos.assign(fold=pos["conditionID"].map(gc))
    seen = {f: set(pos.loc[pos["fold"] != f, "kinase"]) for f in FOLDS}
    f = Pp["conditionID"].map(ctx.fold_of)
    return np.array([k not in seen[ff] for k, ff in zip(Pp["kinase"], f)])


def main():
    ctx = Ctx(ks={"r2np": [100]})
    ev = ctx.ev
    fits = {f: np.load(OUT / f"fit_{f}.npz") for f in FOLDS}
    infos = {f: json.loads((OUT / f"fit_{f}.json").read_text()) for f in FOLDS}
    cf = ctx.cf

    def oof(choice):
        """choice: fold -> key in the fold's npz."""
        A = np.full((len(ctx.conds), len(ctx.kin)), np.nan)
        for f in FOLDS:
            A[cf == f] = fits[f][choice[f]][cf == f]
        return pd.DataFrame(A, index=ctx.conds, columns=ctx.kin)

    def pooled_inner(key, f):
        """M1-sel over the other 4 folds of the fold-f scores (for parameter-free grid selection)."""
        S = pd.DataFrame(fits[f][key], index=ctx.conds, columns=ctx.kin)
        Pp = ev.pairs(S)
        return ev.metric(Pp[Pp["conditionID"].map(ctx.fold_of) != f], "M1-sel")["value"]

    # nested choices
    lam_ch = {}
    for f in FOLDS:
        im = infos[f]["inner_m1"]
        best = max(im, key=im.get)
        lam_ch[f] = "core" if best == "inf" else f"outer_lam{best}"
    beta_inner = {f: {bt: pooled_inner(f"prom_beta{bt}", f) for bt in BETAS} for f in FOLDS}
    beta_ch = {f: "prom_beta" + max(beta_inner[f], key=beta_inner[f].get) for f in FOLDS}
    variants = {
        "core-v1 z_vif (ref)": {f: "core" for f in FOLDS},
        "learned, λ nested (PRIMARY)": lam_ch,
        **{f"learned, λ={l} fixed": {f: f"outer_lam{l}" for f in FOLDS} for l in LAMS},
        "learned logistic loss, λ=0.01": {f: "outer_logistic_lam0.01" for f in FOLDS},
        "promiscuity^-β, β nested": beta_ch,
        **{f"promiscuity^-β, β={bt} fixed": {f: f"prom_beta{bt}" for f in FOLDS} for bt in BETAS[1:]},
    }
    zpost = {"core-v1 z_post (rank output)": {f: "core__zpost" for f in FOLDS},
             "learned λ nested, z_post": {f: lam_ch[f] + "__zpost" for f in FOLDS},
             "promiscuity β nested, z_post": {f: beta_ch[f] + "__zpost" for f in FOLDS}}
    plan = bs.BootstrapPlan.from_evaluator(ev, B=2000, seed=0)
    S_ref = oof(variants["core-v1 z_vif (ref)"])
    Pref = ev.pairs(S_ref)
    um = unseen_mask(ctx, Pref)
    rep_ref_u = metric_replicates(ev, Pref[um], "M1-sel", plan)
    res = {}
    for name, ch in {**variants, **zpost}.items():
        S = oof(ch)
        r = full_eval(ctx, S)
        Pp = ev.pairs(S)
        r["by_sign"] = {s: ev.metric(Pp, f"M1-sel:sign={s}")["value"] for s in ["-1", "1"]}
        if name != "core-v1 z_vif (ref)":
            r["M2_vs_ref"] = ev.m2(S, S_ref)
            d = bs.paired_delta(ev, S, S_ref, "M1-sel", plan)
            r["dM1sel"] = {k: d[k] for k in ["estimate", "ci_low", "ci_high", "p_one_sided"]}
            for ep in ["M1-all", "harmony"]:
                dd = bs.paired_delta(ev, S, S_ref, ep, plan)
                r[f"d{ep}"] = {k: dd[k] for k in ["estimate", "ci_low", "ci_high", "p_one_sided"]}
            ru = metric_replicates(ev, Pp[um], "M1-sel", plan) - rep_ref_u
            pt = ev.metric(Pp[um], "M1-sel")["value"] - ev.metric(Pref[um], "M1-sel")["value"]
            q = np.quantile(ru, [0.025, 0.975])
            r["dT-unseen"] = {"estimate": pt, "ci_low": 2 * pt - q[1], "ci_high": 2 * pt - q[0],
                              "p_one_sided": shifted_null_p(pt, ru)}
        if "__zpost" not in next(iter(ch.values())):
            Zn, Nn, Zp, Np = [], [], [], []
            for f in FOLDS:
                nf = fits[f]["null_fold"] == f
                Zn.append(fits[f][ch[f] + "__null"][nf]); Nn.append(fits[f]["n_null"][nf])
                Zp.append(fits[f][ch[f] + "__perm"][cf == f]); Np.append(fits[f]["n_perm"][cf == f])
            r["calib_null"] = calib_bins(np.vstack(Zn), np.vstack(Nn))
            r["calib_perm"] = calib_bins(np.vstack(Zp), np.vstack(Np))
        r["chosen"] = ch
        res[name] = r
        ledger_row("bayes2_edgeweights", {"variant": name, "chosen": ch}, r)
        print(name, round(r["M1-sel"], 4), r.get("dM1sel"), flush=True)
    feats = infos["fold0"]["features"]
    theta = {l: np.array([infos[f]["theta"][f"outer_lam{l}"] for f in FOLDS]) for l in LAMS}
    theta_summary = {l: {ft: [float(t.mean()), float(t.std()), float(t.min()), float(t.max())]
                         for ft, t in zip(feats, th.T)} for l, th in theta.items()}
    out = {"gt_freeze_sha256": gt_freeze_hash(), "results": res, "lambda_choice": lam_ch,
           "inner_m1": {f: infos[f]["inner_m1"] for f in FOLDS}, "beta_choice": beta_ch, "beta_inner": beta_inner,
           "theta_summary": theta_summary, "rho_e": {f: infos[f]["rho_e"] for f in FOLDS},
           "wall_s_fit": {f: infos[f]["wall_s"] for f in FOLDS}}
    (OUT / "b4b_eval.json").write_text(json.dumps(to_jsonable(out), indent=1, default=str))
    tables(out)


def tables(out):
    f3 = lambda x: "–" if x is None else f"{x:.3f}"
    dfmt = lambda d: "ref" if d is None else f"{d['estimate']:+.3f} [{d['ci_low']:+.3f}, {d['ci_high']:+.3f}] p={d['p_one_sided']:.3f}"
    L = ["| variant | M1-sel | M1-all | harm (0/1-4/5-19/≥20) | H/D/O | sign −/+ | T-uns | nsig | M2 A vs ref | ΔM1-sel | ΔM1-all | ΔT-unseen | null |z|>1.96 all (bins) | perm all (bins) |",
         "|" + "---|" * 15]
    for n, r in out["results"].items():
        hb = r["harmony_per_bin"]
        cal = lambda c: "–" if not c else f"{c['all']:.3f} ({'/'.join(f3(c.get(b)) for b in ['0', '1-4', '5-19', '>=20'])})"
        m2 = r.get("M2_vs_ref")
        L.append(f"| {n} | {r['M1-sel']:.3f} | {r['M1-all']:.3f} | {r['harmony']:.3f} ({'/'.join(f3(hb.get(b)) for b in ['0', '1-4', '5-19', '>=20'])}) | "
                 f"{'/'.join(f3(r['by_source'].get(s)) for s in 'HDO')} | {f3(r['by_sign']['-1'])}/{f3(r['by_sign']['1'])} | "
                 f"{f3(r['T-unseen'])} | {r['n_sig']} | {'–' if not m2 else f'{m2['A']:.3f} vs {m2['B']:.3f}'} | "
                 f"{dfmt(r.get('dM1sel'))} | {dfmt(r.get('dM1-all'))} | {dfmt(r.get('dT-unseen'))} | "
                 f"{cal(r.get('calib_null'))} | {cal(r.get('calib_perm'))} |")
    L += ["", f"λ choice per fold: {out['lambda_choice']}", f"inner M1-sel: {out['inner_m1']}",
          f"β choice: {out['beta_choice']}", "", "θ (outer fits; mean ± sd [min, max] over 5 folds):"]
    for l, d in out["theta_summary"].items():
        L.append(f"- λ={l}: " + "; ".join(f"{k} {v[0]:+.2f}±{v[1]:.2f}" for k, v in d.items()))
    L.append(f"ρ_e refit: {out['rho_e']}")
    (OUT / "b4b_tables.md").write_text("\n".join(L))
    print("\n".join(L))


if __name__ == "__main__":
    main()
