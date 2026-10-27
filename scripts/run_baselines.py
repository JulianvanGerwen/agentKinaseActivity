#!/usr/bin/env python
"""Run all DEV baselines (DESIGN §6) over the method x library grid, evaluate, ledger, benchmark, report.

    /home/plm_env/bin/python scripts/run_baselines.py            # everything (re-uses nothing: always rescored)
    ... --only zscore,ksea   --skip-jumpsem   --skip-bench   --jobs 32

For every config: score ALL DEV-cache conditions (cache/dev, folds 0-4 only), save
cache/scores/<method>__<lib>__<cfg>.parquet (+ .n.parquet with measured-substrate counts), evaluate with
kinact.eval (pooled DEV folds 0-4; M1-sel, M1-all, M1-*-scored + coverage, harmony, by source/sign,
coverage-significant), append to results/dev_ledger.tsv (branch 'baselines'), time with
kinact.runtime.bench_runtime on 200 fixed DEV conditions (seed 0). Writes results/baselines/baselines_dev.{md,json}.
The GT / DEV-cache identity (sha256 of GT_FREEZE.json, gt_pairs.tsv, gt_conditions.tsv, quant_dev.npy) is
recorded in every ledger config and in the report.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from kinact import ledger  # noqa: E402
from kinact.data import DEV_DIR, load_dev, load_dev_meta  # noqa: E402
from kinact.eval import Evaluator, load_gt, read_gt_tsv, to_jsonable  # noqa: E402
from kinact.methods import METHODS, run_jumpsem, study_units  # noqa: E402
from kinact.runtime import bench_runtime  # noqa: E402

FOLDS = [f"fold{i}" for i in range(5)]
SCORES = ROOT / "cache/scores"
OUT = ROOT / "results/baselines"
GT_DIR = ROOT / "results/ground_truth"
L3 = ["psp_human_only", "psp_all", "curated_union"]
SLOW = {"viper", "gsea", "rokai", "phosx_like"}  # process-parallel; n_rep=1 in the benchmark

# (method, library, cfg-label, kwargs, published_default)
#   published defaults: z-score/RoKAI-z, KSEA(KSEAapp), JUMPsem -> PSP human kinase x human substrate
#   (psp_human_only); decoupler methods/GSEA -> OmniPath-style curated KS network (curated_union, the
#   leakage-filtered analogue of dc.op.ksn); RoKAI -> its own network KS layer (rokai_ks); motif methods ->
#   Cantley atlas (KL: top-15 kinases per site; PhosX: top-10).
DEFAULT_LIB = {"zscore": "psp_human_only", "ksea": "psp_human_only", "wmean": "curated_union",
               "norm_wmean": "curated_union", "ulm": "curated_union", "mlm": "curated_union",
               "viper": "curated_union", "gsea": "curated_union", "rokai": "rokai_ks",
               "jumpsem": "psp_human_only", "jumpsem_aligned": "psp_human_only"}


def grid():
    G = []
    for m in ["zscore", "ksea"]:
        for lib in L3:
            for ms in [1, 3, 5]:
                G.append((m, lib, f"min{ms}", {"min_substrates": ms}, lib == DEFAULT_LIB[m] and ms == 3))
    for m in ["wmean", "norm_wmean", "ulm", "mlm", "viper", "gsea"]:
        for lib in L3:
            G.append((m, lib, "min3", {"min_substrates": 3}, lib == DEFAULT_LIB[m]))
    for lib in ["rokai_ks", "psp_all", "curated_union"]:
        G.append(("rokai", lib, "min3", {"min_substrates": 3}, lib == "rokai_ks"))
    for k in [50, 100, 200]:
        G.append(("motif_zscore", f"cantley@{k}", "min3", {"min_substrates": 3}, k == 100))
    for n in [10, 15]:
        G.append(("phosx_like", n, "min3", {"min_substrates": 3}, n == 10))
    for thr in [0.5, 1.0]:
        G.append(("kl_enrich", 15, f"thr{thr}_min3", {"min_substrates": 3, "thr": thr}, thr == 0.5))
    return G


def lib_label(m, lib):
    if m == "phosx_like":
        return f"cantley-site-top{lib}"
    if m == "kl_enrich":
        return f"cantley-site-top{lib}"
    return str(lib).replace("@", "-top")


def sha256(path, chunk=1 << 24):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while b := f.read(chunk):
            h.update(b)
    return h.hexdigest()


def data_identity():
    return {"GT_FREEZE.json": sha256(GT_DIR / "GT_FREEZE.json")[:16],
            "gt_pairs.tsv": sha256(GT_DIR / "gt_pairs.tsv")[:16],
            "gt_conditions.tsv": sha256(GT_DIR / "gt_conditions.tsv")[:16],
            "quant_dev.npy": sha256(DEV_DIR / "quant_dev.npy")[:16],
            "GT_FREEZE_mtime": time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime((GT_DIR / "GT_FREEZE.json").stat().st_mtime))}


def v(x):
    return x["value"] if isinstance(x, dict) and "value" in x else x


def summarise(r):
    return {"M1-sel": v(r["M1-sel"]), "M1-all": v(r["M1-all"]), "M1-sel-scored": v(r["M1-sel-scored"]),
            "M1-all-scored": v(r["M1-all-scored"]), "coverage_M1-sel": r["coverage_M1-sel"],
            "coverage_M1-all": r["coverage_M1-all"], "harmony": v(r.get("harmony")),
            "harmony_per_bin": (r.get("harmony") or {}).get("per_bin"),
            "by_source": {k: v(x) for k, x in r["by_source"].items()},
            "by_sign": {str(k): v(x) for k, x in r["by_sign"].items()},
            "tie_share_M1-sel": r["tie_share_M1-sel"],
            "n_significant": r["coverage_significant"]["n_significant"],
            "n_kinases_M1sel": r["M1-sel"].get("n_kinases")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="")
    ap.add_argument("--skip-jumpsem", action="store_true")
    ap.add_argument("--skip-bench", action="store_true")
    ap.add_argument("--jobs", type=int, default=32)
    ap.add_argument("--no-ledger", action="store_true")
    a = ap.parse_args()
    only = set(a.only.split(",")) - {""}
    t_start = time.time()
    warnings.filterwarnings("ignore", category=FutureWarning)
    ident = data_identity()
    print("data identity", ident, flush=True)

    P, C, U, _ = load_gt()
    bins = read_gt_tsv(GT_DIR / "harmony_bins_group.tsv")
    kin = list(U["kinase"])
    meta = load_dev_meta()
    ev = Evaluator(P, C, U, split_filter=FOLDS, bins=bins)
    missing = [c for c in ev.cond_index if c not in meta.index]
    assert not missing, f"{len(missing)} eligible GT conditions missing from DEV cache: {missing[:5]}"
    assert set(meta["split"]) <= set(FOLDS), set(meta["split"])
    X = load_dev()
    X = X.astype(np.float64)
    all_conds = list(X.columns)
    rng = np.random.default_rng(0)
    bench_conds = sorted(rng.choice(sorted(all_conds), 200, replace=False))
    Xb = X[bench_conds]
    SCORES.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    results_path = OUT / "baselines_dev.json"
    prev = json.loads(results_path.read_text()) if results_path.exists() and only else {}
    results = prev.get("runs", {}) if prev.get("identity", {}).get("gt_pairs.tsv") == ident["gt_pairs.tsv"] else {}

    def record(m, lib, cfg, kw, is_default, S, N, runtime, extra=None):
        name = f"{m}__{lib_label(m, lib)}__{cfg}"
        S.to_parquet(SCORES / f"{name}.parquet")
        if N is not None:
            N.astype(np.float32).to_parquet(SCORES / f"{name}.n.parquet")
        r = ev.evaluate(S)
        s = summarise(r)
        s["coverage_scored_kinase_cond"] = float(S.reindex(index=ev.cond_index).notna().to_numpy().mean())
        config = {"method": m, "library": lib_label(m, lib), "cfg": cfg,
                  "params": {k: kw[k] for k in kw if k != "n_jobs"}, "published_default_library": is_default,
                  "data_identity": {k: ident[k] for k in ["GT_FREEZE.json", "gt_pairs.tsv", "gt_conditions.tsv",
                                                          "quant_dev.npy"]}}
        if not a.no_ledger:
            ledger.append(m, config, "dev", to_jsonable({**s, "runtime": runtime, **(extra or {})}),
                          branch="baselines", folds=range(5))
        results[name] = {"method": m, "library": lib_label(m, lib), "cfg": cfg, "default": is_default,
                         "metrics": s, "runtime": runtime, **(extra or {})}
        print(f"{name}: M1-sel={s['M1-sel']:.4f} M1-all={s['M1-all']:.4f} cov={s['coverage_M1-sel']:.3f} "
              f"harm={s['harmony']:.4f} rt={runtime and runtime.get('ms_per_condition_1t')}", flush=True)
        (OUT / "baselines_dev.json").write_text(json.dumps(to_jsonable(
            {"identity": ident, "runs": results}), indent=1, default=str))

    bench_done = {}
    for m, lib, cfg, kw, is_def in grid():
        if only and m not in only:
            continue
        f = METHODS[m]
        kw2 = dict(kw)
        if m in SLOW:
            kw2["n_jobs"] = a.jobs
        t0 = time.perf_counter()
        S, N = f(X, lib, return_n=True, **kw2)
        t_full = time.perf_counter() - t0
        runtime = {"full_dev_wall_s": t_full, "n_conditions_full": X.shape[1],
                   "full_dev_jobs": a.jobs if m in SLOW else "threads(default)"}
        bkey = (m, str(lib), kw.get("min_substrates"))
        if not a.skip_bench and (m, str(lib)) not in bench_done:
            for nt in (1, 16):
                kwb = {k: v_ for k, v_ in kw.items()}
                if m in SLOW:
                    kwb["n_jobs"] = nt
                fn = (lambda Xx, ll, _f=f, _k=kwb: _f(Xx, ll, **_k))
                b = bench_runtime(fn, Xb, lib, n_rep=1 if m in SLOW else 5, threads=(nt,))
                r_ = b[f"threads_{nt}"]
                runtime[f"ms_per_condition_{nt}t"] = r_["ms_per_condition"]
                runtime[f"compendium_s_{nt}t"] = r_["compendium_extrapolated_s"]
                runtime["bench_n_rep"] = b["n_rep"]
            runtime["bench"] = "bench_runtime, 200 DEV conditions (seed 0), warm, median of n_rep; " + \
                ("process-parallel n_jobs = thread setting" if m in SLOW else "threadpool limits")
            bench_done[(m, str(lib))] = {k: v_ for k, v_ in runtime.items() if k.startswith(("ms_", "compendium", "bench"))}
        elif (m, str(lib)) in bench_done:
            runtime.update(bench_done[(m, str(lib))])
            runtime["bench_note"] = "same method+library benchmarked at another min_substrates"
        record(m, lib, cfg, kw, is_def, S, N, runtime)

    # ------------------------------------------------------------------ JUMPsem (multi-sample, per study)
    if not a.skip_jumpsem and (not only or "jumpsem" in only):
        ok = meta.index[meta["qc_status"] == "ok"]
        units = study_units(meta, ok)
        for lib in L3:
            r = run_jumpsem(X[list(ok)], lib, units=units, min_substrates=3, n_jobs=48,
                            log=lambda *x: None)
            r["units"].to_csv(OUT / f"jumpsem_units__{lib}.tsv", sep="\t", index=False)
            n_scored_c = int(r["native"].notna().any(axis=1).sum())
            runtime = {"full_dev_wall_s": r["wall_s"], "r_cpu_s_sum_units": r["r_cpu_s"],
                       "n_conditions_full": int(units.notna().sum()),
                       "ms_per_condition_1t": 1e3 * r["r_cpu_s"] / max(int(units.notna().sum()), 1),
                       "compendium_s_1t": r["r_cpu_s"] / max(int(units.notna().sum()), 1) * 1937,
                       "bench": "not bench_runtime (transductive per study): summed single-threaded Rscript wall "
                                "time over study units / # conditions in units; full run with 48 parallel units"}
            extra = {"n_units": int(len(r["units"])), "n_units_failed": int((r["units"]["error"] != "").sum()),
                     "n_conditions_in_units": int(units.notna().sum()), "n_conditions_with_scores": n_scored_c}
            for var, mname in [("native", "jumpsem"), ("aligned", "jumpsem_aligned")]:
                S = r[var].reindex(index=all_conds)
                record(mname, lib, "min3", {"min_substrates": 3, "variant": var}, lib == DEFAULT_LIB[mname],
                       S, r["n"].reindex(index=all_conds), runtime, extra)
    write_report(results, ident, time.time() - t_start)


def f3(x):
    return "NA" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.3f}"


def write_report(results, ident, elapsed):
    rows = sorted(results.items(), key=lambda kv: -(kv[1]["metrics"]["M1-sel"] or 0))
    L = ["# DEV baselines (folds 0-4 pooled)", "",
         "Generated by `scripts/run_baselines.py` (do not edit by hand). Data identity (sha256[:16]): "
         + ", ".join(f"{k} `{v}`" for k, v in ident.items()) + ".", "",
         "Evaluation: `kinact.eval.Evaluator` on DEV folds 0-4 (no devconfirm/test), harmony bins from "
         "`harmony_bins_group.tsv`; unscored = NaN (½-rule in M1). `cov` = M1-sel pair coverage (target scored). "
         "`n sig` = coverage-significant kinases (BH q ≤ 0.10). Runtime: `bench_runtime` on 200 fixed DEV conditions "
         "(seed 0), warm, median of 5 (1 for process-parallel methods), ms per condition at 1 / 16 threads "
         "(process-parallel methods: n_jobs = 1 / 16); JUMPsem: summed single-thread Rscript time per condition. "
         "`*` = published-default library of the method.", "",
         "| # | method | library | cfg | M1-sel | M1-all | M1-sel-scored | cov | harmony | H | D | O | sign − | sign + | n sig | ms/cond 1t | ms/cond 16t |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for i, (name, r) in enumerate(rows, 1):
        s, rt = r["metrics"], r["runtime"] or {}
        bs, bg = s["by_source"], s["by_sign"]
        L.append(f"| {i} | {r['method']} | {r['library']}{'*' if r['default'] else ''} | {r['cfg']} | "
                 f"{f3(s['M1-sel'])} | {f3(s['M1-all'])} | {f3(s['M1-sel-scored'])} | {f3(s['coverage_M1-sel'])} | "
                 f"{f3(s['harmony'])} | {f3(bs.get('H'))} | {f3(bs.get('D'))} | {f3(bs.get('O'))} | "
                 f"{f3(bg.get('-1'))} | {f3(bg.get('1'))} | {s['n_significant']} | "
                 f"{f3(rt.get('ms_per_condition_1t'))} | {f3(rt.get('ms_per_condition_16t'))} |")
    L += ["", "## Harmony per bin (M1-sel)", "", "| method | library | cfg | 0 | 1-4 | 5-19 | >=20 |", "|---|---|---|---|---|---|---|"]
    for name, r in rows:
        hb = r["metrics"]["harmony_per_bin"] or {}
        L.append(f"| {r['method']} | {r['library']} | {r['cfg']} | " + " | ".join(f3(hb.get(k)) for k in ["0", "1-4", "5-19", ">=20"]) + " |")
    L += ["", "## Best config per method (M1-sel) and published-default config", "",
          "| method | best library/cfg | best M1-sel | default library/cfg | default M1-sel |", "|---|---|---|---|---|"]
    by_m = {}
    for name, r in rows:
        by_m.setdefault(r["method"], []).append(r)
    for m, rs in by_m.items():
        b = rs[0]
        d = next((x for x in rs if x["default"]), None)
        L.append(f"| {m} | {b['library']} / {b['cfg']} | {f3(b['metrics']['M1-sel'])} | "
                 + (f"{d['library']} / {d['cfg']} | {f3(d['metrics']['M1-sel'])} |" if d else "– | – |"))
    js = [r for r in results.values() if r["method"].startswith("jumpsem")]
    L += ["", "## Notes", "",
          "- All single-sample methods score each condition on its own measured sites; a kinase is scored iff it has "
          "≥ min_substrates measured library substrates (default 3). Formulas and reference checks: "
          "`src/kinact/methods/*.py`, `results/baselines/reference_checks.md`.",
          "- decoupler methods (ulm/mlm/wmean/norm_wmean/zscore/ksea) are vectorised re-implementations matching "
          "decoupler 2.3.1 per condition (Spearman 1.0 except MC noise for norm_wmean/gsea; mlm differs only on "
          "rank-deficient designs where decoupler raises); viper calls decoupler's own `_func_viper` per condition.",
          "- RoKAI: published functional network (KS + STRING PPI + PTMcode co-evolution/structure distance) from "
          "github.com/serhan-yilmaz/RoKAI, propagation re-implemented from the MATLAB source (not run against MATLAB: "
          "no MATLAB licence); inference = RoKAI z-score on refined values. `rokai_ks` = RoKAI's own (2020 PSP) KS layer.",
          "- Motif-only: Cantley atlas library (`cache/libraries/cantley.parquet`). motif_zscore = RoKAI z over each "
          "kinase's top-k sites; phosx_like = per-site top-n kinase assignment + GSEA NES; kl_enrich = KL-style "
          "one-sided hypergeometric up-vs-unchanged / down-vs-unchanged at |log2FC| ≥ thr.",
          "- JUMPsem 1.1.0 run natively in R per study unit (PMID × cell line × control; units < 4 conditions pooled "
          "per PMID; otherwise NaN), complete-case sites within the unit, qc_status == ok conditions only. "
          "`jumpsem` = native output (sign set by lavaan's marker indicator, i.e. arbitrary); `jumpsem_aligned` = "
          "post-hoc sign alignment with the mean centred substrate log2FC (deviation from the package). "
          "Per-unit logs: `results/baselines/jumpsem_units__<lib>.tsv`."]
    for r in js[:1]:
        L.append(f"- JUMPsem coverage: {r.get('n_units')} units, {r.get('n_units_failed')} failed, "
                 f"{r.get('n_conditions_in_units')} conditions in units, {r.get('n_conditions_with_scores')} with ≥1 score.")
    L.append(f"- Total script wall time {elapsed / 60:.1f} min.")
    (OUT / "baselines_dev.md").write_text("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
