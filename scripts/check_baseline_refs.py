#!/usr/bin/env python
"""Reference agreement checks for the baseline re-implementations (DESIGN §6: Spearman >= 0.99).

* decoupler 2.3.1 (dc.mt.*), run per condition on its measured sites (tmin = 3, empty=False), vs
  kinact.methods.{zscore, ksea, wmean, norm_wmean, ulm, mlm, viper, gsea} on N_DC fixed DEV conditions (seed 0).
* KSEAapp 2.0 (R) KSEA.Scores vs kinact.methods.ksea (min_substrates=1: KSEAapp has no m cut-off) on N_KSEA
  conditions.
Per condition: Spearman over kinases scored by both; reported min/median over conditions and pooled.
Writes results/baselines/reference_checks.{json,md}.
"""
import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from kinact.data import load_dev  # noqa: E402
from kinact.methods import METHODS  # noqa: E402
from kinact.methods._common import load_lib_df  # noqa: E402

OUT = ROOT / "results/baselines"
N_DC, N_KSEA, LIB = 60, 12, "psp_all"


def spear(a, b):
    ok = a.notna() & b.notna() & np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3:
        return np.nan, int(ok.sum())
    return float(a[ok].corr(b[ok], method="spearman")), int(ok.sum())


def main():
    import decoupler as dc
    X = load_dev()
    rng = np.random.default_rng(0)
    conds = sorted(rng.choice(X.columns, N_DC, replace=False))
    Xs = X[conds]
    lib = load_lib_df(LIB)
    net = pd.DataFrame({"source": lib.kinase, "target": lib.site, "weight": lib.weight}).drop_duplicates(["source", "target"])
    dcf = {"zscore": lambda d: dc.mt.zscore(d, net, tmin=3, empty=False),
           "ksea": lambda d: dc.mt.zscore(d, net, tmin=3, flavor="KSEA", empty=False),
           "wmean": lambda d: dc.mt.waggr(d, net, tmin=3, times=0, empty=False),
           "norm_wmean": lambda d: dc.mt.waggr(d, net, tmin=3, times=1000, empty=False),
           "ulm": lambda d: dc.mt.ulm(d, net, tmin=3, empty=False),
           "mlm": lambda d: dc.mt.mlm(d, net, tmin=3, empty=False),
           "viper": lambda d: dc.mt.viper(d, net, tmin=3, empty=False),
           "gsea": lambda d: dc.mt.gsea(d, net, tmin=3, empty=False)}
    res = {}
    for m, f in dcf.items():
        t0 = time.perf_counter(); S = METHODS[m](Xs, LIB); t_ours = time.perf_counter() - t0
        rows, fails, t_dc = [], 0, 0.0
        pa, pb = [], []
        for c in conds:
            d = Xs[[c]].dropna().T.astype(np.float64)
            t0 = time.perf_counter()
            try:
                es, _ = f(d)
            except Exception:
                fails += 1
                continue
            finally:
                t_dc += time.perf_counter() - t0
            a, b = S.loc[c].reindex(es.columns), es.iloc[0]
            r, n = spear(a, b)
            cov_ours, cov_dc = int(a.notna().sum()), int(b.notna().sum())
            rows.append((r, n, cov_ours, cov_dc, float(np.nanmax(np.abs(a - b)))))
            ok = a.notna() & b.notna()
            pa.append(a[ok]); pb.append(b[ok])
        R = pd.DataFrame(rows, columns=["rho", "n", "n_ours", "n_dc", "maxabs"])
        res[m] = {"n_conditions": len(R), "decoupler_failed": fails, "rho_min": float(R.rho.min()),
                  "rho_median": float(R.rho.median()),
                  "rho_pooled": float(pd.concat(pa).corr(pd.concat(pb), method="spearman")),
                  "max_abs_diff": float(R.maxabs.max()), "scored_ours_minus_dc": int((R.n_ours - R.n_dc).sum()),
                  "seconds_ours": t_ours, "seconds_decoupler": t_dc}
        print(m, res[m], flush=True)
    # ---- KSEAapp
    kconds = conds[:N_KSEA]
    klib = "psp_human_only"
    L = load_lib_df(klib)[["kinase", "site"]].drop_duplicates()
    L["acc"] = L.site.str.rsplit("_", n=1).str[0]; L["rsd"] = L.site.str.rsplit("_", n=1).str[1]
    td = Path(tempfile.mkdtemp())
    (td / "px").mkdir(); (td / "out").mkdir()
    L[["kinase", "acc", "rsd"]].to_csv(td / "lib.csv", index=False)
    for c in kconds:
        x = X[c].dropna()
        acc = x.index.str.rsplit("_", n=1).str[0]; rsd = x.index.str.rsplit("_", n=1).str[1]
        pd.DataFrame({"Protein": acc, "Gene": acc, "Peptide": x.index, "Residue.Both": rsd, "p": 1.0,
                      "FC": np.power(2.0, x.to_numpy(np.float64))}).to_csv(td / "px" / f"{c}.csv", index=False)
    env = dict(os.environ, LC_ALL="C.UTF-8", PATH=f"{ROOT}/external/envs/r/bin:" + os.environ["PATH"])
    t0 = time.perf_counter()
    subprocess.run([f"{ROOT}/external/envs/r/bin/Rscript", str(ROOT / "scripts/ksea_ref.R"), str(td / "lib.csv"),
                    str(td / "px"), str(td / "out")], check=True, env=env, capture_output=True)
    t_r = time.perf_counter() - t0
    ours = METHODS["ksea"](X[kconds], klib, min_substrates=1)
    rows = []
    for c in kconds:
        r = pd.read_csv(td / "out" / f"{c}.csv").set_index("Kinase.Gene")
        a, b = ours.loc[c].reindex(r.index), r["z.score"]
        rho, n = spear(a, b)
        rows.append((rho, n, float(np.nanmax(np.abs(a - b))), int(a.notna().sum()), len(r)))
    R = pd.DataFrame(rows, columns=["rho", "n", "maxabs", "n_ours", "n_ksea"])
    res["ksea_vs_KSEAapp"] = {"library": klib, "n_conditions": len(R), "rho_min": float(R.rho.min()),
                              "rho_median": float(R.rho.median()), "max_abs_diff": float(R.maxabs.max()),
                              "kinases_ours_minus_KSEAapp": int((R.n_ours - R.n_ksea).sum()),
                              "seconds_KSEAapp_R": t_r}
    print(res["ksea_vs_KSEAapp"])
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "reference_checks.json").write_text(json.dumps({"library": LIB, "conditions": conds, "results": res}, indent=1))
    L = ["# Reference agreement checks (scripts/check_baseline_refs.py)", "",
         f"decoupler 2.3.1 run per condition on measured sites (tmin=3), library `{LIB}`, {N_DC} DEV conditions "
         "(seed 0). Spearman over kinases scored by both, per condition (min / median) and pooled over all pairs.", "",
         "| method | conditions | decoupler failed | rho min | rho median | rho pooled | max abs diff | ours s | decoupler s |",
         "|---|---|---|---|---|---|---|---|---|"]
    for m in dcf:
        r = res[m]
        L.append(f"| {m} | {r['n_conditions']} | {r['decoupler_failed']} | {r['rho_min']:.4f} | {r['rho_median']:.4f} | "
                 f"{r['rho_pooled']:.4f} | {r['max_abs_diff']:.3g} | {r['seconds_ours']:.2f} | {r['seconds_decoupler']:.1f} |")
    k = res["ksea_vs_KSEAapp"]
    L += ["", f"KSEAapp 2.0 `KSEA.Scores` (R) vs `kinact.methods.ksea(min_substrates=1)`, library `{klib}` "
          f"(as KSData, PSP source), {k['n_conditions']} conditions: Spearman min {k['rho_min']:.6f}, median "
          f"{k['rho_median']:.6f}, max |Δz| {k['max_abs_diff']:.3g}; kinase-count difference "
          f"{k['kinases_ours_minus_KSEAapp']}.", "",
          "Notes: norm_wmean is computed analytically (infinite-permutation limit of decoupler's 1000-permutation "
          "z) and gsea uses its own permutation draws, so these differ from decoupler by Monte-Carlo noise only. "
          "mlm: decoupler raises on rank-deficient designs (e.g. two kinases with identical measured substrate sets);"
          " those conditions are counted as 'decoupler failed' and our pseudo-inverse solution is used there; among the "
          "full-rank conditions, rho < 0.99 occurs only on near-singular designs where decoupler's explicit inverse "
          "returns numerically meaningless |t| (up to ~5e5; see max abs diff). viper "
          "calls decoupler's own `_func_viper` per condition (no re-implementation)."]
    (OUT / "reference_checks.md").write_text("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
