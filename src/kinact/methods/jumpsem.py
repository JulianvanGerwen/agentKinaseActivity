"""JUMPsem 1.1.0 (Kong et al. 2025; R package, external/src/JUMPsem) run natively per study via Rscript.

JUMPsem is multi-sample: per study it row-centres each site across the study's conditions and fits, per kinase,
a one-factor CFA (lavaan) over the conditions; activity = lavPredict factor scores (marker-indicator scale,
returned unstandardised). Wrapping choices (documented in results/baselines/baselines_dev.md):
  * study unit = (PMID, Sample, Control) experiment; units with < MIN_CONDITIONS conditions are pooled per PMID;
    if still < MIN_CONDITIONS conditions the conditions stay unscored (NaN). With n samples JUMPsem keeps at most
    n-1 indicators, and a one-factor CFA needs >= 3 indicators -> MIN_CONDITIONS = 4.
  * only sites measured in every condition of the unit are passed (JUMPsem's rowMeans/lavaan do not handle NA;
    a NA row becomes all-NA after centring);
  * input = log2FC, input.log2.norm = FALSE, other arguments default (kmo.off=0, relative.norm.p=TRUE, mdsite=TRUE);
  * database = the chosen library written in PSP format (GENE, ..., SUB_ACC_ID, SUB_MOD_RSD; organisms 'human'),
    restricted to the unit's input sites (does not change results);
  * a kinase is reported only if it has >= min_substrates unit-complete substrates (JUMPsem itself needs >1, and
    lavaan cannot identify a 1-factor model with 2 indicators);
  * 'native' = JUMPsem output as is. Its sign is set by lavaan's marker indicator (first substrate), so it is
    arbitrary w.r.t. activation; 'aligned' additionally flips each (unit, kinase) score vector so that it
    correlates positively across the unit's conditions with the mean centred log2FC of the kinase's substrates
    (a post-hoc deviation, reported separately).
"""
from __future__ import annotations

import os
import subprocess
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from ._common import DEFAULT_MIN_SUB, load_lib_df

ROOT = Path(__file__).resolve().parents[3]
RSCRIPT = ROOT / "external/envs/r/bin/Rscript"
RUNNER = Path(__file__).with_name("jumpsem_run.R")
MIN_CONDITIONS = 4


def study_units(meta: pd.DataFrame, conditions, min_conditions=MIN_CONDITIONS) -> pd.Series:
    """conditionID -> unit label (NaN if the condition cannot be put in a unit of >= min_conditions)."""
    m = meta.loc[list(conditions)]
    key = m["PMID"].astype(str) + "|" + m["Sample"].astype(str) + "|" + m["Control"].astype(str)
    size = key.map(key.value_counts())
    unit = key.where(size >= min_conditions, m["PMID"].astype(str) + "|pooled")
    usize = unit.map(unit.value_counts())
    return unit.where(usize >= min_conditions)


def _split_key(keys):
    s = pd.Series(np.asarray(keys, dtype=object)).astype(str)
    return s.str.rsplit("_", n=1).str[0].to_numpy(), s.str.rsplit("_", n=1).str[1].to_numpy()


def _run_unit(args):
    uid, Xu, lib, min_sub, tmpdir = args
    Xc = Xu.dropna(how="any")
    if Xc.shape[0] < 3:
        return uid, None, None, 0.0, "too few complete sites"
    libu = lib[lib["site"].isin(Xc.index)]
    if libu.empty:
        return uid, None, None, 0.0, "no library sites"
    acc, res = _split_key(Xc.index)
    samples = [f"S{i + 1}" for i in range(Xc.shape[1])]
    inp = pd.DataFrame({"Substrate_ACC": acc, "Substrate_GN": acc, "sites": res, "Peptide": "X"})
    inp = pd.concat([inp, pd.DataFrame(Xc.to_numpy(np.float64), columns=samples)], axis=1)
    la, lr = _split_key(libu["site"])
    db = pd.DataFrame({"GENE": libu["kinase"].to_numpy(), "KINASE": libu["kinase"].to_numpy(), "KIN_ACC_ID": "X",
                       "KIN_ORGANISM": "human", "SUBSTRATE": "s", "SUB_GENE_ID": 1, "SUB_ACC_ID": la,
                       "SUB_GENE": la, "SUB_ORGANISM": "human", "SUB_MOD_RSD": lr})
    for i in range(11, 17):
        db[f"V{i}"] = ""
    d = Path(tempfile.mkdtemp(dir=tmpdir))
    inp.to_csv(d / "in.csv", index=False)
    db.to_csv(d / "db.csv", index=False)
    env = dict(os.environ, LC_ALL="C.UTF-8", PATH=f"{RSCRIPT.parent}:{os.environ.get('PATH', '')}",
               OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    t0 = time.perf_counter()
    p = subprocess.run([str(RSCRIPT), str(RUNNER), str(d / "in.csv"), str(d / "db.csv"), str(d / "act.tsv")],
                       capture_output=True, text=True, env=env, timeout=7200)
    wall = time.perf_counter() - t0
    out = d / "act.tsv"
    if p.returncode != 0 or not out.exists() or out.read_text().startswith("#error"):
        return uid, None, None, wall, f"R failed rc={p.returncode}: {p.stderr[-300:]}"
    act = pd.read_csv(out, sep="\t", keep_default_na=False, na_values=["NA", "NaN", ""])
    act = act.set_index("ENZYME")[samples]
    act.columns = list(Xu.columns)
    # measured-substrate counts and sign alignment reference (mean centred substrate log2FC)
    Xcc = Xc.sub(Xc.mean(axis=1), axis=0)
    g = libu.groupby("kinase")["site"].apply(lambda s: list(pd.unique(s)))
    n = g.apply(len)
    ref = pd.DataFrame({k: Xcc.loc[v].mean(axis=0) for k, v in g.items() if k in act.index}).T
    return uid, act, (n, ref), wall, ""


def run_jumpsem(X: pd.DataFrame, lib, units: pd.Series | None = None, min_substrates=DEFAULT_MIN_SUB,
                n_jobs=None, tmpdir=None, log=print):
    """Returns dict(native, aligned, n) conditions x kinases frames + per-unit log DataFrame."""
    libdf = load_lib_df(lib)[["kinase", "site"]].drop_duplicates()
    if units is None:
        units = pd.Series("all", index=X.columns)
    units = units.reindex(X.columns)
    kin = sorted(libdf["kinase"].unique())
    native = pd.DataFrame(np.nan, index=X.columns, columns=kin)
    aligned = native.copy(); nfr = native.copy()
    own_tmp = tmpdir is None
    tmpdir = tmpdir or tempfile.mkdtemp(prefix="jumpsem_")
    jobs = [(u, X[list(cs)], libdf, min_substrates, tmpdir)
            for u, cs in units.dropna().groupby(units.dropna()).groups.items()]
    jobs.sort(key=lambda j: -j[1].shape[1] * j[1].shape[0])
    rows = []
    n_jobs = n_jobs or min(48, os.cpu_count())
    t0 = time.perf_counter()
    with ThreadPoolExecutor(n_jobs) as ex:
        for uid, act, extra, wall, err in ex.map(_run_unit, jobs):
            cs = list(units.index[units == uid])
            rows.append({"unit": uid, "n_conditions": len(cs), "r_wall_s": wall, "error": err,
                         "n_kinases": 0 if act is None else int(act.shape[0])})
            if act is None:
                continue
            n, ref = extra
            act = act[act.index.isin(n.index[n >= min_substrates])]
            act = act[act.index.isin(native.columns)]
            native.loc[cs, act.index] = act.T.to_numpy()
            nfr.loc[cs, act.index] = np.broadcast_to(n.reindex(act.index).to_numpy(), (len(cs), len(act)))
            sgn = []
            for k in act.index:
                a, r = act.loc[k].to_numpy(float), ref.loc[k].to_numpy(float)
                c = np.corrcoef(a, r)[0, 1] if np.std(a) > 0 and np.std(r) > 0 else 1.0
                sgn.append(-1.0 if c < 0 else 1.0)
            aligned.loc[cs, act.index] = (act.T * np.array(sgn)).to_numpy()
    total = time.perf_counter() - t0
    if own_tmp:
        import shutil
        shutil.rmtree(tmpdir, ignore_errors=True)
    log_df = pd.DataFrame(rows)
    return {"native": native, "aligned": aligned, "n": nfr, "units": log_df, "wall_s": total,
            "r_cpu_s": float(log_df["r_wall_s"].sum()) if len(log_df) else 0.0}


def jumpsem(X, lib, min_substrates=DEFAULT_MIN_SUB, return_n=False, units=None, variant="native", **kw):
    r = run_jumpsem(X, lib, units=units, min_substrates=min_substrates, **kw)
    return (r[variant], r["n"]) if return_n else r[variant]
