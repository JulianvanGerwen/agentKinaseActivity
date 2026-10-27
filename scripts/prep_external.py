"""Prepare the benchmarKIN external sets (TEST-only data; DESIGN §2/§4/§5) into cache/external/.

Reads the already position-mapped benchmarKIN conversions in external/benchmarkin_py (keys `{UniProt}_{res}{pos}`;
mapping by UniProt+position (+15-mer window for CPTAC), logs in *_mapping.json) and writes:
  cache/external/pert_{tyrosine,hernandez_new}.parquet   site x experiment (log2FC as distributed)
  cache/external/pert_targets.tsv       one row per (experiment, target): HGNC symbol, in_U, sign
  cache/external/pert_experiments.tsv   one row per experiment: dataset, pmid, drug, drug keys, drug_in_dev, n_sites
  cache/external/tumour_<COHORT>.parquet site x sample (median-centred, >=30 samples)
  cache/external/tumour_gs_long.tsv     cohort, sample, kinase (HGNC), label (1=GS+,0=GS-), sign (+1/-1), in_U
  cache/external/prep_report.json       sizes, mapping rates, overlap checks
No kinase-activity method is run here. Usage: python scripts/prep_external.py
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BP = ROOT / "external" / "benchmarkin_py"
OUT = ROOT / "cache" / "external"
GT = ROOT / "results" / "ground_truth"
COHORTS = ("BRCA", "CCRCC", "GBM", "HNSCC", "LSCC", "LUAD", "UCEC")
EXTRA_ALIAS = {"sorafinib": "sorafenib", "bi4834": "BI4834", "gsk1120212": "trametinib", "plx4032": "vemurafenib",
               "azd6244": "selumetinib", "gdc0941": "pictilisib", "sti571": "imatinib"}
STIMULI = {  # stimulus label -> regex over (normalised) experiment drug string / DEV Condition+Control text
    "EGF": r"(?<![a-z])egf(?!r)", "HRG": r"hrg|heregulin|nrg1", "VEGF": r"vegf", "TCR": r"cd3|tcr"}


def norm(s) -> str:
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def drug_keys(desc: str, aliases: pd.DataFrame) -> tuple[set, set]:
    """(canonical agent keys, stimulus labels) named in an experiment description."""
    nd = norm(desc)
    agents = {a for al, a in zip(aliases.alias, aliases.agent) if len(norm(al)) >= 4 and norm(al) in nd}
    agents |= {v for k, v in EXTRA_ALIAS.items() if k in nd}
    low = str(desc).lower()
    stim = {k for k, rx in STIMULI.items() if re.search(rx, low)}
    return agents, stim


def dev_drug_sets() -> tuple[set, set]:
    """Agents / stimuli present in any DEV condition (gt_conditions = DEV + DEV-confirm only; never sealed TEST)."""
    gt = pd.read_csv(GT / "gt_conditions.tsv", sep="\t", keep_default_na=False)
    assert "test" not in set(gt["split"]), "gt_conditions contains TEST rows"
    dev = set(gt.conditionID)
    doses = pd.read_csv(GT / "condition_doses.tsv", sep="\t", keep_default_na=False)
    al = pd.read_csv(GT / "drug_aliases.tsv", sep="\t", keep_default_na=False)
    amap = {norm(a): g for a, g in zip(al.alias, al.agent)}
    agents = {amap.get(norm(a), a) for a in doses[doses.conditionID.isin(dev)].agent}
    cond = pd.read_csv(ROOT / "data/01_phosphoproteomics_conditions.csv", dtype=str)
    cond = cond[cond.conditionID.isin(dev)]
    txt = (cond.Condition.fillna("") + " " + cond.Control.fillna("")).str.lower()
    stim = {k for k, rx in STIMULI.items() if txt.str.contains(rx, regex=True).any()}
    return agents, stim


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    U = pd.read_csv(GT / "kinase_universe.tsv", sep="\t", keep_default_na=False)
    acc2sym = dict(zip(U.uniprot, U.kinase))
    rep: dict = {}
    # ---------------------------------------------------------------- perturbation
    meta = pd.read_parquet(BP / "perturbation_meta.parquet")
    data = pd.read_parquet(BP / "perturbation_data.parquet")
    pmap = json.loads((BP / "perturbation_mapping.json").read_text())
    comp = pd.read_csv(ROOT / "data/01_phosphoproteomics_conditions.csv", dtype=str)
    comp_pmids = set(comp.PMID.dropna())
    ex = meta.drop_duplicates("experiment").copy()
    ex["pmid"] = ex.PMID.astype(str)
    ex["pmid_in_compendium"] = ex.pmid.isin(comp_pmids)
    sel = ex[(ex.dataset == "tyrosine") | ((ex.dataset == "hernandez") & ~ex.pmid_in_compendium)]
    assert (sel.dataset == "tyrosine").sum() == 18, "expected 18 tyrosine experiments"
    sets = {"tyrosine": sel[sel.dataset == "tyrosine"], "hernandez_new": sel[sel.dataset == "hernandez"]}
    rep["hernandez_new_pmids"] = sorted(sets["hernandez_new"].pmid.unique())
    rep["hijazi_pmid_in_compendium"] = "31959955" in comp_pmids
    al = pd.read_csv(GT / "drug_aliases.tsv", sep="\t", keep_default_na=False)
    dev_agents, dev_stim = dev_drug_sets()
    rep["dev_agent_count"], rep["dev_stimuli"] = len(dev_agents), sorted(dev_stim)
    exp_rows, tgt_rows = [], []
    for name, e in sets.items():
        d = data[list(e.experiment)].dropna(how="all")
        d.index.name = "site"
        d.to_parquet(OUT / f"pert_{name}.parquet")
        key_ok = d.index.str.fullmatch(r"[A-Z0-9]+(-\d+)?_[STY]\d+")
        rep[f"pert_{name}"] = {"experiments": int(d.shape[1]), "sites": int(d.shape[0]),
                               "median_sites_per_exp": float(d.notna().sum().median()),
                               "key_format_ok_frac": float(key_ok.mean()),
                               "pmids": sorted(e.pmid.unique())}
        src = "tyrosine" if name == "tyrosine" else "hernandez"
        rep[f"pert_{name}"]["source_mapping"] = {k: pmap[src][k] for k in ("source_rows", "mapped_rows", "mapped_frac")}
        for _, r in e.iterrows():
            ag, st = drug_keys(r.drug, al)
            seen = bool(ag & dev_agents) or bool(st & dev_stim)
            t = meta[meta.experiment == r.experiment]
            tt = t.assign(symbol=t.kinase_acc.map(acc2sym)).dropna(subset=["kinase_acc"])
            exp_rows.append(dict(experiment=r.experiment, set=name, pmid=r.pmid, drug=r.drug,
                                 time_min=r.Time_min, sign=int(r.sign), agent_keys=";".join(sorted(ag)),
                                 stimulus_keys=";".join(sorted(st)), drug_in_dev=seen,
                                 n_sites=int(d[r.experiment].notna().sum()),
                                 n_targets_total=int(len(t)), n_targets_in_U=int(tt.symbol.notna().sum())))
            for _, q in tt.iterrows():
                tgt_rows.append(dict(experiment=r.experiment, set=name, kinase_in=q.gene_in, kinase_acc=q.kinase_acc,
                                     kinase=q.symbol, in_U=pd.notna(q.symbol), sign=int(q.sign)))
    E_, T_ = pd.DataFrame(exp_rows), pd.DataFrame(tgt_rows)
    E_.to_csv(OUT / "pert_experiments.tsv", sep="\t", index=False)
    T_.to_csv(OUT / "pert_targets.tsv", sep="\t", index=False)
    rep["pert_targets"] = {"rows": len(T_), "in_U": int(T_.in_U.sum()),
                           "not_in_U": sorted(T_[~T_.in_U].kinase_in.unique())}
    rep["drug_in_dev"] = E_.groupby("set").drug_in_dev.agg(["sum", "size"]).astype(int).to_dict("index")
    # overlap of external site keys with our compendium site universe (informational)
    sk = set(pd.read_parquet(ROOT / "cache/sites.parquet", columns=["site_key"]).site_key)
    for name in sets:
        idx = pd.read_parquet(OUT / f"pert_{name}.parquet").index
        rep[f"pert_{name}"]["sites_in_compendium_universe_frac"] = round(float(idx.isin(sk).mean()), 4)
    # ---------------------------------------------------------------- tumour
    tmap = json.loads((BP / "tumour_mapping.json").read_text())
    rep["tumour"] = {}
    for c in COHORTS:
        d = pd.read_parquet(BP / f"tumour_{c}.parquet")
        d.index.name = "site"
        d.to_parquet(OUT / f"tumour_{c}.parquet")
        rep["tumour"][c] = {"sites": int(d.shape[0]), "samples": int(d.shape[1]),
                            "source_mapped_frac": tmap[c]["mapped_frac"]}
    gs = pd.read_parquet(BP / "tumour_gs.parquet")
    gs = gs.assign(kinase=gs.kinase_acc.map(acc2sym))
    gs["in_U"] = gs.kinase.notna()
    gs["sign"] = gs.label.map({1: 1, 0: -1})
    gl = gs.rename(columns={"gene_in": "kinase_in"})[["cohort", "sample", "kinase", "kinase_acc", "label", "sign",
                                                      "in_U"]]
    gl.to_csv(OUT / "tumour_gs_long.tsv", sep="\t", index=False)
    g = gl[gl.in_U]
    rep["tumour_gs"] = {"pairs_total": len(gl), "pairs_in_U": len(g), "kinases_in_U": int(g.kinase.nunique()),
                        "pos": int((g.label == 1).sum()), "neg": int((g.label == 0).sum()),
                        "note": "benchmarKIN GS = top/bottom 5% kinase PROTEIN level per cohort (not activating-site "
                                "based); label 1 -> +1, 0 -> -1"}
    (OUT / "prep_report.json").write_text(json.dumps(rep, indent=1, default=str))
    print(json.dumps(rep, indent=1, default=str))


if __name__ == "__main__":
    sys.exit(main())
