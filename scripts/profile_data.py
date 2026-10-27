#!/usr/bin/env python
"""Profile the phosphoproteomics compendium + PSP coverage. Writes docs/data_profile.md and
results/data_profile/*.tsv|json. Requires scripts/prep_data.py and scripts/build_candidate_gt.py to have run.
LEGACY: the candidate GT annotation sections read the retired pre-G1 draft files, now moved to
cache/sealed/draft_candidate_v0/ (sealed; not for development use). They are no longer in results/ground_truth/, so this script
fails at the annotation step unless re-pointed deliberately (not part of the pipeline). The kinact harness does not depend on this script."""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "cache"
OUT = ROOT / "results" / "data_profile"
GT = ROOT / "results" / "ground_truth"
DOC = ROOT / "docs" / "data_profile.md"
HIJAZI = "31959955"
THR = [1, 3, 5, 10, 20]


def q(x, ps=(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)):
    return {f"q{int(p * 100)}": round(float(np.nanquantile(x, p)), 3) for p in ps}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    conds = pd.read_parquet(CACHE / "conditions.parquet")
    sites = pd.read_parquet(CACHE / "sites.parquet")
    W = np.load(CACHE / "quant_wide_f32.npy", mmap_mode="r")
    W = np.asarray(W)
    obs = ~np.isnan(W)
    rep = {}
    nsite_cond = obs.sum(0)
    ncond_site = obs.sum(1)
    rep["n_conditions"] = int(len(conds))
    rep["n_studies"] = int(conds["PMID"].nunique())
    rep["n_samples_celltypes_raw"] = int(conds["Sample"].nunique())
    rep["n_sites"] = int(len(sites))
    rep["n_values"] = int(obs.sum())
    rep["density"] = round(float(obs.mean()), 4)
    rep["sites_per_condition"] = q(nsite_cond) | {"min": int(nsite_cond.min()), "max": int(nsite_cond.max()),
                                                   "mean": round(float(nsite_cond.mean()), 1)}
    rep["conditions_per_site"] = q(ncond_site) | {"mean": round(float(ncond_site.mean()), 1),
                                                   "sites_in_1_condition": int((ncond_site == 1).sum()),
                                                   "sites_in_>=100": int((ncond_site >= 100).sum()),
                                                   "sites_in_>=50pct": int((ncond_site >= 0.5 * len(conds)).sum())}
    rep["source_counts"] = conds["source"].value_counts().to_dict()
    rep["labeling"] = conds["Labeling.method"].value_counts().head(8).to_dict()
    vals = W[obs]
    rep["log2FC_overall"] = q(vals) | {"mean": round(float(vals.mean()), 4), "sd": round(float(vals.std()), 3),
                                       "min": round(float(vals.min()), 2), "max": round(float(vals.max()), 2),
                                       "frac_exact_zero": round(float((vals == 0).mean()), 4),
                                       "frac_abs_gt_1": round(float((np.abs(vals) > 1).mean()), 4)}
    med = np.nanmedian(W, axis=0)
    sd = np.nanstd(W, axis=0)
    rep["per_condition_median"] = q(med) | {"frac_abs_median_lt_0.1": round(float((np.abs(med) < 0.1).mean()), 3),
                                            "frac_abs_median_gt_0.5": round(float((np.abs(med) > 0.5).mean()), 3)}
    rep["per_condition_sd"] = q(sd)
    pc = conds[["conditionID", "PMID", "Condition", "Control", "Sample", "Timepoint", "source"]].copy()
    pc["n_sites"] = nsite_cond
    pc["median"] = med.round(4)
    pc["sd"] = sd.round(4)
    pc["frac_zero"] = ((W == 0).sum(0) / np.maximum(nsite_cond, 1)).round(4)
    # study-level
    st = pc.groupby("PMID").agg(n_conditions=("conditionID", "size"), source=("source", "first"),
                                n_cells=("Sample", "nunique"), median_sites=("n_sites", "median"),
                                median_of_medians=("median", "median"), median_sd=("sd", "median"))
    rep["studies_by_source"] = st["source"].value_counts().to_dict()
    rep["dup_condition_names_within_study"] = int(conds.duplicated(["PMID", "Condition", "Control", "Sample", "Timepoint"]).sum())
    # GT annotation
    gtc = pd.read_csv(GT / "candidate_gt_conditions.tsv", sep="\t", dtype=str).fillna("")
    pc = pc.merge(gtc[["conditionID", "perturbation_class", "kinases_direct", "kinases_downstream",
                       "expected_direction", "best_confidence"]], on="conditionID", how="left").fillna("")
    # QC flags (heuristic): SD > 3 suggests values are not log2 ratios or contain extreme outliers;
    # |median| > 0.5 suggests no centring; < 2500 sites = small condition
    pc["max_abs"] = np.nanmax(np.abs(W), axis=0).round(2)
    pc["qc_flag"] = ""
    pc.loc[pc["sd"] > 3, "qc_flag"] += "sd>3;"
    pc.loc[pc["max_abs"] > 20, "qc_flag"] += "max|log2FC|>20;"
    pc.loc[pc["median"].abs() > 0.5, "qc_flag"] += "|median|>0.5;"
    pc.loc[pc["n_sites"] < 2500, "qc_flag"] += "n_sites<2500;"
    rep["qc_flag_counts"] = {"sd>3": int((pc["sd"] > 3).sum()), "max|log2FC|>20": int((pc["max_abs"] > 20).sum()),
                             "|median|>0.5": int((pc["median"].abs() > 0.5).sum()),
                             "n_sites<2500": int((pc["n_sites"] < 2500).sum()),
                             "any": int((pc["qc_flag"] != "").sum()),
                             "gt_conditions_flagged_sd_or_max": sorted(pc.loc[(pc["best_confidence"] != "") & pc["qc_flag"].str.contains("sd>3|max"), "conditionID"])}
    pc.to_csv(OUT / "conditions_profile.tsv", sep="\t", index=False)
    st.to_csv(OUT / "studies_profile.tsv", sep="\t")

    # ---- overlap check: Wilkes 2015 (26060313, Cutillas lab) MCF7 vs Hijazi MCF7 for shared drugs
    cidx = pd.Series(conds.index, index=conds["conditionID"])
    pairs = [("26060313_Ochoa_1289_272_log2", "31959955_179_log2", "MK-2206"),
             ("26060313_Ochoa_1299_272_log2", "31959955_121_log2", "Torin1"),
             ("26060313_Ochoa_1300_272_log2", "31959955_13_log2", "PF4708671"),
             ("26060313_Ochoa_1302_272_log2", "31959955_162_log2", "GDC0941"),
             ("26060313_Ochoa_1298_272_log2", "31959955_118_log2", "KU0063794"),
             ("26060313_Ochoa_1293_272_log2", "31959955_119_log2", "PD153035"),
             ("26060313_Ochoa_1290_272_log2", "31959955_178_log2", "KN93")]
    ov = []
    rng = np.random.default_rng(0)
    for a, b, d in pairs:
        if a not in cidx or b not in cidx:
            continue
        x, y = W[:, cidx[a]], W[:, cidx[b]]
        m = ~np.isnan(x) & ~np.isnan(y)
        r = float(np.corrcoef(x[m], y[m])[0, 1]) if m.sum() > 10 else np.nan
        # reference: Wilkes condition vs a random other Hijazi MCF7 condition
        refs = [c for c in conds.loc[conds["Condition"].astype(str).str.startswith("MCF7_"), "conditionID"] if c != b]
        rr = []
        for c in rng.choice(refs, 10, replace=False):
            z = W[:, cidx[c]]
            mm = ~np.isnan(x) & ~np.isnan(z)
            rr.append(np.corrcoef(x[mm], z[mm])[0, 1])
        ov.append({"drug": d, "wilkes": a, "hijazi": b, "n_shared_sites": int(m.sum()), "pearson": round(r, 3),
                   "ref_median_pearson_other_hijazi_mcf7": round(float(np.median(rr)), 3),
                   "identical_values_frac": round(float((x[m] == y[m]).mean()), 3) if m.any() else np.nan})
    rep["wilkes_vs_hijazi_mcf7"] = ov

    # ---- PSP coverage
    ks = pd.read_parquet(CACHE / "psp_ks_human.parquet")
    ks = ks[ks["in_compendium"]]
    sidx = pd.Series(sites["site_idx"].values, index=sites["peptide_UID"].values)
    ks = ks.assign(site_idx=sidx.reindex(ks["peptide_UID"].astype(str)).values)
    kss = ks.drop_duplicates(["kinase", "site_idx"])[["kinase", "site_idx"]]
    hz_cols = conds.index[conds["PMID"] == HIJAZI].values
    any_all = ncond_site > 0
    any_hz = obs[:, hz_cols].any(1)
    per_cond_hz = obs[:, hz_cols]
    cov = []
    allks = pd.read_parquet(CACHE / "psp_ks_human.parquet")
    n_psp = allks.drop_duplicates(["kinase", "hs_acc", "hs_residue", "hs_pos"]).groupby("kinase").size()
    for k, g in kss.groupby("kinase"):
        si = g["site_idx"].values.astype(int)
        cov.append({"kinase": k, "n_psp_sites": int(n_psp.get(k, 0)), "n_sites_compendium": int(any_all[si].sum()),
                    "n_sites_hijazi_union": int(any_hz[si].sum()),
                    "median_sites_per_hijazi_condition": float(np.median(per_cond_hz[si].sum(0))),
                    "median_sites_per_condition_all": float(np.median(obs[si].sum(0)))})
    cov = pd.DataFrame(cov)
    missing = sorted(set(n_psp.index) - set(cov["kinase"]))
    cov = pd.concat([cov, pd.DataFrame({"kinase": missing, "n_psp_sites": [int(n_psp[m]) for m in missing]})]).fillna(0)
    cov = cov.sort_values("n_sites_compendium", ascending=False)
    hzt = pd.read_csv(GT / "hijazi_conditions.tsv", sep="\t", keep_default_na=False)
    hz_targets = sorted(set(";".join(hzt["targets"]).split(";")) - {""})
    cov["is_hijazi_target"] = cov["kinase"].isin(hz_targets)
    cov.to_csv(OUT / "psp_kinase_coverage.tsv", sep="\t", index=False)
    rep["psp_n_human_kinases"] = int(len(n_psp))
    rep["psp_kinases_by_n_substrate_sites_quantified"] = {
        col: {f">={t}": int((cov[col] >= t).sum()) for t in THR}
        for col in ["n_psp_sites", "n_sites_compendium", "n_sites_hijazi_union", "median_sites_per_hijazi_condition",
                    "median_sites_per_condition_all"]}
    ht = pd.DataFrame({"kinase": hz_targets}).merge(cov, on="kinase", how="left").fillna(0)
    rep["hijazi_targets"] = {"n": len(hz_targets),
                             "not_in_psp_as_kinase": int((ht["n_psp_sites"] == 0).sum()),
                             **{f"union_>={t}": int((ht["n_sites_hijazi_union"] >= t).sum()) for t in THR},
                             **{f"median_per_condition_>={t}": int((ht["median_sites_per_hijazi_condition"] >= t).sum()) for t in THR}}
    # evaluable (condition, target) pairs in Hijazi: >=3 PSP substrates quantified in that condition
    pair_rows = []
    cid2col = pd.Series(conds.index, index=conds["conditionID"])
    ksite = kss.groupby("kinase")["site_idx"].apply(lambda s: s.values.astype(int)).to_dict()
    gt = pd.read_csv(GT / "candidate_gt.tsv", sep="\t", dtype=str).fillna("")
    gt = gt[(gt["expected_direction"] != "0") & (gt["relation"] != "")]
    for _, r in gt.iterrows():
        si = ksite.get(r["kinase"])
        n = int(obs[si, cid2col[r["conditionID"]]].sum()) if si is not None else 0
        pair_rows.append(n)
    gt["n_psp_sites_in_condition"] = pair_rows
    gt[["conditionID", "kinase", "relation", "confidence", "PMID", "n_psp_sites_in_condition"]].to_csv(
        OUT / "gt_pair_psp_coverage.tsv", sep="\t", index=False)
    ev = {}
    for name, sub in [("hijazi", gt[gt["PMID"] == HIJAZI]), ("non_hijazi", gt[gt["PMID"] != HIJAZI]),
                      ("non_hijazi_high", gt[(gt["PMID"] != HIJAZI) & (gt["confidence"] == "high")])]:
        ev[name] = {"pairs": int(len(sub)), "kinases": int(sub["kinase"].nunique()),
                    **{f"pairs_>={t}": int((sub["n_psp_sites_in_condition"] >= t).sum()) for t in (1, 3, 5)},
                    **{f"kinases_>=3_in_some_condition": int(sub.loc[sub["n_psp_sites_in_condition"] >= 3, "kinase"].nunique())},
                    **{f"conditions_with_>=1_pair_>=3": int(sub.loc[sub["n_psp_sites_in_condition"] >= 3, "conditionID"].nunique())}}
    rep["gt_pairs_evaluable"] = ev
    top = cov.head(10)[["kinase", "n_sites_compendium"]].values.tolist()
    rep["top10_kinases_by_sites"] = top
    tot = cov["n_sites_compendium"].sum()
    rep["share_of_kinase_site_annotations_top10_pct"] = round(100 * cov.head(10)["n_sites_compendium"].sum() / tot, 1)
    (OUT / "profile_summary.json").write_text(json.dumps(rep, indent=1, default=str))
    write_md(rep, pc, st, conds)
    print(json.dumps(rep, indent=1, default=str)[:6000])


def write_md(rep, pc, st, conds):
    L = ["# Data profile — phosphoproteomics compendium", "",
         "Generated by `scripts/profile_data.py` (inputs: caches from `scripts/prep_data.py`, GT from "
         "`scripts/build_candidate_gt.py`). Machine-readable numbers: `results/data_profile/profile_summary.json`; "
         "per-condition table: `results/data_profile/conditions_profile.tsv`; per-study: `studies_profile.tsv`; "
         "site mapping: `site_mapping_report.json`; PSP coverage: `psp_kinase_coverage.tsv`.", "",
         "## Overview", ""]
    L += [f"- Conditions: **{rep['n_conditions']}**; studies (PMIDs): **{rep['n_studies']}**; distinct `Sample` "
          f"strings (cell lines/tissues, not harmonised): **{rep['n_samples_celltypes_raw']}**.",
          f"- Source column: {rep['source_counts']} (studies: {rep['studies_by_source']}).",
          f"- Sites: **{rep['n_sites']}** (all `_M1`), values: **{rep['n_values']:,}**, matrix density "
          f"**{100 * rep['density']:.2f}%** (95.4% missing).",
          f"- Sites per condition: {rep['sites_per_condition']}",
          f"- Conditions per site: {rep['conditions_per_site']}",
          f"- Labeling: {rep['labeling']}", "",
          "## log2FC distribution / normalisation", "",
          f"- Overall: {rep['log2FC_overall']}",
          f"- Per-condition median: {rep['per_condition_median']}",
          f"- Per-condition SD: {rep['per_condition_sd']}",
          "- Interpretation: see *Findings and caveats* at the end (from docs/data_profile_notes.md).", ""]
    L += ["## Overlap check: Wilkes 2015 (26060313, Ochoa2016 source) vs Hijazi MCF7", "",
          "| drug | n shared sites | Pearson | ref. median r (other Hijazi MCF7 conds) | identical values |",
          "|---|---|---|---|---|"]
    for o in rep["wilkes_vs_hijazi_mcf7"]:
        L.append(f"| {o['drug']} | {o['n_shared_sites']} | {o['pearson']} | {o['ref_median_pearson_other_hijazi_mcf7']} | {o['identical_values_frac']} |")
    L += ["", "## PSP coverage", "",
          f"Human kinases in PSP KS (human kinase; human or SITE_GRP-mapped substrate): {rep['psp_n_human_kinases']}.", "",
          "| measure | " + " | ".join(f"≥{t}" for t in THR) + " |", "|---|" + "---|" * len(THR)]
    for k, v in rep["psp_kinases_by_n_substrate_sites_quantified"].items():
        L.append(f"| {k} | " + " | ".join(str(v[f'>={t}']) for t in THR) + " |")
    L += ["", f"Hijazi target kinases (kinhib matrix, HGNC-normalised): {rep['hijazi_targets']}", "",
          f"GT (condition, kinase) pairs with PSP substrates quantified in that condition: {rep['gt_pairs_evaluable']}", "",
          f"Top-10 kinases by quantified substrate sites: {rep['top10_kinases_by_sites']} "
          f"({rep['share_of_kinase_site_annotations_top10_pct']}% of all kinase–site annotations in the compendium).", ""]
    L += ["## Studies and conditions (grouped by perturbation class from candidate GT; 'unannotated' = no kinase GT)", ""]
    for pmid, g in pc.groupby("PMID", sort=False):
        s = st.loc[pmid]
        cells = ", ".join(map(str, g["Sample"].unique()[:8])) + (" …" if g["Sample"].nunique() > 8 else "")
        tps = ", ".join(map(str, g["Timepoint"].dropna().unique()[:8]))
        L.append(f"### PMID {pmid} — {len(g)} conditions ({s['source']}); samples: {cells}; time: {tps}")
        g = g.assign(cls=g["perturbation_class"].replace("", "unannotated"))
        for cls, gg in g.groupby("cls"):
            names = pd.unique(gg["Condition"].astype(str) + " vs " + gg["Control"].astype(str))
            if pmid == "36926954":
                names = pd.unique(gg["Condition"].astype(str).str.rsplit("_", n=1).str[0] + " [" + gg["Sample"].astype(str) + "]")
            txt = "; ".join(names[:40]) + (f"; … (+{len(names) - 40})" if len(names) > 40 else "")
            L.append(f"- **{cls}** ({len(gg)}): {txt}")
        L.append("")
    notes = ROOT / "docs" / "data_profile_notes.md"
    if notes.exists():
        L += ["", notes.read_text()]
    DOC.write_text("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
