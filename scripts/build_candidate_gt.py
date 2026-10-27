#!/usr/bin/env python
"""Build the draft ground-truth table of conditions with plausibly known kinase activity changes.

Inputs
  cache/conditions.parquet, cache/hgnc_symbols.parquet        (from scripts/prep_data.py)
  data/kinhib_inhib_comb__1.csv                                (Hijazi inhibitor x kinase matrix)
  data/Klaeger/aan4368_Table_S1.xlsx, Table_S3.xlsx            (Klaeger 2017 kinobeads)
  results/ground_truth/drug_targets_manual.tsv                 (manual: agent -> direct targets)
  results/ground_truth/condition_rules_manual.tsv              (manual: regex rules over condition metadata)

Outputs (results/ground_truth/)
  candidate_gt.tsv             long: one row per (conditionID, kinase); expected_direction -1/+1/0
  candidate_gt_conditions.tsv  one row per condition (kinases joined)
  hijazi_conditions.tsv        Hijazi condition -> cell line, inhibitor, kinhib column, Klaeger name, targets
  klaeger_hijazi_targets.tsv   per Hijazi drug x kinase: in kinhib matrix? Klaeger Kd_app (nM)
  gt_summary.json              counts

Dose rule (decryptM 36926954 only, where every condition is one dose of a series): if a Klaeger Kd_app
for (drug, kinase) exists: dose < Kd/10 -> confidence 'low'; Kd/10 <= dose < Kd -> downgrade one level.
Without Kd: dose < 10 nM -> 'low'; 10 <= dose < 100 nM -> downgrade one level. 0 nM conditions are
emitted as null controls (expected_direction 0, all kinases).
"""
import json
import warnings
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from prep_data import normalise_symbols  # noqa: E402

DATA = ROOT / "data"
CACHE = ROOT / "cache"
GT = ROOT / "results" / "ground_truth"
HIJAZI_PMID = "31959955"
DECRYPTM_PMID = "36926954"
CONF = ["low", "medium", "high"]


def norm(s):
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def load_klaeger(hg):
    ann = pd.read_excel(DATA / "Klaeger" / "aan4368_Table_S1.xlsx", sheet_name="Inhibitor annotation")
    syn = {}
    for _, r in ann.iterrows():
        for s in [r["Drug"]] + str(r["Synonyms"]).split(","):
            if s and s != "nan":
                syn.setdefault(norm(s), set()).add(r["Drug"])
    kd = pd.read_excel(DATA / "Klaeger" / "Table_S3.xlsx", sheet_name="Kinobeads Drugmatrix detailed")
    kd = kd.rename(columns={"Name": "gene"})
    long = kd.melt(id_vars=["gene", "Kinase", "Direct binder"], var_name="drug", value_name="kd")
    long["kd_nM"] = pd.to_numeric(long["kd"], errors="coerce")
    long["status"] = np.where(long["kd_nM"].notna(), "kd", long["kd"].astype(str))
    long["gene1"] = long["gene"].astype(str).str.split(";").str[0]
    sym = normalise_symbols(long["gene1"], hg).set_index("name")["hgnc_symbol"]
    long["kinase"] = sym.reindex(long["gene1"]).values
    return ann, syn, long


def perturbation_class(t):
    if t.startswith("null_control") or t == "negative_control":
        return "negative_control"
    if t.startswith("inhibitor") or t in ("antibody", "chemical_genetic"):
        return "inhibitor"
    if t in ("knockout", "knockdown", "overexpression", "overexpression_or_rescue", "mutation_knockin") or t.startswith("genetic"):
        return "genetic"
    if t.startswith("cell_cycle"):
        return "cell_cycle"
    return "stimulus"


def adjust(conf, dose, kd):
    i = CONF.index(conf)
    if kd is not None and np.isfinite(kd):
        if dose < kd / 10:
            return "low"
        if dose < kd:
            return CONF[max(0, i - 1)]
        return conf
    if dose < 10:
        return "low"
    if dose < 100:
        return CONF[max(0, i - 1)]
    return conf


def main():
    warnings.filterwarnings("ignore", message="This pattern is interpreted")
    conds = pd.read_parquet(CACHE / "conditions.parquet")
    hg = pd.read_parquet(CACHE / "hgnc_symbols.parquet")
    ann, syn, kdl = load_klaeger(hg)
    kd_lookup = kdl.dropna(subset=["kd_nM"]).groupby(["drug", "kinase"])["kd_nM"].min()

    drugs = pd.read_csv(GT / "drug_targets_manual.tsv", sep="\t", dtype=str).fillna("")
    drugs = drugs.set_index("agent")
    rules = pd.read_csv(GT / "condition_rules_manual.tsv", sep="\t", dtype=str).fillna("")
    missing = sorted(set(rules["agent"]) - set(drugs.index) - {""})
    assert not missing, f"agents missing from drug_targets_manual.tsv: {missing}"
    assert not drugs.index.duplicated().any() and not rules["rule_id"].duplicated().any()

    rows = []
    # ------------------------------------------------------------------ Hijazi
    kin = pd.read_csv(DATA / "kinhib_inhib_comb__1.csv", index_col=0)
    ksym = normalise_symbols(pd.Series(kin.index), hg).set_index("name")
    kin.index = [ksym.loc[k, "hgnc_symbol"] for k in kin.index]
    kin = kin.groupby(level=0).max()
    hz = conds[conds["PMID"] == HIJAZI_PMID].copy()
    hz["cell_line"] = hz["Condition"].str.split("_", n=1).str[0]
    hz["drug"] = hz["Condition"].str.split("_", n=1).str[1]
    hrows = []
    for _, r in hz.iterrows():
        col = r["drug"] if r["drug"] in kin.columns else None
        targets = sorted(kin.index[kin[col] == 1]) if col else []
        kl = sorted(syn.get(norm(r["drug"]), []))
        hrows.append({"conditionID": r["conditionID"], "cell_line": r["cell_line"], "drug": r["drug"],
                      "kinhib_column": col or "", "klaeger_name": ";".join(kl), "n_targets": len(targets),
                      "targets": ";".join(targets), "control": r["Control"], "timepoint": r["Timepoint"]})
        if not targets:
            rows.append(dict(conditionID=r["conditionID"], PMID=HIJAZI_PMID, treatment=r["drug"],
                             perturbation_type="inhibitor_no_kinase_target", kinase="", relation="",
                             expected_direction=0, confidence="", inferred="no", rule_id="HIJAZI",
                             evidence=f"{r['drug']} not a kinase inhibitor / no targets in kinhib matrix",
                             source="kinhib_inhib_comb__1.csv", cell_line=r["cell_line"]))
        for k in targets:
            rows.append(dict(conditionID=r["conditionID"], PMID=HIJAZI_PMID, treatment=r["drug"],
                             perturbation_type="inhibitor", kinase=k, relation="direct", expected_direction=-1,
                             confidence="high", inferred="no", rule_id="HIJAZI",
                             evidence=f"Hijazi kinhib matrix target of {r['drug']} ({len(targets)} targets); 1h, 1-10 uM",
                             source="kinhib_inhib_comb__1.csv", cell_line=r["cell_line"]))
    hzt = pd.DataFrame(hrows)
    hzt.to_csv(GT / "hijazi_conditions.tsv", sep="\t", index=False)

    # Klaeger cross-check for Hijazi drugs
    krows = []
    for d in sorted(hzt["drug"].unique()):
        col = d if d in kin.columns else None
        names = sorted(syn.get(norm(d), []))
        mat = set(kin.index[kin[col] == 1]) if col else set()
        kl = kdl[kdl["drug"].isin(names)]
        klk = kl.dropna(subset=["kd_nM"]).groupby("kinase")["kd_nM"].min()
        for k in sorted(mat | set(klk.index)):
            krows.append({"drug": d, "klaeger_name": ";".join(names), "kinase": k, "in_kinhib": k in mat,
                          "klaeger_kd_nM": klk.get(k, np.nan),
                          "klaeger_status": ("kd" if k in klk.index else
                                             (kl.loc[kl["kinase"] == k, "status"].iloc[0]
                                              if (kl["kinase"] == k).any() else ("drug_not_profiled" if not names else "not_in_kinobeads")))})
    kh = pd.DataFrame(krows)
    kh.to_csv(GT / "klaeger_hijazi_targets.tsv", sep="\t", index=False)

    # ------------------------------------------------------------------ rule-based
    dose_re = re.compile(r"_([0-9.]+(?:e\+?\d+)?)nM$")
    unmatched_rules = []
    for _, ru in rules.iterrows():
        c = conds[conds["PMID"] == ru["PMID"]]
        m = c["Condition"].astype(str).str.contains(ru["cond_regex"], regex=True)
        if ru["ctl_regex"]:
            m &= c["Control"].astype(str).str.contains(ru["ctl_regex"], regex=True)
        if ru["sample_regex"]:
            m &= c["Sample"].astype(str).str.contains(ru["sample_regex"], regex=True)
        hit = c[m]
        if hit.empty:
            unmatched_rules.append(ru["rule_id"])
            continue
        agent = ru["agent"]
        if ru["kinases"]:
            kins = ru["kinases"].split(";")
        else:
            kins = drugs.loc[agent, "direct_targets"].split(";")
        if kins == ["PAN"]:
            continue
        treat = drugs.loc[agent, "canonical_name"] if agent else ""
        klname = drugs.loc[agent, "klaeger_name"] if agent else ""
        for _, r in hit.iterrows():
            dose = None
            if ru["PMID"] == DECRYPTM_PMID:
                dm = dose_re.search(str(r["Condition"]))
                dose = float(dm.group(1)) if dm else None
            for k in kins:
                conf = ru["confidence"]
                ptype = ru["perturbation_type"]
                direction = int(ru["expected_direction"])
                ev = ru["note"]
                if dose is not None:
                    if dose == 0:
                        continue  # handled by null-control block below
                    kd = kd_lookup.get((klname, k), None) if klname else None
                    conf = adjust(conf, dose, kd)
                    ev = f"{ev}; dose {dose:g} nM" + (f", Klaeger Kd_app {kd:g} nM" if kd is not None else "")
                rows.append(dict(conditionID=r["conditionID"], PMID=ru["PMID"],
                                 treatment=treat or str(r["Condition"]), perturbation_type=ptype, kinase=k,
                                 relation=ru["relation"], expected_direction=direction, confidence=conf,
                                 inferred=ru["inferred"], rule_id=ru["rule_id"], evidence=ev,
                                 source="condition_rules_manual.tsv", cell_line=r["Sample"]))
    # decryptM 0 nM null controls (condition and control both untreated)
    dz = conds[(conds["PMID"] == DECRYPTM_PMID) & conds["Condition"].astype(str).str.contains(r"_0nM$")]
    for _, r in dz.iterrows():
        rows.append(dict(conditionID=r["conditionID"], PMID=DECRYPTM_PMID, treatment=r["Condition"],
                         perturbation_type="null_control_0nM", kinase="ALL", relation="", expected_direction=0,
                         confidence="high", inferred="no", rule_id="DECRYPTM_0nM",
                         evidence="0 nM vs 0 nM (vehicle vs vehicle); no kinase change expected",
                         source="conditions metadata", cell_line=r["Sample"]))

    gt = pd.DataFrame(rows)
    # downstream rows without own perturbation_type inherit the condition's direct-row type
    pt = gt[gt["perturbation_type"] != ""].groupby("conditionID")["perturbation_type"].first()
    empty = gt["perturbation_type"] == ""
    gt.loc[empty, "perturbation_type"] = gt.loc[empty, "conditionID"].map(pt).fillna("stimulus_or_inhibitor")
    gt["perturbation_class"] = gt["perturbation_type"].map(perturbation_class)
    # merge duplicates (same condition x kinase): keep max confidence, join rule ids
    gt["conf_rank"] = gt["confidence"].map({c: i for i, c in enumerate(CONF)}).fillna(-1)
    agg = []
    for (cid, k), g in gt.groupby(["conditionID", "kinase"], sort=False):
        g = g.sort_values("conf_rank", ascending=False)
        top = g.iloc[0].copy()
        dirs = set(g["expected_direction"])
        if len(dirs) > 1:
            top["evidence"] += " | CONFLICT directions " + str(sorted(dirs))
            top["confidence"] = "low"
        top["rule_id"] = ";".join(dict.fromkeys(g["rule_id"]))
        if (g["relation"] == "direct").any():
            top["relation"] = "direct"
        agg.append(top)
    gt = pd.DataFrame(agg).drop(columns="conf_rank")
    gt = gt.merge(conds[["conditionID", "Condition", "Control", "Timepoint"]], on="conditionID", how="left")
    cols = ["conditionID", "PMID", "treatment", "perturbation_type", "perturbation_class", "kinase", "relation", "expected_direction",
            "confidence", "inferred", "evidence", "rule_id", "source", "cell_line", "Condition", "Control", "Timepoint"]
    gt = gt[cols].sort_values(["PMID", "conditionID", "relation", "kinase"])
    gt.to_csv(GT / "candidate_gt.tsv", sep="\t", index=False)

    # condition-level summary
    def best(s):
        r = [CONF.index(x) for x in s if x in CONF]
        return CONF[max(r)] if r else ""
    cs = gt.groupby("conditionID").agg(
        PMID=("PMID", "first"), treatment=("treatment", lambda s: ";".join(dict.fromkeys(s))),
        perturbation_type=("perturbation_type", lambda s: ";".join(dict.fromkeys(s))),
        perturbation_class=("perturbation_class", lambda s: ";".join(dict.fromkeys(s))),
        kinases_direct=("kinase", lambda s: ";".join(s[gt.loc[s.index, "relation"] == "direct"])),
        kinases_downstream=("kinase", lambda s: ";".join(s[gt.loc[s.index, "relation"] == "downstream"])),
        expected_direction=("expected_direction", lambda s: ";".join(map(str, sorted(set(s))))),
        best_confidence=("confidence", best), inferred=("inferred", lambda s: "yes" if (s == "yes").any() else "no"),
        cell_line=("cell_line", "first"), Condition=("Condition", "first"), Control=("Control", "first"))
    cs.reset_index().to_csv(GT / "candidate_gt_conditions.tsv", sep="\t", index=False)

    # summary
    g_nonnull = gt[gt["expected_direction"] != 0]
    ptype_group = g_nonnull["perturbation_class"]
    cond_level = g_nonnull.assign(pg=ptype_group).groupby("conditionID").agg(
        pg=("pg", lambda s: ";".join(sorted(set(s)))), conf=("confidence", best), PMID=("PMID", "first"))
    summ = {
        "n_conditions_with_kinase_GT": int(cond_level.shape[0]),
        "n_studies_with_kinase_GT": int(cond_level["PMID"].nunique()),
        "n_condition_kinase_pairs": int(len(g_nonnull)),
        "n_condition_kinase_pairs_direct": int((g_nonnull["relation"] == "direct").sum()),
        "conditions_by_best_confidence": cond_level["conf"].value_counts().to_dict(),
        "conditions_by_type": cond_level["pg"].value_counts().to_dict(),
        "conditions_by_PMID": cond_level["PMID"].value_counts().to_dict(),
        "conditions_hijazi": int((cond_level["PMID"] == HIJAZI_PMID).sum()),
        "conditions_nonhijazi": int((cond_level["PMID"] != HIJAZI_PMID).sum()),
        "conditions_by_class_and_confidence": cond_level.groupby(["pg", "conf"]).size().unstack(fill_value=0).to_dict(orient="index"),
        "nonhijazi_by_best_confidence": cond_level[cond_level["PMID"] != HIJAZI_PMID]["conf"].value_counts().to_dict(),
        "unique_kinases_direct_any": int(g_nonnull.loc[g_nonnull["relation"] == "direct", "kinase"].nunique()),
        "unique_kinases_direct_highconf_nonhijazi": int(g_nonnull.loc[(g_nonnull["relation"] == "direct")
                                                                       & (g_nonnull["confidence"] == "high")
                                                                       & (g_nonnull["PMID"] != HIJAZI_PMID), "kinase"].nunique()),
        "null_control_conditions": int(gt.loc[gt["expected_direction"] == 0, "conditionID"].nunique()),
        "rules_without_match": unmatched_rules,
        "direction_conflicts": int(gt["evidence"].str.contains("CONFLICT").sum()),
        "hijazi": {
            "n_conditions": int(len(hzt)), "n_drugs": int(hzt["drug"].nunique()),
            "cell_lines": hzt["cell_line"].value_counts().to_dict(),
            "drugs_without_kinhib_column": sorted(hzt.loc[hzt["kinhib_column"] == "", "drug"].unique()),
            "drugs_with_zero_targets": sorted(hzt.loc[(hzt["kinhib_column"] != "") & (hzt["n_targets"] == 0), "drug"].unique()),
            "kinhib_columns_without_conditions": sorted(set(kin.columns) - set(hzt["drug"])),
            "n_unique_target_kinases": int(len(set(";".join(hzt["targets"]).split(";")) - {""})),
            "targets_per_drug": hzt.drop_duplicates("drug").set_index("drug")["n_targets"].describe().round(2).to_dict(),
            "drugs_matched_to_klaeger": int((hzt.drop_duplicates("drug")["klaeger_name"] != "").sum()),
        },
    }
    khm = kh[kh["klaeger_status"] != "drug_not_profiled"]
    summ["klaeger_vs_kinhib"] = {
        "drugs_profiled": int(khm["drug"].nunique()),
        "kinhib_pairs_for_profiled_drugs": int(khm["in_kinhib"].sum()),
        "kinhib_pairs_with_klaeger_kd": int((khm["in_kinhib"] & khm["klaeger_kd_nM"].notna()).sum()),
        "kinhib_pairs_with_kd_le_1uM": int((khm["in_kinhib"] & (khm["klaeger_kd_nM"] <= 1000)).sum()),
        "kinhib_pairs_not_in_kinobeads_or_ni": int((khm["in_kinhib"] & khm["klaeger_kd_nM"].isna()).sum()),
        "klaeger_kd_le_1uM_not_in_kinhib": int((~khm["in_kinhib"] & (khm["klaeger_kd_nM"] <= 1000)).sum()),
        "klaeger_kd_le_100nM_not_in_kinhib": int((~khm["in_kinhib"] & (khm["klaeger_kd_nM"] <= 100)).sum()),
    }
    (GT / "gt_summary.json").write_text(json.dumps(summ, indent=1, default=str))
    print(json.dumps(summ, indent=1, default=str))


if __name__ == "__main__":
    main()
