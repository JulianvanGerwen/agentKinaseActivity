"""Tests for the ground truth (scripts/build_gt.py) and split (scripts/make_split.py). DESIGN §3/§4, B5, B9, B10, B12.

Reads the agent-visible DEV files (results/ground_truth) and the sealed TEST files (cache/sealed) - this test module is
PI infrastructure, not a development agent."""
import hashlib
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
GT = ROOT / "results" / "ground_truth"
SEALED = ROOT / "cache" / "sealed"
PY = sys.executable
LIPID = {"PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG", "PIK3C2A", "PIK3C2B", "PIK3C2G", "PIK3C3", "PI4KA", "PI4KB",
         "PIP4K2A", "PIP4K2B", "PIP4K2C", "PIP5K1A", "PIP5K1B", "PIP5K1C", "PIKFYVE"}


def rd(path, **kw):
    return pd.read_csv(path, sep="\t", keep_default_na=False, **kw)


@pytest.fixture(scope="module")
def gt():
    pairs = pd.concat([rd(GT / "gt_pairs.tsv"), rd(SEALED / "gt_pairs.tsv")], ignore_index=True)
    conds = pd.concat([rd(GT / "gt_conditions.tsv", dtype={"PMID": str}),
                       rd(SEALED / "gt_conditions.tsv", dtype={"PMID": str})], ignore_index=True)
    groups = pd.concat([rd(GT / "condition_groups.tsv", dtype=str), rd(SEALED / "condition_groups.tsv", dtype=str)],
                       ignore_index=True)
    U = set(rd(GT / "kinase_universe.tsv")["kinase"])
    return pairs, conds, groups, U


def test_universe():
    u = rd(GT / "kinase_universe.tsv")
    assert len(u) == 503 and u["kinase"].is_unique
    assert list(u.columns) == ["kinase", "uniprot", "manning_group", "manning_family", "is_pseudokinase"]
    assert not (set(u["kinase"]) & LIPID)


def test_columns_contract(gt):
    pairs, conds, _, _ = gt
    assert list(rd(GT / "gt_pairs.tsv").columns) == ["conditionID", "kinase", "cls", "sign", "confidence", "rule_id", "note"]
    assert list(rd(GT / "gt_conditions.tsv").columns) == ["conditionID", "PMID", "source", "group_id", "n_targets",
                                                        "in_M1sel", "role", "split"]
    assert set(pairs["cls"]) <= {"direct", "downstream", "masked"}
    assert set(conds["source"]) <= {"H", "D", "O"}
    assert set(conds["role"]) <= {"gt", "null", "specificity", "excluded"}
    assert set(conds["split"]) <= {f"fold{i}" for i in range(5)} | {"devconfirm", "test"}
    assert not pairs.duplicated(["conditionID", "kinase"]).any()


def test_gt_kinases_in_universe(gt):
    pairs, _, _, U = gt
    assert set(pairs["kinase"]) <= U


def test_no_bad_tokens(gt):
    pairs, _, _, _ = gt
    k = set(pairs["kinase"])
    assert not (k & LIPID)
    for bad in ("BCR/ABL", "Q6ZSR9", "ALL", "", "PRKY", "PDK1"):
        assert bad not in k
    assert not any("/" in x or ";" in x or " " in x for x in k)


def test_signs(gt):
    pairs, _, _, _ = gt
    pos = pairs[pairs["cls"].isin(["direct", "downstream"])]
    assert set(pos["sign"].astype(int)) <= {-1, 1}
    assert (pairs.loc[pairs["cls"] == "masked", "sign"].astype(int) == 0).all()
    assert set(pos["confidence"]) <= {"high", "medium", "low"}


def test_n_targets_and_m1sel(gt):
    pairs, conds, _, _ = gt
    d = pairs[pairs["cls"] == "direct"].groupby("conditionID").size()
    hm = pairs[(pairs["cls"] == "direct") & pairs["confidence"].isin(["high", "medium"])].groupby("conditionID").size()
    c = conds.set_index("conditionID")
    assert (c["n_targets"].astype(int) == d.reindex(c.index, fill_value=0)).all()
    sel = c["in_M1sel"].astype(str) == "True"
    exp = (c["role"] == "gt") & c["n_targets"].astype(int).between(1, 10) & (hm.reindex(c.index, fill_value=0) >= 1)
    assert (sel == exp).all()


def test_conditions_exist_in_compendium(gt):
    _, conds, groups, _ = gt
    comp = pd.read_parquet(ROOT / "cache" / "conditions.parquet")
    ids = set(comp["conditionID"])
    assert set(conds["conditionID"]) <= ids
    assert set(groups["conditionID"]) == ids and groups["conditionID"].is_unique and len(groups) == 1937
    pmid = comp.set_index("conditionID")["PMID"].astype(str)
    assert (conds.set_index("conditionID")["PMID"] == pmid.reindex(conds["conditionID"]).values).all()


def test_group_disjointness(gt):
    _, conds, groups, _ = gt
    s = groups.groupby("group_id")["split"].nunique()
    assert (s == 1).all(), "a group spans several splits"
    test_g = set(groups.loc[groups["split"] == "test", "group_id"])
    dev_g = set(groups.loc[groups["split"] != "test", "group_id"])
    assert not (test_g & dev_g)
    m = conds.merge(groups[["conditionID", "group_id", "split"]], on="conditionID", suffixes=("", "_g"))
    assert (m["group_id"] == m["group_id_g"]).all() and (m["split"] == m["split_g"]).all()


def test_no_drug_shared_across_test_dev(gt):
    _, _, groups, _ = gt
    g = groups.assign(drug=groups["drugs"].str.split(";")).explode("drug")
    g = g[g["drug"] != ""]
    t = set(g.loc[g["split"] == "test", "drug"])
    d = set(g.loc[g["split"] != "test", "drug"])
    assert not (t & d), f"drugs in both TEST and DEV: {sorted(t & d)[:10]}"
    # Hijazi <-> Wilkes linked through shared drugs
    for drug in ("torin1", "MK2206", "BIM1", "trametinib"):
        assert g.loc[g["drug"] == drug, "split"].nunique() == 1


def test_null_controls_follow_drug_groups(gt):
    _, conds, groups, _ = gt
    comp = pd.read_parquet(ROOT / "cache" / "conditions.parquet")
    nulls = conds[conds["role"] == "null"]
    assert len(nulls) == int(((comp.PMID.astype(str) == "36926954") & (comp.Condition == comp.Control)).sum()) == 16
    gm = groups.set_index("conditionID")
    for cid in nulls["conditionID"]:
        drugs = gm.at[cid, "drugs"].split(";")
        assert drugs and drugs[0]
        same = groups[groups["drugs"].str.split(";").map(lambda x: bool(set(x) & set(drugs)))]
        assert same["split"].nunique() == 1 and same["group_id"].nunique() == 1


def test_results_contain_no_test_rows(gt):
    salt = (SEALED / "salt.txt").read_text().strip()
    hashes = set((GT / "test_condition_hashes.txt").read_text().split())
    for f in ("gt_pairs.tsv", "gt_conditions.tsv", "harmony_bins.tsv", "condition_groups.tsv"):
        ids = set(rd(GT / f)["conditionID"])
        assert not {hashlib.sha256((salt + c).encode()).hexdigest() for c in ids} & hashes, f
    sealed = set(rd(SEALED / "gt_conditions.tsv")["split"])
    assert sealed == {"test"}
    for old in ("candidate_gt.tsv", "candidate_gt_conditions.tsv", "hijazi_conditions.tsv"):
        assert not (GT / old).exists()


def test_dev_cache_sealed(gt):
    _, _, groups, _ = gt
    dev = ROOT / "cache" / "dev"
    ids = [l.strip() for l in open(dev / "conditions_dev.txt") if l.strip()]
    sp = groups.set_index("conditionID")["split"]
    assert set(sp.reindex(ids)) <= {f"fold{i}" for i in range(5)}
    assert set(ids) == set(sp[sp.str.startswith("fold")].index)
    X = np.load(dev / "quant_dev.npy", mmap_mode="r")
    assert X.shape[1] == len(ids)


def test_paradox_rule(gt):
    """Type I/I½ RAF inhibitor in BRAF-wt (or unknown) cells -> no RAF-family positives (r3 M2)."""
    pairs, conds, groups, _ = gt
    full = rd(SEALED / "full" / "gt_conditions_full.tsv", dtype={"PMID": str})
    para = full[full["agents"].str.contains("vemurafenib|dabrafenib|SB590885|BRAFi") & (full["braf"] != "V600")]
    assert len(para) > 0
    p = pairs[pairs["conditionID"].isin(para["conditionID"]) & pairs["cls"].isin(["direct", "downstream"])]
    assert not set(p["kinase"]) & {"BRAF", "ARAF", "RAF1", "MAPK1", "MAPK3", "MAP2K1", "MAP2K2"}
    v600 = full[full["agents"].str.contains("vemurafenib") & (full["braf"] == "V600") & (full["role"] == "gt")]
    pv = pairs[pairs["conditionID"].isin(v600["conditionID"]) & (pairs["cls"] == "direct")]
    assert "BRAF" in set(pv["kinase"])


def test_precedence_downstream_not_masked(gt):
    pairs, _, _, _ = gt
    # EGF conditions: MAPK1/3 direct +1 (stimulus table), RSK downstream (kept despite the EGF pathway mask)
    comp = pd.read_parquet(ROOT / "cache" / "conditions.parquet")
    egf = comp[(comp.PMID.astype(str) == "17081983")]["conditionID"]
    p = pairs[pairs.conditionID.isin(egf)].set_index(["conditionID", "kinase"])
    for c in egf:
        assert p.loc[(c, "MAPK1"), "cls"] == "direct" and int(p.loc[(c, "MAPK1"), "sign"]) == 1
        assert p.loc[(c, "RPS6KA1"), "cls"] == "downstream"


def test_stimulus_time_window(gt):
    pairs, conds, _, _ = gt
    full = rd(SEALED / "full" / "gt_conditions_full.tsv", dtype={"PMID": str})
    late = full[(pd.to_numeric(full["minutes"], errors="coerce") > 240)]
    rules = rd(SEALED / "full" / "gt_pairs_full.tsv")
    stim_late = rules[rules.conditionID.isin(late.conditionID) & rules.rule_id.str.contains("TIME_GT4H")]
    assert len(stim_late) > 0 and (stim_late["cls"] == "masked").all()


def _sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def test_deterministic_rebuild(tmp_path):
    env = {"PYTHONHASHSEED": "0"}
    import os
    env = {**os.environ, **env}
    subprocess.run([PY, str(ROOT / "scripts" / "build_gt.py"), "--outdir", str(tmp_path)], check=True, env=env,
                   cwd=ROOT / "scripts", capture_output=True)
    subprocess.run([PY, str(ROOT / "scripts" / "make_split.py"), "--outdir", str(tmp_path), "--fulldir",
                    str(tmp_path / "sealed_full")], check=True, env=env, cwd=ROOT / "scripts", capture_output=True)
    for f in ("gt_pairs_full.tsv", "gt_conditions_full.tsv", "harmony_full.tsv"):
        assert _sha(tmp_path / "sealed_full" / f) == _sha(SEALED / "full" / f), f
    for f in ("gt_pairs.tsv", "gt_conditions.tsv", "harmony_bins.tsv", "harmony_bins_group.tsv", "condition_groups.tsv",
              "kinase_universe.tsv", "condition_doses.tsv", "condition_qc.tsv", "split_candidates.tsv",
              "test_condition_hashes.txt"):
        assert _sha(tmp_path / "results" / f) == _sha(GT / f), f
    for f in ("gt_pairs.tsv", "gt_conditions.tsv", "harmony_bins.tsv", "condition_groups.tsv"):
        assert _sha(tmp_path / "sealed" / f) == _sha(SEALED / f), f


def test_combination_ambiguous_low_and_cellular_gate(gt):
    """Audit 2026-10-27: decryptM combinations whose per-drug doses are ambiguous carry only low-confidence labels;
    cellular-potency overrides are used as the dose gate (gefitinib/EGFR 54 nM)."""
    pairs, conds = gt[0], gt[1]
    doses = rd(GT / "condition_doses.tsv")
    amb = set(doses.loc[doses.dose_source.str.contains("ambiguous"), "conditionID"])
    assert amb, "no ambiguous combination conditions found"
    pos = pairs[pairs.conditionID.isin(amb) & pairs.cls.isin(["direct", "downstream"])]
    assert (pos.confidence == "low").all()
    assert not conds.loc[conds.conditionID.isin(amb), "in_M1sel"].astype(str).eq("True").any()
    g = pairs[pairs.note.str.contains("gefitinib", case=False) & (pairs.kinase == "EGFR") & (pairs.cls == "direct")]
    assert g.note.str.contains("cellular IC50 54").any()
