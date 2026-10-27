#!/usr/bin/env python
"""Build the ground truth (GT) of kinase activity changes for the compendium (DESIGN v2 §3 + B5, B9, B11, B12, B16;
PI decisions of r3_final: precedence, stimulus table, sealing, harmony bins, decryptM per-dose M1-sel).

Deterministic, rule-file driven. Every biological decision is in results/ground_truth/*.tsv (note columns):
  drug_targets_manual.tsv    agent -> designated direct targets, target_class (pathway-mask keys), mechanism, raf_type
  drug_aliases.tsv           text token -> agent (also used by make_split.py for drug grouping)
  condition_rules_manual.tsv regex rules over condition metadata (non-Hijazi studies incl. decryptM)
  stimulus_kinases.tsv       stimulus class -> direct (receptor + first-tier) / downstream kinases
  pathway_masks.tsv          target/stimulus class -> masked downstream/feedback kinases
  offtarget_literature.tsv   literature off-targets (masked at any dose)
  drug_potency_manual.tsv    IC50/Kd for the decryptM dose rule (allosteric/covalent/PIKK/antibody or missing Kd_app)
  cell_line_genotypes.tsv    BRAF/RAS status for the paradoxical-RAF rule
  condition_doses_manual.tsv study-specific dose notations

Label classes per (condition c, kinase k in U): direct (s = -1 inhibited / +1 activated), downstream (secondary), masked
(neither positive nor negative); every other k in U is a negative. Precedence (r3 M2 / PI 2026-10-27):
  paradox-drop > direct > downstream > masked > negative.
Masks: pathway class (downstream/feedback), paradoxical RAF (BRAF-non-V600 lines + type I/I½ RAF inhibitor: RAF
positives dropped, MEK/ERK/RSK masked), off-target (Kd_app <= 3 x dose, profile precedence Klaeger > Reinecke > Davis >
Hijazi kinhib matrix; kinoplex holds no drug-kinase affinities), literature off-targets, PSP-substrate paralogs (a
kinase whose PSP substrates are >= 80 % substrates of a positive), decryptM sub-potency targets, stimulus labels > 4 h.

n_targets (B5): number of direct positives in U after the dose/Kd and paradox rules, before masking (all confidence
levels). in_M1sel: role == gt and 1 <= n_targets <= 10 and >= 1 high/medium direct positive. decryptM: each dose is its
own condition (r3 S12). Roles: gt (>= 1 direct or downstream positive), null (decryptM Condition == Control,
"0 nM vs 0 nM"), specificity (inhibitor / negative-control conditions with no positive), excluded (QC: scale error or
duplicate vector; B-§2).

Outputs (FULL tables incl. future TEST rows -> sealed area; make_split.py splits them):
  cache/sealed/full/gt_pairs_full.tsv, gt_conditions_full.tsv, harmony_full.tsv
Label-free outputs (results/ground_truth/): kinase_universe.tsv, condition_doses.tsv, drug_profile_sources.tsv,
  condition_qc.tsv, paralog_pairs.tsv, build_gt_log.json

Usage: python scripts/build_gt.py [--outdir DIR]   (DIR replaces cache/sealed/full and results/ground_truth outputs;
       used by tests for the determinism check)
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from gtlib import (CACHE, CONF, DATA, DECRYPTM_PMID, EXT, GT, HIJAZI_PMID, SEALED, DrugMatcher, SymbolMapper,
                   load_universe, norm, parse_doses_nM, parse_minutes, read_rules, write_tsv)

DEFAULT_DOSE_NM = 10000.0      # pre-registered default (B11)
OFFTARGET_FACTOR = 3.0         # mask if Kd_app <= 3 x dose (DESIGN §3(iii))
PARALOG_SHARE = 0.8            # paralog mask threshold (DESIGN §3(iv))
STIM_MAX_MIN = 240.0           # stimulus labels only <= 4 h (B12)
QC_SD, QC_MAXABS = 3.0, 20.0   # scale-error QC (data_profile_notes)
ANTIBODY_NM_PER_NGML = 1 / 148.0   # 148 kDa IgG: 1 ng/mL = 1/148 nM
RAF_FAMILY = ["BRAF", "ARAF", "RAF1"]
PARADOX_MASK = ["MAP2K1", "MAP2K2", "MAPK1", "MAPK3", "RPS6KA1", "RPS6KA2", "RPS6KA3", "RPS6KA6"]
PRIO = {"paradox": 0, "direct": 1, "downstream": 2, "masked": 3}
warnings.filterwarnings("ignore", message="This pattern is interpreted")
warnings.filterwarnings("ignore", category=UserWarning, module="openpyxl")


# ============================================================================ drug-kinase affinity profiles
def _xlsx(path, sheet):
    import openpyxl
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    return list(wb[sheet].iter_rows(values_only=True))


def _num(v):
    if v is None:
        return None
    if isinstance(v, (int, float)):
        return float(v)
    s = str(v).strip()
    if s == "n.i.":
        return math.inf
    try:
        return float(s)
    except ValueError:
        return None


def load_profiles(mapper):
    """Long table (source, src_drug, kinase, kd_nM) with protein groups expanded to all genes, and a synonym table
    (source, src_drug, name)."""
    recs, syn = [], []

    def add(source, drug, genes, kd):
        for g in str(genes).split(";"):
            k = mapper(g.strip())
            if k is not None:
                recs.append((source, str(drug), k, kd))

    # Klaeger 2017 (Kd_app nM; 'n.i.' = measured, not inhibited -> inf). Detailed matrix + cell-line re-profiling sheets.
    for sheet in ("Kinobeads Drugmatrix detailed", "in BT-474 Drugmatrix by Drug", "with Caki-1 Drugmatrix by Drug"):
        rows = _xlsx(DATA / "Klaeger" / "Table_S3.xlsx", sheet)
        hdr = rows[0]
        for r in rows[1:]:
            for d, v in zip(hdr[3:], r[3:]):
                kd = _num(v)
                if d and kd is not None and not (sheet != "Kinobeads Drugmatrix detailed" and math.isinf(kd)):
                    add("klaeger2017", d, r[0], kd)
    for r in _xlsx(DATA / "Klaeger" / "aan4368_Table_S1.xlsx", "Inhibitor annotation")[1:]:
        if r[0]:
            syn += [("klaeger2017", r[0], x.strip()) for x in re.split(r"[,;]", str(r[3] or "")) + [r[0]] if x.strip()]
    # Reinecke 2024 (kinobeads, numeric EC50/Kd_app only; empty = not reported)
    rows = _xlsx(EXT / "reinecke2024" / "41589_2023_1459_MOESM3_ESM.xlsx", "Kinobeads Drugmatrix - all")
    hdr = rows[0]
    for r in rows[1:]:
        for d, v in zip(hdr[1:], r[1:]):
            if isinstance(v, (int, float)):
                add("reinecke2024", d, r[0], float(v))
    for sheet in ("Tool compounds", "Clinical KIs"):
        for r in _xlsx(EXT / "reinecke2024" / "41589_2023_1459_MOESM2_ESM.xlsx", sheet)[1:]:
            if r[0]:
                syn += [("reinecke2024", r[0], x.strip()) for x in re.split(r"[,;]", str(r[1] or "")) + [r[0]] if x.strip()]
    # Davis 2011 (KINOMEscan Kd nM; blank = tested, Kd > 10 uM -> inf; mutant / phospho constructs skipped)
    import xlrd
    sh = xlrd.open_workbook(EXT / "davis2011" / "41587_2011_BFnbt1990_MOESM5_ESM.xls").sheet_by_index(0)
    hdr = sh.row_values(0)
    m3 = xlrd.open_workbook(EXT / "davis2011" / "41587_2011_BFnbt1990_MOESM3_ESM.xls").sheet_by_index(0)
    mutant = {m3.row_values(i)[2]: m3.row_values(i)[3] == "YES" for i in range(1, m3.nrows)}
    for i in range(1, sh.nrows):
        r = sh.row_values(i)
        if mutant.get(r[2], False) or ("(" in str(r[2]) and "phosph" in str(r[2]).lower()):
            continue
        for d, v in zip(hdr[3:], r[3:]):
            kd = math.inf if v in ("", None) else _num(v)
            if kd is not None:
                add("davis2011", d, r[1], kd)
    m4 = xlrd.open_workbook(EXT / "davis2011" / "41587_2011_BFnbt1990_MOESM4_ESM.xls").sheet_by_index(0)
    for i in range(1, m4.nrows):
        r = m4.row_values(i)
        if r[0]:
            syn += [("davis2011", r[0], x.strip()) for x in re.split(r"[,;/]", str(r[1] or "")) + [r[0]] if x.strip()]
    for s in ("klaeger2017", "reinecke2024", "davis2011"):
        pass
    prof = pd.DataFrame(recs, columns=["source", "src_drug", "kinase", "kd_nM"])
    prof = prof.groupby(["source", "src_drug", "kinase"], as_index=False)["kd_nM"].min()
    syn = pd.DataFrame(syn, columns=["source", "src_drug", "name"]).drop_duplicates()
    # Davis column headers like 'AZD-6244/ARRY-886' contain several names
    extra = [("davis2011", d, x) for d in prof.loc[prof.source == "davis2011", "src_drug"].unique()
             for x in re.split(r"[/]", d)]
    extra += [(s, d, d) for s, d in prof[["source", "src_drug"]].drop_duplicates().itertuples(index=False)]
    syn = pd.concat([syn, pd.DataFrame(extra, columns=syn.columns)]).drop_duplicates()
    syn["key"] = syn["name"].map(norm)
    syn = syn[syn["key"].str.len() >= 4]
    return prof, syn


def agent_profiles(drugs, aliases, prof, syn, kinhib):
    """agent -> {kinase: (kd_nM, source)}. Highest-priority source that profiled the drug is primary; lower sources fill
    kinases it did not measure. Hijazi kinhib matrix (Kd := 0 for listed targets) is the last source."""
    order = ["klaeger2017", "reinecke2024", "davis2011"]
    names = defaultdict(set)
    for r in drugs.itertuples():
        for x in [r.agent, r.klaeger_name] + re.split(r"[;()]", r.canonical_name) + r.profile_names.split(";"):
            if x and x.strip():
                names[r.agent].add(norm(x.strip()))
    for r in aliases.itertuples():
        if str(r.is_class) != "1" and not r.pmid:
            names[r.agent].add(norm(r.alias))
    by_key = {(s, k): d for s, d, k in syn[["source", "src_drug", "key"]].itertuples(index=False)}
    psrc = {s: {d: g for d, g in prof[prof.source == s].groupby("src_drug")} for s in order}
    kh_cols = {c: c for c in kinhib.columns}
    out, rep = {}, []
    for agent in sorted(names):
        ks = sorted(k for k in names[agent] if len(k) >= 4)
        prof_a = {}
        for s in order:
            hit = sorted({by_key[(s, k)] for k in ks if (s, k) in by_key})
            for d in hit:
                g = psrc[s].get(d)
                if g is None:
                    continue
                n_new = 0
                for k, kd in zip(g.kinase, g.kd_nM):
                    if k not in prof_a:
                        prof_a[k] = (kd, s)
                        n_new += 1
                rep.append(dict(agent=agent, source=s, matched_name=d, n_kinases=len(g),
                                n_finite=int(np.isfinite(g.kd_nM).sum()), n_used=n_new))
        # Hijazi kinhib matrix
        for col in kh_cols:
            if norm(col) in names[agent]:
                tg = kinhib.index[kinhib[col] == 1]
                n_new = 0
                for k in tg:
                    if k not in prof_a:
                        prof_a[k] = (0.0, "hijazi_kinhib")
                        n_new += 1
                rep.append(dict(agent=agent, source="hijazi_kinhib", matched_name=col, n_kinases=len(tg),
                                n_finite=len(tg), n_used=n_new))
        out[agent] = prof_a
    return out, pd.DataFrame(rep, columns=["agent", "source", "matched_name", "n_kinases", "n_finite", "n_used"])


# ============================================================================ helpers
def split_tokens(s):
    return [x for x in str(s).split(";") if x.strip()]


class Labels:
    """Collects candidate labels per (condition, kinase) and resolves them by precedence."""

    def __init__(self):
        self.rows = []

    def add(self, cid, k, kind, sign, conf, rule, note=""):
        self.rows.append((cid, k, kind, int(sign), conf, rule, note))

    def resolve(self, U):
        cand = pd.DataFrame(self.rows, columns=["conditionID", "kinase", "kind", "sign", "confidence", "rule_id", "note"])
        cand = cand[cand["kinase"].isin(U)]
        out = []
        for (cid, k), g in cand.groupby(["conditionID", "kinase"], sort=True):
            best = min(PRIO[x] for x in g["kind"])
            top = g[g["kind"].map(PRIO) == best]
            rules = ";".join(dict.fromkeys(g["rule_id"] if best == 0 else top["rule_id"]))
            if best == 0:
                out.append((cid, k, "masked", 0, "", rules, "paradoxical RAF rule: " + ";".join(dict.fromkeys(top["note"]))))
                continue
            if best == 3:
                out.append((cid, k, "masked", 0, "", rules, ";".join(x for x in dict.fromkeys(top["note"]) if x)))
                continue
            cls = "direct" if best == 1 else "downstream"
            signs = set(top["sign"])
            if len(signs) > 1:
                out.append((cid, k, "masked", 0, "", rules, f"CONFLICT: {cls} signs {sorted(signs)}"))
                continue
            conf = max(top["confidence"], key=lambda c: CONF.index(c))
            out.append((cid, k, cls, signs.pop(), conf, rules, ";".join(x for x in dict.fromkeys(top["note"]) if x)))
        return pd.DataFrame(out, columns=["conditionID", "kinase", "cls", "sign", "confidence", "rule_id", "note"])


def min_conf(a, b):
    return CONF[min(CONF.index(a), CONF.index(b))]


# ============================================================================ main
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=None, help="write all outputs below this dir (tests); default = repo locations")
    args = ap.parse_args(argv)
    out_sealed = Path(args.outdir) / "sealed_full" if args.outdir else SEALED / "full"
    out_res = Path(args.outdir) / "results" if args.outdir else GT
    log = {}

    conds = pd.read_parquet(CACHE / "conditions.parquet")
    conds["PMID"] = conds["PMID"].astype(str)
    hg = pd.read_parquet(CACHE / "hgnc_symbols.parquet")
    Udf = load_universe(hg)
    U = set(Udf["kinase"])
    mapper = SymbolMapper(Udf, hg)
    write_tsv(Udf.drop(columns=["ksm_gene"]), out_res / "kinase_universe.tsv")
    log["U_size"] = len(U)

    drugs = read_rules("drug_targets_manual.tsv")
    assert drugs["agent"].is_unique
    drugs = drugs.set_index("agent", drop=False)
    aliases = read_rules("drug_aliases.tsv")
    matcher = DrugMatcher(aliases)
    rules = read_rules("condition_rules_manual.tsv")
    assert rules["rule_id"].is_unique, "duplicate rule ids"
    stim = read_rules("stimulus_kinases.tsv")
    pmask = read_rules("pathway_masks.tsv")
    lit = read_rules("offtarget_literature.tsv")
    pot = read_rules("drug_potency_manual.tsv")
    geno = read_rules("cell_line_genotypes.tsv")
    dose_man = read_rules("condition_doses_manual.tsv")
    unknown_agents = sorted(set(rules["agent"]) - set(drugs.index) - {""})
    assert not unknown_agents, f"rule agents missing from drug_targets_manual.tsv: {unknown_agents}"

    # rule-file symbols outside U (reported, ignored)
    notU = {}
    for name, df in (("stimulus_kinases", stim), ("pathway_masks", pmask), ("offtarget_literature", lit)):
        bad = sorted(set(df["kinase"]) - U)
        if bad:
            notU[name] = bad
    for r in drugs.itertuples():
        bad = [k for k in split_tokens(r.direct_targets) if mapper(k) is None]
        if bad:
            notU.setdefault("drug_targets_manual", []).extend(f"{r.agent}:{k}" for k in bad)
    log["rule_symbols_not_in_U"] = notU
    mask_of = pmask.groupby("mask_class")["kinase"].apply(list).to_dict()
    lit_of = lit.groupby("agent")["kinase"].apply(list).to_dict()
    lit_note = {(a, k): n for a, k, n in lit[["agent", "kinase", "note"]].itertuples(index=False)}
    pot_of = {(a, k): float(v) for a, k, v in pot[["agent", "kinase", "potency_nM"]].itertuples(index=False)}

    # Hijazi kinhib matrix
    kin = pd.read_csv(DATA / "kinhib_inhib_comb__1.csv", index_col=0)
    tok_map, dropped_tok = {}, defaultdict(int)
    for t in kin.index:
        tt = "ABL1" if t == "BCR/ABL" else t
        k = mapper(tt)
        tok_map[t] = k
        if k is None:
            dropped_tok[t] += int(kin.loc[t].sum())
    kin_u = kin[[tok_map[t] is not None for t in kin.index]].copy()
    kin_u.index = [tok_map[t] for t in kin_u.index]
    kin_u = kin_u.groupby(level=0).max()
    log["kinhib_tokens_dropped_not_in_U"] = dict(sorted(dropped_tok.items()))

    # affinity profiles
    prof, syn = load_profiles(mapper)
    aprof, prep = agent_profiles(drugs, aliases, prof, syn, kin_u)
    write_tsv(prep.sort_values(["agent", "source", "matched_name"]), out_res / "drug_profile_sources.tsv")
    log["profile_sources"] = {"kinoplex": "not used: kinaseSubMain kinoplex files are a structural S/T/Y atlas, no drug-kinase affinities"}

    def kd_of(agent, k):
        v = aprof.get(agent, {}).get(k)
        return v if v else (None, None)

    # cell-line genotype
    def braf_status(sample):
        for r in geno.itertuples():
            if re.search(r.sample_regex, str(sample)):
                return r.BRAF
        return "unknown"

    # ---------------------------------------------------------------- per-condition drugs, doses, time
    crec = {}
    dose_rows = []
    dm_re = re.compile(r"^(.+?)_([0-9.]+(?:e\+?\d+)?)(nM|ng/mL)$")
    for r in conds.itertuples():
        cid = r.conditionID
        if r.PMID == HIJAZI_PMID:
            drug = str(r.Condition).split("_", 1)[1]
            cdrugs = matcher.find(drug, r.PMID)
            kdrugs = []
        else:
            cdrugs = matcher.find(r.Condition, r.PMID)
            kdrugs = matcher.find(r.Control, r.PMID)
        mins, how = parse_minutes(r.Condition, r.Timepoint)
        doses = {}
        if cdrugs or kdrugs:
            if r.PMID == HIJAZI_PMID:
                doses = {a: (DEFAULT_DOSE_NM, "hijazi_default_10uM") for a in cdrugs}
            elif r.PMID == DECRYPTM_PMID:
                m = re.match(r"^(.+?)(?:_(\d+)to(\d+))?_([0-9.]+(?:e\+?\d+)?)(nM|ng/mL)$", str(r.Condition))
                if m:
                    names = m.group(1).split("+")
                    ratios = [float(m.group(2)), float(m.group(3))] if m.group(2) else [1.0] * len(names)
                    base = float(m.group(4)) * (ANTIBODY_NM_PER_NGML if m.group(5) == "ng/mL" else 1.0)
                    for nm_, ra in zip(names, ratios):
                        for a in matcher.find(nm_, r.PMID):
                            doses[a] = (base * ra, "decryptM_label" + ("_ratio" if m.group(2) else "")
                                        + ("_ngmL_to_nM" if m.group(5) == "ng/mL" else ""))
            else:
                for dm in dose_man[dose_man["PMID"] == r.PMID].itertuples():
                    if re.search(dm.cond_regex, str(r.Condition)):
                        for a in cdrugs:
                            if dm.agent in ("*", a):
                                doses[a] = (float(dm.dose_nM), "condition_doses_manual")
                segs = re.split(r"\+|,", str(r.Condition))
                for sgm in segs:
                    ds = parse_doses_nM(sgm)
                    for a in matcher.find(sgm, r.PMID):
                        if a not in doses and len(ds) == 1:
                            doses[a] = (ds[0], "parsed_text")
                allds = parse_doses_nM(r.Condition)
                for a in cdrugs:
                    if a not in doses and len(allds) == 1:
                        doses[a] = (allds[0], "parsed_text")
            for a in sorted(set(cdrugs) | set(kdrugs)):
                if a not in doses:
                    doses[a] = (DEFAULT_DOSE_NM, "default_10uM")
                dose_rows.append(dict(conditionID=cid, PMID=r.PMID, side="condition" if a in cdrugs else "control",
                                      agent=a, dose_nM=doses[a][0], dose_source=doses[a][1]))
        crec[cid] = dict(PMID=r.PMID, Condition=r.Condition, Control=r.Control, Sample=r.Sample, minutes=mins,
                         time_how=how, cdrugs=cdrugs, kdrugs=kdrugs, doses=doses, braf=braf_status(r.Sample))
    dose_df = pd.DataFrame(dose_rows)
    write_tsv(dose_df, out_res / "condition_doses.tsv")

    L = Labels()
    role_hint = {}       # cid -> set of hints
    src_of = {}
    agents_in = defaultdict(set)   # cid -> agents to mask off-targets for
    mask_classes = defaultdict(set)
    def_targets = defaultdict(set)  # cid -> designated targets (exempt from off-target masks)

    # ---------------------------------------------------------------- GT-H Hijazi
    for r in conds[conds["PMID"] == HIJAZI_PMID].itertuples():
        cid = r.conditionID
        src_of[cid] = "H"
        drug = str(r.Condition).split("_", 1)[1]
        ags = crec[cid]["cdrugs"]
        assert len(ags) == 1, f"Hijazi drug {drug} -> {ags}"
        a = ags[0]
        agents_in[cid].add(a)
        info = drugs.loc[a]
        col = drug if drug in kin_u.columns else None
        tg = sorted(kin_u.index[kin_u[col] == 1]) if col else []
        for k in tg:
            L.add(cid, k, "direct", -1, "high", "HIJAZI", f"kinhib target of {drug}")
            def_targets[cid].add(k)
        for c in split_tokens(info.target_class):
            mask_classes[cid].add(c)
        role_hint[cid] = {"inhibitor"}

    # ---------------------------------------------------------------- rule-based (GT-D, GT-O)
    unmatched = []
    null_cids = set(conds.loc[(conds["PMID"] == DECRYPTM_PMID) & (conds["Condition"] == conds["Control"]), "conditionID"])
    for ru in rules.itertuples():
        c = conds[conds["PMID"] == ru.PMID]
        m = c["Condition"].astype(str).str.contains(ru.cond_regex, regex=True)
        if ru.ctl_regex:
            m &= c["Control"].astype(str).str.contains(ru.ctl_regex, regex=True)
        if ru.sample_regex:
            m &= c["Sample"].astype(str).str.contains(ru.sample_regex, regex=True)
        hit = c[m]
        hit = hit[~hit["conditionID"].isin(null_cids)]
        if hit.empty:
            unmatched.append(ru.rule_id)
            continue
        ptype = ru.perturbation_type
        direction = int(ru.expected_direction or 0)
        is_stim = ptype.startswith("stimulus") or bool(ru.stim_class)
        for h in hit.itertuples():
            cid = h.conditionID
            rec = crec[cid]
            src_of[cid] = "D" if ru.PMID == DECRYPTM_PMID else "O"
            role_hint.setdefault(cid, set())
            if ptype == "negative_control":
                role_hint[cid].add("negative_control")
                for a in rec["cdrugs"] + rec["kdrugs"]:
                    agents_in[cid].add(a)
                continue
            role_hint[cid].add("stimulus" if is_stim else ("inhibitor" if ru.agent else ptype))
            items = []   # (kinase, kind, sign, conf, note)
            if ru.agent:
                info = drugs.loc[ru.agent]
                agents_in[cid].add(ru.agent)
                for cc in split_tokens(info.target_class):
                    mask_classes[cid].add(cc)
                if info.mechanism == "pan":
                    role_hint[cid].add("pan")
                    continue
                kins = split_tokens(ru.kinases) if ru.kinases else split_tokens(info.direct_targets)
                rel = ru.relation or "direct"
                for t in kins:
                    k = mapper(t)
                    if k is None:
                        continue
                    items.append((k, rel, direction, ru.confidence, ru.note))
                    def_targets[cid].add(k)
            elif ru.kinases:
                for t in split_tokens(ru.kinases):
                    k = mapper(t)
                    if k is not None:
                        items.append((k, ru.relation or "direct", direction, ru.confidence, ru.note))
            if ru.stim_class:
                mask_classes[cid].add(ru.stim_class)
                for s in stim[stim["stim_class"] == ru.stim_class].itertuples():
                    items.append((s.kinase, s.cls, direction * int(s.sign), min_conf(ru.confidence, s.confidence),
                                  f"{ru.stim_class}: {s.note}"))
            if ru.mask_class:
                mask_classes[cid].add(ru.mask_class)
            for a in rec["cdrugs"] + rec["kdrugs"]:
                agents_in[cid].add(a)
            # stimulus time window (B12)
            tmin = rec["minutes"]
            for k, kind, sgn, conf, note in items:
                rid = ru.rule_id
                if is_stim and np.isfinite(tmin) and tmin > STIM_MAX_MIN:
                    L.add(cid, k, "masked", 0, "", rid + ";TIME_GT4H", f"stimulus label at {tmin:g} min > 4 h -> masked")
                    continue
                if is_stim and not np.isfinite(tmin):
                    note = (note + "; " if note else "") + "treatment time unknown"
                # decryptM dose/potency rule (GT-D)
                if ru.PMID == DECRYPTM_PMID and ru.agent:
                    a = ru.agent
                    dose = rec["doses"].get(a, (None,))[0]
                    if dose is None:
                        raise RuntimeError(f"no decryptM dose for {cid} {a}")
                    info = drugs.loc[a]
                    thr, how = None, None
                    if kind == "downstream" and ru.relation == "downstream":
                        thr, how = pot_of.get((a, "*")), "potency(*)"
                    elif kind == "direct":
                        kd, ks = kd_of(a, k)
                        if info.mechanism == "atp" and kd is not None and np.isfinite(kd) and ks != "hijazi_kinhib":
                            thr, how = kd, f"Kd_app {ks}"
                        elif (a, k) in pot_of:
                            thr, how = pot_of[(a, k)], "potency table"
                    if kind == "direct" or (kind == "downstream" and ru.relation == "downstream"):
                        if thr is None:
                            raise RuntimeError(f"decryptM {a}/{k}: no Kd_app and no potency entry")
                        note = f"{note}; dose {dose:g} nM vs {how} {thr:g} nM"
                        if dose < thr:
                            L.add(cid, k, "masked", 0, "", rid + ";SUBKD",
                                  f"decryptM dose {dose:g} nM < {how} {thr:g} nM: target masked (B11)")
                            continue
                L.add(cid, k, kind, sgn, conf, rid, note)
    log["rules_without_match"] = unmatched

    # decryptM curated downstream rows without own agent: keep only if a target of a condition drug passes the dose rule
    lab = pd.DataFrame(L.rows, columns=["conditionID", "kinase", "kind", "sign", "confidence", "rule_id", "note"])
    dm_c = set(conds.loc[conds["PMID"] == DECRYPTM_PMID, "conditionID"])
    has_direct = set(lab.loc[(lab.kind == "direct") & lab.conditionID.isin(dm_c), "conditionID"])
    for cid in dm_c:   # agents without targets in U (e.g. lipid-kinase-only PI3Ki): agent-level potency ('*') rule
        for a, (dose, _) in crec[cid]["doses"].items():
            if (a, "*") in pot_of and dose >= pot_of[(a, "*")]:
                has_direct.add(cid)
    new_rows = []
    for row in L.rows:
        cid, k, kind, sign, conf, rid, note = row
        if cid in dm_c and kind == "downstream" and cid not in has_direct and "potency(*)" not in note:
            new_rows.append((cid, k, "masked", 0, "", rid + ";SUBKD", "decryptM: no drug target above potency -> downstream masked"))
        else:
            new_rows.append(row)
    L.rows = new_rows

    # null controls
    for cid in sorted(null_cids):
        src_of[cid] = "D"
        role_hint[cid] = {"null"}

    # ---------------------------------------------------------------- paradoxical RAF (B12)
    for cid, ags in agents_in.items():
        rec = crec[cid]
        par = [a for a in ags if a in drugs.index and drugs.loc[a, "raf_type"] == "paradox"]
        if par and rec["braf"] != "V600":
            for k in RAF_FAMILY:
                L.add(cid, k, "paradox", 0, "", "PARADOX_RAF", f"{','.join(par)} in BRAF-{rec['braf']} cells")
            for k in PARADOX_MASK:
                L.add(cid, k, "paradox", 0, "", "PARADOX_RAF", f"{','.join(par)} in BRAF-{rec['braf']} cells")

    # ---------------------------------------------------------------- masks: pathway, off-target, literature
    for cid in sorted(set(src_of) - null_cids):
        for c in sorted(mask_classes[cid]):
            for k in mask_of.get(c, []):
                L.add(cid, k, "masked", 0, "", f"PATHWAY_{c}", f"pathway mask {c}")
        rec = crec[cid]
        for a in sorted(agents_in[cid]):
            dose = rec["doses"].get(a, (DEFAULT_DOSE_NM, "default"))[0]
            for k, (kd, s) in sorted(aprof.get(a, {}).items()):
                if np.isfinite(kd) and kd <= OFFTARGET_FACTOR * dose and k not in def_targets[cid]:
                    L.add(cid, k, "masked", 0, "", f"OFFTARGET_{s}", f"{a} Kd_app {kd:g} nM <= 3x dose {dose:g} nM ({s})")
            for k in lit_of.get(a, []):
                L.add(cid, k, "masked", 0, "", "OFFTARGET_LIT", f"{a}: {lit_note[(a, k)]}")

    # ---------------------------------------------------------------- paralog mask (PSP substrate overlap)
    psp = pd.read_parquet(CACHE / "psp_ks_human.parquet")
    psp = psp[psp["kinase"].isin(U) & psp["site_key_psp"].notna()]
    subs = psp.groupby("kinase")["site_key_psp"].apply(set).to_dict()
    par_rows = []
    for n, Sn in subs.items():
        for p_, Sp in subs.items():
            if n != p_ and len(Sn & Sp) / len(Sn) >= PARALOG_SHARE:
                par_rows.append((n, p_, len(Sn), len(Sn & Sp), round(len(Sn & Sp) / len(Sn), 4)))
    par = pd.DataFrame(par_rows, columns=["kinase_masked", "positive_kinase", "n_psp_substrates", "n_shared", "share"])
    write_tsv(par.sort_values(["positive_kinase", "kinase_masked"]), out_res / "paralog_pairs.tsv")
    par_of = par.groupby("positive_kinase")["kinase_masked"].apply(list).to_dict()

    res = L.resolve(U)
    extra = []
    for (cid,), g in res[res.cls.isin(["direct", "downstream"])].groupby(["conditionID"]):
        labelled = set(res.loc[res.conditionID == cid, "kinase"])
        for p_ in sorted(set(g.kinase)):
            for n in par_of.get(p_, []):
                if n not in labelled:
                    extra.append((cid, n, "masked", 0, "", "PARALOG_PSP80", f">=80% of PSP substrates shared with positive {p_}"))
    if extra:
        ex = pd.DataFrame(extra, columns=res.columns)
        ex = ex.groupby(["conditionID", "kinase"], as_index=False).agg(
            cls=("cls", "first"), sign=("sign", "first"), confidence=("confidence", "first"),
            rule_id=("rule_id", "first"), note=("note", lambda s: "; ".join(sorted(set(s)))))
        res = pd.concat([res, ex[res.columns]], ignore_index=True)
    res = res.sort_values(["conditionID", "kinase"]).reset_index(drop=True)
    assert not res.duplicated(["conditionID", "kinase"]).any()
    assert (res.loc[res.cls != "masked", "sign"] != 0).all(), "sign 0 for a positive"

    # ---------------------------------------------------------------- QC (scale errors, duplicate vectors)
    X = np.load(CACHE / "quant_wide_f32.npy", mmap_mode="r")
    cols = [l.strip() for l in open(CACHE / "quant_wide_conditions.txt")]
    assert cols == list(conds["conditionID"])
    X = np.asarray(X)
    qrows = []
    hashes = {}
    notnan = ~np.isnan(X)
    for j, cid in enumerate(cols):
        v = X[:, j]
        ok = notnan[:, j]
        vv = v[ok]
        sd = float(np.std(vv)) if vv.size else float("nan")
        mx = float(np.max(np.abs(vv))) if vv.size else float("nan")
        h = hashlib.sha256(np.where(ok, v, np.float32(-9999.0)).tobytes()).hexdigest()
        dup = hashes.get(h, "")
        hashes.setdefault(h, cid)
        status = "scale_error" if (sd > QC_SD or mx > QC_MAXABS) else ("duplicate" if dup else "ok")
        qrows.append(dict(conditionID=cid, n_sites=int(ok.sum()), sd=round(sd, 4), max_abs=round(mx, 3),
                          duplicate_of=dup, qc_status=status))
    qc = pd.DataFrame(qrows)
    write_tsv(qc, out_res / "condition_qc.tsv")
    qc_of = qc.set_index("conditionID")["qc_status"].to_dict()

    # ---------------------------------------------------------------- conditions table
    crow = []
    for cid in sorted(src_of):
        g = res[res.conditionID == cid]
        nd = int((g.cls == "direct").sum())
        nds = int((g.cls == "downstream").sum())
        hm = int(((g.cls == "direct") & g.confidence.isin(["high", "medium"])).sum())
        hints = role_hint.get(cid, set())
        time_masked = g.rule_id.str.contains("TIME_GT4H").any()
        if "null" in hints:
            role = "null"
        elif nd + nds > 0:
            role = "gt"
        elif time_masked and not (hints & {"inhibitor", "negative_control"}):
            role = "excluded"    # stimulus condition whose labels all fall outside the <=4 h window
        else:
            role = "specificity"
        if qc_of[cid] != "ok" and role != "null":
            role = "excluded"
        crow.append(dict(conditionID=cid, PMID=crec[cid]["PMID"], source=src_of[cid], n_targets=nd,
                         in_M1sel=bool(role == "gt" and 1 <= nd <= 10 and hm >= 1), role=role,
                         n_direct_highmed=hm, n_downstream=nds, n_masked=int((g.cls == "masked").sum()),
                         qc_status=qc_of[cid], minutes=crec[cid]["minutes"], time_source=crec[cid]["time_how"],
                         agents=";".join(sorted(agents_in.get(cid, []))), braf=crec[cid]["braf"]))
    gc = pd.DataFrame(crow)
    # pairs of excluded / null / specificity conditions are kept (masks only matter for gt); positives of excluded
    # conditions stay in the table but the harness never evaluates role != gt.

    # ---------------------------------------------------------------- harmony (PSP substrates quantified per condition)
    import scipy.sparse as sp
    sites = [l.strip() for l in open(CACHE / "quant_wide_sites.txt")]
    spos = {s: i for i, s in enumerate(sites)}
    q = psp[psp["peptide_UID"].notna() & psp["peptide_UID"].isin(spos)][["kinase", "peptide_UID"]].drop_duplicates()
    kl = sorted(U)
    kpos = {k: i for i, k in enumerate(kl)}
    A = sp.csr_matrix((np.ones(len(q), np.float32), ([kpos[k] for k in q.kinase], [spos[s] for s in q.peptide_UID])),
                      shape=(len(kl), len(sites)))
    cpos = {c: j for j, c in enumerate(cols)}
    gcols = [cpos[c] for c in gc.conditionID]
    cnt = np.asarray(A @ notnan[:, gcols].astype(np.float32))
    cnt_df = pd.DataFrame(cnt.astype(int), index=kl, columns=list(gc.conditionID))
    hp = res[res.cls.isin(["direct", "downstream"])][["conditionID", "kinase", "cls"]].copy()
    hp["n_psp_quantified"] = [int(cnt_df.at[k, c]) for c, k in zip(hp.conditionID, hp.kinase)]
    hp["bin"] = pd.cut(hp["n_psp_quantified"], [-1, 0, 4, 19, 10 ** 9], labels=["0", "1-4", "5-19", ">=20"]).astype(str)

    # ---------------------------------------------------------------- write
    pairs_cols = ["conditionID", "kinase", "cls", "sign", "confidence", "rule_id", "note"]
    write_tsv(res[pairs_cols], out_sealed / "gt_pairs_full.tsv")
    write_tsv(gc, out_sealed / "gt_conditions_full.tsv")
    write_tsv(hp[["conditionID", "kinase", "cls", "n_psp_quantified", "bin"]], out_sealed / "harmony_full.tsv")
    log.update({
        "n_gt_conditions_by_role": gc.role.value_counts().to_dict(),
        "n_pairs_by_cls": res.cls.value_counts().to_dict(),
        "n_conditions_default_dose": int(dose_df.loc[dose_df.dose_source.str.startswith("default"), "conditionID"].nunique()) if len(dose_df) else 0,
        "dose_sources": dose_df.dose_source.value_counts().to_dict() if len(dose_df) else {},
        "qc_status_all_conditions": qc.qc_status.value_counts().to_dict(),
    })
    (out_res / "build_gt_log.json").write_text(json.dumps(log, indent=1, sort_keys=True, default=str) + "\n")
    print(json.dumps({k: log[k] for k in ("n_gt_conditions_by_role", "n_pairs_by_cls", "rules_without_match")}, indent=1))


if __name__ == "__main__":
    main()
