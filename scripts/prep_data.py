#!/usr/bin/env python
"""Build reusable caches from the raw files in data/ (never committed).

Outputs (cache/, git-ignored):
  conditions.parquet          condition metadata (01_*.csv) + integer cond_idx
  sites.parquet               compendium sites: peptide_UID, acc, pos, residue, site_key, window, site_idx
  quant_long.parquet          long table: site_idx, cond_idx (int32), log2FC (float32), peptide_UID, conditionID
  quant_wide_f32.npy          dense float32 matrix [n_sites x n_conditions], NaN = not quantified
  quant_wide_sites.txt        row order (peptide_UID), == sites.parquet order
  quant_wide_conditions.txt   column order (conditionID), == conditions.parquet order
  psp_sites.parquet           Phosphorylation_site_dataset (phospho, all organisms) parsed
  psp_ks_all.parquet          Kinase_Substrate_Dataset, every row, + human-mapped substrate site + compendium match
  psp_ks_human.parquet        subset: human kinase, substrate human or mouse/rat mapped to human
  psp_regulatory.parquet      Regulatory_sites, phospho only, human + mouse/rat mapped to human, + compendium match
  hgnc_symbols.parquet        HGNC symbol / previous / alias -> approved symbol, UniProt (if external/hgnc present)
Small, committable report:
  results/data_profile/site_mapping_report.json, site_mapping_report.md

Canonical site key: f"{UniProtAcc}_{residue}{position}" e.g. "P28482_T185" (isoform suffix kept if present).
peptide_UID "Q9NRG9_495_M1" = accession _ position _ multiplicity (all M1 = singly phosphorylated site).

Usage: python scripts/prep_data.py [--skip-wide] [--skip-quant]
"""
import argparse
import json
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.csv as pacsv
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
CACHE = ROOT / "cache"
REPORT_DIR = ROOT / "results" / "data_profile"
HGNC = ROOT / "external" / "hgnc" / "hgnc_complete_set.txt"
PSP_SKIPROWS = 3  # date code, license, blank line


def log(msg, t0=[time.time()]):
    print(f"[{time.time() - t0[0]:7.1f}s] {msg}", flush=True)


def read_psp(name):
    return pd.read_csv(DATA / name, sep="\t", skiprows=PSP_SKIPROWS, dtype=str, encoding="latin-1")


def psp_date_code(name):
    with open(DATA / name, encoding="latin-1") as f:
        return f.readline().strip()


def split_rsd(s):
    """'S336' or 'S336-p' -> ('S', 336). Non-parsable -> (None, None)."""
    s = s.str.replace(r"-p$", "", regex=True)
    m = s.str.extract(r"^([A-Z])(\d+)$")
    return m[0], pd.to_numeric(m[1], errors="coerce").astype("Int64")


def win_upper(w):
    return w.str.upper()


# --------------------------------------------------------------------------- compendium
def prep_conditions():
    c = pd.read_csv(DATA / "01_phosphoproteomics_conditions.csv")
    c = c.drop(columns=[c.columns[0]])
    c["cond_idx"] = np.arange(len(c), dtype=np.int32)
    c["PMID"] = c["PMID"].astype(str)
    c.to_parquet(CACHE / "conditions.parquet", index=False)
    return c


def prep_sites():
    s = pd.read_csv(DATA / "02_phosphosite_information.csv").drop(columns=[0], errors="ignore")
    s = s[["peptide_UID", "sequence_window"]].rename(columns={"sequence_window": "window"})
    m = s["peptide_UID"].str.extract(r"^(?P<acc>.+)_(?P<pos>\d+)_M(?P<mult>\d+)$")
    assert m.notna().all().all(), "unparsable peptide_UID"
    s["acc"] = m["acc"]
    s["acc_base"] = s["acc"].str.replace(r"-\d+$", "", regex=True)
    s["pos"] = m["pos"].astype(int)
    s["multiplicity"] = m["mult"].astype(int)
    s["residue"] = s["window"].str[7]
    s["site_key"] = s["acc"] + "_" + s["residue"] + s["pos"].astype(str)
    s["site_idx"] = np.arange(len(s), dtype=np.int32)
    assert s["site_key"].is_unique
    s.to_parquet(CACHE / "sites.parquet", index=False)
    return s


def prep_quant(conds, sites, wide=True):
    t = pacsv.read_csv(
        DATA / "03_phosphosite_quantifications.csv",
        convert_options=pacsv.ConvertOptions(
            include_columns=["peptide_UID", "conditionID", "log2FC"],
            column_types={"log2FC": pa.float64()},
        ),
    )
    log(f"read quant csv: {t.num_rows} rows")
    df = t.to_pandas()
    sidx = pd.Series(sites["site_idx"].values, index=sites["peptide_UID"].values)
    cidx = pd.Series(conds["cond_idx"].values, index=conds["conditionID"].values)
    df["site_idx"] = sidx.reindex(df["peptide_UID"].values).values
    df["cond_idx"] = cidx.reindex(df["conditionID"].values).values
    assert df["site_idx"].notna().all() and df["cond_idx"].notna().all(), "quant IDs missing from metadata"
    df["site_idx"] = df["site_idx"].astype(np.int32)
    df["cond_idx"] = df["cond_idx"].astype(np.int32)
    n_nonfinite = int((~np.isfinite(df["log2FC"])).sum())
    dup = int(df.duplicated(["site_idx", "cond_idx"]).sum())
    log(f"non-finite log2FC: {n_nonfinite}; duplicated site x condition: {dup}")
    df["log2FC"] = df["log2FC"].astype(np.float32)
    df = df.sort_values(["cond_idx", "site_idx"])
    out = pa.Table.from_pandas(
        df[["site_idx", "cond_idx", "log2FC", "peptide_UID", "conditionID"]], preserve_index=False
    )
    out = out.set_column(3, "peptide_UID", out["peptide_UID"].dictionary_encode())
    out = out.set_column(4, "conditionID", out["conditionID"].dictionary_encode())
    pq.write_table(out, CACHE / "quant_long.parquet", compression="zstd")
    log("wrote quant_long.parquet")
    if wide:
        W = np.full((len(sites), len(conds)), np.nan, dtype=np.float32)
        W[df["site_idx"].values, df["cond_idx"].values] = df["log2FC"].values
        np.save(CACHE / "quant_wide_f32.npy", W)
        (CACHE / "quant_wide_sites.txt").write_text("\n".join(sites["peptide_UID"]) + "\n")
        (CACHE / "quant_wide_conditions.txt").write_text("\n".join(conds["conditionID"]) + "\n")
        log(f"wrote wide matrix {W.shape}, {W.nbytes / 1e9:.2f} GB")
    return {"n_rows": int(len(df)), "non_finite": n_nonfinite, "duplicates": dup}


# --------------------------------------------------------------------------- HGNC
def prep_hgnc():
    if not HGNC.exists():
        log("HGNC file missing; skipping symbol normalisation")
        return None
    h = pd.read_csv(HGNC, sep="\t", dtype=str,
                    usecols=["symbol", "alias_symbol", "prev_symbol", "uniprot_ids", "locus_group"])
    rows = []
    for kind, col in [("symbol", "symbol"), ("prev_symbol", "prev_symbol"), ("alias_symbol", "alias_symbol")]:
        x = pd.DataFrame({"symbol": h["symbol"], "q": h[col], "uniprot_ids": h["uniprot_ids"]}).dropna(subset=["q"])
        x["query"] = x["q"].str.split("|")
        x = x.explode("query")
        x["kind"] = kind
        rows.append(x[["query", "symbol", "uniprot_ids", "kind"]])
    m = pd.concat(rows, ignore_index=True)
    m["query_upper"] = m["query"].str.upper()
    m.to_parquet(CACHE / "hgnc_symbols.parquet", index=False)
    return m


def normalise_symbols(names, hg):
    """Map gene-like names to HGNC approved symbols. Priority: approved symbol > previous > alias;
    ambiguous (several approved symbols at the same priority) -> left unmapped and flagged."""
    out = []
    pri = {"symbol": 0, "prev_symbol": 1, "alias_symbol": 2}
    by_q = {q: g for q, g in hg.groupby("query_upper")} if hg is not None else {}
    for n in pd.unique(pd.Series(names).dropna()):
        g = by_q.get(str(n).upper())
        if g is None:
            out.append((n, n, None, "unmapped_kept_as_is"))
            continue
        g = g.assign(p=g["kind"].map(pri))
        best = g[g["p"] == g["p"].min()]
        syms = best["symbol"].unique()
        if len(syms) == 1:
            out.append((n, syms[0], best["uniprot_ids"].iloc[0], best["kind"].iloc[0]))
        else:
            out.append((n, n, None, "ambiguous:" + "|".join(syms)))
    return pd.DataFrame(out, columns=["name", "hgnc_symbol", "uniprot", "match_kind"])


# --------------------------------------------------------------------------- PSP
def prep_psp_sites():
    ps = read_psp("Phosphorylation_site_dataset")
    ps["residue"], ps["pos"] = split_rsd(ps["MOD_RSD"])
    ps["site_key"] = ps["ACC_ID"] + "_" + ps["residue"] + ps["pos"].astype(str)
    ps["window_u"] = win_upper(ps["SITE_+/-7_AA"])
    ps.to_parquet(CACHE / "psp_sites.parquet", index=False)
    return ps


def human_site_for_group(ps):
    """SITE_GRP_ID -> one human site (prefer canonical accession without isoform suffix)."""
    h = ps[ps["ORGANISM"] == "human"].copy()
    h["is_iso"] = h["ACC_ID"].str.contains("-")
    h = h.sort_values(["SITE_GRP_ID", "is_iso", "ACC_ID"])
    n_per_grp = h.groupby("SITE_GRP_ID").size()
    h = h.drop_duplicates("SITE_GRP_ID")
    h = h.set_index("SITE_GRP_ID")[["ACC_ID", "GENE", "residue", "pos", "site_key", "window_u"]]
    h["n_human_rows_in_group"] = n_per_grp.reindex(h.index).values
    return h


def match_to_compendium(df, sites, acc_col, res_col, pos_col, win_col, prefix=""):
    """Add compendium match columns. Strategy:
      1. exact key (acc+residue+pos). window_agree = PSP 15-mer == compendium 15-mer.
      2. if no exact key: unique 15-mer match on same base accession (isoform/position shift).
      3. else: unique 15-mer match anywhere (accession change) -> recorded but flagged.
    Final peptide_UID uses 1 (always) or 2; 3 is reported only (map_method='window_other_acc')."""
    key = df[acc_col] + "_" + df[res_col] + df[pos_col].astype(str)
    s_by_key = sites.set_index("site_key")
    hit = key.isin(s_by_key.index)
    res = pd.DataFrame(index=df.index)
    res["site_key_psp"] = key
    res["peptide_UID"] = pd.Series(pd.NA, index=df.index, dtype="string")
    res["map_method"] = "none"
    res.loc[hit, "peptide_UID"] = s_by_key.loc[key[hit], "peptide_UID"].values
    comp_win = s_by_key.loc[key[hit], "window"].values
    res["window_agree"] = pd.NA
    res.loc[hit, "window_agree"] = (df.loc[hit, win_col].values == comp_win)
    res.loc[hit, "map_method"] = "exact_key"
    # strategy 2 / 3: window lookup
    acc_base = df[acc_col].str.replace(r"-\d+$", "", regex=True)
    sw = sites[["acc_base", "window", "peptide_UID"]]
    win_counts = sw.groupby("window").size()
    uniq_win = sw[sw["window"].map(win_counts) == 1].set_index("window")
    pair = sw.drop_duplicates(["acc_base", "window"], keep=False).set_index(["acc_base", "window"])["peptide_UID"]
    miss = ~hit
    if miss.any():
        idx = pd.MultiIndex.from_arrays([acc_base[miss].values, df.loc[miss, win_col].values])
        p2 = pair.reindex(idx).values
        ok2 = pd.notna(p2)
        mi = df.index[miss]
        res.loc[mi[ok2], "peptide_UID"] = p2[ok2]
        res.loc[mi[ok2], "map_method"] = "window_same_acc"
        rest = mi[~ok2]
        p3 = uniq_win["peptide_UID"].reindex(df.loc[rest, win_col].values).values
        ok3 = pd.notna(p3)
        res.loc[rest[ok3], "map_method"] = "window_other_acc(not_used)"
        res.loc[rest[ok3], "window_other_acc_UID"] = p3[ok3]
    res["in_compendium"] = res["peptide_UID"].notna()
    return res.add_prefix(prefix)


def prep_psp_ks(sites, ps, hg):
    ks = read_psp("Kinase_Substrate_Dataset")
    ks["residue"], ks["pos"] = split_rsd(ks["SUB_MOD_RSD"])
    ks["window_u"] = win_upper(ks["SITE_+/-7_AA"])
    grp = human_site_for_group(ps)
    # human substrate coordinates
    is_h = ks["SUB_ORGANISM"] == "human"
    ks["hs_acc"] = np.where(is_h, ks["SUB_ACC_ID"], None)
    ks["hs_gene"] = np.where(is_h, ks["SUB_GENE"], None)
    ks["hs_residue"] = np.where(is_h, ks["residue"], None)
    ks["hs_pos"] = ks["pos"].where(is_h)
    ks["hs_window"] = np.where(is_h, ks["window_u"], None)
    ks["sub_mapping"] = np.where(is_h, "human", "unmapped")
    nh = ~is_h & ks["SITE_GRP_ID"].isin(grp.index)
    g = grp.reindex(ks.loc[nh, "SITE_GRP_ID"])
    ks.loc[nh, "hs_acc"] = g["ACC_ID"].values
    ks.loc[nh, "hs_gene"] = g["GENE"].values
    ks.loc[nh, "hs_residue"] = g["residue"].values
    ks.loc[nh, "hs_pos"] = g["pos"].values
    ks.loc[nh, "hs_window"] = g["window_u"].values
    ks.loc[nh, "sub_mapping"] = "via_SITE_GRP_ID_from_" + ks.loc[nh, "SUB_ORGANISM"]
    # kinase symbol normalisation (human kinases)
    sym = normalise_symbols(ks["GENE"], hg).set_index("name")
    ks["kinase"] = sym["hgnc_symbol"].reindex(ks["GENE"]).values
    ks["kinase_symbol_match"] = sym["match_kind"].reindex(ks["GENE"]).values
    ok = ks["hs_acc"].notna() & ks["hs_pos"].notna()
    m = match_to_compendium(ks[ok].assign(hs_pos=ks.loc[ok, "hs_pos"].astype(int)), sites,
                            "hs_acc", "hs_residue", "hs_pos", "hs_window")
    ks = ks.join(m)
    ks["in_compendium"] = ks["in_compendium"].fillna(False).astype(bool)
    ks["hs_pos"] = ks["hs_pos"].astype("Int64")
    for col in ks.columns:
        if ks[col].dtype == object:
            ks[col] = ks[col].astype("string")
    ks.to_parquet(CACHE / "psp_ks_all.parquet", index=False)
    hum = ks[(ks["KIN_ORGANISM"] == "human") & ks["hs_acc"].notna()
             & ks["hs_residue"].isin(["S", "T", "Y"])].copy()
    hum.to_parquet(CACHE / "psp_ks_human.parquet", index=False)
    return ks, hum


def prep_regulatory(sites, ps):
    r = read_psp("Regulatory_sites")
    r = r[r["MOD_RSD"].str.endswith("-p", na=False)].copy()
    r["residue"], r["pos"] = split_rsd(r["MOD_RSD"])
    r["window_u"] = win_upper(r["SITE_+/-7_AA"])
    grp = human_site_for_group(ps)
    is_h = r["ORGANISM"] == "human"
    r["hs_acc"] = np.where(is_h, r["ACC_ID"], None)
    r["hs_residue"] = np.where(is_h, r["residue"], None)
    r["hs_pos"] = r["pos"].where(is_h)
    r["hs_window"] = np.where(is_h, r["window_u"], None)
    r["sub_mapping"] = np.where(is_h, "human", "unmapped")
    nh = ~is_h & r["SITE_GRP_ID"].isin(grp.index)
    g = grp.reindex(r.loc[nh, "SITE_GRP_ID"])
    r.loc[nh, "hs_acc"] = g["ACC_ID"].values
    r.loc[nh, "hs_residue"] = g["residue"].values
    r.loc[nh, "hs_pos"] = g["pos"].values
    r.loc[nh, "hs_window"] = g["window_u"].values
    r.loc[nh, "sub_mapping"] = "via_SITE_GRP_ID_from_" + r.loc[nh, "ORGANISM"]
    ok = r["hs_acc"].notna() & r["hs_pos"].notna()
    m = match_to_compendium(r[ok].assign(hs_pos=r.loc[ok, "hs_pos"].astype(int)), sites,
                            "hs_acc", "hs_residue", "hs_pos", "hs_window")
    r = r.join(m)
    r["in_compendium"] = r["in_compendium"].fillna(False).astype(bool)
    r["hs_pos"] = r["hs_pos"].astype("Int64")
    for col in r.columns:
        if r[col].dtype == object:
            r[col] = r[col].astype("string")
    r.to_parquet(CACHE / "psp_regulatory.parquet", index=False)
    return r


# --------------------------------------------------------------------------- report
def mapping_report(sites, ps, ks_h, reg):
    rep = {}
    # (A) compendium sites vs PSP human site catalogue
    h = ps[ps["ORGANISM"] == "human"]
    hk = h.drop_duplicates("site_key").set_index("site_key")
    in_psp = sites["site_key"].isin(hk.index)
    agree = sites.loc[in_psp, "window"].values == hk.loc[sites.loc[in_psp, "site_key"], "window_u"].values
    win_any = sites["window"].isin(set(h["window_u"]))
    rep["compendium_sites"] = int(len(sites))
    rep["compendium_sites_with_isoform_acc"] = int(sites["acc"].str.contains("-").sum())
    rep["compendium_multiplicity_values"] = sites["multiplicity"].value_counts().to_dict()
    rep["compendium_residues"] = sites["residue"].value_counts().to_dict()
    rep["compendium_windows_with_terminal_padding"] = int(sites["window"].str.contains("_").sum())
    rep["compendium_sites_key_in_PSP_human"] = int(in_psp.sum())
    rep["compendium_sites_key_in_PSP_human_frac"] = round(float(in_psp.mean()), 4)
    rep["of_those_window_agrees"] = round(float(agree.mean()), 4)
    rep["compendium_window_found_anywhere_in_PSP_human"] = round(float(win_any.mean()), 4)
    rep["compendium_key_not_in_PSP_but_window_is"] = int((~in_psp & win_any).sum())
    # (B) PSP human-kinase KS sites (unique human substrate sites) -> compendium
    u = ks_h.drop_duplicates(["hs_acc", "hs_residue", "hs_pos"])
    rep["ks_human_kinase_rows"] = int(len(ks_h))
    rep["ks_rows_by_substrate_mapping"] = ks_h["sub_mapping"].value_counts().to_dict()
    rep["ks_unique_human_sites"] = int(len(u))
    rep["ks_unique_sites_isoform_acc"] = int(u["hs_acc"].str.contains("-").sum())
    rep["ks_unique_sites_by_map_method"] = u["map_method"].value_counts().to_dict()
    ex = u[u["map_method"] == "exact_key"]
    rep["ks_exact_key_window_agree_frac"] = round(float(ex["window_agree"].astype(bool).mean()), 4)
    rep["ks_exact_key_window_disagree_n"] = int((~ex["window_agree"].astype(bool)).sum())
    rep["ks_unique_sites_in_compendium"] = int(u["in_compendium"].sum())
    rep["ks_unique_sites_in_compendium_frac"] = round(float(u["in_compendium"].mean()), 4)
    rep["ks_rows_in_compendium"] = int(ks_h["in_compendium"].sum())
    # (C) regulatory
    ru = reg.drop_duplicates(["hs_acc", "hs_residue", "hs_pos"])
    rep["regulatory_phospho_unique_human_sites"] = int(ru["hs_acc"].notna().sum())
    rep["regulatory_sites_in_compendium"] = int(ru["in_compendium"].sum())
    rep["regulatory_by_map_method"] = ru["map_method"].value_counts().to_dict()
    # examples of disagreement
    bad = ex[~ex["window_agree"].astype(bool)].head(10)
    rep["examples_exact_key_window_disagree"] = [
        {"key": k, "psp": a, "compendium": sites.set_index("site_key").loc[k, "window"]}
        for k, a in zip(bad["site_key_psp"], bad["hs_window"])]
    return rep


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-wide", action="store_true")
    ap.add_argument("--skip-quant", action="store_true")
    a = ap.parse_args()
    CACHE.mkdir(exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    conds = prep_conditions(); log(f"conditions {len(conds)}")
    sites = prep_sites(); log(f"sites {len(sites)}")
    qinfo = {} if a.skip_quant else prep_quant(conds, sites, wide=not a.skip_wide)
    hg = prep_hgnc(); log("hgnc done")
    ps = prep_psp_sites(); log(f"psp sites {len(ps)}")
    ks, ks_h = prep_psp_ks(sites, ps, hg); log(f"psp ks all {len(ks)}, human-kinase {len(ks_h)}")
    reg = prep_regulatory(sites, ps); log(f"regulatory phospho rows {len(reg)}")
    rep = mapping_report(sites, ps, ks_h, reg)
    rep["quant"] = qinfo
    rep["psp_date_codes"] = {n: psp_date_code(n) for n in
                             ["Kinase_Substrate_Dataset", "Phosphorylation_site_dataset", "Regulatory_sites"]}
    unm = ks_h.loc[ks_h["kinase_symbol_match"] != "symbol", ["GENE", "kinase", "kinase_symbol_match"]]
    rep["psp_kinase_symbols_renamed_or_unmapped"] = unm.drop_duplicates().astype(str).values.tolist()
    (REPORT_DIR / "site_mapping_report.json").write_text(json.dumps(rep, indent=1, default=str))
    log(json.dumps({k: v for k, v in rep.items() if not k.startswith("examples")}, indent=1, default=str))


if __name__ == "__main__":
    main()
