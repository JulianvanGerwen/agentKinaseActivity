#!/usr/bin/env python
"""Build all kinase–substrate libraries (deterministic) -> cache/libraries/*.parquet, plus committable
results/libraries/{excluded_pmids.tsv, filter_counts.tsv, build_manifest.json, kinase_universe_used.tsv,
LIBRARY_REPORT.md, overlap.tsv, predicted_vs_psp.tsv}.

Usage:
  python scripts/build_libraries.py                 # build everything + report
  python scripts/build_libraries.py --report-only   # recompute the report from cache/libraries
  python scripts/build_libraries.py --no-fetch      # never call PubMed (uses cache/pubmed only)
Needs: cache/ from scripts/prep_data.py; external/{omnipath,iptmnet,sugiyama2019,ksm_index,ksm_pred,uniprot,
benchmarkin_py,hgnc} (see DATA.md). PubMed metadata are fetched once and cached in cache/pubmed/.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from kinact import library as L  # noqa: E402

WELL_STUDIED = ["MAPK1", "MAPK3", "CDK1", "CDK2", "AKT1", "PRKACA", "CSNK2A1", "GSK3B", "PRKCA", "MTOR",
                "ATM", "SRC", "EGFR", "ABL1", "PLK1", "AURKB", "CHEK1", "MAPK14", "RPS6KB1", "CSNK1D"]
SANITY_KINASES = ["MAPK1", "CDK1", "AKT1", "PRKACA"]


def validate(name, lib, U):
    assert list(lib.columns[:5]) == L.LIB_COLUMNS, (name, lib.columns)
    assert lib["kinase"].isin(set(U["kinase"])).all(), f"{name}: kinase not in U"
    assert L.valid_site_keys(lib["site"]).all(), f"{name}: bad site key"
    assert (lib["weight"] > 0).all() and np.isfinite(lib["weight"]).all(), name
    assert not lib.duplicated(["kinase", "site"]).any(), f"{name}: duplicate edges"


def build(fetch=True):
    L.LIBDIR.mkdir(parents=True, exist_ok=True)
    L.RESLIB.mkdir(parents=True, exist_ok=True)
    U = L.load_universe()
    U.to_csv(L.RESLIB / "kinase_universe_used.tsv", sep="\t", index=False)
    L.log(f"U: {len(U)} kinases ({U['universe_source'].iloc[0]})")
    sm = L.SiteMapper()
    km = L.KinaseMapper(U)
    counts = {}
    libs = {}

    psp, c = L.build_psp(U, sm, km)
    counts.update(c)
    libs.update(psp)
    L.log(f"PSP: {[(k, len(v)) for k, v in psp.items()]}")

    # ---- excluded PMIDs (compendium + HT-MS)
    omni = L.load_omnipath(U, sm, km)
    ipt = L.load_iptmnet(U, sm, km)
    allp = pd.concat([omni["pmids"], ipt["pmids"]], ignore_index=True)
    cited = pd.Series([p for s in allp for p in s]).value_counts()
    cited = pd.DataFrame({"pmid": cited.index.astype(str), "n_edges": cited.values})
    cited["pmid_int"] = cited["pmid"].astype(int)
    cited = cited.sort_values("pmid_int").drop(columns="pmid_int")
    excluded = L.build_excluded_pmids(cited, fetch=fetch)
    excluded.to_csv(L.RESLIB / "excluded_pmids.tsv", sep="\t", index=False)
    counts["cited_pmids_total"] = int(len(cited))
    counts["cited_pmids_with_pubmed_record"] = int(len(set(cited["pmid"]) & set(L.pubmed_metadata([], fetch=False))))
    for r in ["compendium", "ht_ms_text", "ht_ms_edgecount"]:
        counts[f"excluded_pmids_{r}"] = int(excluded["reason"].str.contains(r).sum())
    counts["excluded_pmids_total"] = int(len(excluded))
    L.log(f"excluded PMIDs: {len(excluded)}")

    cu, audit, _, _ = L.build_curated_union(libs["psp_all"], U, sm, km, excluded, counts)
    libs["curated_union"] = cu
    audit_s = audit.assign(kept=audit["drop_reason"] == "").groupby(["db", "drop_reason"]).size()
    L.log(f"curated_union {len(cu)}\n{audit_s}")
    libs["invitro_sugiyama"] = L.build_sugiyama(U, sm, km, counts)
    L.log(f"sugiyama {len(libs['invitro_sugiyama'])}")

    # ---- predicted
    pred_info = {}
    ksm_sites = L.candidate_sites(sm)
    counts["predicted_candidate_sites_ksm"] = int(ksm_sites["cand"].sum())
    counts["predicted_candidate_sites_in_compendium"] = int(ksm_sites["in_comp"].sum())
    counts["predicted_candidate_sites_benchmarkin_only"] = int((ksm_sites["in_bench"] & ~ksm_sites["in_comp"]).sum())
    for name in L.PREDICTED_LIBS:
        pred_info[name] = L.verify_predicted_file(name)
        libs[name] = L.build_predicted(name, U, sm, ksm_sites=ksm_sites)
        L.log(f"{name}: {len(libs[name])} edges, {libs[name]['kinase'].nunique()} kinases")

    # ---- validate + write
    out_hash = {}
    for name in L.ALL_LIBS:
        lib = libs[name].sort_values(["kinase", "site"] if name not in L.PREDICTED_LIBS else ["kinase", "rank"])
        lib = lib.reset_index(drop=True)
        lib["weight"] = lib["weight"].astype(np.float64)
        validate(name, lib, U)
        lib.to_parquet(L.LIBDIR / f"{name}.parquet", index=False)
        out_hash[name] = L.library_content_hash(lib)
    counts_df = pd.DataFrame([(k, json.dumps(v) if isinstance(v, list) else v) for k, v in counts.items()],
                             columns=["key", "value"])
    counts_df.to_csv(L.RESLIB / "filter_counts.tsv", sep="\t", index=False)
    inputs = {}
    for rel in ["cache/psp_ks_all.parquet", "cache/sites.parquet", "cache/hgnc_symbols.parquet",
                "external/omnipath/omnipath_enzsub.tsv", "external/iptmnet/ptm.txt",
                "external/sugiyama2019/supp_PMC6642169.zip", "external/ksm_index/kinase_list.tsv",
                "external/ksm_index/kinase_synonyms.tsv", "external/ksm_index/kinases.parquet",
                "external/ksm_index/sites.parquet", "external/ksm_index/universe.parquet",
                "external/uniprot/UP000005640_annotations.tsv.gz",
                "external/uniprot/UP000005640_reviewed_isoforms.fasta.gz", "external/hgnc/hgnc_complete_set.txt",
                "data/01_phosphoproteomics_conditions.csv"]:
        inputs[rel] = L.sha256_file(ROOT / rel)
    man = {"inputs_sha256": inputs, "predicted_files": pred_info, "library_content_sha256": out_hash,
           "k_max": L.K_MAX, "n_edges": {k: int(len(v)) for k, v in libs.items()}}
    (L.RESLIB / "build_manifest.json").write_text(json.dumps(man, indent=1, sort_keys=True))
    audit_s.rename("n").reset_index().to_csv(L.RESLIB / "curated_filter_audit.tsv", sep="\t", index=False)
    return libs, counts


# ------------------------------------------------------------------------------------------- report
def _md(df, floatfmt="{:.3g}"):
    cols = list(df.columns)
    lines = ["| " + " | ".join(map(str, cols)) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        cells = []
        for v in r:
            if isinstance(v, (float, np.floating)):
                cells.append(floatfmt.format(v) if np.isfinite(v) else "")
            else:
                cells.append(str(v))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def report():
    U = pd.read_csv(L.RESLIB / "kinase_universe_used.tsv", sep="\t", dtype=str)
    comp = set(pd.read_parquet(L.CACHE / "sites.parquet", columns=["site_key"])["site_key"])
    libs = {n: pd.read_parquet(L.LIBDIR / f"{n}.parquet") for n in L.ALL_LIBS}
    counts = pd.read_csv(L.RESLIB / "filter_counts.tsv", sep="\t", dtype=str).set_index("key")["value"]
    man = json.loads((L.RESLIB / "build_manifest.json").read_text())
    views = {}
    for n in L.CURATED_LIBS:
        views[n] = libs[n]
    for n in L.PREDICTED_LIBS:
        for k in (50, 100, 200, 500):
            views[f"{n}@{k}"] = libs[n][libs[n]["rank"] <= k]

    # 1. per-library summary
    rows = []
    for n, lib in views.items():
        inc = lib[lib["site"].isin(comp)]
        per_k = inc.groupby("kinase").size()
        rows.append({"library": n, "kinases": lib["kinase"].nunique(), "edges": len(lib),
                     "sites": lib["site"].nunique(), "edges_site_in_compendium": len(inc),
                     "sites_in_compendium": inc["site"].nunique(),
                     "kinases_>=1_quant": int((per_k >= 1).sum()), "kinases_>=5_quant": int((per_k >= 5).sum()),
                     "kinases_>=20_quant": int((per_k >= 20).sum())})
    summ = pd.DataFrame(rows)

    # 2. overlap (edge level) between libraries
    ov_names = L.CURATED_LIBS + [f"{n}@100" for n in L.PREDICTED_LIBS]
    sets = {n: set(zip(views[n]["kinase"], views[n]["site"])) for n in ov_names}
    ov = []
    for a in ov_names:
        r = {"library": a}
        for b in ov_names:
            r[b] = len(sets[a] & sets[b]) / max(1, len(sets[a]))
        ov.append(r)
    ov = pd.DataFrame(ov)
    jac = pd.DataFrame([{"A": a, "B": b, "n_A": len(sets[a]), "n_B": len(sets[b]), "n_AB": len(sets[a] & sets[b]),
                         "jaccard": len(sets[a] & sets[b]) / max(1, len(sets[a] | sets[b]))}
                        for i, a in enumerate(ov_names) for b in ov_names[i + 1:]])
    ov.to_csv(L.RESLIB / "overlap.tsv", sep="\t", index=False)
    jac.to_csv(L.RESLIB / "overlap_jaccard.tsv", sep="\t", index=False)

    # 3. predicted top-k vs PSP
    sm = L.SiteMapper()
    cand = L.candidate_sites(sm)
    cand_keys = set(cand.loc[cand["cand"], "site"])
    res_of = {s: s.rsplit("_", 1)[1][0] for s in cand_keys}
    n_st = sum(1 for r in res_of.values() if r in "ST")
    n_y = len(res_of) - n_st
    rc = dict(zip(U["kinase"], U["residue_class"])) if "residue_class" in U else {}
    psp = libs["psp_all"]
    psp = psp[psp["site"].isin(cand_keys)]
    psp_sets = psp.groupby("kinase")["site"].apply(set).to_dict()
    pv = []
    for n in L.PREDICTED_LIBS:
        lib = libs[n]
        for k in (100, 200):
            top = lib[lib["rank"] <= k].groupby("kinase")["site"].apply(set).to_dict()
            for kin in sorted(set(top) & set(psp_sets)):
                P = psp_sets[kin]
                if len(P) < 20:
                    continue
                N = n_y if rc.get(kin) == "Y" else n_st
                P_cls = {s for s in P if (res_of[s] == "Y") == (rc.get(kin) == "Y")}
                T = top[kin]
                x = len(T & P_cls)
                exp = len(T) * len(P_cls) / N
                jac_obs = x / len(T | P_cls)
                jac_rand = exp / (len(T) + len(P_cls) - exp)
                p = hypergeom.sf(x - 1, N, len(P_cls), len(T))
                pv.append({"library": n, "k": k, "kinase": kin, "n_psp_cand": len(P_cls), "n_top": len(T),
                           "overlap": x, "expected_random": exp, "fold": x / exp if exp > 0 else np.nan,
                           "jaccard": jac_obs, "jaccard_random": jac_rand, "p_hypergeom": p})
    pv = pd.DataFrame(pv)
    pv.to_csv(L.RESLIB / "predicted_vs_psp.tsv", sep="\t", index=False)
    pv_sum = pv.groupby(["library", "k"]).agg(kinases=("kinase", "size"), median_jaccard=("jaccard", "median"),
                                             median_jaccard_random=("jaccard_random", "median"),
                                             median_fold=("fold", "median"),
                                             frac_p_lt_0_01=("p_hypergeom", lambda s: float((s < 0.01).mean()))
                                             ).reset_index()
    sanity = pv[(pv["k"] == 200) & pv["kinase"].isin(SANITY_KINASES)][
        ["library", "kinase", "n_psp_cand", "overlap", "expected_random", "fold", "jaccard", "jaccard_random",
         "p_hypergeom"]]
    ws = pv[(pv["k"] == 100) & pv["kinase"].isin(WELL_STUDIED)].pivot(index="kinase", columns="library",
                                                                     values="jaccard").reset_index()

    # 4. leakage filter counts
    audit = pd.read_csv(L.RESLIB / "curated_filter_audit.tsv", sep="\t").fillna({"drop_reason": ""})
    audit["drop_reason"] = audit["drop_reason"].replace("", "KEPT")
    audit = audit.pivot(index="drop_reason", columns="db", values="n").fillna(0).astype(int).reset_index()
    ex = pd.read_csv(L.RESLIB / "excluded_pmids.tsv", sep="\t", dtype=str)
    top_ht = ex[ex["reason"].str.contains("ht_ms")].assign(n=lambda d: d["n_edges_citing"].astype(int)) \
        .sort_values("n", ascending=False).head(15)[["pmid", "reason", "n_edges_citing", "year", "title"]]

    def cnt(k):
        return counts.get(k, "n/a")

    pf = man["predicted_files"]
    lines = [
        "# KS library report", "",
        "Generated by `scripts/build_libraries.py` (do not edit by hand). Libraries: `cache/libraries/<name>.parquet` "
        "(kinase, site, weight, source, evidence [+ rank, score for predicted]); load with "
        "`kinact.library.load_library(name, k=..., gamma=...)`. Kinase universe U: "
        f"{len(U)} kinases ({U['universe_source'].iloc[0]}). Build manifest with input/output sha256: "
        "`results/libraries/build_manifest.json`.", "",
        "## 1. Library sizes", "",
        "`*_quant` = number of kinases with ≥1/5/20 library substrates quantified in the compendium "
        "(cache/sites.parquet; any condition — the DEV-only restriction of N10 is applied at scoring time).", "",
        _md(summ), "",
        "Notes:",
        f"- psp_all: PSP rows with kinase and substrate organism ∈ {{human, mouse, rat}}, non-human substrates mapped "
        f"to the human site via SITE_GRP_ID; {cnt('psp_orthologue_duplicate_rows_collapsed')} orthologue/duplicate rows "
        f"collapsed into {cnt('psp_all_edges_after_orthologue_collapse')} (kinase, human site) edges; "
        f"{cnt('psp_rows_kinase_not_in_U')} rows ({cnt('psp_kinase_names_not_in_U')} kinase names) dropped as kinase ∉ U.",
        "- psp_invivo: ≥1 contributing PSP row flagged IN_VIVO_RXN. psp_human_only: human kinase × human substrate rows only.",
        f"- Predicted libraries: candidate sites = kinaseSubMain site index ∩ (compendium-quantified ∪ benchmarKIN "
        f"external data) = {cnt('predicted_candidate_sites_ksm')} sites ({cnt('predicted_candidate_sites_in_compendium')} "
        f"in compendium, {cnt('predicted_candidate_sites_benchmarkin_only')} benchmarKIN-only). Per kinase the top 500 "
        "(rank by raw score, ties by site key) are stored; weight = within-kinase percentile among candidates "
        "= 1 − (rank−1)/n_cand; raw score + rank + n_cand in `evidence`. Cantley = unprimed published-definition "
        "atlas score (kinases without an own matrix have no rows; no imputation).", "",
        "## 2. Edge overlap between libraries", "",
        "Row library A, column B: fraction of A's edges also in B (predicted at k=100).", "",
        _md(ov, "{:.2f}"), "",
        "Pairwise Jaccard: `results/libraries/overlap_jaccard.tsv`.", "",
        "## 3. Predicted top-k vs PSP (sanity)", "",
        "Kinases with ≥20 psp_all substrates among candidate sites; random expectation from the hypergeometric "
        "(population = candidate sites of the kinase's residue class). R2/R3 were trained on PSP (incl. these pairs), "
        "so their agreement partly reflects training fit; Cantley is label-free and is the cleanest check.", "",
        _md(pv_sum), "",
        "Pre-registered sanity kinases at k=200:", "",
        _md(sanity), "",
        "Jaccard (k=100) for well-studied kinases:", "",
        _md(ws, "{:.3f}"), "",
        "## 4. Leakage filters (curated_union, DESIGN B8)", "",
        f"Excluded PMIDs: {cnt('excluded_pmids_total')} (`results/libraries/excluded_pmids.tsv`): "
        f"{cnt('excluded_pmids_compendium')} compendium PMIDs (178 from 01_*.csv + qPTM 36165955 + Ochoa 2016 "
        f"27909043), {cnt('excluded_pmids_ht_ms_text')} flagged high-throughput MS by the PubMed title/abstract/MeSH "
        f"heuristic (MS evidence AND HT evidence; `kinact.library.classify_ht_ms`), "
        f"{cnt('excluded_pmids_ht_ms_edgecount')} by the edge-count rule (cited by ≥{L.HT_EDGECOUNT_THRESHOLD} "
        f"edges). {cnt('cited_pmids_total')} PMIDs are cited by OmniPath/iPTMnet human phospho edges; "
        f"{cnt('cited_pmids_with_pubmed_record')} had a PubMed record.", "",
        "Edges per filter reason (first applicable reason, in this order: kinase ∉ U → only MIMP/ProtMapper/"
        "text-mining support → cites a compendium PMID → cites an HT-MS PMID → no reference and no whitelisted "
        "curated source):", "",
        _md(audit), "",
        f"- Edges dropped for citing an HT-MS PMID that also cite another (non-excluded) PMID: OmniPath "
        f"{cnt('omnipath_cites_ht_ms_pmid_but_has_other_refs')}, iPTMnet {cnt('iptmnet_cites_ht_ms_pmid_but_has_other_refs')} "
        "(dropped anyway: strict B8 rule).",
        f"- Reference-less edges kept because a whitelisted curated source supports them: OmniPath "
        f"{cnt('omnipath_no_reference_kept_whitelisted')} (whitelist {sorted(L.OMNI_NOREF_WHITELIST)}), iPTMnet "
        f"{cnt('iptmnet_no_reference_kept_whitelisted')} (whitelist {sorted(L.IPTM_NOREF_WHITELIST)}).",
        "- MIMP, ProtMapper (and every `*_MIMP`/`*_ProtMapper` derived source; ProtMapper includes REACH/Sparser/"
        "RLIMS-P text mining) and iPTMnet RLIMS-P (`rlim+`) do not count as support. iPTMnet PhosphoGRID (`pgrd`, yeast) "
        "is ignored. iPTMnet restricted to human substrates.",
        "- PSP edges carry no PMIDs and cannot be checked against compendium studies (r2_data N9); they enter "
        "curated_union unfiltered with evidence `PSP`. No curated_union edge lists a compendium or HT-MS PMID "
        "(tests/test_libraries.py).",
        f"- curated_union = psp_all ∪ cleaned OmniPath ∪ cleaned iPTMnet: {cnt('curated_union_edges')} edges, "
        f"{cnt('curated_union_edges_not_in_psp_all')} not in psp_all.",
        "- Dropped sources: kinaseSubMain `ikipdb_zenodo_*` (\"iKiP-DB\"; mislabelled cell-line atlas files, Table7 "
        "is iKiP∪PSP — r2_data N3); SELPHI2, PTMsigDB, kinaseSubMain release/clean_bench (banned, DESIGN §2); "
        "Phosformer (B17); NetworKIN (prediction). SIGNOR enters via OmniPath/iPTMnet (no separate download).", "",
        "Top HT-MS PMIDs by number of citing edges:", "",
        _md(top_ht), "",
        "## 5. In vitro (Sugiyama 2019) — separate library `invitro_sugiyama`", "",
        f"Sci Rep 9:10503 (PMID {L.SUGIYAMA_PMID}), Table S2: {cnt('sugiyama_rows_total')} rows → "
        f"{cnt('sugiyama_rows_wildtype')} wild-type kinase → {cnt('sugiyama_rows_protein_kinase')} protein kinase → "
        f"{cnt('sugiyama_rows_sidic_confirmed')} SIDIC-confirmed; residue mismatches vs UniProt canonical dropped: "
        f"{cnt('sugiyama_rows_residue_mismatch_dropped')}; unmapped entry names {cnt('sugiyama_rows_entry_name_unmapped')}; "
        f"Carna kinase names not mapped to U: {cnt('sugiyama_kinase_names_unmapped_or_not_in_U')}. "
        "Provenance clean: in vitro kinase reactions on HeLa-lysate-derived proteins/peptides, not a compendium "
        "study, CC BY 4.0. Kept out of curated_union (in vitro substrates hurt as labels per Müller-Dott 2025).", "",
        "## 6. Predicted-file provenance (kinaseSubMain, docs/kinasesub_audit.md)", "",
        _md(pd.DataFrame([{"library": k, "origin": v["origin"], "sha256": v["sha256"][:16] + "…",
                           "check": v["manifest_check"]} for k, v in pf.items()])), "",
        "Each file's sha256 is asserted equal to the audited value and absent from the banned list (kinaseSubMain "
        "release_scores, r2_deploy, r2r3_release outputs, clean_bench), and paths are checked for banned tokens "
        "(`kinact.library.verify_predicted_file`). R2-nophospho: feature groups F1,F0,F2a,F3,F4,F5,F7 (no F6 "
        "phospho/inhibitor features), DEV-only PSP labels. R3 deploy: ESM-2 windows + protein token + U_ms bit, "
        "DEV+TEST PSP labels. Residual leakage for both: site existence (U_ms includes compendium sites) — "
        "carries no perturbation sign/magnitude. Model selection in kinaseSubMain used PSP-label CV (macro "
        "hard-negative AUROC), not Hijazi/benchmarKIN perturbation results (audit §1–2).", "",
        "## 7. Caveats", "",
        "- 'Quantified in the compendium' uses all compendium conditions (incl. TEST) to define the candidate site "
        "set of predicted libraries; this is site existence only (r2_data N10). Methods must additionally restrict "
        "to sites measured in the scored condition / DEV cache (`load_library(..., sites=...)`).",
        "- iPTMnet licence: its readme states both CC BY-NC-SA 4.0 (download terms) and CC BY 4.0; we apply the "
        "more restrictive CC BY-NC-SA 4.0 (compatible with our non-commercial use; same terms as PSP). Derived "
        "libraries are not committed.",
    ]
    (L.RESLIB / "LIBRARY_REPORT.md").write_text("\n".join(lines) + "\n")
    L.log("report written")
    return summ, pv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--report-only", action="store_true")
    ap.add_argument("--no-fetch", action="store_true")
    a = ap.parse_args()
    if not a.report_only:
        build(fetch=not a.no_fetch)
    report()


if __name__ == "__main__":
    main()
