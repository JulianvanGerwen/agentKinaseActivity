"""Kinase–substrate (KS) library builder and loader — the ONE shared library builder (DESIGN §2, §7A, B8, B17).

Library files: ``cache/libraries/<name>.parquet`` with columns
    kinase   HGNC symbol in U (kinase universe)
    site     canonical site key ``{UniProt}_{res}{pos}`` (compendium key from cache/sites.parquet when the
             site is in the compendium, otherwise the canonical UniProt key)
    weight   float > 0 (1.0 for curated; within-kinase percentile in (0, 1] for predicted)
    source   ';'-joined contributing sources
    evidence free text (PMIDs, organisms, in vivo/vitro flags, raw predicted score)
  predicted libraries additionally: rank (1 = best, within kinase), score (raw model score).

Libraries
    psp_all          PSP KS, kinase organism human/mouse/rat, substrate human or mouse/rat mapped to human via
                     SITE_GRP_ID; orthologue duplicates collapsed to one (kinase, human site) edge
    psp_invivo       psp_all edges with >=1 in-vivo PSP row
    psp_human_only   human kinase + human substrate rows only
    curated_union    psp_all ∪ leakage-filtered OmniPath ∪ iPTMnet (see ``filter_*``)
    invitro_sugiyama Sugiyama 2019 (Sci Rep, PMID 31324866) in vitro HeLa-lysate kinase reactions
    r2np, r3, cantley top-500 predicted sites per kinase among sites quantified in the compendium or in the
                     benchmarKIN external data (external/benchmarkin_py); use ``load_library(name, k=...)``.

Usage: ``from kinact.library import load_library; lib = load_library('psp_all')``.
Build: ``python scripts/build_libraries.py``.
"""
from __future__ import annotations

import gzip
import os
import hashlib
import json
import re
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
CACHE = ROOT / "cache"
LIBDIR = CACHE / "libraries"
EXT = ROOT / "external"
RESLIB = ROOT / "results" / "libraries"
PUBMED_CACHE = CACHE / "pubmed" / "pubmed_meta.jsonl"

LIB_COLUMNS = ["kinase", "site", "weight", "source", "evidence"]
CURATED_LIBS = ["psp_all", "psp_invivo", "psp_human_only", "curated_union", "invitro_sugiyama"]
PREDICTED_LIBS = ["r2np", "r3", "cantley"]
ALL_LIBS = CURATED_LIBS + PREDICTED_LIBS
K_MAX = 500
SITE_KEY_RE = re.compile(r"^(?P<acc>[A-Z0-9]+(?:-\d+)?)_(?P<res>[STY])(?P<pos>\d+)$")

# Papers that are the source of compendium conditions but are not in its PMID column (01_*.csv `source`):
# qPTM (Yu et al. 2023 NAR) and Ochoa et al. 2016 MSB ("An atlas of human kinase regulation");
# resolved with PubMed esearch on 2026-10-27.
EXTRA_COMPENDIUM_PMIDS = {"36165955": "qPTM database (source of 461 compendium conditions)",
                          "27909043": "Ochoa 2016 MSB atlas (source of 39 compendium conditions)"}
SUGIYAMA_PMID = "31324866"

# Predicted score files (copied from kinaseSubMain, READ-ONLY origin) and the banned kinaseSubMain variants.
KSM = Path(os.environ.get("KINASESUB_DIR", "kinaseSubMain"))  # sibling kinaseSubMain checkout
PREDICTED_FILES = {
    "r2np": {"path": EXT / "ksm_pred/r2np/uobs_predictions.parquet",
             "origin": "cache/runs/r2_final_nophospho_uobs__20261026-213950__s0/uobs_predictions.parquet",
             "sha256": "1edbd2691e3c2d638f449e1a423e97088c6d1299fc791a1eaf10f5d2aca2c052",
             "manifest": EXT / "ksm_pred/r2np/manifest.json"},
    "r3": {"path": EXT / "ksm_pred/r3/uobs_predictions.parquet",
           "origin": "cache/release/r3_deploy_20261027-005816/uobs_predictions.parquet",
           "sha256": "2f62988cf798f15ace6cd9925dce77e45ca2c6cf297549203734d256e9d002ba",
           "manifest": EXT / "ksm_pred/r3/README.json"},
    "cantley": {"path": EXT / "ksm_pred/cantley/cantley_unprimed_score.npy",
                "origin": "cache/features/motif/cantley_unprimed_score.npy",
                "sha256": "e5597327729ab3605ca1649d510d89fbeffb0af1c416cd1fe96c6527bd4da2dc",
                "manifest": EXT / "ksm_pred/cantley/cantley_unprimed_score.json"},
}
# sha256 of kinaseSubMain outputs that must never be used (DESIGN §2 banned list; from their manifests).
BANNED_SHA256 = {
    "cceb32c73be9c092923e8d123bf8dd74a676afe2bea845e002f33f296de839f5": "release r2r3_release release_scores.parquet",
    "10e25be1127e4d97d1a2b215e9525fc6c1c08a4f3b97f69ebeddd11b4a5b6e76": "release r2_deploy r2_uobs_scores.parquet (full F6)",
    "ff5bfe80c034ffe87a04adc7af89962f0826fb2cf0116fe6ac7d3dc57e716073": "r2_final_cleanbench uobs_predictions (F6 minus benchmarKIN)",
    "98f348ed330784c8dbe52f9d327a14226975a74ab07b1e7612ebb6cf58d64567": "r2r3_release output",
    "61e20f51b4eb043365ef492a0f20c31168f30bc6bacd637d40154fd0a8c3729f": "r2r3_release output",
}
BANNED_PATH_TOKENS = ["r2r3_release", "r2_deploy", "cleanbench", "clean_bench", "release_scores",
                      "kinase_top200", "site_top5", "ikipdb", "selphi", "f6clean"]

# OmniPath: sources that are not curation (prediction / text mining) — edges supported only by these are dropped.
OMNI_NONCURATED = {"MIMP", "ProtMapper"}  # plus every '<X>_MIMP' / '<X>_ProtMapper' derived source
# Sources whose edges are accepted WITHOUT a reference (manually curated or in vitro, B8/N2).
OMNI_NOREF_WHITELIST = {"PhosphoSite", "SIGNOR", "phosphoELM", "HPRD", "PhosphoNetworks", "Li2012"}
# iPTMnet: rlim+ = RLIMS-P text mining (dropped); pgrd = PhosphoGRID (yeast). Whitelist for no-reference edges.
IPTM_TEXTMINING = {"rlim+", "rlimsp", "efip"}
IPTM_NOREF_WHITELIST = {"sign", "hprd", "unip", "pro", "npro", "psp", "pelm"}
HT_EDGECOUNT_THRESHOLD = 150  # r2_data N1 operational rule: a paper cited by >=150 edges is high-throughput


def log(msg, t0=[time.time()]):
    print(f"[library {time.time() - t0[0]:7.1f}s] {msg}", flush=True)


def sha256_file(path, chunk=1 << 24):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


# --------------------------------------------------------------------------------------- site keys
def parse_site_key(key: str):
    """'P28482_T185' -> ('P28482', 'T', 185); raises ValueError if malformed."""
    m = SITE_KEY_RE.match(str(key))
    if not m:
        raise ValueError(f"bad site key {key!r}")
    return m["acc"], m["res"], int(m["pos"])


def valid_site_keys(keys) -> np.ndarray:
    return pd.Series(keys, dtype="string").str.fullmatch(SITE_KEY_RE.pattern).fillna(False).to_numpy()


class SiteMapper:
    """(acc, res, pos) -> canonical site key. Compendium key if the site is in the compendium (exact key, or via
    kinaseSubMain's universe comp_uid isoform/shift mapping), else the UniProt key."""

    def __init__(self):
        s = pd.read_parquet(CACHE / "sites.parquet", columns=["peptide_UID", "site_key"])
        self.comp_keys = set(s["site_key"])
        self.uid2key = dict(zip(s["peptide_UID"], s["site_key"]))
        u = pd.read_parquet(EXT / "ksm_index/universe.parquet", columns=["site_id", "comp_uid"])
        u = u[u["comp_uid"].notna()]
        ck = u["comp_uid"].map(self.uid2key)
        self.ksm2comp = dict(zip(u["site_id"], ck))

    def key(self, acc, res, pos):
        acc = pd.Series(acc, dtype="string").reset_index(drop=True)
        res = pd.Series(res, dtype="string").reset_index(drop=True)
        pos = pd.Series(pos).reset_index(drop=True).astype("Int64").astype("string")
        k = acc + "_" + res + pos
        out = k.copy()
        not_comp = ~k.isin(self.comp_keys).fillna(False)
        mapped = k[not_comp].map(self.ksm2comp)
        out.loc[not_comp & mapped.reindex(k.index).notna()] = mapped.dropna()
        return out.to_numpy(dtype=object)

    def from_keys(self, keys):
        keys = pd.Series(keys, dtype="string")
        m = keys.map(self.ksm2comp)
        out = keys.where(keys.isin(self.comp_keys) | m.isna(), m)
        return out.to_numpy(dtype=object)

    def in_compendium(self, keys):
        return pd.Series(keys).isin(self.comp_keys).to_numpy()


# --------------------------------------------------------------------------------------- kinases
def load_universe() -> pd.DataFrame:
    """U: results/ground_truth/kinase_universe.tsv if present, else kinaseSubMain include==True list (503),
    HGNC-normalised with external/hgnc (approved symbol of the UniProt accession; fallback symbol lookup)."""
    f = ROOT / "results/ground_truth/kinase_universe.tsv"
    if f.exists():
        u = pd.read_csv(f, sep="\t", dtype=str)
        u["universe_source"] = "results/ground_truth/kinase_universe.tsv"
        if "residue_class" not in u.columns:
            k = pd.read_csv(EXT / "ksm_index/kinase_list.tsv", sep="\t", dtype=str)
            u["residue_class"] = u["uniprot"].map(dict(zip(k["uniprot_acc"], k["residue_class"])))
        return u[["kinase", "uniprot"] + [c for c in u.columns if c not in ("kinase", "uniprot")]]
    k = pd.read_csv(EXT / "ksm_index/kinase_list.tsv", sep="\t", dtype=str)
    k = k[k["include"] == "True"].copy()
    h = pd.read_csv(EXT / "hgnc/hgnc_complete_set.txt", sep="\t", dtype=str,
                    usecols=["symbol", "uniprot_ids", "prev_symbol", "alias_symbol"])
    acc2sym = {}
    for sym, ids in zip(h["symbol"], h["uniprot_ids"].fillna("")):
        for a in ids.split("|"):
            if a:
                acc2sym.setdefault(a, set()).add(sym)
    approved = set(h["symbol"])
    syms, how = [], []
    for g, a in zip(k["gene"], k["uniprot_acc"]):
        s = acc2sym.get(a, set())
        if len(s) == 1:
            syms.append(next(iter(s))); how.append("hgnc_by_uniprot")
        elif g in approved:
            syms.append(g); how.append("hgnc_symbol")
        else:
            syms.append(g); how.append("kept_unmapped")
    u = pd.DataFrame({"kinase": syms, "uniprot": k["uniprot_acc"].values, "ksm_gene": k["gene"].values,
                      "hgnc_match": how, "manning_group": k["manning_group"].values,
                      "manning_family": k["manning_family"].values, "residue_class": k["residue_class"].values,
                      "is_pseudokinase": k["is_pseudokinase"].values})
    u["universe_source"] = "kinaseSubMain results/labels/kinase_list.tsv include==True (HGNC-normalised)"
    assert u["kinase"].is_unique and u["uniprot"].is_unique and len(u) == 503
    return u.sort_values("kinase").reset_index(drop=True)


class KinaseMapper:
    """Map kinase accessions / names (any organism, PSP/Carna/Manning names) to U symbols."""

    def __init__(self, U: pd.DataFrame):
        self.U = U
        self.symbols = set(U["kinase"])
        self.acc2k = dict(zip(U["uniprot"], U["kinase"]))
        h = pd.read_parquet(CACHE / "hgnc_symbols.parquet")
        pri = {"symbol": 0, "prev_symbol": 1, "alias_symbol": 2}
        h = h.assign(p=h["kind"].map(pri))
        self.name2k = {}
        for q, g in h.groupby("query_upper"):
            best = g[g["p"] == g["p"].min()]["symbol"].unique()
            inU = [s for s in best if s in self.symbols]
            if len(best) == 1 and best[0] in self.symbols:
                self.name2k[q] = best[0]
            elif len(inU) == 1 and len(best) == 1:
                self.name2k[q] = inU[0]
        # kinaseSubMain synonym table (PSP / Cantley / Manning / KinHub / UniProt names) keyed by accession
        syn = pd.read_csv(EXT / "ksm_index/kinase_synonyms.tsv", sep="\t", dtype=str)
        syn = syn[syn["use_for_join"] == "True"]
        self.syn2k = {}
        for s, a in zip(syn["synonym_uc"], syn["uniprot_acc"]):
            if a in self.acc2k:
                self.syn2k.setdefault(s, set()).add(self.acc2k[a])
        self.syn2k = {s: next(iter(v)) for s, v in self.syn2k.items() if len(v) == 1}
        for s in self.symbols:  # U symbols always win
            self.name2k[s.upper()] = s

    def by_acc(self, acc):
        if acc is None or (isinstance(acc, float) and np.isnan(acc)):
            return None
        acc = str(acc)
        return self.acc2k.get(acc) or self.acc2k.get(acc.split("-")[0])

    def by_name(self, name):
        if name is None or (isinstance(name, float) and np.isnan(name)):
            return None
        q = str(name).upper().strip()
        return self.name2k.get(q) or self.syn2k.get(q)


# --------------------------------------------------------------------------------------- PMIDs
def compendium_pmids() -> set:
    c = pd.read_csv(ROOT / "data/01_phosphoproteomics_conditions.csv", dtype=str)
    p = set()
    for v in c["PMID"].dropna():
        p |= set(re.findall(r"\d{6,9}", v))
    return p | set(EXTRA_COMPENDIUM_PMIDS)


def _efetch(pmids, tries=4):
    data = urllib.parse.urlencode({"db": "pubmed", "id": ",".join(pmids), "retmode": "xml",
                                   "tool": "kinact_library"}).encode()
    for t in range(tries):
        try:
            with urllib.request.urlopen("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                                        data=data, timeout=120) as r:
                return r.read()
        except Exception as e:  # noqa: BLE001
            log(f"efetch retry {t}: {e}")
            time.sleep(2 + 3 * t)
    raise RuntimeError("efetch failed")


def _parse_pubmed_xml(xml_bytes):
    out = []
    root = ET.fromstring(xml_bytes)
    for art in root.iter("PubmedArticle"):
        pmid = art.findtext(".//MedlineCitation/PMID")
        title = "".join(art.find(".//ArticleTitle").itertext()) if art.find(".//ArticleTitle") is not None else ""
        abstract = " ".join("".join(a.itertext()) for a in art.findall(".//Abstract/AbstractText"))
        mesh = [m.text for m in art.findall(".//MeshHeading/DescriptorName") if m.text]
        ptypes = [p.text for p in art.findall(".//PublicationType") if p.text]
        year = art.findtext(".//JournalIssue/PubDate/Year") or art.findtext(".//JournalIssue/PubDate/MedlineDate") or ""
        journal = art.findtext(".//Journal/ISOAbbreviation") or ""
        out.append({"pmid": pmid, "title": title, "abstract": abstract, "mesh": mesh, "ptypes": ptypes,
                    "year": year[:4], "journal": journal})
    return out


def pubmed_metadata(pmids, fetch=True) -> dict:
    """PubMed title/abstract/MeSH per PMID, cached in cache/pubmed/pubmed_meta.jsonl (fetched once)."""
    meta = {}
    if PUBMED_CACHE.exists():
        with open(PUBMED_CACHE) as f:
            for line in f:
                d = json.loads(line)
                meta[d["pmid"]] = d
    missing = sorted(set(pmids) - set(meta) - _known_missing())
    if missing and fetch:
        log(f"fetching PubMed metadata for {len(missing)} PMIDs")
        PUBMED_CACHE.parent.mkdir(parents=True, exist_ok=True)
        got_all = set()
        with open(PUBMED_CACHE, "a") as f:
            for i in range(0, len(missing), 200):
                batch = missing[i:i + 200]
                recs = _parse_pubmed_xml(_efetch(batch))
                for d in recs:
                    if d["pmid"] in set(batch) and d["pmid"] not in meta:
                        meta[d["pmid"]] = d
                        f.write(json.dumps(d) + "\n")
                        got_all.add(d["pmid"])
                time.sleep(0.4)
        nf = sorted(set(missing) - got_all)
        (PUBMED_CACHE.parent / "not_found.txt").write_text("\n".join(sorted(_known_missing() | set(nf))) + "\n")
    return meta


def _known_missing():
    f = PUBMED_CACHE.parent / "not_found.txt"
    return set(f.read_text().split()) if f.exists() else set()


_HT_RE = re.compile(
    r"phosphoproteom|large[- ]scale|proteome[- ]wide|kinome[- ]wide|high[- ]throughput|"
    r"systems?[- ]wide|global(?:ly)?\s+(?:phospho|analysis|quantitative|profil|survey|identification|mapping)|"
    r"in[- ]depth\s+(?:phospho|proteom|analysis of the)|comprehensive\s+(?:phospho|proteom|catalog|survey|identification of phospho)|"
    r"quantitative\s+(?:phospho)?proteomic|thousands?\s+of\s+(?:phospho|sites|proteins|peptides)|"
    r"(?<![\w-])\d[\d,]{2,}\s+(?:(?:unique|distinct|different|novel|non-redundant|phosphorylation|phosphorylated|in vivo|class I)\s+){0,2}"
    r"(?:sites|phosphosites|phosphopeptides|phospho-sites|proteins|phosphoproteins|peptides)\b|"
    r"phosphopeptide enrichment|SILAC|iTRAQ|\bTMT\b|tandem mass tag|label[- ]free quantif", re.I)
_MS_RE = re.compile(r"mass spectromet|\bLC-?MS|MS/MS|\bMS\b|proteomic|SILAC|iTRAQ|\bTMT\b|tandem mass tag|orbitrap|"
                    r"\bIMAC\b|\bTiO2\b|titanium dioxide", re.I)
_MS_MESH = {"Mass Spectrometry", "Tandem Mass Spectrometry", "Proteomics", "Spectrometry, Mass, Electrospray Ionization",
            "Spectrometry, Mass, Matrix-Assisted Laser Desorption-Ionization", "Chromatography, Liquid",
            "Isotope Labeling", "Proteome", "Phosphoproteins"}
_HT_MESH = {"Proteomics", "Proteome", "High-Throughput Screening Assays"}


def classify_ht_ms(meta: dict) -> tuple[bool, str]:
    """Pre-declared heuristic: a paper is high-throughput MS if (MS evidence) AND (HT evidence), where MS evidence =
    MS terms in title/abstract or MS MeSH, HT evidence = HT terms in title/abstract (phosphoproteome, large-scale,
    global/quantitative proteomic, SILAC/iTRAQ/TMT, '>=100 sites/proteins/phosphopeptides', ...)."""
    text = f"{meta.get('title', '')} {meta.get('abstract', '')}"
    mesh = set(meta.get("mesh", []))
    ms_t = _MS_RE.search(text)
    ms_m = mesh & _MS_MESH
    ht_t = _HT_RE.search(text)
    ht_m = mesh & _HT_MESH
    is_ms = bool(ms_t) or bool(ms_m - {"Phosphoproteins", "Proteome"})
    is_ht = bool(ht_t)  # MeSH 'Proteomics' alone is too broad (assigned to many single-protein MS studies)
    why = []
    if ms_t:
        why.append(f"ms_text:{ms_t.group(0)}")
    if ms_m:
        why.append("ms_mesh:" + "|".join(sorted(ms_m)))
    if ht_t:
        why.append(f"ht_text:{ht_t.group(0)}")
    if ht_m:
        why.append("ht_mesh:" + "|".join(sorted(ht_m)))
    return bool(is_ms and is_ht), ";".join(why)


def build_excluded_pmids(cited: pd.DataFrame, fetch=True) -> pd.DataFrame:
    """cited: pmid, n_edges (number of curated non-PSP phospho edges citing it, OmniPath+iPTMnet, human, any kinase).
    Returns excluded PMIDs with reason ∈ {compendium, ht_ms_text, ht_ms_edgecount}."""
    comp = compendium_pmids()
    meta = pubmed_metadata(cited["pmid"].tolist(), fetch=fetch)
    rows = []
    for pm, n in zip(cited["pmid"], cited["n_edges"]):
        m = meta.get(pm)
        reasons = []
        if pm in comp:
            reasons.append("compendium")
        why = ""
        if m is not None:
            is_ht, why = classify_ht_ms(m)
            if is_ht:
                reasons.append("ht_ms_text")
        # edge-count rule (r2_data N1) only for papers with MS evidence (or no PubMed record): many-edge
        # low-throughput papers (e.g. tau/GSK3 biochemistry) are not HT-MS
        if n >= HT_EDGECOUNT_THRESHOLD and (m is None or "ms_" in why):
            reasons.append("ht_ms_edgecount")
        if reasons:
            rows.append({"pmid": pm, "reason": ";".join(reasons), "n_edges_citing": int(n),
                         "year": (m or {}).get("year", ""), "journal": (m or {}).get("journal", ""),
                         "title": ((m or {}).get("title", "") or "")[:200].replace("\t", " "),
                         "evidence": why if m is not None else "no_pubmed_record"})
    # compendium PMIDs not cited by any edge are listed too (complete exclusion list)
    for pm in sorted(comp - set(cited["pmid"])):
        rows.append({"pmid": pm, "reason": "compendium", "n_edges_citing": 0, "year": "", "journal": "",
                     "title": EXTRA_COMPENDIUM_PMIDS.get(pm, ""), "evidence": "compendium PMID"})
    ex = pd.DataFrame(rows)
    ex["pmid_int"] = ex["pmid"].astype(int)
    return ex.sort_values("pmid_int").drop(columns="pmid_int").reset_index(drop=True)


# --------------------------------------------------------------------------------------- PSP
def _collapse(df, extra_cols=()):
    """Collapse rows to one edge per (kinase, site)."""
    agg = {"source": lambda s: ";".join(sorted(set(";".join(s).split(";")))),
           "evidence": lambda s: " | ".join(sorted(set(s)))}
    for c in extra_cols:
        agg[c] = "max"
    out = df.groupby(["kinase", "site"], sort=True).agg(agg).reset_index()
    out["weight"] = 1.0
    return out


def build_psp(U, sm: SiteMapper, km: KinaseMapper):
    ks = pd.read_parquet(CACHE / "psp_ks_all.parquet")
    n0 = len(ks)
    counts = {"psp_rows_total": n0}
    orgs = {"human", "mouse", "rat"}
    ks = ks[ks["KIN_ORGANISM"].isin(orgs) & ks["SUB_ORGANISM"].isin(orgs)].copy()
    counts["psp_rows_kin_and_sub_human_mouse_rat"] = len(ks)
    ks = ks[ks["hs_acc"].notna() & ks["hs_pos"].notna() & ks["hs_residue"].isin(["S", "T", "Y"])]
    counts["psp_rows_with_human_site"] = len(ks)
    # kinase -> U: human kinase accession first, otherwise symbol (mouse/rat genes upper-cased)
    kin = [km.by_acc(a) if o == "human" else None for a, o in zip(ks["KIN_ACC_ID"], ks["KIN_ORGANISM"])]
    kin = [k if k is not None else km.by_name(g) for k, g in zip(kin, ks["GENE"])]
    ks["kinase"] = kin
    counts["psp_rows_kinase_not_in_U"] = int(ks["kinase"].isna().sum())
    counts["psp_kinase_names_not_in_U"] = int(ks.loc[ks["kinase"].isna(), "GENE"].str.upper().nunique())
    ks = ks[ks["kinase"].notna()].copy()
    site = sm.key(ks["hs_acc"], ks["hs_residue"], ks["hs_pos"])
    comp_uid = ks["peptide_UID"].astype("string").to_numpy()
    site_from_uid = pd.Series(comp_uid).map(sm.uid2key).to_numpy(dtype=object)
    ks["site"] = np.where(pd.notna(site_from_uid), site_from_uid, site)
    bad = ~valid_site_keys(ks["site"])  # non-UniProt substrate ids (RefSeq NP_, VAR_ variants)
    counts["psp_rows_non_uniprot_substrate_dropped"] = int(bad.sum())
    ks = ks[~bad].copy()
    ks["vivo"] = ks["IN_VIVO_RXN"].fillna("").str.strip().eq("X")
    ks["vitro"] = ks["IN_VITRO_RXN"].fillna("").str.strip().eq("X")
    ks["human_pair"] = (ks["KIN_ORGANISM"] == "human") & (ks["SUB_ORGANISM"] == "human")
    ks["source"] = "PSP"
    ks["evidence"] = ("PSP:" + ks["KIN_ORGANISM"] + "_kin/" + ks["SUB_ORGANISM"] + "_sub;"
                      + np.where(ks["vivo"], "vivo", "") + np.where(ks["vitro"], "+vitro", ""))
    counts["psp_rows_used"] = len(ks)
    lib_all = _collapse(ks, extra_cols=("vivo", "vitro", "human_pair"))
    counts["psp_all_edges_after_orthologue_collapse"] = len(lib_all)
    counts["psp_orthologue_duplicate_rows_collapsed"] = int(len(ks) - len(lib_all))
    libs = {
        "psp_all": lib_all,
        "psp_invivo": lib_all[lib_all["vivo"]].copy(),
        "psp_human_only": _collapse(ks[ks["human_pair"]]),
    }
    for k in libs:
        libs[k] = libs[k][LIB_COLUMNS]
    return libs, counts


# --------------------------------------------------------------------------------------- OmniPath / iPTMnet
def load_omnipath(U, sm, km):
    o = pd.read_csv(EXT / "omnipath/omnipath_enzsub.tsv", sep="\t", dtype=str)
    o = o[o["modification"] == "phosphorylation"].copy()
    o["kinase"] = [km.by_acc(a) or km.by_name(g) for a, g in zip(o["enzyme"], o["enzyme_genesymbol"])]
    o = o[o["residue_type"].isin(["S", "T", "Y"])]
    o["site"] = sm.key(o["substrate"], o["residue_type"], o["residue_offset"].astype(int))
    o["src_set"] = o["sources"].fillna("").str.split(";").map(lambda x: frozenset(s for s in x if s))
    o["pmids"] = o["references"].fillna("").map(lambda x: frozenset(re.findall(r":(\d+)", x)))
    o["ref_text"] = o["references"].fillna("")
    return o[["kinase", "enzyme_genesymbol", "site", "src_set", "pmids", "ref_text"]].rename(
        columns={"enzyme_genesymbol": "enzyme_name"}).reset_index(drop=True)


def load_iptmnet(U, sm, km):
    p = pd.read_csv(EXT / "iptmnet/ptm.txt", sep="\t", header=None, dtype=str,
                    names=["type", "source", "sub", "subg", "org", "site", "enz", "enzg", "note", "pmid"])
    p = p[(p["type"] == "PHOSPHORYLATION") & p["enz"].notna() & p["org"].str.startswith("Homo sapiens", na=False)].copy()
    m = p["site"].str.extract(r"^([STY])(\d+)$")
    p = p[m[0].notna()].copy()
    m = m[m[0].notna()]
    p["kinase"] = [km.by_acc(a) or km.by_name(g) for a, g in zip(p["enz"], p["enzg"])]
    acc = p["sub"].str.replace(r"-1$", "", regex=True)  # '-1' = canonical isoform
    p["site_key"] = sm.key(acc, m[0], m[1].astype(int))
    p["pmids"] = p["pmid"].fillna("").map(lambda x: frozenset(re.findall(r"\d{5,9}", x)))
    g = p.groupby(["kinase", "enzg", "site_key"], dropna=False).agg(
        src_set=("source", lambda s: frozenset(s)), pmids=("pmids", lambda s: frozenset().union(*s)),
        ref_text=("pmid", lambda s: ",".join(sorted(set(",".join(s.dropna()).split(",")) - {""})))).reset_index()
    return g.rename(columns={"site_key": "site", "enzg": "enzyme_name"})


def filter_curated(df: pd.DataFrame, db: str, excluded: pd.DataFrame, counts: dict) -> pd.DataFrame:
    """Leakage filter for OmniPath / iPTMnet edges (B8). Order of reasons (first that applies is counted):
    kinase_not_in_U, noncurated_only (MIMP/ProtMapper/text mining), cites_compendium_pmid, cites_ht_ms_pmid,
    no_reference_not_whitelisted. Mutates `counts`."""
    comp = set(excluded.loc[excluded["reason"].str.contains("compendium"), "pmid"])
    htms = set(excluded.loc[excluded["reason"].str.contains("ht_ms"), "pmid"])
    if db == "omnipath":
        def curated(srcs):
            return frozenset(s for s in srcs if s not in OMNI_NONCURATED
                             and not s.endswith("_MIMP") and not s.endswith("_ProtMapper"))
        white = OMNI_NOREF_WHITELIST
    else:
        def curated(srcs):
            return frozenset(s for s in srcs if s not in IPTM_TEXTMINING and s != "pgrd")
        white = IPTM_NOREF_WHITELIST
    df = df.copy()
    df["cur_src"] = df["src_set"].map(curated)
    reason = pd.Series("", index=df.index, dtype=object)
    reason[df["kinase"].isna()] = "kinase_not_in_U"
    r = reason == ""
    reason[r & (df["cur_src"].map(len) == 0)] = "noncurated_only(MIMP/ProtMapper/text-mining)"
    r = reason == ""
    reason[r & df["pmids"].map(lambda s: len(s & comp) > 0)] = "cites_compendium_pmid"
    r = reason == ""
    reason[r & df["pmids"].map(lambda s: len(s & htms) > 0)] = "cites_ht_ms_pmid"
    r = reason == ""
    reason[r & (df["pmids"].map(len) == 0) & df["cur_src"].map(lambda s: len(s & white) == 0)] = \
        "no_reference_not_whitelisted"
    df["drop_reason"] = reason
    vc = reason.replace("", "kept").value_counts()
    counts[f"{db}_edges_total"] = int(len(df))
    for k, v in vc.items():
        counts[f"{db}_{k}"] = int(v)
    # informative: HT-MS edges that also had a non-excluded reference
    ht = df["drop_reason"] == "cites_ht_ms_pmid"
    counts[f"{db}_cites_ht_ms_pmid_but_has_other_refs"] = int(df.loc[ht, "pmids"].map(lambda s: len(s - htms - comp) > 0).sum())
    counts[f"{db}_no_reference_kept_whitelisted"] = int(((df["drop_reason"] == "") & (df["pmids"].map(len) == 0)).sum())
    kept = df[df["drop_reason"] == ""].copy()
    tag = "OmniPath" if db == "omnipath" else "iPTMnet"
    kept["source"] = kept["cur_src"].map(lambda s: ";".join(f"{tag}:{x}" for x in sorted(s)))
    kept["evidence"] = kept["pmids"].map(lambda s: "PMID:" + ",".join(sorted(s, key=int)) if s else "PMID:none")
    return kept, df


def build_curated_union(psp_all, U, sm, km, excluded, counts):
    omni = load_omnipath(U, sm, km)
    ipt = load_iptmnet(U, sm, km)
    o_kept, o_all = filter_curated(omni, "omnipath", excluded, counts)
    i_kept, i_all = filter_curated(ipt, "iptmnet", excluded, counts)
    parts = [psp_all.assign(evidence="PSP")[LIB_COLUMNS[:2] + ["source", "evidence"]],
             o_kept[["kinase", "site", "source", "evidence"]], i_kept[["kinase", "site", "source", "evidence"]]]
    u = pd.concat(parts, ignore_index=True)
    u = u[valid_site_keys(u["site"])]

    def merge_ev(s):
        pm = set()
        for e in s:
            if e.startswith("PMID:") and e != "PMID:none":
                pm |= set(e[5:].split(","))
        has_psp = any(e == "PSP" for e in s)
        return ";".join((["PSP"] if has_psp else []) + (["PMID:" + ",".join(sorted(pm, key=int))] if pm else []))\
            or "no_reference(whitelisted source)"
    lib = u.groupby(["kinase", "site"], sort=True).agg(
        source=("source", lambda s: ";".join(sorted(set(";".join(s).split(";"))))),
        evidence=("evidence", merge_ev)).reset_index()
    lib["weight"] = 1.0
    counts["curated_union_edges"] = len(lib)
    counts["curated_union_edges_not_in_psp_all"] = int(len(lib) - len(psp_all))
    # leakage audit table (edge-level reasons) for the report
    audit = pd.concat([o_all.assign(db="omnipath"), i_all.assign(db="iptmnet")], ignore_index=True)
    return lib[LIB_COLUMNS], audit, omni, ipt


# --------------------------------------------------------------------------------------- Sugiyama 2019
SUGIYAMA_MANUAL = {  # Carna Biosciences names not resolved by synonym tables (checked against HGNC / Carna)
    "ACK": "TNK2", "ALK1": "ACVRL1", "ALK2": "ACVR1", "ALK4": "ACVR1B", "AMPKA1": "PRKAA1", "AMPKA2": "PRKAA2",
    "AURA": "AURKA", "AURB": "AURKB", "AURC": "AURKC", "BARK1": "GRK2", "BARK2": "GRK3", "BRK": "PTK6",
    "BUBR1": "BUB1B", "CDC2": "CDK1", "CDC7": "CDC7", "CDK11": "CDK19", "CGK2": "PRKG2", "CHK1": "CHEK1",
    "CHK2": "CHEK2", "CK1A": "CSNK1A1", "CK1D": "CSNK1D", "CK1E": "CSNK1E", "CK1G1": "CSNK1G1",
    "CK1G2": "CSNK1G2", "CK1G3": "CSNK1G3", "CK2A1": "CSNK2A1", "CK2A2": "CSNK2A2", "COT": "MAP3K8",
    "CRIK": "CIT", "CRK7": "CDK12", "CTK": "MATK", "CAMK1A": "CAMK1", "CAMK1B": "PNCK", "CAMK1D": "CAMK1D",
    "CAMK2A": "CAMK2A", "CAMK2B": "CAMK2B", "CAMK2D": "CAMK2D", "CAMK2G": "CAMK2G", "CAMK4": "CAMK4",
    "DLK": "MAP3K12", "DMPK1": "DMPK", "DRAK2": "STK17B", "EPHA7": "EPHA7", "ERBB2": "ERBB2", "ERBB4": "ERBB4",
    "ERK1": "MAPK3", "ERK2": "MAPK1", "ERK5": "MAPK7", "ERK7": "MAPK15", "FAK": "PTK2", "FMS": "CSF1R",
    "FRAP": "MTOR", "GAK": "GAK", "GCK": "MAP4K2", "GCN2": "EIF2AK4", "GPRK5": "GRK5", "GPRK7": "GRK7",
    "GSK3A": "GSK3A", "GSK3B": "GSK3B", "HGK": "MAP4K4", "HRI": "EIF2AK1", "HASPIN": "HASPIN", "ICK": "CILK1",
    "IKKA": "CHUK", "IKKB": "IKBKB", "IKKE": "IKBKE", "IRR": "INSRR", "JNK1": "MAPK8", "JNK2": "MAPK9",
    "JNK3": "MAPK10", "KHS1": "MAP4K5", "LKB1": "STK11", "LOK": "STK10", "LYNA": "LYN", "LYNB": "LYN",
    "LZK": "MAP3K13", "MER": "MERTK", "MINK": "MINK1", "MLK1": "MAP3K9", "MLK2": "MAP3K10", "MLK3": "MAP3K11",
    "MNK1": "MKNK1", "MNK2": "MKNK2", "MPSK1": "STK16", "MRCKA": "CDC42BPA", "MRCKB": "CDC42BPB",
    "MSK1": "RPS6KA5", "MSK2": "RPS6KA4", "MSSK1": "SRPK3", "MST1": "STK4", "MST2": "STK3", "MST3": "STK24",
    "MST4": "STK26", "MYT1": "PKMYT1", "NDR1": "STK38", "NDR2": "STK38L", "NIK": "MAP3K14", "NIM1": "NIM1K",
    "NUAK1": "NUAK1", "NUAK2": "NUAK2", "OSR1": "OXSR1", "PCTAIRE1": "CDK16", "PDGFRA": "PDGFRA",
    "PDGFRB": "PDGFRB", "PDHK1": "PDK1", "PDHK2": "PDK2", "PDHK3": "PDK3", "PDHK4": "PDK4", "PEK": "EIF2AK3",
    "PFTAIRE1": "CDK14", "PKACA": "PRKACA", "PKACB": "PRKACB", "PKACG": "PRKACG", "PKCA": "PRKCA",
    "PKCB1": "PRKCB", "PKCB2": "PRKCB", "PKCD": "PRKCD", "PKCE": "PRKCE", "PKCG": "PRKCG", "PKCH": "PRKCH",
    "PKCI": "PRKCI", "PKCT": "PRKCQ", "PKCZ": "PRKCZ", "PKD1": "PRKD1", "PKD2": "PRKD2", "PKD3": "PRKD3",
    "PKG1A": "PRKG1", "PKG1B": "PRKG1", "PKN1": "PKN1", "PKR": "EIF2AK2", "PRPK": "TP53RK", "PYK2": "PTK2B",
    "QIK": "SIK2", "RHOK": "GRK1", "RON": "MST1R", "ROS": "ROS1", "RSK1": "RPS6KA1", "RSK2": "RPS6KA3",
    "RSK3": "RPS6KA2", "RSK4": "RPS6KA6", "SGK496": "PKDCC", "SIK": "SIK1", "SRM": "SRMS", "STLK3": "STK39",
    "TAK1": "MAP3K7", "TAO1": "TAOK1", "TAO2": "TAOK2", "TGFBR1": "TGFBR1", "TIE2": "TEK", "TRKA": "NTRK1",
    "TRKB": "NTRK2", "TRKC": "NTRK3", "TSSK1": "TSSK1B", "TSSK2": "TSSK2", "YES": "YES1", "YSK1": "STK25",
    "ZAK": "MAP3K20", "P38A": "MAPK14", "P38B": "MAPK11", "P38D": "MAPK13", "P38G": "MAPK12",
    "P70S6K": "RPS6KB1", "P70S6KB": "RPS6KB2", "SKMLCK": "MYLK2", "SMMLCK": "MYLK", "ACTR2": "ACVR2A",
    "ACTR2B": "ACVR2B", "ULK3": "ULK3", "VRK1": "VRK1", "VRK2": "VRK2",
}


def build_sugiyama(U, sm, km, counts):
    import io
    import zipfile
    z = zipfile.ZipFile(EXT / "sugiyama2019/supp_PMC6642169.zip")
    x = pd.read_excel(io.BytesIO(z.read("41598_2019_46385_MOESM3_ESM.xlsx")), header=None, skiprows=2, dtype=str)
    x = x.iloc[:, :7]
    x.columns = ["type", "kname", "entry", "desc", "pos", "sidic", "ptmscore"]
    counts["sugiyama_rows_total"] = len(x)
    x = x[~x["type"].str.contains("mutant", na=False) & ~x["kname"].str.contains(r"\[", na=False)]
    counts["sugiyama_rows_wildtype"] = len(x)
    x = x[x["type"] != "LK"]
    counts["sugiyama_rows_protein_kinase"] = len(x)
    x = x[x["sidic"].notna()]
    counts["sugiyama_rows_sidic_confirmed"] = len(x)
    base = x["kname"].str.split("/").str[0].str.strip().str.upper()
    kin = []
    for b in base:
        k = km.by_name(SUGIYAMA_MANUAL.get(b, b)) if b in SUGIYAMA_MANUAL else km.by_name(b)
        kin.append(k)
    x["kinase"] = kin
    counts["sugiyama_kinase_names_unmapped_or_not_in_U"] = sorted(set(x.loc[x["kinase"].isna(), "kname"]))
    x = x[x["kinase"].notna()]
    an = pd.read_csv(EXT / "uniprot/UP000005640_annotations.tsv.gz", sep="\t", dtype=str,
                     usecols=["Entry", "Entry Name"])
    e2a = dict(zip(an["Entry Name"], an["Entry"]))
    x["acc"] = x["entry"].map(e2a)
    counts["sugiyama_rows_entry_name_unmapped"] = int(x["acc"].isna().sum())
    x = x[x["acc"].notna()]
    m = x["pos"].str.extract(r"^([STY])(\d+)$")
    x = x[m[0].notna()]
    m = m[m[0].notna()]
    # verify residue against UniProt canonical sequence
    seqs = _uniprot_canonical_seqs()
    ok = [(a in seqs and 0 < int(p) <= len(seqs[a]) and seqs[a][int(p) - 1] == r)
          for a, r, p in zip(x["acc"], m[0], m[1])]
    counts["sugiyama_rows_residue_mismatch_dropped"] = int(len(ok) - sum(ok))
    x, m = x[ok], m[ok]
    x["site"] = sm.key(x["acc"], m[0], m[1].astype(int))
    x["source"] = "Sugiyama2019_invitro"
    x["evidence"] = f"PMID:{SUGIYAMA_PMID};carna:" + x["kname"]
    lib = _collapse(x)
    counts["invitro_sugiyama_edges"] = len(lib)
    return lib[LIB_COLUMNS]


def _uniprot_canonical_seqs():
    seqs, cur, buf = {}, None, []
    with gzip.open(EXT / "uniprot/UP000005640_reviewed_isoforms.fasta.gz", "rt") as f:
        for line in f:
            if line.startswith(">"):
                if cur:
                    seqs[cur] = "".join(buf)
                cur = line[1:].split("|")[1]
                buf = []
            else:
                buf.append(line.strip())
    if cur:
        seqs[cur] = "".join(buf)
    return seqs


# --------------------------------------------------------------------------------------- predicted
def candidate_sites(sm: SiteMapper):
    """KSM site index -> canonical key, restricted to sites quantified in the compendium (KSM comp_uid) or present
    in the benchmarKIN external data. Returns (ksm_sites DataFrame with columns site_idx, site, cand, in_comp)."""
    ks = pd.read_parquet(EXT / "ksm_index/sites.parquet", columns=["site_id"])
    ks["site_idx"] = np.arange(len(ks), dtype=np.int64)
    bench = set()
    bp = EXT / "benchmarkin_py"
    bench |= set(pd.read_parquet(bp / "perturbation_data.parquet", columns=[]).index)
    for f in sorted(bp.glob("tumour_[A-Z]*.parquet")):
        if f.name == "tumour_gs.parquet":
            continue
        bench |= set(pd.read_parquet(f, columns=[]).index)
    ks["site"] = sm.from_keys(ks["site_id"])
    ks["in_comp"] = sm.in_compendium(ks["site"])
    ks["in_bench"] = ks["site_id"].isin(bench).to_numpy()
    ks["cand"] = ks["in_comp"] | ks["in_bench"]
    return ks


def _rank_topk(kin_idx, site_idx, score, ksm_sites, kinases, U, name, kmax=K_MAX):
    """Per kinase: dedupe by canonical site (max score), rank among candidates, keep top kmax."""
    codes, uniq = pd.factorize(ksm_sites["site"], sort=True)  # code order == lexicographic site key order
    ok = np.isfinite(score)
    kin_idx, site_idx, score = kin_idx[ok].astype(np.int64), site_idx[ok], score[ok].astype(np.float64)
    code = codes[site_idx].astype(np.int64)
    order = np.lexsort((code, -score, kin_idx))  # kinase, score desc, site key asc (deterministic ties)
    key = kin_idx[order] * (len(uniq) + 1) + code[order]
    _, first = np.unique(key, return_index=True)  # first = best-scoring row of each (kinase, canonical site)
    sel = order[np.sort(first)]
    k_s, c_s, v_s = kin_idx[sel], code[sel], score[sel]
    starts = np.r_[0, np.flatnonzero(np.diff(k_s)) + 1]
    lens = np.diff(np.r_[starts, len(k_s)])
    rank = np.arange(len(k_s)) - np.repeat(starts, lens) + 1
    n_cand = np.repeat(lens, lens)
    top = rank <= kmax
    df = pd.DataFrame({"kidx": k_s[top], "site": np.asarray(uniq, dtype=object)[c_s[top]], "score": v_s[top],
                       "rank": rank[top], "n_cand": n_cand[top]})
    df["weight"] = 1.0 - (df["rank"] - 1) / df["n_cand"]
    acc = kinases["uniprot_acc"].to_numpy()
    acc2k = dict(zip(U["uniprot"], U["kinase"]))
    df["kinase"] = [acc2k.get(acc[i]) for i in df["kidx"]]
    df = df[df["kinase"].notna()]
    df["source"] = name
    df["evidence"] = ("score=" + df["score"].map(lambda v: f"{v:.6g}") + ";rank=" + df["rank"].astype(str)
                      + ";n_cand=" + df["n_cand"].astype(str))
    out = df[LIB_COLUMNS + ["rank", "score"]].sort_values(["kinase", "rank"]).reset_index(drop=True)
    out["rank"] = out["rank"].astype(np.int32)
    return out


def verify_predicted_file(name, check_sha=True) -> dict:
    spec = PREDICTED_FILES[name]
    info = {"name": name, "path": str(spec["path"].relative_to(ROOT)), "origin": f"kinaseSubMain/{spec['origin']}"}
    for tok in BANNED_PATH_TOKENS:
        assert tok not in spec["origin"].lower() and tok not in str(spec["path"]).lower(), (name, tok)
    if check_sha:
        sha = sha256_file(spec["path"])
        assert sha == spec["sha256"], f"{name}: sha256 mismatch {sha}"
        assert sha not in BANNED_SHA256, f"{name} is a banned kinaseSubMain variant: {BANNED_SHA256[sha]}"
        info["sha256"] = sha
    man = json.loads(spec["manifest"].read_text())
    if name == "r2np":
        groups = man["config"]["groups"]
        assert not any(g.startswith("F6") for g in groups), groups  # no phospho features
        assert man["config"].get("comp_know") is False and man["split"] == "final_dev"
        assert man["outputs"]["uobs_predictions"]["sha256"] == spec["sha256"]
        info["manifest_check"] = f"groups={groups} (no F6), split={man['split']}, run_id={man['run_id']}"
    elif name == "r3":
        assert man["sha256"] == spec["sha256"] and man["path"] == spec["origin"]
        assert "r3_deploy" in man["run_id"]
        info["manifest_check"] = (f"run_id={man['run_id']}, ESM windows + U_ms phospho-channel bit only "
                                  f"(no perturbation features; docs/kinasesub_audit.md)")
    elif name == "cantley":
        assert man["meta"]["label_derived"] is False and man["meta"]["imputed"] is False
        assert man["meta"]["variant"] == "unprimed"
        info["manifest_check"] = "unprimed published-definition atlas score, label_derived=False, imputed=False"
    return info


def build_predicted(name, U, sm, ksm_sites=None, kinases_subset=None, kmax=K_MAX):
    """Top-kmax candidate sites per kinase for one predicted source. kinases_subset (U symbols) for tests."""
    if ksm_sites is None:
        ksm_sites = candidate_sites(sm)
    kin = pd.read_parquet(EXT / "ksm_index/kinases.parquet")
    acc2k = dict(zip(U["uniprot"], U["kinase"]))
    k_in_U = np.array([acc2k.get(a) is not None and (kinases_subset is None or acc2k.get(a) in kinases_subset)
                       for a in kin["uniprot_acc"]])
    cand = ksm_sites["cand"].to_numpy()
    spec = PREDICTED_FILES[name]
    if name == "cantley":
        M = np.load(spec["path"], mmap_mode="r")
        assert M.shape == (len(kin), len(ksm_sites))
        cidx = np.flatnonzero(cand)
        ks_, ss_, sc_ = [], [], []
        for ki in np.flatnonzero(k_in_U):
            row = np.asarray(M[ki, cidx], dtype=np.float32)
            ok = np.isfinite(row)
            if not ok.any():
                continue
            ks_.append(np.full(ok.sum(), ki, dtype=np.int32)); ss_.append(cidx[ok]); sc_.append(row[ok])
        kin_idx, site_idx, score = np.concatenate(ks_), np.concatenate(ss_), np.concatenate(sc_)
    else:
        import pyarrow.parquet as pq
        pf = pq.ParquetFile(spec["path"])
        kk, ss, sc = [], [], []
        for b in pf.iter_batches(batch_size=4_000_000, columns=["kinase_idx", "site_idx", "score"]):
            k = b.column(0).to_numpy(); s = b.column(1).to_numpy(); v = b.column(2).to_numpy()
            keep = cand[s] & k_in_U[k]
            kk.append(k[keep]); ss.append(s[keep]); sc.append(v[keep])
        kin_idx, site_idx, score = np.concatenate(kk), np.concatenate(ss), np.concatenate(sc)
    return _rank_topk(kin_idx, site_idx, score, ksm_sites, kin, U, name, kmax=kmax)


# --------------------------------------------------------------------------------------- loader
def load_library(name: str, k: int | None = None, gamma: float = 1.0, compendium_only: bool = False,
                 sites=None) -> pd.DataFrame:
    """Load a library. Predicted libraries: keep rank <= k (default all stored, k<=500) and weight = pct**gamma.
    `sites`: optional iterable of measured site keys to restrict to (e.g. per-condition / DEV-cache measured)."""
    lib = pd.read_parquet(LIBDIR / f"{name}.parquet")
    if k is not None:
        if "rank" not in lib.columns:
            raise ValueError(f"{name} is not a ranked (predicted) library")
        if k > K_MAX:
            raise ValueError(f"k={k} > stored K_MAX={K_MAX}")
        lib = lib[lib["rank"] <= k]
    if gamma != 1.0:
        lib = lib.assign(weight=lib["weight"] ** gamma)
    if compendium_only:
        comp = set(pd.read_parquet(CACHE / "sites.parquet", columns=["site_key"])["site_key"])
        lib = lib[lib["site"].isin(comp)]
    if sites is not None:
        lib = lib[lib["site"].isin(set(sites))]
    return lib.reset_index(drop=True)


def library_content_hash(lib: pd.DataFrame) -> str:
    cols = [c for c in LIB_COLUMNS + ["rank", "score"] if c in lib.columns]
    x = lib[cols].sort_values(["kinase", "site"]).reset_index(drop=True)
    return hashlib.sha256(pd.util.hash_pandas_object(x, index=False).values.tobytes()).hexdigest()
