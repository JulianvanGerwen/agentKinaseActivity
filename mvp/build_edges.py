"""Build mvp/model/edges.tsv.gz (the kinase-substrate library of core-v3) from your own PhosphoSitePlus download.

PhosphoSitePlus (PSP) is licensed CC BY-NC-SA 3.0 and requires registration, so PSP-derived edges are not shipped
with this repository. Download `Kinase_Substrate_Dataset.gz` (and, recommended, `Phosphorylation_site_dataset.gz`
for mapping mouse/rat sites to human) from https://www.phosphosite.org/staticDownloads, then run:

    python mvp/build_edges.py --psp-ks Kinase_Substrate_Dataset.gz [--psp-sites Phosphorylation_site_dataset.gz] \
        [--omnipath omnipath_enzsub.tsv]

Optional OmniPath edges (recommended; adds the non-PSP curated resources) are downloaded with
    curl -o omnipath_enzsub.tsv "https://omnipathdb.org/enzsub?genesymbols=yes&fields=sources,references&format=tsv&license=academic"

Edges written (one per kinase x site, weight 1):
  cur   PSP kinase-substrate pairs (human, mouse, rat), kinase mapped to the 503-kinase universe (model/kinases.tsv)
        by human UniProt accession or gene symbol; substrate site mapped to the human site key `{UniProt}_{res}{pos}`
        directly (human rows) or through the PSP site group (mouse/rat rows; needs --psp-sites).
  r2np  predicted edges shipped in model/predicted_edges.tsv.gz (top-100 per kinase), minus pairs already in `cur`.

  cur   (--omnipath) curated OmniPath phosphorylation edges: MIMP/ProtMapper predictions dropped; edges need a
        reference or a manually curated source (PhosphoSite, SIGNOR, phospho.ELM, HPRD, PhosphoNetworks, Li2012).

Approximation vs the evaluated model: the research library `curated_union` also used iPTMnet and dropped edges that
cite the benchmark studies (needed for a fair evaluation, not for use). Versions of PSP/OmniPath also drift. Scores are
therefore close to, but not identical with, the TEST results in REPORT.md.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
MODEL = HERE / "model"
ORGS = {"human", "mouse", "rat"}
OMNI_NOREF_OK = {"PhosphoSite", "SIGNOR", "phosphoELM", "HPRD", "PhosphoNetworks", "Li2012"}


def read_psp(path: str) -> pd.DataFrame:
    """PSP flat files: 3 header lines (date code, licence, blank), latin-1, tab-separated."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="latin-1") as fh:
        lines = fh.readlines()
    start = next(i for i, l in enumerate(lines) if l.startswith("GENE\t"))
    from io import StringIO
    return pd.read_csv(StringIO("".join(lines[start:])), sep="\t", dtype=str, keep_default_na=False)


def split_rsd(rsd: pd.Series):
    """'S473' or 'S473-p' -> ('S', 473); anything else -> (None, None)."""
    m = rsd.str.extract(r"^([STY])(\d+)")
    return m[0], pd.to_numeric(m[1], errors="coerce")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--psp-ks", required=True, help="PSP Kinase_Substrate_Dataset(.gz)")
    ap.add_argument("--psp-sites", default=None, help="PSP Phosphorylation_site_dataset(.gz), maps mouse/rat sites")
    ap.add_argument("--omnipath", default=None, help="OmniPath enzsub TSV (see module docstring)")
    ap.add_argument("--out", default=str(MODEL / "edges.tsv.gz"))
    a = ap.parse_args(argv)

    kin = pd.read_csv(MODEL / "kinases.tsv", sep="\t", keep_default_na=False)
    by_acc = dict(zip(kin["uniprot"], kin["kinase"]))
    by_name = {k.upper(): k for k in kin["kinase"]}

    ks = read_psp(a.psp_ks)
    ks = ks[ks["KIN_ORGANISM"].isin(ORGS) & ks["SUB_ORGANISM"].isin(ORGS)].copy()
    acc = ks["KIN_ACC_ID"].str.split("-").str[0]
    k_acc = acc.map(by_acc).where(ks["KIN_ORGANISM"].eq("human"))
    ks["kinase"] = k_acc.fillna(ks["GENE"].str.upper().map(by_name))
    ks = ks[ks["kinase"].notna()]

    res, pos = split_rsd(ks["SUB_MOD_RSD"])
    human = ks["SUB_ORGANISM"].eq("human") & ~ks["SUB_ACC_ID"].str.contains("-") & res.notna()
    ks["site"] = (ks["SUB_ACC_ID"] + "_" + res + pos.astype("Int64").astype(str)).where(human)
    if a.psp_sites:                                   # mouse/rat substrate sites -> human site in the same group
        ps = read_psp(a.psp_sites)
        ps = ps[ps["ORGANISM"].eq("human") & ~ps["ACC_ID"].str.contains("-")]
        r2, p2 = split_rsd(ps["MOD_RSD"])
        ps = ps.assign(site=ps["ACC_ID"] + "_" + r2 + p2.astype("Int64").astype(str))[r2.notna()]
        grp = ps.drop_duplicates("SITE_GRP_ID").set_index("SITE_GRP_ID")["site"]
        ks["site"] = ks["site"].fillna(ks["SITE_GRP_ID"].map(grp))
    n_unmapped = int(ks["site"].isna().sum())
    cur = ks.dropna(subset=["site"])[["site", "kinase"]]
    n_psp = len(cur.drop_duplicates())
    if a.omnipath:
        o = pd.read_csv(a.omnipath, sep="\t", dtype=str, keep_default_na=False)
        o = o[o["modification"].eq("phosphorylation") & o["residue_type"].isin(["S", "T", "Y"])
              & ~o["substrate"].str.contains("-")].copy()
        src = o["sources"].str.split(";").map(lambda x: {v for v in x if v and v not in ("MIMP", "ProtMapper")
                                                        and not v.endswith(("_MIMP", "_ProtMapper"))})
        o = o[src.map(len).gt(0) & (o["references"].str.len().gt(0) | src.map(lambda x: bool(x & OMNI_NOREF_OK)))]
        o["kinase"] = o["enzyme"].map(by_acc).fillna(o["enzyme_genesymbol"].str.upper().map(by_name))
        o = o[o["kinase"].notna()]
        o["site"] = o["substrate"] + "_" + o["residue_type"] + o["residue_offset"]
        cur = pd.concat([cur, o[["site", "kinase"]]])
    cur = cur.drop_duplicates()
    cur = cur.assign(weight=1.0, source="cur")

    pred = pd.read_csv(MODEL / "predicted_edges.tsv.gz", sep="\t")
    key = set(zip(cur["site"], cur["kinase"]))
    pred = pred[[(s, k) not in key for s, k in zip(pred["site"], pred["kinase"])]]
    edges = pd.concat([cur, pred[["site", "kinase", "weight", "source"]]], ignore_index=True)
    edges = edges.sort_values(["kinase", "source", "site"]).reset_index(drop=True)
    edges.to_csv(a.out, sep="\t", index=False, compression="gzip")

    sha = hashlib.sha256(Path(a.out).read_bytes()).hexdigest()
    print(f"wrote {a.out}: {len(cur)} curated edges ({n_psp} from PSP; {cur['kinase'].nunique()} kinases; {n_unmapped} PSP rows "
          f"without a human site{'' if a.psp_sites else ' - pass --psp-sites to map mouse/rat'}), "
          f"{len(pred)} predicted edges; sha256 {sha[:12]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
