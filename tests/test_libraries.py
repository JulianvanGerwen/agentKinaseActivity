"""Tests for the shared KS library builder (src/kinact/library.py, scripts/build_libraries.py).

Require the built libraries in cache/libraries (python scripts/build_libraries.py) and the external inputs;
the module is skipped otherwise.
"""
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import hypergeom

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from kinact import library as L  # noqa: E402

pytestmark = pytest.mark.skipif(not all((L.LIBDIR / f"{n}.parquet").exists() for n in L.ALL_LIBS),
                                reason="libraries not built")
SANITY = ["MAPK1", "CDK1", "AKT1", "PRKACA"]


@pytest.fixture(scope="module")
def libs():
    return {n: pd.read_parquet(L.LIBDIR / f"{n}.parquet") for n in L.ALL_LIBS}


@pytest.fixture(scope="module")
def U():
    return L.load_universe()


@pytest.fixture(scope="module")
def manifest():
    return json.loads((L.RESLIB / "build_manifest.json").read_text())


def _pmids(evidence):
    out = set()
    for e in evidence:
        for m in re.finditer(r"PMID:([\d,]+)", e):
            out |= set(m.group(1).split(","))
    return out


def test_no_edge_cites_compendium_pmid(libs):
    comp = L.compendium_pmids()
    assert len(comp) >= 178
    for name, lib in libs.items():
        assert not (_pmids(lib["evidence"]) & comp), name


def test_no_curated_edge_cites_excluded_pmid(libs):
    ex = pd.read_csv(L.RESLIB / "excluded_pmids.tsv", sep="\t", dtype=str)
    assert set(L.compendium_pmids()) <= set(ex["pmid"])
    bad = set(ex["pmid"])
    assert not (_pmids(libs["curated_union"]["evidence"]) & bad)


def test_curated_union_has_no_noncurated_only_sources(libs):
    src = libs["curated_union"]["source"].str.split(";")
    for s in src:
        assert any(not (x.endswith("MIMP") or x.endswith("ProtMapper") or x == "iPTMnet:rlim+") for x in s)
    assert not libs["curated_union"]["source"].str.contains("ikip", case=False).any()


def test_all_kinases_in_universe(libs, U):
    Uset = set(U["kinase"])
    assert len(Uset) == len(U)
    for name, lib in libs.items():
        assert lib["kinase"].isin(Uset).all(), name


def test_site_keys_parse(libs):
    for name, lib in libs.items():
        for s in lib["site"].unique():
            acc, res, pos = L.parse_site_key(s)
            assert res in "STY" and pos > 0
    with pytest.raises(ValueError):
        L.parse_site_key("P28482-T185")


def test_contract_columns_and_weights(libs):
    for name, lib in libs.items():
        assert list(lib.columns[:5]) == L.LIB_COLUMNS
        assert (lib["weight"] > 0).all() and (lib["weight"] <= 1).all(), name
        assert not lib.duplicated(["kinase", "site"]).any(), name
    for name in L.PREDICTED_LIBS:
        lib = libs[name]
        assert lib["rank"].max() <= L.K_MAX
        g = lib.groupby("kinase")
        assert (g["rank"].min() == 1).all()
        assert (g["rank"].max() == g.size()).all()  # contiguous ranks
        d = lib.sort_values(["kinase", "rank"])
        same = d["kinase"].values[1:] == d["kinase"].values[:-1]
        assert (np.diff(d["weight"].values)[same] < 0).all()
        assert (np.diff(d["score"].values)[same] <= 0).all()


def test_psp_libraries_nested(libs):
    a = set(zip(libs["psp_all"]["kinase"], libs["psp_all"]["site"]))
    for n in ["psp_invivo", "psp_human_only"]:
        assert set(zip(libs[n]["kinase"], libs[n]["site"])) <= a, n
    cu = set(zip(libs["curated_union"]["kinase"], libs["curated_union"]["site"]))
    assert a <= cu


def test_deterministic_rebuild_curated(libs, U, manifest):
    sm, km = L.SiteMapper(), L.KinaseMapper(U)
    psp, _ = L.build_psp(U, sm, km)
    for n, lib in psp.items():
        assert L.library_content_hash(lib) == manifest["library_content_sha256"][n] == \
            L.library_content_hash(libs[n]), n
    ex = pd.read_csv(L.RESLIB / "excluded_pmids.tsv", sep="\t", dtype=str)
    cu, *_ = L.build_curated_union(psp["psp_all"], U, sm, km, ex, {})
    assert L.library_content_hash(cu) == manifest["library_content_sha256"]["curated_union"]
    sg = L.build_sugiyama(U, sm, km, {})
    assert L.library_content_hash(sg) == manifest["library_content_sha256"]["invitro_sugiyama"]


@pytest.mark.parametrize("name", L.PREDICTED_LIBS)
def test_deterministic_rebuild_predicted_subset(libs, U, name):
    sm = L.SiteMapper()
    sub = L.build_predicted(name, U, sm, kinases_subset=set(SANITY))
    stored = libs[name][libs[name]["kinase"].isin(SANITY)].reset_index(drop=True)
    assert len(sub) == len(stored) > 0
    assert L.library_content_hash(sub) == L.library_content_hash(stored)


@pytest.mark.parametrize("name", L.PREDICTED_LIBS)
def test_predicted_files_are_audited_variants(name):
    info = L.verify_predicted_file(name)
    assert info["sha256"] not in L.BANNED_SHA256


@pytest.fixture(scope="module")
def candidates():
    c = L.candidate_sites(L.SiteMapper())
    return c[c["cand"]]


@pytest.mark.parametrize("name", L.PREDICTED_LIBS)
@pytest.mark.parametrize("kinase", SANITY)
def test_predicted_top200_enriched_for_psp(libs, U, candidates, name, kinase):
    rc = dict(zip(U["kinase"], U["residue_class"]))[kinase]
    keys = candidates["site"].drop_duplicates()
    is_y = keys.str.contains(r"_Y\d+$")
    pop = set(keys[is_y] if rc == "Y" else keys[~is_y])
    psp = libs["psp_all"]
    P = set(psp.loc[psp["kinase"] == kinase, "site"]) & pop
    T = set(libs[name].loc[(libs[name]["kinase"] == kinase) & (libs[name]["rank"] <= 200), "site"])
    assert len(T) == 200 and len(P) >= 20
    x = len(T & P)
    p = hypergeom.sf(x - 1, len(pop), len(P), len(T))
    assert x > 200 * len(P) / len(pop) and p < 0.01, (name, kinase, x, p)


def test_load_library_topk():
    lib = L.load_library("r3", k=50, gamma=2.0)
    assert lib["rank"].max() == 50
    assert np.allclose(lib["weight"], (1 - (lib["rank"] - 1) / lib["evidence"].str.extract(r"n_cand=(\d+)")[0]
                                       .astype(float)) ** 2)
    with pytest.raises(ValueError):
        L.load_library("psp_all", k=10)
