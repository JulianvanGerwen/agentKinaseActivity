"""Tests for core-v3 in the user-facing scorer (kinact.core with config/core_v3.json):
late-fused perturbation-response signatures (B7b) + E1 per-kinase permutation-null scale (B8).

* toy: fusion algebra, sign symmetry, non-signature kinases == core-v1, min_sites, NaN contract, core-v2 fallback;
* DEV (skipped without cache/dev + models/core_v3): equivalence with the research code path on >= 200 DEV conditions
  (bayes2 core z_vif + explore.signatures fit + signatures_late.sig_project/rho_sig + final_box.perm_scales /
  fuse_calibrated, every quantity recomputed independently of the artefacts).
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from kinact import core
from kinact.core import KinaseActivityModel

ROOT = Path(__file__).resolve().parents[1]
CFG3 = ROOT / "config/core_v3.json"
ART = ROOT / "models/core_v3"
DEV_OK = (ROOT / "cache/dev/quant_dev.npy").exists() and (ART / "edges.tsv.gz").exists() \
    and (ART / "signatures.tsv.gz").exists()
try:
    import torch
    HAS_CUDA = torch.cuda.is_available()
except ImportError:                                                        # pragma: no cover
    torch, HAS_CUDA = None, False

from test_core import toy_data, toy_model  # noqa: E402


# ----------------------------------------------------------------------------- toy
@pytest.fixture(scope="module")
def toy3():
    m, e, kin, sites, unknown = toy_model()
    rng = np.random.default_rng(3)
    rows = []
    # signatures for K2 (planted, sites in the library) and K5 (partly on non-library sites)
    for k, ss in (("K2", list(e.loc[e.kinase == "K2", "site"])[:8]),
                  ("K5", list(rng.choice(sites, 10, replace=False)) + unknown[:6])):
        for s in ss:
            rows.append((k, s, float(rng.normal(1.0, 0.3))))
    sig = pd.DataFrame(rows, columns=["kinase", "site", "weight"])
    scales = pd.DataFrame({"kinase": ["K2", "K5"], "v_both": [2.5, 1.7], "v_sig": [1.3, 0.9]})
    m3, *_ = toy_model()
    m3.attach_signatures(sig, scales, rho_sig=0.8, alpha=1.0)
    X = toy_data(m, e, sites, unknown)
    return m, m3, sig, scales, X


def ref_zsig(X, sig, rho, clip=4.0, min_sites=5):
    V = X.to_numpy(float)
    med = np.nanmedian(V, 0); R = V - med
    s = 1.4826 * np.nanmedian(np.abs(R), 0); Y = R / s
    Y = Y / np.sqrt(np.nanmean(Y * Y, 0))
    out = {}
    for k, g in sig.groupby("kinase"):
        ri = X.index.get_indexer(g.site)
        ok = ri >= 0
        y = Y[ri[ok]]; w = g.weight.to_numpy()[ok]
        m = np.isfinite(y)
        num = np.nansum(np.clip(np.where(m, y, 0), -clip, clip) * w[:, None], 0)
        den = np.sqrt((m * (w ** 2)[:, None]).sum(0))
        out[k] = np.where(m.sum(0) >= min_sites, num / den, np.nan) / rho
    return out


def test_toy_fusion_formula(toy3):
    m1, m3, sig, scales, X = toy3
    r1, r3 = m1.score(X), m3.score(X)
    zs = ref_zsig(X, sig, 0.8)
    sc = scales.set_index("kinase")
    for k in ("K2", "K5"):
        np.testing.assert_allclose(r3.z_sig[k].to_numpy(), zs[k], atol=1e-12, equal_nan=True)
        zc = r1.z[k].to_numpy()
        want = np.where(np.isfinite(zc) & np.isfinite(zs[k]), (zc + zs[k]) / np.sqrt(sc.loc[k, "v_both"]),
                        np.where(np.isfinite(zc), zc, zs[k] / np.sqrt(sc.loc[k, "v_sig"])))
        np.testing.assert_allclose(r3.z[k].to_numpy(), want, atol=1e-12, equal_nan=True)
    other = [k for k in r1.z.columns if k not in ("K2", "K5")]
    pd.testing.assert_frame_equal(r3.z[other], r1.z[other])
    pd.testing.assert_frame_equal(r3.z_core, r1.z)
    pd.testing.assert_frame_equal(r3.mean, r1.mean)
    assert r3.info["n_signature_sites_matched"] > 0


def test_toy_sign_flip_and_long(toy3):
    _, m3, _, _, X = toy3
    a, b = m3.score(X), m3.score(-X)
    np.testing.assert_allclose(a.z.to_numpy(), -b.z.to_numpy(), atol=1e-10, equal_nan=True)
    lg = a.long()
    assert {"z_core", "z_sig"} <= set(lg.columns) and len(lg) == a.z.size


def test_toy_min_sites_and_core_v2_fallback(toy3):
    m1, _, sig, _, X = toy3
    m2, *_ = toy_model()
    m2.attach_signatures(sig, None, rho_sig=1.0, alpha=0.5)            # core-v2: fixed sqrt(1 + a^2)
    r1, r2 = m1.score(X), m2.score(X)
    zs = ref_zsig(X, sig, 1.0)["K2"]
    zc = r1.z["K2"].to_numpy()
    np.testing.assert_allclose(r2.z["K2"].to_numpy(), np.where(np.isfinite(zs), (zc + .5 * zs) / np.sqrt(1.25), zc),
                               atol=1e-12, equal_nan=True)
    # fewer than min_sites measured signature sites -> z_sig NaN -> z = z_core
    Xm = X.copy()
    Xm.loc[sig.loc[sig.kinase == "K2", "site"][:4], :] = np.nan
    r = m2.score(Xm)
    assert r.z_sig["K2"].isna().all()
    np.testing.assert_allclose(r.z["K2"].to_numpy(), r.z_core["K2"].to_numpy(), equal_nan=True)


# ----------------------------------------------------------------------------- DEV equivalence
@pytest.fixture(scope="module")
def dev3():
    from kinact.data import load_dev
    X = load_dev().astype(np.float64)
    rng = np.random.default_rng(0)
    conds = sorted(rng.choice(sorted(X.columns), 200, replace=False))      # the fixed bench set (seed 0)
    return X, conds


@pytest.mark.skipif(not DEV_OK or torch is None, reason="DEV cache / core-v3 artefacts / torch missing")
def test_config_and_manifest():
    m = KinaseActivityModel.from_config(CFG3)
    cfg = json.loads(CFG3.read_text())
    assert m.sig is not None and m.cfg["name"] == "core-v3"
    assert m.sig["alpha"] == 1.0 and m.sig["rho"] == cfg["signatures"]["rho_sig"]
    man = m.manifest
    from kinact.library import sha256_file
    for f, h in man["artifact_sha256"].items():
        assert sha256_file(ART / f) == h, f
    assert man["signatures"]["n_top"] == cfg["signatures"]["n_top"]
    assert man["signatures"]["training_splits"] == [f"fold{i}" for i in range(5)]
    # core part identical to core-v1
    m1 = KinaseActivityModel.from_config()
    assert m1.sites.equals(m.sites) and (m1._W != m._W).nnz == 0


@pytest.mark.skipif(not DEV_OK or torch is None, reason="DEV cache / core-v3 artefacts / torch missing")
def test_equivalence_research_dev(dev3):
    """>= 200 DEV conditions: deployed core-v3 z == research path recomputed from scratch."""
    from kinact import controls as ctl
    from kinact.data import load_dev, load_dev_meta
    from kinact.eval import Evaluator, load_gt, read_gt_tsv
    from kinact.explore import bayes2 as B2, final_box as FB, signatures as SG, signatures_late as SL
    from kinact.library import load_library
    X, conds = dev3
    cfg = json.loads(CFG3.read_text())
    sg = cfg["signatures"]
    m = KinaseActivityModel.from_config(CFG3)
    FOLDS = [f"fold{i}" for i in range(5)]
    P, C, U, _ = load_gt()
    Udf = read_gt_tsv(ROOT / "results/ground_truth/kinase_universe.tsv")
    kin = list(Udf["kinase"])
    meta = load_dev_meta()
    X32 = load_dev()
    # ---- research signatures (B7 fit, clean+contrast, all DEV-CV folds)
    Yd, Md, _ = B2.standardise(X32.to_numpy(np.float64), "wsd", None)
    Yfit = np.where(Md > 0, Yd, np.nan).astype(np.float32)
    Y32, M32 = Yd.astype(np.float32), Md.astype(np.float32)
    tr = SG.training_conditions(C, meta, FOLDS, [])
    pos = SG.positives(P, C, tr, 10)
    sig = SG.fit_signatures(Yfit, X32.columns, tr, pos, meta["group_id"], kappa=float(sg["kappa"]),
                            n_max=int(sg["n_top"]), contrast=True)
    art = pd.read_csv(ART / "signatures.tsv.gz", sep="\t")
    assert set(sig.sites) == set(art.kinase)
    for k, s in sig.sites.items():
        a = art[art.kinase == k].sort_values("rank")
        assert list(X32.index[s]) == list(a.site)
        np.testing.assert_allclose(sig.weights[k], a.weight.to_numpy(), rtol=1e-7)
    Ws = sig.matrix(len(X32.index), kin, int(sg["n_top"]))
    allc = pd.Index(X32.columns)
    nulls = [c for c in C.loc[(C.role == "null") & C.split.isin(FOLDS), "conditionID"] if c in allc]
    rho_ref = SL.rho_sig(SL.sig_project(Ws, Y32[:, allc.get_indexer(nulls)], M32[:, allc.get_indexer(nulls)]))
    assert abs(rho_ref - m.sig["rho"]) < 1e-5 * rho_ref
    # ---- research core z_vif (bayes2 batched scorer)
    md = B2.build_multi_design(X32.index, kin, load_library("curated_union"), {"r2np": load_library("r2np", k=100)})
    bcfg = B2.Bayes2Config(tau=cfg["tau"], rho_f=cfg["rho_f"], rho_g=cfg["rho_g"], alphas={"r2np": 1.0},
                           output="z_vif", rho_e=cfg["rho_e"], scale="wsd", winsor=None, borrow=True)
    dev_ = "cuda" if HAS_CUDA else "cpu"
    # ---- research E1 scales on site-permuted DEV-CV GT conditions
    bins = read_gt_tsv(ROOT / "results/ground_truth/harmony_bins_group.tsv")
    gconds = list(Evaluator(P, C, U, split_filter=FOLDS, bins=bins, zero_policy="ignore").cond_index)
    Xp = ctl.site_permutation(X32[gconds], seed=0, mode="within_condition")
    zc_p = B2.score(Xp, md, Udf, bcfg, device=dev_).reindex(columns=kin).to_numpy()
    Yp, Mp, _ = B2.standardise(Xp.to_numpy(np.float64), "wsd", None)
    zs_p = SL.sig_project(Ws, Yp.astype(np.float32), Mp.astype(np.float32)) / rho_ref
    has = np.isin(kin, list(sig.sites))
    v, vs = FB.perm_scales(zc_p, zs_p, float(sg["alpha"]), has)
    sc = pd.read_csv(ART / "e1_scales.tsv", sep="\t").set_index("kinase").reindex(kin)
    np.testing.assert_allclose(sc.v_both.to_numpy()[has], v[has], rtol=1e-4)
    np.testing.assert_allclose(sc.v_sig.to_numpy()[has], vs[has], rtol=1e-4)
    # ---- 200 DEV bench conditions + all 372 DEV-CV GT conditions: fused z
    for cs in (conds, gconds):
        R = m.score(X[cs])
        zc = B2.score(X32[cs], md, Udf, bcfg, device=dev_).reindex(columns=kin).to_numpy()
        ci = allc.get_indexer(cs)
        zs = SL.sig_project(Ws, Y32[:, ci], M32[:, ci]) / rho_ref
        ref = FB.fuse_calibrated(zc, zs, float(sg["alpha"]), v, vs)
        got = R.z.to_numpy()
        assert len(cs) >= 200
        assert np.array_equal(np.isnan(got), np.isnan(ref))
        assert np.nanmax(np.abs(got - ref)) < 1e-4, np.nanmax(np.abs(got - ref))
        np.testing.assert_allclose(R.z_sig.to_numpy(), zs, atol=1e-5, equal_nan=True)
        # non-signature kinases: unchanged core-v1 z
        np.testing.assert_array_equal(R.z.to_numpy()[:, ~has], R.z_core.to_numpy()[:, ~has])
