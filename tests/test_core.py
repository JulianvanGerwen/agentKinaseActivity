"""Tests for the user-facing core-v1 scorer (src/kinact/core.py) and its CLI."""
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from kinact import core
from kinact.core import KinaseActivityModel

ROOT = Path(__file__).resolve().parents[1]
CFG = json.loads((ROOT / "config/core_v1.json").read_text())
HAS_CUDA = False
try:
    import torch
    HAS_CUDA = torch.cuda.is_available()
except ImportError:                                                        # pragma: no cover
    torch = None


# ----------------------------------------------------------------------------- toy model
def toy_model(seed=0, n_kin=30, n_sites=1200, rho_e=0.01, **over):
    rng = np.random.default_rng(seed)
    kin = pd.DataFrame({"kinase": [f"K{i}" for i in range(n_kin)], "uniprot": [f"Q{i:05d}" for i in range(n_kin)],
                        "manning_group": [["AGC", "CMGC", "Other"][i % 3] for i in range(n_kin)],
                        "manning_family": [f"F{i // 4}" if i < n_kin - 3 else None for i in range(n_kin)]})
    sites = [f"P{i:05d}_{'STY'[i % 3]}{10 + i}" for i in range(n_sites)]
    rows = []
    for k in range(n_kin - 1):                                            # last kinase: no substrates at all
        for s in rng.choice(n_sites // 2, rng.integers(4, 25), replace=False):
            rows.append((sites[s], f"K{k}", 1.0, "cur"))
        for s in rng.choice(n_sites, rng.integers(5, 30), replace=False):
            rows.append((sites[s], f"K{k}", 1.0, "pred"))
    e = pd.DataFrame(rows, columns=["site", "kinase", "weight", "source"]).drop_duplicates(["site", "kinase", "source"])
    cur = set(zip(e.loc[e.source == "cur", "site"], e.loc[e.source == "cur", "kinase"]))
    e = e[(e.source == "cur") | ~pd.Series(list(zip(e.site, e.kinase)), index=e.index).isin(cur)]
    cfg = {**CFG, "library": {"curated": "toy", "predicted": {"pred": {"k": 100, "alpha": 1.0}}}, "rho_e": rho_e, **over}
    # unknown (non-library) sites
    unknown = [f"O{i:05d}_S{i + 1}" for i in range(800)]
    return KinaseActivityModel(cfg, e, kin), e, kin, sites, unknown


def toy_data(model, e, sites, unknown, n_cond=6, planted=(("K2", 2.0), ("K9", -2.0)), seed=1, frac=0.85):
    rng = np.random.default_rng(seed)
    idx = sites + unknown
    X = rng.normal(0, 1, (len(idx), n_cond))
    for k, v in planted:
        sk = e.loc[e.kinase == k, "site"]
        X[pd.Index(idx).get_indexer(sk), 1:] += v
    X[rng.random(X.shape) > frac] = np.nan
    return pd.DataFrame(X, index=idx, columns=[f"c{j}" for j in range(n_cond)])


def dense_reference(X, e, kin, cfg):
    """Independent dense implementation of the core-v1 posterior (bayes2 algebra) for toy checks."""
    V = X.to_numpy(float)
    med = np.nanmedian(V, 0); R = V - med
    s = 1.4826 * np.nanmedian(np.abs(R), 0); Y = R / s
    Y = Y / np.sqrt(np.nanmean(Y * Y, 0))
    kidx = pd.Index(kin.kinase)
    alpha = {"cur": 1.0, **{n: v["alpha"] for n, v in cfg["library"]["predicted"].items()}}
    W = np.zeros((len(X), len(kidx)))
    si = X.index.get_indexer(e.site); ki = kidx.get_indexer(e.kinase)
    ok = si >= 0
    np.add.at(W, (si[ok], ki[ok]), e.weight.to_numpy()[ok] * e.source.map(alpha).to_numpy()[ok])
    fam = kin.manning_family.fillna("").to_numpy(); grp = kin.manning_group.fillna("").to_numpy()
    K = core.kinase_kernel(fam, grp, cfg["rho_f"], cfg["rho_g"])
    P0 = np.linalg.inv(K) / cfg["tau"] ** 2
    rel = (K > 1e-12) & ~np.eye(len(K), dtype=bool)
    Z = np.full((X.shape[1], len(kidx)), np.nan)
    for j in range(X.shape[1]):
        m = np.isfinite(V[:, j])
        Wm = W[m]
        G = Wm.T @ Wm
        Li = np.linalg.inv(G + P0)
        mu = Li @ (Wm.T @ Y[m, j]); q = np.diag(Li @ G @ Li); suma = Li @ Wm.sum(0)
        with np.errstate(all="ignore"):
            me = np.where(q > 0, suma ** 2 / np.where(q > 0, q, 1), 0)
            z = mu / np.sqrt(q * (1 + cfg["rho_e"] * np.maximum(me - 1, 0)))
        has = (Wm != 0).sum(0) >= 1
        rep = has | (has.astype(float) @ rel > 0)
        Z[j] = np.where(rep, z, np.nan)
    return pd.DataFrame(Z, index=X.columns, columns=kidx)


@pytest.fixture(scope="module")
def toy():
    m, e, kin, sites, unknown = toy_model()
    X = toy_data(m, e, sites, unknown)
    return m, e, kin, sites, unknown, X


# ----------------------------------------------------------------------------- toy tests
def test_toy_matches_dense_reference(toy):
    m, e, kin, _, _, X = toy
    R = m.score(X, n_jobs=1)
    ref = dense_reference(X, e, kin, m.cfg)
    assert np.array_equal(np.isnan(R.z.to_numpy()), np.isnan(ref.to_numpy()))
    np.testing.assert_allclose(R.z.to_numpy(), ref.to_numpy(), atol=1e-9, equal_nan=True)


def test_planted_kinases_recovered(toy):
    m, *_, X = toy
    Z = m.score(X).z
    for c in X.columns[1:]:
        z = Z.loc[c]
        assert z.idxmax() == "K2" and z.idxmin() == "K9"
        assert z["K2"] > 3 and z["K9"] < -3
    assert Z.loc["c0"].abs().max() < 4                                     # no planted signal in c0


def test_sign_flip(toy):
    m, *_, X = toy
    a, b = m.score(X), m.score(-X)
    np.testing.assert_allclose(b.z.to_numpy(), -a.z.to_numpy(), atol=1e-10, equal_nan=True)
    np.testing.assert_allclose(b.mean.to_numpy(), -a.mean.to_numpy(), atol=1e-12, equal_nan=True)
    np.testing.assert_allclose(b.sd.to_numpy(), a.sd.to_numpy(), atol=1e-12, equal_nan=True)


def test_nan_contract(toy):
    m, e, kin, sites, unknown, X = toy
    R = m.score(X)
    # K29 has no substrates and only a group relation inside 'Other' (no coupling) and no family -> never scored
    assert R.z["K29"].isna().all() and R.mean["K29"].isna().all() and R.sd["K29"].isna().all()
    assert (R.n_substrates["K29"] == 0).all()
    scored = R.z.notna()
    assert np.array_equal(scored.to_numpy(), R.mean.notna().to_numpy())
    assert np.array_equal(scored.to_numpy(), R.sd.notna().to_numpy())
    assert np.isfinite(R.z.to_numpy()[scored.to_numpy()]).all()
    assert (R.n_substrates.to_numpy() >= 0).all() and R.n_substrates.dtypes.eq(np.int64).all()
    # kinases with measured substrates are always scored; borrowed ones have n_substrates == 0
    assert scored.to_numpy()[R.n_substrates.to_numpy() >= 1].all()
    # an all-missing condition -> every kinase NaN
    Xn = X.copy(); Xn["empty"] = np.nan
    assert m.score(Xn).z.loc["empty"].isna().all()


def test_unknown_and_malformed_sites(toy):
    m, e, kin, sites, unknown, X = toy
    lib_only = X.loc[sites]
    R_all, R_lib = m.score(X), m.score(lib_only)
    n_unk = int((~X.index.isin(m.sites)).sum())
    assert R_all.info["n_unknown_sites"] == n_unk and R_lib.info["n_unknown_sites"] == n_unk - len(unknown)
    assert R_all.info["n_library_sites_matched"] == len(X) - n_unk
    # unknown sites only enter centre/scale: same ranking per condition, different scale
    for c in X.columns:
        a, b = R_all.z.loc[c].dropna(), R_lib.z.loc[c].dropna()
        assert a.index.equals(b.index)
    # malformed keys are kept for centre/scale and ignored by the model (no crash)
    Xb = pd.concat([X, pd.DataFrame(np.ones((3, X.shape[1])), index=["junk", "AKT1", "P1_X5"], columns=X.columns)])
    Rb = m.score(Xb)
    assert Rb.info["n_unknown_sites"] == n_unk + 3
    # lower-case / padded keys are normalised
    Xl = X.copy(); Xl.index = [f" {k.lower()} " for k in X.index]
    np.testing.assert_allclose(m.score(Xl).z.to_numpy(), R_all.z.to_numpy(), atol=1e-12, equal_nan=True)
    # nothing matches the library -> warning and all NaN
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        Ru = m.score(X.loc[unknown])
    assert any("no input site matches" in str(x.message) for x in w)
    assert Ru.z.isna().all().all()


def test_single_column_series_long_inputs(toy):
    m, *_, X = toy
    R = m.score(X)
    s = X["c3"]
    for inp in (s, X[["c3"]], s.rename(None)):
        r1 = m.score(inp)
        assert r1.z.shape == (1, len(m.kinases))
        np.testing.assert_allclose(r1.z.to_numpy()[0], R.z.loc["c3"].to_numpy(), atol=1e-12, equal_nan=True)
    # long table
    acc, rest = zip(*[k.rsplit("_", 1) for k in X.index])
    L = pd.DataFrame({"uniprot": acc, "residue": [r[0] for r in rest], "position": [int(r[1:]) for r in rest],
                      "value": s.to_numpy()}).dropna()
    rl = m.score(L, condition="c3")
    np.testing.assert_allclose(rl.z.to_numpy()[0], R.z.loc["c3"].to_numpy(), atol=1e-12, equal_nan=True)
    Lc = pd.concat([L.assign(condition="a"), L.assign(condition="b", value=-L["value"])])
    r2 = m.score(Lc)
    np.testing.assert_allclose(r2.z.loc["b"].to_numpy(), -r2.z.loc["a"].to_numpy(), atol=1e-10, equal_nan=True)
    # duplicated site keys are averaged
    Xd = pd.concat([X[["c3"]], X[["c3"]]])
    np.testing.assert_allclose(m.score(Xd).z.to_numpy()[0], R.z.loc["c3"].to_numpy(), atol=1e-12, equal_nan=True)
    # output shapes
    assert set(r1.wide().columns) == {"z", "mean", "sd", "n_substrates"}
    lg = R.long()
    assert len(lg) == X.shape[1] * len(m.kinases) and len(R.long(dropna=True)) == int(R.z.notna().sum().sum())


def test_shared_patterns_and_woodbury_equal_full():
    m, e, kin, sites, unknown = toy_model(seed=3)
    X = toy_data(m, e, sites, unknown, n_cond=10, seed=4)
    # near-duplicate measured sets (a few sites dropped / added) and exact duplicates
    X["c1"] = X["c0"].where(np.arange(len(X)) % 97 != 0)
    X["c2"] = X["c0"] * 2 + 0.3
    X.iloc[:15, 3] = 1.0
    st = m.prepare(X)
    m.cfg["woodbury_r_max"] = 0
    full = m.score(st, n_jobs=1)
    m.cfg["woodbury_r_max"] = 10_000
    wb = m.score(st, n_jobs=1)
    m.cfg.pop("woodbury_r_max")
    assert wb.info["n_woodbury"] > 0 and full.info["n_woodbury"] == 0
    assert wb.info["n_patterns"] < X.shape[1]
    for a in ("z", "mean", "sd"):
        np.testing.assert_allclose(getattr(wb, a).to_numpy(), getattr(full, a).to_numpy(), atol=1e-9, equal_nan=True)
    thr = m.score(st, n_jobs=4)
    np.testing.assert_allclose(thr.z.to_numpy(), full.z.to_numpy(), atol=1e-9, equal_nan=True)
    ref = dense_reference(X, e, kin, m.cfg)
    np.testing.assert_allclose(full.z.to_numpy(), ref.to_numpy(), atol=1e-9, equal_nan=True)


def test_numpy_fallback_equal(toy):
    m, *_, X = toy
    st = m.prepare(X)
    po = m.posterior(st)
    out = {k: np.zeros_like(po["mu"]) for k in ("mu", "sd", "q", "suma")}
    Mf = st["M"].astype(float)
    H = np.asarray(m._WT @ (st["Y"] * Mf)); G0 = np.asarray(m._WT @ Mf)
    _, first, pat = np.unique(np.packbits(st["M"], axis=0).T, axis=0, return_index=True, return_inverse=True)
    m._run_cpu_numpy(st["M"], H, G0, first, pat.ravel(), out)
    for k in out:
        np.testing.assert_allclose(out[k], po[k], atol=1e-10)


@pytest.mark.skipif(not HAS_CUDA, reason="no CUDA")
def test_gpu_equals_cpu(toy):
    m, *_, X = toy
    a, b = m.score(X), m.score(X, device="cuda")
    for k in ("z", "mean", "sd"):
        np.testing.assert_allclose(getattr(b, k).to_numpy(), getattr(a, k).to_numpy(), atol=1e-9, equal_nan=True)
    assert b.n_substrates.equals(a.n_substrates)


def test_cli_roundtrip(tmp_path, toy, monkeypatch):
    m, *_, X = toy
    from kinact import cli
    monkeypatch.setattr(core.KinaseActivityModel, "from_config", classmethod(lambda cls, path=None: m))
    inp = tmp_path / "in.tsv"
    X.rename_axis("site").to_csv(inp, sep="\t")
    out = tmp_path / "out.tsv"
    assert cli.main(["score", str(inp), "-o", str(out), "-q"]) == 0
    L = pd.read_csv(out, sep="\t")
    assert list(L.columns) == ["condition", "kinase", "z", "mean", "sd", "n_substrates"]
    R = m.score(X)
    piv = L.pivot(index="condition", columns="kinase", values="z").loc[R.z.index, R.z.columns]
    np.testing.assert_allclose(piv.to_numpy(), R.z.to_numpy(), rtol=1e-5, atol=1e-5, equal_nan=True)
    outw = tmp_path / "w.tsv"
    assert cli.main(["score", str(inp), "-o", str(outw), "--format", "wide", "-q"]) == 0
    Wd = pd.read_csv(outw, sep="\t", index_col=0)
    assert "z:c0" in Wd.columns and len(Wd) == len(m.kinases)


# ----------------------------------------------------------------------------- DEV equivalence with bayes2
DEV_OK = (ROOT / "cache/dev/quant_dev.npy").exists() and (ROOT / "models/core_v1/edges.tsv.gz").exists()


@pytest.fixture(scope="module")
def dev():
    from kinact.data import load_dev
    X = load_dev().astype(np.float64)
    rng = np.random.default_rng(0)
    conds = sorted(rng.choice(sorted(X.columns), 200, replace=False))     # bench set (scripts/run_baselines.py)
    return X[conds]


@pytest.mark.skipif(not DEV_OK or torch is None, reason="DEV cache / artefacts / torch missing")
def test_config_and_artifacts():
    m = KinaseActivityModel.from_config()
    assert (m.cfg["tau"], m.cfg["rho_f"], m.cfg["rho_g"], m.cfg["rho_e"]) == (0.125, 3.0, 0.25, 0.01)
    assert m.manifest["rho_e_refit"] == m.cfg["rho_e"]
    assert len(m.kinases) == 503


@pytest.mark.skipif(not DEV_OK or torch is None, reason="DEV cache / artefacts / torch missing")
def test_equivalence_bayes2_dev(dev):
    from kinact.eval import read_gt_tsv
    from kinact.explore import bayes2 as B2
    from kinact.library import load_library
    m = KinaseActivityModel.from_config()
    U = read_gt_tsv(ROOT / "results/ground_truth/kinase_universe.tsv")
    kin = list(U["kinase"])
    md = B2.build_multi_design(dev.index, kin, load_library("curated_union"), {"r2np": load_library("r2np", k=100)})
    cfg = B2.Bayes2Config(tau=m.cfg["tau"], rho_f=m.cfg["rho_f"], rho_g=m.cfg["rho_g"], alphas={"r2np": 1.0},
                          output="z_vif", rho_e=m.cfg["rho_e"], scale="wsd", winsor=None, borrow=True)
    R = m.score(dev)
    # (a) bayes2 batched scorer (float32 sufficient statistics) on all 200 bench conditions
    S = B2.score(dev, md, U, cfg, device="cuda" if HAS_CUDA else "cpu").reindex(columns=m.kinases)
    assert np.array_equal(np.isnan(S.to_numpy()), R.z.isna().to_numpy())
    assert np.nanmax(np.abs(S.to_numpy() - R.z.to_numpy())) < 1e-5
    # (b) bayes2 float64 single-sample path: exact up to rounding
    ss = B2.SingleSample2(md, U, cfg)
    for c in dev.columns[:25]:
        z = ss.score_one(dev[c].to_numpy())[0]
        np.testing.assert_allclose(R.z.loc[c].to_numpy(), z, atol=1e-9, rtol=0, equal_nan=True)
