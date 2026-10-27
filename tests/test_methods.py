"""Toy tests for the baseline methods (DESIGN §9): planted kinase effects are ranked top, sign-flip
antisymmetry, NaN (never 0) for kinases without measured substrates / below min_substrates."""
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from kinact.controls import sign_flip_check  # noqa: E402
from kinact.methods import (gsea, kl_enrich, ksea, mlm, norm_wmean, ulm, viper, wmean,  # noqa: E402
                            zscore)
from kinact.methods.jumpsem import RSCRIPT, run_jumpsem  # noqa: E402

N_SITES, N_K, N_SUB = 600, 12, 8


def toy(seed=0, n_cond=4, effect=2.5, missing=0.1):
    rng = np.random.default_rng(seed)
    sites = [f"P{i:05d}_S{10 + i}" for i in range(N_SITES)]
    kin = [f"K{k}" for k in range(N_K)]
    rows = []
    for k in range(N_K):
        for s in rng.choice(N_SITES, N_SUB, replace=False):
            rows.append((kin[k], sites[s], 1.0))
    rows.append(("KNOSUB", "Q99999_S1", 1.0))          # substrate not in data -> must be NaN
    rows.append(("KTWO", sites[0], 1.0)); rows.append(("KTWO", sites[1], 1.0))  # 2 < min_substrates
    lib = pd.DataFrame(rows, columns=["kinase", "site", "weight"]).drop_duplicates(["kinase", "site"])
    X = pd.DataFrame(rng.normal(0, 0.5, (N_SITES, n_cond)), index=sites, columns=[f"c{j}" for j in range(n_cond)])
    planted = {}
    for j, c in enumerate(X.columns):
        k = kin[j % N_K]
        sgn = 1 if j % 2 == 0 else -1
        subs = lib.loc[lib.kinase == k, "site"]
        X.loc[subs, c] += sgn * effect
        planted[c] = (k, sgn)
    X = X.mask(rng.random(X.shape) < missing)
    # keep planted kinase's substrates measured
    for c, (k, _) in planted.items():
        subs = lib.loc[lib.kinase == k, "site"]
        X.loc[subs, c] = X.loc[subs, c].fillna(effect)
    return X, lib, planted


SINGLE = {"zscore": zscore, "ksea": ksea, "wmean": wmean, "norm_wmean": norm_wmean, "ulm": ulm, "mlm": mlm,
          "viper": lambda X, l, **k: viper(X, l, n_jobs=1, **k), "gsea": lambda X, l, **k: gsea(X, l, n_jobs=1, **k),
          "kl_enrich": lambda X, l, **k: kl_enrich(X, l, thr=1.0, **k)}


@pytest.mark.parametrize("name", list(SINGLE))
def test_planted_kinase_ranked_top(name):
    X, lib, planted = toy()
    S = SINGLE[name](X, lib)
    for c, (k, sgn) in planted.items():
        t = (sgn * S.loc[c]).dropna()
        assert t.idxmax() == k, (name, c, t.sort_values().tail(3))


@pytest.mark.parametrize("name", list(SINGLE))
def test_nan_for_unscorable_kinases(name):
    X, lib, _ = toy()
    S = SINGLE[name](X, lib)
    assert S["KNOSUB"].isna().all() and S["KTWO"].isna().all()
    assert S.drop(columns=["KNOSUB", "KTWO"]).notna().any().all()
    S1 = SINGLE[name](X, lib, min_substrates=1)
    assert S1["KNOSUB"].isna().all()
    assert not (S.fillna(1) == 0).all().any()


@pytest.mark.parametrize("name", ["zscore", "ksea", "wmean", "norm_wmean", "ulm", "mlm", "viper", "kl_enrich"])
def test_sign_flip_exact(name):
    X, lib, _ = toy(seed=1)
    r = sign_flip_check(lambda X_, l: SINGLE[name](X_, l), X, lib, atol=1e-8, rtol=1e-6)
    assert r["passed"], r


def test_sign_flip_gsea_approx():
    # permutation-normalised: antisymmetric up to Monte-Carlo noise of the null
    X, lib, _ = toy(seed=1)
    a, b = SINGLE["gsea"](X, lib), SINGLE["gsea"](-X, lib)
    assert np.array_equal(a.isna().values, b.isna().values)
    ok = a.notna().values
    assert np.corrcoef(a.values[ok], -b.values[ok])[0, 1] > 0.98


def test_zscore_formula():
    X, lib, _ = toy(seed=2)
    S = zscore(X, lib)
    c, k = "c0", "K3"
    x = X[c].dropna()
    s = x.reindex(lib.loc[lib.kinase == k, "site"]).dropna()
    assert np.isclose(S.loc[c, k], s.sum() / (x.std(ddof=1) * np.sqrt(len(s))))
    K = ksea(X, lib)
    assert np.isclose(K.loc[c, k], (s.mean() - x.mean()) * np.sqrt(len(s)) / x.std(ddof=1))


@pytest.mark.skipif(not Path(RSCRIPT).exists(), reason="R env not installed")
def test_jumpsem_planted_time_course():
    """Multi-sample: one kinase's substrates follow a latent activity across 6 conditions."""
    rng = np.random.default_rng(3)
    sites = [f"P{i:05d}_S{10 + i}" for i in range(80)]
    act = np.array([-2, -1, 0, 1, 2, 3.0])
    X = pd.DataFrame(rng.normal(0, 0.2, (80, 6)), index=sites, columns=[f"c{j}" for j in range(6)])
    X.iloc[:6] += np.outer(np.linspace(0.6, 1.2, 6), act)
    lib = pd.DataFrame({"kinase": ["KA"] * 6 + ["KB"] * 5 + ["KN"],
                        "site": sites[:6] + sites[20:25] + ["Q99999_S1"], "weight": 1.0})
    r = run_jumpsem(X, lib, n_jobs=1)
    a = r["aligned"]["KA"].to_numpy()
    assert np.corrcoef(a, act)[0, 1] > 0.95
    assert r["aligned"]["KN"].isna().all()
    assert abs(np.corrcoef(r["native"]["KA"].to_numpy(), act)[0, 1]) > 0.95
