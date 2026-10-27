"""Toy tests for the branch-2 Bayesian linear-Gaussian core (src/kinact/explore/bayes.py)."""
import numpy as np
import pandas as pd
import pytest

torch = pytest.importorskip("torch")
from kinact.explore import bayes as B  # noqa: E402

DEV = "cuda" if torch.cuda.is_available() else "cpu"


def toy(seed=0, n_sites=800, n_kin=30, n_cond=6, planted=((0, 3.0), (5, -3.0)), frac_meas=0.8, sub_per_kin=25):
    rng = np.random.default_rng(seed)
    kin = [f"K{i}" for i in range(n_kin)]
    sites = [f"P{i:05d}_S{i}" for i in range(n_sites)]
    edges = [(kin[k], sites[s], 1.0) for k in range(n_kin) for s in rng.choice(n_sites, sub_per_kin, replace=False)]
    lib = pd.DataFrame(edges, columns=["kinase", "site", "weight"]).drop_duplicates(["kinase", "site"])
    U = pd.DataFrame({"kinase": kin, "manning_family": [f"F{i // 3}" for i in range(n_kin)],
                      "manning_group": ["G1"] * n_kin})
    d = B.build_design(sites, kin, lib)
    W = d.Wc.toarray()  # sites x kinases (lib rows only)
    Wfull = np.zeros((n_sites, n_kin)); Wfull[d.lib_rows] = W
    a = np.zeros((n_kin, n_cond))
    for k, val in planted:
        a[k, 1:] = val            # condition 0 = null
    X = Wfull @ a * 0.5 + rng.normal(0, 1, (n_sites, n_cond))
    X[rng.random(X.shape) > frac_meas] = np.nan
    return pd.DataFrame(X, index=sites, columns=[f"c{j}" for j in range(n_cond)]), d, U, kin


def test_planted_recovered():
    X, d, U, kin = toy()
    cfg = B.BayesConfig(tau=1.0, rho_f=0.0, rho_g=0.0, intercept="median")
    Z = B.score(X, d, U, cfg, device=DEV)
    for c in X.columns[1:]:
        top = Z.loc[c].sort_values()
        assert top.index[-1] == "K0" and top.index[0] == "K5"
        assert top["K0"] > 2.5 and top["K5"] < -2.5
    assert Z.loc["c0"].abs().max() < 4


def test_closed_form_matches_numpy_and_cpu_path():
    X, d, U, kin = toy(seed=1)
    cfg = B.BayesConfig(tau=0.7, rho_f=1.0, rho_g=0.5, intercept="median")
    S, parts = B.score(X, d, U, cfg, device=DEV, return_all=True)
    K = B.kinase_kernel(U, kin, cfg.rho_f, cfg.rho_g)
    x = X["c2"].to_numpy()
    med = np.nanmedian(x); sig = 1.4826 * np.nanmedian(np.abs(x - med))
    xl = x[d.lib_rows]; m = ~np.isnan(xl)
    Wm = d.Wc.toarray()[m]
    Lam = Wm.T @ Wm + np.linalg.inv(K) / cfg.tau ** 2
    mu = np.linalg.solve(Lam, Wm.T @ ((xl[m] - med) / sig))
    sd = np.sqrt(np.diag(np.linalg.inv(Lam)))
    np.testing.assert_allclose(parts["mean"].loc["c2"].to_numpy(), mu, rtol=1e-3, atol=1e-4)
    np.testing.assert_allclose(parts["sd"].loc["c2"].to_numpy(), sd, rtol=1e-3, atol=1e-5)
    ss = B.SingleSample(d, U, cfg, None)
    z1, mu1, sd1, n1 = ss.score_one(xl, med, sig)
    np.testing.assert_allclose(mu1, mu, rtol=1e-8, atol=1e-10)
    np.testing.assert_allclose(z1, S.loc["c2"].to_numpy(), rtol=1e-3, atol=1e-4)


def test_family_borrowing_scores_unlinked_kinase():
    X, d, U, kin = toy(seed=2)
    lib = pd.DataFrame({"kinase": [], "site": [], "weight": []})
    # drop K1's edges (same family F0 as planted K0): it can only be scored via K
    from kinact.explore.bayes import Design
    Wc = d.Wc.tolil(); Wc[:, 1] = 0; d2 = Design(d.kinases, d.lib_rows, Wc.tocsr(), d.Wp)
    cfg = B.BayesConfig(tau=1.0, rho_f=2.0, rho_g=0.0, borrow=False)
    assert B.score(X, d2, U, cfg, device=DEV)["K1"].isna().all()
    cfg.borrow = True
    Z = B.score(X, d2, U, cfg, device=DEV)
    assert (Z.loc["c1":, "K1"] > 0.2).all()         # borrowed from planted K0 (same family)
    assert Z["K29"].notna().all()


def test_sign_flip_and_shift_invariance():
    X, d, U, kin = toy(seed=3)
    cfg = B.BayesConfig(tau=1.0, rho_f=1.0, rho_g=0.0, intercept="joint")
    Z = B.score(X, d, U, cfg, device=DEV)
    np.testing.assert_allclose(B.score(-X, d, U, cfg, device=DEV).to_numpy(), -Z.to_numpy(), atol=1e-4)
    np.testing.assert_allclose(B.score(X + 5.0, d, U, cfg, device=DEV).to_numpy(), Z.to_numpy(), atol=1e-3)


def test_null_calibration_toy():
    X, d, U, kin = toy(seed=4, planted=(), n_cond=40, n_sites=2000, sub_per_kin=40)
    Z = B.score(X, d, U, B.BayesConfig(tau=3.0, rho_f=0.0, rho_g=0.0), device=DEV)
    frac = (Z.abs() > 1.96).to_numpy().mean()
    assert 0.02 < frac < 0.09


def test_noise_variance_downweights_noisy_sites():
    rng = np.random.default_rng(5)
    V = rng.normal(0, 1, (1000, 200)); V[:100] *= 3
    v = B.site_noise_variance(V)
    assert np.median(v[:100]) > 5 and abs(np.median(v[100:]) - 1) < 0.15
