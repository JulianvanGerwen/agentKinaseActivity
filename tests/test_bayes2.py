"""Toy tests for branch B2.1/B3 (src/kinact/explore/bayes2.py)."""
import numpy as np
import pandas as pd
import pytest

torch = pytest.importorskip("torch")
from kinact.explore import bayes2 as B2  # noqa: E402

DEV = "cuda" if torch.cuda.is_available() else "cpu"


def toy(seed=0, n_sites=1500, n_kin=40, n_cond=8, planted=((0, 3.0), (5, -3.0)), frac_meas=0.8, heavy=False):
    rng = np.random.default_rng(seed)
    kin = [f"K{i}" for i in range(n_kin)]
    sites = [f"P{i:05d}_S{i}" for i in range(n_sites)]
    sizes = rng.integers(3, 60, n_kin)
    cur = pd.DataFrame([(kin[k], sites[s], 1.0) for k in range(n_kin) for s in rng.choice(n_sites, sizes[k] // 2, replace=False)],
                       columns=["kinase", "site", "weight"]).drop_duplicates(["kinase", "site"])
    pred = {}
    for src in ("pa", "pb"):
        e = [(kin[k], sites[s], 1.0, r + 1) for k in range(n_kin)
             for r, s in enumerate(rng.choice(n_sites, sizes[k], replace=False))]
        pred[src] = pd.DataFrame(e, columns=["kinase", "site", "weight", "rank"])
    U = pd.DataFrame({"kinase": kin, "manning_family": [f"F{i // 4}" for i in range(n_kin)],
                      "manning_group": ["G1"] * n_kin})
    md = B2.build_multi_design(sites, kin, cur, pred)
    W = md.W({"pa": 1.0, "pb": 0.0}).toarray()
    Wf = np.zeros((n_sites, n_kin)); Wf[md.lib_rows] = W
    a = np.zeros((n_kin, n_cond))
    for k, v in planted:
        a[k, 1:] = v
    noise = rng.standard_t(3, (n_sites, n_cond)) if heavy else rng.normal(0, 1, (n_sites, n_cond))
    X = Wf @ a * 0.3 + noise
    X[rng.random(X.shape) > frac_meas] = np.nan
    return pd.DataFrame(X, index=sites, columns=[f"c{j}" for j in range(n_cond)]), md, U, kin


def test_blocks_combine_equal_direct():
    X, md, U, kin = toy()
    Y, M, _ = B2.standardise(X.to_numpy(), "wsd", 3.0)
    bl = B2.block_stats(md, Y, M, device=DEV)
    al = {"pa": 0.7, "pb": 0.3}
    G, h, g0, n = B2.combine(bl, al)
    W = md.W(al).toarray()
    m = M[md.lib_rows, 2]; y = Y[md.lib_rows, 2]
    np.testing.assert_allclose(G[2].cpu().numpy(), (W * m[:, None]).T @ W, rtol=1e-4, atol=1e-3)
    np.testing.assert_allclose(h[2].cpu().numpy(), W.T @ (m * y), rtol=1e-4, atol=1e-3)
    np.testing.assert_allclose(g0[2].cpu().numpy(), W.T @ m, rtol=1e-4, atol=1e-3)


def test_predicted_edges_exclude_curated():
    X, md, U, kin = toy()
    C = md.mats["cur"].toarray() != 0
    for n in ("pa", "pb"):
        assert not np.any(C & (md.mats[n].toarray() != 0))


def test_posterior_quantities_and_cpu_path():
    X, md, U, kin = toy(seed=1)
    cfg = B2.Bayes2Config(tau=0.5, rho_f=1.0, rho_g=0.25, alphas={"pa": 1.0, "pb": 0.5}, output="z_freq",
                          scale="wsd", winsor=4.0)
    S, parts = B2.score(X, md, U, cfg, device=DEV, return_all=True)
    post = parts["post"]
    Y, M, _ = B2.standardise(X.to_numpy(), "wsd", 4.0)
    W = md.W(cfg.alphas).toarray(); m = M[md.lib_rows, 3] > 0; Wm = W[m]
    L = Wm.T @ Wm + np.linalg.inv(parts["K"]) / cfg.tau ** 2
    A = np.linalg.solve(L, Wm.T)
    np.testing.assert_allclose(post["mu"][3], A @ Y[md.lib_rows, 3][m], rtol=1e-3, atol=1e-4)
    np.testing.assert_allclose(post["q"][3], (A * A).sum(1), rtol=1e-3, atol=1e-6)
    np.testing.assert_allclose(post["suma"][3], A.sum(1), rtol=1e-3, atol=1e-5)
    ss = B2.SingleSample2(md, U, cfg)
    sc, zc, mu, q, n = ss.score_one(X["c3"].to_numpy())
    np.testing.assert_allclose(sc, S.loc["c3"].to_numpy(), rtol=1e-3, atol=1e-4)


def test_planted_recovered_all_outputs():
    X, md, U, kin = toy(seed=2, planted=((0, 8.0), (5, -8.0)))
    for out in ["mean", "z_post", "z_freq", "hybrid", "empnull"]:
        cfg = B2.Bayes2Config(tau=1.0, rho_f=0.0, rho_g=0.0, alphas={"pa": 1.0}, output=out, borrow=False)
        Z = B2.score(X, md, U, cfg, device=DEV)
        for c in X.columns[1:]:
            s = Z.loc[c].sort_values()
            assert s.index[-1] == "K0" and s.index[0] == "K5", out


def test_z_freq_calibrated_under_site_permutation():
    """Heavy-tailed noise, no signal: z_freq with the winsorised-SD scale has ~5% |z|>1.96 in small and large
    substrate sets (posterior z would not)."""
    X, md, U, kin = toy(seed=3, planted=(), n_cond=60, heavy=True, n_sites=3000)
    cfg = B2.Bayes2Config(tau=0.125, rho_f=1.0, rho_g=0.0, alphas={"pa": 1.0, "pb": 1.0}, output="z_freq",
                          scale="wsd", winsor=3.0, borrow=False)
    Z, parts = B2.score(X, md, U, cfg, device=DEV, return_all=True)
    z = Z.to_numpy(); n = parts["n"]
    ok = ~np.isnan(z)
    for sel in (ok & (n < 40), ok & (n >= 40)):
        assert 0.02 < np.mean(np.abs(z[sel]) > 1.96) < 0.09


def test_empirical_null_removes_m_dependent_inflation():
    rng = np.random.default_rng(0)
    me = rng.uniform(1, 80, (50, 400))
    Z = rng.normal(0, 1, me.shape) * 1.3 * np.sqrt(1 + 0.05 * (me - 1))
    zc, s0, rho = B2.empirical_null(Z, me, np.ones_like(Z, bool))
    assert abs(np.median(rho) - 0.05) < 0.02 and abs(np.median(s0) - 1.3) < 0.1
    for sel in (me < 20, me > 60):
        assert 0.03 < np.mean(np.abs(zc[sel]) > 1.96) < 0.075
    r, _ = B2.fit_rho_null(Z, me, np.ones_like(Z, bool))
    assert 0.03 <= r <= 0.075


def test_similarity_kernel_psd_and_sign_flip():
    X, md, U, kin = toy(seed=4)
    lib = pd.concat([md_lib for md_lib in [pd.DataFrame({"kinase": kin * 3, "site": [f"s{i % 7}" for i in range(len(kin) * 3)],
                                                         "rank": 1})]])
    S = B2.topk_similarity(lib, kin, 500, power=2)
    K = B2.kernel2(U, kin, 1.0, 0.25, 3.0, S)
    assert np.allclose(np.diag(K), 1) and np.linalg.eigvalsh(K).min() > -1e-9
    cfg = B2.Bayes2Config(tau=0.25, rho_f=1.0, rho_s=1.0, alphas={"pa": 1.0}, output="z_freq", scale="wsd", winsor=3.0)
    Z1 = B2.score(X, md, U, cfg, S=S, device=DEV); Z2 = B2.score(-X, md, U, cfg, S=S, device=DEV)
    np.testing.assert_allclose(Z2.to_numpy(), -Z1.to_numpy(), atol=1e-4)
