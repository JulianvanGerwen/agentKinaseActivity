"""Toy tests for branch B4b (src/kinact/explore/edgeweights.py)."""
import numpy as np
import pandas as pd
import pytest

torch = pytest.importorskip("torch")
from kinact.explore import bayes2 as B2  # noqa: E402
from kinact.explore import edgeweights as EW  # noqa: E402
from kinact.explore.bayes import report_mask  # noqa: E402

DEV = "cuda" if torch.cuda.is_available() else "cpu"


def toy(seed=0, n_sites=800, n_kin=30, n_cond=6):
    rng = np.random.default_rng(seed)
    kin = [f"K{i}" for i in range(n_kin)]
    sites = [f"P{i:05d}_{'Y' if i % 7 == 0 else 'S'}{i}" for i in range(n_sites)]
    cur = pd.DataFrame([(kin[k], sites[s], 1.0, "PSP;OmniPath:SIGNOR") for k in range(n_kin)
                        for s in rng.choice(n_sites, 10, replace=False)],
                       columns=["kinase", "site", "weight", "source"]).drop_duplicates(["kinase", "site"])
    pred = pd.DataFrame([(kin[k], sites[s], 1 - r / 40, r + 1) for k in range(n_kin)
                         for r, s in enumerate(rng.choice(n_sites, 40, replace=False))],
                        columns=["kinase", "site", "weight", "rank"])
    U = pd.DataFrame({"kinase": kin, "manning_family": [f"F{i // 3}" for i in range(n_kin)],
                      "manning_group": ["TK" if i < 5 else "AGC" for i in range(n_kin)]})
    md = B2.build_multi_design(sites, kin, cur, {"r2np@20": pred[pred["rank"] <= 20]})
    X = rng.normal(0, 1, (n_sites, n_cond))
    X[rng.random(X.shape) > 0.7] = np.nan
    X = pd.DataFrame(X, index=sites, columns=[f"c{j}" for j in range(n_cond)])
    E = EW.edge_table(md, sites, U, cur, cur.iloc[:150], cur.iloc[:80], {"r2np": pred, "r3": pred, "cantley": pred})
    return X, md, U, kin, E


def setup(theta_dim=None):
    X, md, U, kin, E = toy()
    Y, M, _ = B2.standardise(X.to_numpy(), "wsd", None)
    K = B2.kernel2(U, kin, 3.0, 0.25)
    core = EW.DiffCore(E, len(md.lib_rows), len(kin), K, 0.125, DEV)
    Yt = torch.as_tensor(Y[md.lib_rows].astype(np.float32), device=DEV)
    Mt = torch.as_tensor(M[md.lib_rows].astype(np.float32), device=DEV)
    return X, md, U, kin, E, K, core, Yt, Mt


def test_edge_table_features():
    X, md, U, kin, E = toy()
    assert len(E) == sum(m.nnz for m in md.mats.values())
    assert set(EW.STATIC) <= set(E.columns)
    assert ((E["psp_vivo"] + E["psp_vitro"] + E["is_r2np"]) <= 1).all()
    assert (E.loc[E["is_r2np"] == 1, "n_src"] == 0).all()
    # residue mismatch: TK & S-site or non-TK & Y-site
    tk = E["k"] < 5
    y = E["site"].str.contains("_Y")
    assert (E["res_mismatch"].astype(bool) == (tk ^ y)).all()
    ls, mf = EW.site_stats(np.zeros((len(md.lib_rows), 3)), np.ones((len(md.lib_rows), 3)), np.arange(3))
    Phi = EW.feature_matrix(E, ls, mf)
    assert Phi.shape == (len(E), len(EW.FEATURES))
    assert np.allclose(Phi.mean(0), 0, atol=1e-5)


def test_theta_zero_reproduces_bayes2():
    X, md, U, kin, E, K, core, Yt, Mt = setup()
    cfg = B2.Bayes2Config(tau=0.125, rho_f=3.0, rho_g=0.25, alphas={"r2np@20": 1.0}, output="z_vif", rho_e=0.01,
                          scale="wsd")
    ref = B2.score(X, md, U, cfg, device=DEV).to_numpy()
    Phi = torch.zeros(len(E), 3, device=DEV)
    Y, M, _ = B2.standardise(X.to_numpy(), "wsd", None)
    n = M[md.lib_rows].T @ (md.W({"r2np@20": 1.0}).toarray() != 0)
    rep = report_mask(n, K, 1, True)
    for th in (None, np.zeros(3, np.float32)):
        S = EW.score_theta(core, th, Phi, Yt, Mt, rep, out="z_vif", rho_e=0.01)
        assert np.array_equal(np.isfinite(S), np.isfinite(ref))
        assert np.nanmax(np.abs(S - ref)) < 1e-4


def test_weights_change_scores_and_gradient_flows():
    X, md, U, kin, E, K, core, Yt, Mt = setup()
    Phi = torch.as_tensor(EW.feature_matrix(E, np.zeros(len(md.lib_rows)), np.zeros(len(md.lib_rows)),
                                            EW.STATIC), device=DEV)
    th = torch.zeros(Phi.shape[1], device=DEV, requires_grad=True)
    z = core.z(core.W(th, Phi), Yt, Mt, "z_vif", 0.01)
    z[:, 0].sum().backward()
    assert th.grad is not None and torch.isfinite(th.grad).all() and th.grad.abs().sum() > 0
    th2 = torch.full((Phi.shape[1],), 0.5, device=DEV)
    z2 = core.z(core.W(th2, Phi), Yt, Mt, "z_vif", 0.01)
    assert not torch.allclose(z.detach(), z2)


def test_site_stats_exclude_columns():
    rng = np.random.default_rng(1)
    Y = rng.normal(size=(50, 10)); M = (rng.random((50, 10)) > 0.2).astype(float); Y = Y * M
    a = EW.site_stats(Y, M, np.arange(6))
    Y2, M2 = Y.copy(), M.copy()
    Y2[:, 6:] = rng.normal(size=(50, 4)) * 100; M2[:, 6:] = 1
    b = EW.site_stats(Y2, M2, np.arange(6))
    assert np.allclose(a[0], b[0]) and np.allclose(a[1], b[1])


def test_sigmoid_loss_matches_auroc_at_low_temperature():
    rng = np.random.default_rng(2)
    z = torch.as_tensor(rng.normal(size=(3, 20)), dtype=torch.float64)
    ki = torch.tensor([0, 4, 7]); sign = torch.tensor([1.0, -1.0, 1.0], dtype=torch.float64)
    neg = torch.ones(3, 20, dtype=torch.float64)
    neg[0, 0] = neg[1, 4] = neg[2, 7] = 0
    nN = neg.sum(1); w = torch.full((3,), 1 / 3, dtype=torch.float64)
    loss = EW.pairwise_loss(z, torch.arange(3), ki, sign, w, neg, nN, T=1e-4, kind="sigmoid")
    zz = z.numpy() * sign.numpy()[:, None]
    au = [np.mean(zz[i, ki[i]] > np.delete(zz[i], ki[i])) for i in range(3)]
    assert abs(float(loss) - (1 - np.mean(au))) < 1e-6
