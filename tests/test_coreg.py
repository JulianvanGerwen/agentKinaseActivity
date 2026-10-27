"""Leakage and correctness tests for branch B4 co-regulation (src/kinact/explore/coreg.py).

The key test fails if any held-out condition (or a condition sharing a control with one) influences the
co-regulation fit used to score that held-out fold.
"""
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

torch = pytest.importorskip("torch")
from kinact.explore import coreg as CR  # noqa: E402

DEV = "cuda" if torch.cuda.is_available() else "cpu"


def toy_meta(n_cond=60, seed=0):
    rng = np.random.default_rng(seed)
    conds = [f"c{i}" for i in range(n_cond)]
    split = [f"fold{i % 3}" for i in range(n_cond)]
    pmid = [str(100 + i // 6) for i in range(n_cond)]
    ctrl = ["DMSO" if i % 2 else "vehicle" for i in range(n_cond)]
    meta = pd.DataFrame({"conditionID": conds, "split": split, "PMID": pmid, "Control": ctrl,
                         "Sample": ["HeLa"] * n_cond, "qc_status": ["ok"] * n_cond}).set_index("conditionID", drop=False)
    meta.loc["c5", "qc_status"] = "scale_error"
    return meta


def fit(Y, meta, conds, fold, W, q=0):
    """The B4 fitting path (as in scripts/exp_b4_precompute.py) for one fold."""
    al = CR.allowed_mask(meta, conds, CR.fold_heldout(meta, fold))
    Ya = Y[:, al]
    P = CR.seed_profiles(Ya, W)
    F = CR.global_factors(P, q) if q else None
    Yr = CR.residualize(Ya, F, device=DEV) if q else Ya
    Pr = CR.residualize(P, F, device=DEV) if q else P
    R, N = CR.profile_site_corr(Pr, Yr, device=DEV)
    nb = CR.site_site_topk(Yr, np.arange(10), k=5, nmin=5, device=DEV)
    return al, R, N, nb


def test_allowed_mask_excludes_heldout_and_shared_controls():
    meta = toy_meta()
    conds = list(meta.index)
    held = CR.fold_heldout(meta, "fold0")
    al = CR.allowed_mask(meta, conds, held)
    keys = CR.control_key(meta)
    held_keys = set(keys[held])
    for c, a in zip(conds, al):
        if c in held or keys[c] in held_keys or meta.loc[c, "qc_status"] != "ok":
            assert not a, c
    assert al.sum() > 0 or len(held_keys) == len(set(keys))  # toy: every key shared -> may be empty
    # a TEST / DEV-confirm condition may never be passed in
    m2 = meta.copy(); m2.loc["c7", "split"] = "test"
    with pytest.raises(ValueError):
        CR.allowed_mask(m2, conds, [])


def test_heldout_values_do_not_enter_fit():
    """Perturbing every held-out condition AND every shared-control condition leaves the fold's fit unchanged;
    perturbing an allowed condition changes it (so the test is sensitive)."""
    rng = np.random.default_rng(1)
    meta = toy_meta()
    # make control keys distinct per pair of conditions so that some non-held conditions survive
    meta["Control"] = [f"ctl{i // 2}" for i in range(len(meta))]
    conds = list(meta.index)
    S, C, K = 80, len(conds), 4
    Y = rng.normal(size=(S, C)); Y[rng.random(Y.shape) < 0.2] = np.nan
    W = sp.csr_matrix((np.ones(K * 5), (np.repeat(np.arange(K), 5), rng.choice(S, K * 5, replace=False))), shape=(K, S))
    for q in (0, 2):
        al, R, N, nb = fit(Y, meta, conds, "fold1", W, q)
        assert (~al).sum() > (meta["split"] == "fold1").sum()  # shared-control conditions also excluded
        Y2 = Y.copy()
        Y2[:, ~al] = rng.normal(size=(S, (~al).sum())) * 10      # held-out + shared-control + qc-failed columns
        _, R2, N2, nb2 = fit(Y2, meta, conds, "fold1", W, q)
        np.testing.assert_allclose(R, R2, atol=1e-6)
        np.testing.assert_array_equal(N, N2)
        np.testing.assert_array_equal(nb[0], nb2[0])
        Y3 = Y.copy()
        j = np.flatnonzero(al)[0]
        Y3[:, j] = rng.normal(size=S) * 10
        _, R3, _, _ = fit(Y3, meta, conds, "fold1", W, q)
        assert np.abs(R - R3).max() > 1e-3


def test_exp_b4_cache_is_fold_safe():
    """On the real DEV cache: the allowed lists written by exp_b4_precompute exclude every condition of the fold
    and every condition sharing (PMID, Control, Sample) with it."""
    import json
    from pathlib import Path
    p = Path(__file__).resolve().parents[1] / "cache/b4/allowed.json"
    if not p.exists():
        pytest.skip("B4 precompute not run")
    from kinact.data import load_dev_meta
    meta = load_dev_meta()
    keys = CR.control_key(meta)
    allowed = json.load(open(p))["allowed"]
    for f, al in allowed.items():
        held = set(meta.index[meta["split"] == f])
        hk = set(keys[list(held)])
        assert not (set(al) & held), f
        assert not any(keys[c] in hk for c in al), f
        assert all(meta.loc[c, "qc_status"] == "ok" for c in al)
        assert set(meta.loc[al, "split"]) <= {"fold0", "fold1", "fold2", "fold3", "fold4"}


def test_pairwise_corr_matches_pandas():
    rng = np.random.default_rng(2)
    Y = rng.normal(size=(30, 50)); Y[rng.random(Y.shape) < 0.3] = np.nan
    P = rng.normal(size=(3, 50)); P[rng.random(P.shape) < 0.2] = np.nan
    R, N = CR.profile_site_corr(P, Y, device=DEV)
    for k in range(3):
        for j in range(30):
            m = ~np.isnan(P[k]) & ~np.isnan(Y[j])
            assert N[k, j] == m.sum()
            assert R[k, j] == pytest.approx(np.corrcoef(P[k, m], Y[j, m])[0, 1], abs=1e-5)


def test_expansion_edges_rules():
    top = {"idx": np.array([[3, 1, 2, -1]]), "r": np.array([[0.9, -0.8, 0.2, 0.0]]),
           "rs": np.array([[0.8, -0.7, 0.1, 0.0]]), "n": np.array([[50, 50, 50, 0]])}
    e = CR.expansion_edges(top, ["K"], ["a", "b", "c", "d"], m=2)
    assert list(e["site"]) == ["d", "b"] and list(e["weight"]) == [0.8, -0.7]
    e = CR.expansion_edges(top, ["K"], ["a", "b", "c", "d"], m=2, sign="pos")
    assert list(e["site"]) == ["d", "c"]
