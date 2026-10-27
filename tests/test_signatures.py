"""Tests for branch B7 perturbation-response signatures (kinact.explore.signatures): leakage, sign symmetry,
no-intercept / alpha=0 identity, controls."""
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from kinact.explore import signatures as SG


def toy(seed=0, n_sites=400, n_groups=12, per=3):
    rng = np.random.default_rng(seed)
    conds = [f"c{g}_{i}" for g in range(n_groups) for i in range(per)]
    group = pd.Series({c: f"G{c.split('_')[0][1:]}" for c in conds})
    Y = rng.standard_normal((n_sites, len(conds))).astype(np.float32)
    Y[rng.random(Y.shape) < 0.1] = np.nan
    # kinase K1 is a positive (inhibited, s=-1) in groups 0..5; its response = sites 0..19 go down
    rows = []
    for j, c in enumerate(conds):
        g = int(c.split("_")[0][1:])
        if g < 6:
            Y[:20, j] -= 2.0
            rows.append((c, "K1", -1, group[c]))
        else:
            rows.append((c, "K2", 1, group[c]))
    pos = pd.DataFrame(rows, columns=["conditionID", "kinase", "sign", "group"])
    return Y, conds, group, pos


def test_signature_recovers_planted_sites_and_is_signed():
    Y, conds, group, pos = toy()
    s = SG.fit_signatures(Y, conds, conds, pos, group, kappa=0.0, n_max=20, contrast=False)
    # K1 inhibited -> sites drop; sign-aligned (s=-1) response is positive at the planted sites
    assert set(s.sites["K1"]) == set(range(20))
    assert (s.weights["K1"] > 0).all()
    assert "K2" in s.sites and s.n_groups["K1"] == 6


def test_heldout_conditions_never_enter_fit():
    Y, conds, group, pos = toy()
    held = [c for c in conds if group[c] in ("G0", "G7")]
    train = [c for c in conds if c not in held]
    ptr = pos[pos["conditionID"].isin(train)]
    a = SG.fit_signatures(Y, conds, train, ptr, group, eval_conds=held, kappa=1.0)
    Y2 = Y.copy()
    ci = pd.Index(conds).get_indexer(held)
    Y2[:, ci] = np.random.default_rng(5).standard_normal((Y.shape[0], len(ci))) * 50
    b = SG.fit_signatures(Y2, conds, train, ptr, group, eval_conds=held, kappa=1.0)
    for k in a.sites:
        np.testing.assert_array_equal(a.sites[k], b.sites[k])
        np.testing.assert_allclose(a.weights[k], b.weights[k])


def test_leak_guards_raise():
    Y, conds, group, pos = toy()
    held = [c for c in conds if group[c] == "G0"]
    # evaluated condition in training set
    with pytest.raises(SG.SignatureLeak):
        SG.fit_signatures(Y, conds, conds, pos, group, eval_conds=held)
    # a training condition sharing the drug group of an evaluated condition
    train = [c for c in conds if c != held[0]]
    with pytest.raises(SG.SignatureLeak):
        SG.fit_signatures(Y, conds, train, pos[pos["conditionID"].isin(train)], group, eval_conds=[held[0]])
    # positives table containing a non-training condition
    train = [c for c in conds if c not in held]
    with pytest.raises(SG.SignatureLeak):
        SG.fit_signatures(Y, conds, train, pos, group, eval_conds=held)


def test_training_conditions_real_dev_folds_disjoint():
    """On the real DEV GT: for every outer and inner (f, g) split, no condition of a held-out fold, and no condition
    sharing a drug group with one, is a signature training condition."""
    from kinact.eval import load_gt
    from kinact.data import load_dev_meta
    _, C, _, _ = load_gt()
    meta = load_dev_meta()
    folds = [f"fold{i}" for i in range(5)]
    for f in folds:
        for g in folds:
            held = {f, g}
            ev = list(meta.index[meta["split"].isin(held)])
            tr = SG.training_conditions(C, meta, [x for x in folds if x not in held], ev)
            assert len(tr) > 100
            assert not set(tr) & set(ev)
            assert not set(meta.loc[tr, "group_id"]) & set(meta.loc[ev, "group_id"])
            assert set(meta.loc[tr, "split"]) <= set(folds) - held


def test_sign_flip_invariance_of_fit_and_odd_scores():
    Y, conds, group, pos = toy()
    a = SG.fit_signatures(Y, conds, conds, pos, group, kappa=1.0, contrast=True)
    pf = pos.assign(sign=-pos["sign"])
    b = SG.fit_signatures(-Y, conds, conds, pf, group, kappa=1.0, contrast=True)
    for k in a.sites:
        np.testing.assert_array_equal(a.sites[k], b.sites[k])
        np.testing.assert_allclose(a.weights[k], b.weights[k], rtol=1e-5)


@pytest.mark.skipif(not __import__("torch").cuda.is_available(), reason="needs CUDA")
def test_blocks_alpha0_identity_and_oddness():
    import torch
    rng = np.random.default_rng(1)
    S, p, B = 60, 5, 4
    Wb = sp.csr_matrix((rng.random((S, p)) < 0.2).astype(float))
    Ws = sp.csr_matrix(np.where(rng.random((S, p)) < 0.1, rng.standard_normal((S, p)), 0.0))
    Y = rng.standard_normal((S, B)).astype(np.float32); M = (rng.random((S, B)) < 0.8).astype(np.float32)
    Y = Y * M
    Wd = Wb.toarray()
    base = dict(G=torch.as_tensor(np.stack([Wd.T @ np.diag(M[:, b]) @ Wd for b in range(B)]), dtype=torch.float32, device="cuda"),
                h=torch.as_tensor((Wd.T @ (M * Y)).T, dtype=torch.float32, device="cuda"),
                g0=torch.as_tensor((Wd.T @ M).T, dtype=torch.float32, device="cuda"), n=np.zeros((B, p)))
    sb = SG.sig_blocks(Wb, Ws, Y, M)
    G0, h0, *_ = SG.combine_sig(base, sb, 0.0)
    assert torch.equal(G0, base["G"]) and torch.equal(h0, base["h"])
    a = 0.7
    G, h, g0, _ = SG.combine_sig(base, sb, a)
    Wt = Wd + a * Ws.toarray()
    for b in range(B):
        np.testing.assert_allclose(G[b].cpu().numpy(), Wt.T @ np.diag(M[:, b]) @ Wt, atol=1e-4)
        np.testing.assert_allclose(h[b].cpu().numpy(), Wt.T @ (M[:, b] * Y[:, b]), atol=1e-4)
    # oddness: -y -> -h, G unchanged (no intercept: score is linear in the condition's data)
    sbn = SG.sig_blocks(Wb, Ws, -Y, M)
    torch.testing.assert_close(sbn["h_s"], -sb["h_s"]); torch.testing.assert_close(sbn["G_ss"], sb["G_ss"])


def test_controls_preserve_size_and_weights():
    Y, conds, group, pos = toy()
    s = SG.fit_signatures(Y, conds, conds, pos, group, kappa=1.0, n_max=30)
    sh = s.shuffled(0)
    gc = s.generic_control()
    for k in s.sites:
        assert len(sh.sites[k]) == len(s.sites[k]) and len(set(sh.sites[k])) == len(sh.sites[k])
        np.testing.assert_allclose(sh.weights[k], s.weights[k])
        assert set(sh.sites[k]) <= set(s.eligible_sites)
        np.testing.assert_allclose(np.linalg.norm(gc.weights[k]), np.linalg.norm(s.weights[k]))
    W = s.matrix(Y.shape[0], ["K1", "K2", "K3"], 10, alpha=2.0)
    assert W.shape == (Y.shape[0], 3) and W[:, 2].nnz == 0 and W[:, 0].nnz == 10
