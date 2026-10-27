"""Toy-data tests with hand-computed answers for the frozen harness (docs/harness.md)."""
import json
import time

import numpy as np
import pandas as pd
import pytest

from kinact import bootstrap as bs
from kinact import controls as ctl
from kinact import ledger
from kinact.eval import (Evaluator, TestSplitLocked, bh, calibration, check_split_access,
                         evaluate, load_sealed_test, psp_bin, sealed_gt_hash)
from kinact.runtime import bench_runtime

K5 = ["A", "B", "C", "D", "E"]


def make_gt(spec, split="fold0"):
    """spec: list of (cond, group, source, [(kinase, cls, sign, conf)], n_targets, in_sel)."""
    P, C = [], []
    for cond, grp, src, rows, nt, sel in spec:
        for k, cls, s, conf in rows:
            P.append(dict(conditionID=cond, kinase=k, cls=cls, sign=s, confidence=conf,
                          rule_id="t", note=""))
        C.append(dict(conditionID=cond, PMID=1, source=src, group_id=grp, n_targets=nt,
                      in_M1sel=sel, role="gt", split=split))
    return pd.DataFrame(P), pd.DataFrame(C)


def ev_of(P, C, uni=K5, **kw):
    kw.setdefault("zero_policy", "ignore")
    return Evaluator(P, C, uni, ["fold0", "fold1"], **kw)


def top_or_bottom(conds, target_rank):
    """scores where target A (sign +1) ranks top (1) or bottom (0) among 5 kinases."""
    rows = {}
    for c, r in zip(conds, target_rank):
        rows[c] = [5.0, 1, 2, 3, 4] if r == 1 else [-5.0, 1, 2, 3, 4]
    return pd.DataFrame.from_dict(rows, orient="index", columns=K5)


# ------------------------------------------------------------------------------- basics
def test_pair_auroc_half_rule_by_hand():
    P, C = make_gt([("c1", "g1", "H", [("A", "direct", -1, "high")], 1, True)])
    S = pd.DataFrame({"A": [-3.0], "B": [-1.0], "C": [0.0], "D": [1.0], "E": [np.nan]},
                     index=["c1"])
    pr = ev_of(P, C).pairs(S)
    # t_A=3 beats t_B=1,t_C=0,t_D=-1 ; E unscored -> ½ : (3+0.5)/4
    assert pr["auroc"].iloc[0] == pytest.approx(0.875)
    assert pr["auroc_scored"].iloc[0] == pytest.approx(1.0)
    # ties count ½: B ties with A
    S2 = S.assign(B=-3.0)
    assert ev_of(P, C).pairs(S2)["auroc"].iloc[0] == pytest.approx((2 + 0.5 + 0.5) / 4)
    # unscored target -> ½
    assert ev_of(P, C).pairs(S.assign(A=np.nan))["auroc"].iloc[0] == 0.5


def test_oracle_antioracle_unscored():
    spec = [(f"c{i}", f"g{i % 3}", "HDO"[i % 3], [("A", "direct", (-1) ** i, "high"),
                                                   ("B", "direct", 1, "medium")], 2, True)
            for i in range(9)]
    P, C = make_gt(spec)
    ev = ev_of(P, C)
    conds = C["conditionID"]
    assert ev.evaluate(ctl.oracle(P, conds, K5), n_perm=50)["M1-sel"]["value"] == 1.0
    assert ev.evaluate(ctl.oracle(P, conds, K5, anti=True), n_perm=50)["M1-sel"]["value"] == 0.0
    nan = pd.DataFrame(np.nan, index=conds, columns=K5)
    r = ev.evaluate(nan, n_perm=50)
    assert r["M1-sel"]["value"] == 0.5 and r["M1-all"]["value"] == 0.5
    assert r["coverage_M1-sel"] == 0.0 and r["coverage_significant"]["n_significant"] == 0
    # low-confidence pairs are not targets
    P2 = P.assign(confidence="low")
    with pytest.raises(ValueError, match="no evaluable"):
        ev_of(P2, C).evaluate(nan)  # no evaluable conditions left


def big_gt(n_cond=300, n_k=200, seed=1, n_groups=60):
    rng = np.random.default_rng(seed)
    ks = [f"K{i}" for i in range(n_k)]
    spec = []
    for i in range(n_cond):
        nt = int(rng.integers(1, 4))
        tg = rng.choice(n_k, nt + 2, replace=False)
        rows = [(ks[j], "direct", int(rng.choice([-1, 1])), "high") for j in tg[:nt]]
        rows += [(ks[tg[nt]], "masked", 0, "high"), (ks[tg[nt + 1]], "downstream", 1, "high")]
        spec.append((f"c{i}", f"g{i % n_groups}", "HDO"[i % 3], rows, nt, nt <= 2))
    P, C = make_gt(spec)
    C.loc[C.index % 2 == 1, "split"] = "fold1"
    return P, C, ks


def test_random_near_half_and_shift_scale_invariance():
    P, C, ks = big_gt()
    ev = ev_of(P, C, ks)
    R = ctl.random_scores(C["conditionID"], ks, seed=3, coverage=0.6)
    r = ev.evaluate(R, n_perm=200)
    assert abs(r["M1-sel"]["value"] - 0.5) < 0.02
    assert abs(r["M1-all"]["value"] - 0.5) < 0.02
    assert r["coverage_M1-sel"] == pytest.approx(0.6, abs=0.05)
    r2 = ev.evaluate(R * 2.5 - 3.0, n_perm=200)
    for ep in ["M1-sel", "M1-all", "M1-sel-scored"]:
        assert r2[ep]["value"] == pytest.approx(r[ep]["value"], abs=1e-12)
    # false discoveries under the null stay small
    assert r["coverage_significant"]["n_significant"] <= 3


def test_averaging_order_differs_from_pooling():
    # kinase A: 3 conditions in g1 (AUROC 1) + 1 in g2 (AUROC 0); kinase B: 1 cond AUROC 1
    spec = [(c, g, "H", [("A", "direct", 1, "high")], 1, True)
            for c, g in [("c1", "g1"), ("c2", "g1"), ("c3", "g1"), ("c4", "g2")]]
    spec.append(("c5", "g3", "H", [("B", "direct", 1, "high")], 1, True))
    P, C = make_gt(spec)
    S = top_or_bottom(["c1", "c2", "c3", "c4"], [1, 1, 1, 0])
    S.loc["c5"] = [1, 9, 2, 3, 4]  # B top
    r = ev_of(P, C).evaluate(S, n_perm=50)
    # hierarchy: A = mean(g1=1, g2=0) = .5 ; B = 1 ; macro = .75
    assert r["M1-sel"]["value"] == pytest.approx(0.75)
    pooled_a = 0.75  # naive pooling over A's 4 pairs
    assert (pooled_a + 1) / 2 == pytest.approx(0.875) != r["M1-sel"]["value"]
    pk = r["per_kinase"].set_index("kinase")
    assert pk.loc["A", "M1-sel"] == pytest.approx(0.5)


def test_masked_and_downstream_excluded_from_negatives():
    P, C = make_gt([("c1", "g1", "H", [("A", "direct", 1, "high"), ("B", "masked", 0, "high"),
                                       ("C", "downstream", 1, "high")], 1, True)])
    S = pd.DataFrame([[1.0, 9.0, 8.0, 0.0, -1.0]], index=["c1"], columns=K5)
    pr = ev_of(P, C).pairs(S)
    assert pr["n_neg"].iloc[0] == 2 and pr["auroc"].iloc[0] == 1.0
    # downstream secondary positive set
    ds = ev_of(P, C, positive_cls="downstream").pairs(S)
    assert list(ds["kinase"]) == ["C"] and ds["auroc"].iloc[0] == 1.0


def test_one_over_n_weighting():
    P, C = make_gt([
        ("c1", "g1", "H", [("A", "direct", 1, "high")], 1, True),
        ("c2", "g1", "H", [("A", "direct", 1, "high"), ("B", "direct", 1, "high")], 2, True)])
    S = pd.DataFrame([[5.0, 1, 2, 3, 4], [-5.0, 9, 2, 3, 4]], index=["c1", "c2"], columns=K5)
    r = ev_of(P, C).evaluate(S, n_perm=50)
    pk = r["per_kinase"].set_index("kinase")
    assert pk.loc["A", "M1-sel"] == pytest.approx(0.5)
    assert pk.loc["A", "M1-all"] == pytest.approx((1 * 1 + 0.5 * 0) / 1.5)  # 2/3
    # B in c2: negatives C,D,E (A is a positive) -> top -> 1
    assert r["M1-all"]["value"] == pytest.approx((2 / 3 + 1) / 2)


def test_m1sel_by_source_sign_unseen_harmony():
    spec = [("c1", "g1", "H", [("A", "direct", 1, "high")], 1, True),
            ("c2", "g2", "D", [("B", "direct", -1, "high")], 1, True),
            ("c3", "g3", "O", [("C", "direct", 1, "high")], 11, False)]
    P, C = make_gt(spec)
    S = pd.DataFrame([[5.0, 1, 2, 3, 4], [1.0, 9, 2, 3, 4], [1.0, 2, 9, 3, 4]],
                     index=["c1", "c2", "c3"], columns=K5)
    bins = pd.DataFrame({"conditionID": ["c1", "c2", "c3"], "kinase": ["A", "B", "C"],
                         "n_psp_quantified": [0, 7, 30]})
    ev = ev_of(P, C, bins=bins, train_pos_kinases={"A"})
    r = ev.evaluate(S, n_perm=50)
    assert r["M1-sel"]["value"] == pytest.approx(0.5)  # A=1, B=0 (wrong sign)
    assert r["by_source"]["H"]["value"] == 1 and r["by_source"]["D"]["value"] == 0
    assert "O" not in r["by_source"] or np.isnan(r["by_source"]["O"]["value"])
    assert r["by_sign"][-1]["value"] == 0 and r["by_sign"][1]["value"] == 1
    assert r["T-unseen"]["value"] == 0.0  # only B unseen in M1-sel
    assert r["harmony"]["per_bin"] == {"0": 1.0, "5-19": 0.0}
    assert r["harmony"]["value"] == pytest.approx(0.5)
    assert list(psp_bin([0, 1, 4, 4.5, 5, 19, 20])) == ["0", "1-4", "1-4", "1-4", "5-19", "5-19", ">=20"]
    # frozen (kinase, group_id, bin) table form
    b2 = pd.DataFrame({"kinase": ["A", "B"], "group_id": ["g1", "g2"], "bin": [">=20", ">=20"]})
    assert ev_of(P, C, bins=b2).evaluate(S, n_perm=10)["harmony"]["per_bin"] == {">=20": 0.5}


def test_m2_matched_universe():
    P, C = make_gt([("c1", "g1", "H", [("A", "direct", 1, "high")], 1, True),
                    ("c2", "g2", "H", [("B", "direct", 1, "high")], 1, True)])
    SA = pd.DataFrame([[5.0, 1, 2, 3, 4], [0.0, 2.5, 2, 3, 4]], index=["c1", "c2"], columns=K5)
    SB = pd.DataFrame([[3.5, np.nan, np.nan, 3, 4], [np.nan] * 5], index=["c1", "c2"], columns=K5)
    m = ev_of(P, C).m2(SA, SB)
    # only c1 survives; universe {A,D,E}: A (5) beats D,E in A; 3.5 < 4 in B -> ½
    assert m["n_pairs"] == 1 and m["A"] == 1.0 and m["B"] == 0.5 and m["delta"] == 0.5


def test_holm_and_bh_examples():
    np.testing.assert_allclose(bs.holm([0.01, 0.04, 0.03, 0.005]), [0.03, 0.06, 0.06, 0.02])
    np.testing.assert_allclose(bh([0.01, 0.04, 0.03, 0.005]), [0.02, 0.04, 0.04, 0.02])


def test_bootstrap_reproducible_and_matches_bruteforce():
    P, C, ks = big_gt(n_cond=120, n_k=40, n_groups=30)
    ev = ev_of(P, C, ks)
    A = ctl.random_scores(C["conditionID"], ks, seed=1, coverage=0.7)
    O = ctl.oracle(P, C["conditionID"], ks)
    Bs = O + ctl.random_scores(C["conditionID"], ks, seed=2) * 0.8
    p1 = bs.BootstrapPlan.from_evaluator(ev, B=200, seed=0)
    p2 = bs.BootstrapPlan.from_evaluator(ev, B=200, seed=0)
    p3 = bs.BootstrapPlan.from_evaluator(ev, B=200, seed=1)
    assert np.array_equal(p1.W, p2.W) and not np.array_equal(p1.W, p3.W)
    # strata respected: each resample has as many groups per stratum as the original
    for s in p1.strata.unique():
        cols = (p1.strata == s).values
        assert np.all(p1.W[:, cols].sum(1) == cols.sum())
    d1 = bs.paired_delta(ev, Bs, A, "M1-sel", p1)
    d2 = bs.paired_delta(ev, Bs, A, "M1-sel", p2)
    assert d1["estimate"] == d2["estimate"] and np.array_equal(d1["replicates"], d2["replicates"])
    assert d1["estimate"] > 0.1 and d1["p_one_sided"] == pytest.approx(1 / 201)
    # brute force: materialise resample 7 by duplicating groups, recompute from scratch
    pa = ev.pairs(Bs)
    for ep in ["M1-sel", "M1-all"]:
        reps = bs.metric_replicates(ev, pa, ep, p1)
        b = 7
        parts = []
        for gi, g in enumerate(p1.groups):
            for r in range(int(p1.W[b, gi])):
                parts.append(pa[pa["group"] == g].assign(group=f"{g}#{r}"))
        assert ev.metric(pd.concat(parts), ep)["value"] == pytest.approx(reps[b], abs=1e-9)


def test_harmony_bootstrap_bruteforce():
    P, C, ks = big_gt(n_cond=90, n_k=30, n_groups=20)
    rng = np.random.default_rng(0)
    bins = P[P["cls"] == "direct"][["conditionID", "kinase"]].assign(
        n_psp_quantified=rng.integers(0, 30, (P["cls"] == "direct").sum()))
    ev = ev_of(P, C, ks, bins=bins)
    S = ctl.random_scores(C["conditionID"], ks, seed=4, coverage=0.8)
    plan = bs.BootstrapPlan.from_evaluator(ev, B=50)
    pa = ev.pairs(S)
    reps = bs.metric_replicates(ev, pa, "harmony", plan)
    b = 3
    parts = [pa[pa["group"] == g].assign(group=f"{g}#{r}")
             for gi, g in enumerate(plan.groups) for r in range(int(plan.W[b, gi]))]
    assert ev.metric(pd.concat(parts), "harmony")["value"] == pytest.approx(reps[b], abs=1e-9)


def test_sign_flip_check():
    X = pd.DataFrame(np.random.default_rng(0).standard_normal((50, 4)),
                     index=[f"s{i}" for i in range(50)], columns=list("wxyz"))
    lib = {"A": [f"s{i}" for i in range(10)], "B": [f"s{i}" for i in range(10, 30)]}

    def zscore(X, lib):
        return pd.DataFrame({k: X.loc[v].mean() * np.sqrt(len(v)) for k, v in lib.items()})

    def absmethod(X, lib):
        return zscore(X, lib).abs()
    assert ctl.sign_flip_check(zscore, X, lib)["passed"]
    assert not ctl.sign_flip_check(absmethod, X, lib)["passed"]


def test_controls_prior_shift_permutation_shuffle():
    P, C = make_gt([("c1", "g1", "H", [("A", "direct", -1, "high")], 1, True),
                    ("c2", "g2", "H", [("A", "direct", -1, "high"), ("B", "direct", -1, "high")], 2, True),
                    ("c3", "g3", "O", [("C", "direct", 1, "high")], 1, True),
                    ("c4", "g4", "H", [("B", "direct", -1, "high")], 1, True)])
    C.loc[C["conditionID"] == "c4", "split"] = "fold1"
    f = ctl.training_frequency(P, C, ["c1", "c2", "c3"], K5)
    assert f.to_dict() == pytest.approx({"A": 2 / 3, "B": 1 / 3, "C": 1 / 3, "D": 0, "E": 0})
    po = ctl.prior_only(P, C, ["c1", "c2", "c3"], ["c3", "c4"], K5)
    assert po.loc["c4", "A"] == pytest.approx(-2 / 3) and po.loc["c3", "C"] == pytest.approx(1 / 3)
    Ct = C.assign(split=np.where(C["conditionID"] == "c3", "test", C["split"]))
    with pytest.raises(ValueError):
        ctl.training_frequency(P, Ct, ["c1", "c3"], K5)
    X = pd.DataFrame({"c1": [1.0, 3.0, np.nan], "c2": [-2.0, 0.0, 2.0]})
    sp = ctl.shift_prior(X, f)
    assert sp.loc["c1", "A"] == pytest.approx(2 * 2 / 3) and sp.loc["c2", "A"] == 0
    Xp = ctl.site_permutation(X, seed=0)
    assert Xp.isna().equals(X.isna()) and sorted(Xp["c2"]) == sorted(X["c2"])
    lib = pd.DataFrame({"kinase": list("AAABBC"), "site": ["s1", "s2", "s3", "s1", "s4", "s5"],
                        "weight": 1.0})
    sh = ctl.shuffled_library(lib, seed=3)
    assert sh.groupby("kinase").size().equals(lib.groupby("kinase").size())
    assert not sh.duplicated(["kinase", "site"]).any() and set(sh["site"]) == set(lib["site"])


def test_zero_encoding_flag():
    P, C = make_gt([("c1", "g1", "H", [("A", "direct", 1, "high")], 1, True)])
    S = pd.DataFrame([[1.0, 0, 0, 0, 2]], index=["c1"], columns=K5)
    with pytest.raises(ValueError):
        Evaluator(P, C, K5, ["fold0"]).pairs(S)
    with pytest.warns(UserWarning):
        Evaluator(P, C, K5, ["fold0"], zero_policy="warn").pairs(S)


def test_coverage_significant_detects_signal():
    P, C, ks = big_gt(n_cond=300, n_k=60, n_groups=100)
    O = ctl.oracle(P, C["conditionID"], ks)
    S = O * 3 + ctl.random_scores(C["conditionID"], ks, seed=5)
    r = ev_of(P, C, ks).evaluate(S, n_perm=500)
    cs = r["coverage_significant"]
    assert cs["n_significant"] >= 0.8 * cs["n_kinases"]


def test_calibration():
    Z = ctl.random_scores([f"n{i}" for i in range(200)], [f"K{i}" for i in range(100)], seed=0)
    ns = pd.DataFrame(np.tile(np.arange(100) % 30, (200, 1)), index=Z.index, columns=Z.columns)
    cal = calibration(Z, list(Z.index), n_substrates=ns)
    assert cal["frac_exceed"] == pytest.approx(0.05, abs=0.005)
    assert set(cal["per_bin"]) == {"0", "1-4", "5-19", ">=20"}
    assert cal["ks_uniform"] < 0.02
    assert calibration(Z * 3, list(Z.index))["frac_exceed"] > 0.4


# ------------------------------------------------------------------------------- guards
def test_test_split_guard(tmp_path):
    P, C = make_gt([("c1", "g1", "H", [("A", "direct", 1, "high")], 1, True)], split="test")
    fz = tmp_path / "FREEZE.json"
    with pytest.raises(TestSplitLocked):
        Evaluator(P, C, K5, None, freeze=fz)
    with pytest.raises(TestSplitLocked):
        Evaluator(P, C, K5, ["test"], freeze=fz)
    with pytest.raises(TestSplitLocked):
        evaluate(pd.DataFrame(), P, C, ["fold0", "test"], universe=K5, freeze=fz)
    with pytest.raises(TestSplitLocked):
        check_split_access(["test"], fz)
    Evaluator(P, C, K5, ["fold0"], freeze=fz)  # no test rows selected -> fine
    fz.write_text("{}")
    Evaluator(P, C, K5, ["test"], freeze=fz)
    # sealed loader: needs FREEZE + matching hash
    sd = tmp_path / "sealed"
    sd.mkdir()
    P.to_csv(sd / "gt_pairs.tsv", sep="\t", index=False)
    C.to_csv(sd / "gt_conditions.tsv", sep="\t", index=False)
    with pytest.raises(TestSplitLocked):
        load_sealed_test(sd, tmp_path / "nofreeze.json")
    with pytest.raises(TestSplitLocked):
        load_sealed_test(sd, fz)
    fz.write_text(json.dumps({"gt_hash": sealed_gt_hash(sd)}))
    assert len(load_sealed_test(sd, fz)["gt_pairs"]) == 1


def test_ledger_append_only(tmp_path):
    lp = tmp_path / "dev_ledger.tsv"
    r1 = ledger.append("zscore", {"lib": "psp"}, "dev", {"M1-sel": {"value": 0.51}},
                       branch="b1", folds=[0, 1, 2, 3, 4], path=lp, commit="abc")
    first = lp.read_bytes()
    r2 = ledger.append("zscore", {"lib": "psp+"}, "dev", {"M1-sel": {"value": 0.52}},
                       branch="b1", path=lp, commit="abc")
    r3 = ledger.append("zscore", {"lib": "psp"}, "dev", {}, branch="b1", path=lp, commit="abc")
    r4 = ledger.append("mlm", {"lib": "psp"}, "dev", {}, branch="b2", path=lp, commit="abc")
    assert lp.read_bytes().startswith(first)
    assert [r["n_configs_tried"] for r in (r1, r2, r3, r4)] == [1, 2, 2, 1]
    rows = ledger.read_ledger(lp)
    assert len(rows) == 4 and rows[0]["folds"] == "0,1,2,3,4"
    assert json.loads(rows[0]["metrics_json"])["M1-sel"]["value"] == 0.51
    assert lp.read_text().count("timestamp") == 1
    with pytest.raises((PermissionError, TestSplitLocked)):
        ledger.append("x", {}, "test", {}, path=lp)


def test_runtime_bench():
    X = pd.DataFrame(np.random.default_rng(0).standard_normal((2000, 20)))

    class M:
        def prepare(self, X, lib):
            return X.to_numpy()

        def score(self, st):
            return st.T @ np.ones((st.shape[0], 5))
    r = bench_runtime(M(), X, n_rep=3, threads=(1, 2))
    assert {"load_mapping_s", "scoring_s", "ms_per_condition"} <= set(r["threads_1"])
    r2 = bench_runtime(lambda X, lib: X.mean(), X, n_rep=3, threads=(1,))
    assert r2["threads_1"]["total_s"] >= 0


def test_speed_paired_bootstrap_1000x500():
    P, C, ks = big_gt(n_cond=1000, n_k=500, n_groups=250, seed=7)
    t0 = time.perf_counter()
    ev = ev_of(P, C, ks)
    A = ctl.random_scores(C["conditionID"], ks, seed=1, coverage=0.5)
    B = ctl.random_scores(C["conditionID"], ks, seed=2, coverage=0.9)
    plan = bs.BootstrapPlan.from_evaluator(ev, B=2000)
    d = bs.paired_delta(ev, A, B, "M1-sel", plan)
    dt = time.perf_counter() - t0
    print(f"\n[timing] index+plan+paired ΔM1-sel B=2000 1000x500: {dt:.2f}s; SE={d['se']:.4f}")
    assert dt < 60
