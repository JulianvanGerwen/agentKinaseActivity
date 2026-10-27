"""Harness sanity for the external benchmarKIN sets (synthetic scores only; no real method is run)."""
import numpy as np
import pandas as pd
import pytest

from kinact import external_eval as X

pytestmark = pytest.mark.skipif(not (X.EXT / "pert_experiments.tsv").exists(), reason="run scripts/prep_external.py")


def test_site_keys_and_targets():
    for n in X.PERT_SETS:
        d, e, t = X.load_pert(n)
        assert d.index.str.fullmatch(r"[A-Z0-9]+(-\d+)?_[STY]\d+").all()
        assert set(t.kinase) <= set(X.load_universe())
        assert d.shape[1] == len(e) and t.experiment.isin(e.experiment).all()
    assert len(X.load_pert("tyrosine")[1]) == 18
    # key convention matches the compendium site universe
    sk = set(pd.read_parquet(X.ROOT / "cache/sites.parquet", columns=["site_key"]).site_key)
    idx = X.load_pert("tyrosine")[0].index
    assert idx.isin(sk).mean() > 0.7


def test_hernandez_new_not_in_compendium():
    comp = set(pd.read_csv(X.ROOT / "data/01_phosphoproteomics_conditions.csv", dtype=str).PMID)
    _, e, _ = X.load_pert("hernandez_new")
    assert len(e) == 15 and not (set(e.pmid) & comp)


@pytest.mark.parametrize("name", X.PERT_SETS)
def test_perturbation_random_and_oracle(name):
    _, e, tg = X.load_pert(name)
    R = np.mean([X.m1_perturbation(X.random_scores(e.experiment, seed=s), name)["M1-all"]["value"] for s in range(20)])
    assert abs(R - 0.5) < 0.05
    O = X.oracle_pert_scores(name)
    assert X.m1_perturbation(O, name)["M1-all"]["value"] > 0.99
    assert X.benchmarkin_perturbation(O, name, n_iter=200)["auroc_full"] > 0.99
    rb = np.mean([X.benchmarkin_perturbation(X.random_scores(e.experiment, seed=s), name, n_iter=100)["auroc_mean"]
                  for s in range(10)])
    assert abs(rb - 0.5) < 0.08
    # sign flip of the oracle -> ~0
    assert X.m1_perturbation(-O, name)["M1-all"]["value"] < 0.01


def test_unscored_is_half():
    _, e, _ = X.load_pert("tyrosine")
    S = pd.DataFrame(np.nan, index=e.experiment, columns=X.load_universe())
    assert X.m1_perturbation(S, "tyrosine")["M1-all"]["value"] == pytest.approx(0.5)


def test_tumour_random_and_oracle():
    data, gs = X.load_tumour()
    assert set(gs.cohort) == set(X.COHORTS)
    O = X.oracle_tumour_scores(gs, data)
    assert X.m1_tumour(O, gs)["M1-kinase"] > 0.99
    assert X.benchmarkin_tumour(O, gs, n_iter=100)["auroc_mean"] > 0.99
    r1 = [X.m1_tumour({c: X.random_scores(d.columns, seed=s) for c, d in data.items()}, gs)["M1-kinase"] for s in range(3)]
    assert abs(np.mean(r1) - 0.5) < 0.03
    R = {c: X.random_scores(d.columns, seed=5) for c, d in data.items()}
    assert abs(X.benchmarkin_tumour(R, gs, n_iter=100)["auroc_mean"] - 0.5) < 0.03


def test_metric_plumbing():
    rng = np.random.default_rng(0)
    p, n = rng.normal(1, 1, 300), rng.normal(0, 1, 500)
    from scipy.stats import mannwhitneyu
    assert X.auroc(p, n) == pytest.approx(mannwhitneyu(p, n).statistic / (300 * 500))
    assert abs(X.mc_auroc(p, n, 500).mean() - X.auroc(p, n)) < 0.02
    A = pd.DataFrame(rng.normal(size=(50, 3)))
    assert np.allclose(np.sqrt((X.rms_scale(A) ** 2).sum(axis=0) / 49), 1)
    with pytest.raises(RuntimeError):
        X.run_external(lambda d: d)
