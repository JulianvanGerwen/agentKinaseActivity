"""Guards of the single TEST script (scripts/run_test.py). Never touches cache/sealed/ or FREEZE.json."""
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def rt():
    spec = importlib.util.spec_from_file_location("run_test", ROOT / "scripts/run_test.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_refuses_without_freeze(rt, monkeypatch, tmp_path):
    monkeypatch.setattr(rt, "FREEZE", tmp_path / "FREEZE.json")
    monkeypatch.delenv("KINACT_FREEZE_PATH", raising=False)
    with pytest.raises(SystemExit, match="does not exist"):
        rt.check_freeze()


def test_refuses_on_env_redirect(rt, monkeypatch, tmp_path):
    monkeypatch.setenv("KINACT_FREEZE_PATH", str(tmp_path / "other.json"))
    with pytest.raises(SystemExit, match="points away"):
        rt.check_freeze()


def test_verify_freeze_rules(rt):
    want = {"gt_hash": "a", "gt_freeze_sha256": "b", "config_sha256": "c", "model_manifest_sha256": "d",
            "code_commit": "e"}
    ok = {"mismatch": {}}
    rt.verify_freeze(dict(want), want, [], ok)                                   # passes
    for k in want:
        with pytest.raises(SystemExit, match="hash_mismatch"):
            rt.verify_freeze({**want, k: "x"}, want, [], ok)
        F = dict(want); F.pop(k)
        with pytest.raises(SystemExit):
            rt.verify_freeze(F, want, [], ok)
    with pytest.raises(SystemExit, match="modified_tracked_files"):
        rt.verify_freeze(dict(want), want, [" M src/kinact/core.py"], ok)
    with pytest.raises(SystemExit, match="artifact_mismatch"):
        rt.verify_freeze(dict(want), want, [], {"mismatch": {"edges.tsv.gz": "h"}})


def test_manifest_artifacts_match(rt):
    if not (ROOT / "models/core_v3/edges.tsv.gz").exists():
        pytest.skip("core-v3 artefacts not built")
    assert rt.check_manifest_artifacts()["mismatch"] == {}


def test_frozen_comparators(rt):
    labs = [x[0] for x in rt.HOLM_COMPARATORS]
    assert labs == ["z-score", "KSEA", "JUMPsem", "best-other-baseline", "prior-only"]
    assert rt.B == 2000 and rt.SEED == 0 and rt.FINAL == "core-v3"


def test_external_block_synthetic(rt, monkeypatch):
    """external_block end-to-end on tiny SYNTHETIC sets (loaders patched; no external data is read)."""
    import numpy as np
    import pandas as pd
    from kinact import external_eval as XE
    if not (ROOT / "models/core_v3/edges.tsv.gz").exists():
        pytest.skip("core-v3 artefacts not built")
    rng = np.random.default_rng(0)
    edges = pd.read_csv(ROOT / "models/core_v1/edges.tsv.gz", sep="\t")
    sites = list(pd.unique(edges["site"]))[:3000]
    U = XE.load_universe()
    kins = ["EGFR", "ABL1", "MAPK1", "AKT1", "CDK1", "MTOR"]

    def load_pert(name, ext=None):
        ex = [f"{name}_e{i}" for i in range(8)]
        d = pd.DataFrame(rng.normal(size=(len(sites), 8)), index=sites, columns=ex)
        e = pd.DataFrame({"experiment": ex, "pmid": [1, 1, 2, 2, 3, 3, 4, 4], "drug": list("aabbccdd"),
                          "sign": -1, "drug_in_dev": ["True", "False"] * 4})
        t = pd.DataFrame({"experiment": ex, "kinase": [kins[i % 6] for i in range(8)], "sign": -1})
        return d, e, t

    def load_tumour(ext=None):
        data, rows = {}, []
        for c in XE.COHORTS:
            sm = [f"{c}_s{i}" for i in range(40)]
            data[c] = pd.DataFrame(rng.normal(size=(len(sites), 40)), index=sites, columns=sm)
            for k in kins:
                for i, s in enumerate(sm[:10]):
                    rows.append((c, s, k, int(i < 5)))
        return data, pd.DataFrame(rows, columns=["cohort", "sample", "kinase", "label"])
    monkeypatch.setattr(XE, "load_pert", load_pert)
    monkeypatch.setattr(XE, "load_tumour", load_tumour)
    monkeypatch.setattr(rt, "phosx_like", lambda X, *a, **k: rt.zscore(X, "curated_union", 1))
    out = rt.external_block(lambda *a: None)
    assert set(out["paired_vs_final"]) == {"tyrosine", "hernandez_new", "tumour"}
    assert "random" in out["results"] and "core-v3" in out["results"]
    t = out["paired_vs_final"]["tumour"]["core-v1"]["M1-kinase"]
    assert np.isfinite(t["estimate"]) and t["n_kinases"] == len(kins)
    import json
    dry = ROOT / "results/test_dryrun/test_results.json"
    if dry.exists():
        from kinact.eval import to_jsonable
        R = json.loads(dry.read_text())
        R["external"] = json.loads(json.dumps(to_jsonable(out), default=str))
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            rt.write_ledger(Path(d), R, dry=False)
            txt = (Path(d) / "LEDGER.md").read_text()
        assert "## External sets" in txt and "tumour | core-v1 | M1-kinase" in txt


def test_verify_freeze_refuses_on_input_mismatch(rt):
    want = {"gt_hash": "a"}
    rt.verify_freeze(dict(want), want, [], {"mismatch": {}}, {"mismatch": {}})
    with pytest.raises(SystemExit, match="input_mismatch"):
        rt.verify_freeze(dict(want), want, [], {"mismatch": {}},
                         {"mismatch": {"cache/quant_wide_f32.npy": {"frozen": "x", "disk": "y"}}})


def test_frozen_inputs_match(rt):
    if not (ROOT / "cache/quant_wide_f32.npy").exists():
        pytest.skip("compendium cache not present")
    r = rt.check_frozen_inputs()
    assert r["n_checked"] > 20 and r["mismatch"] == {}


def _toy_ev():
    import pandas as pd
    from kinact.eval import Evaluator
    kins = [f"K{i}" for i in range(8)]
    conds = [f"c{i}" for i in range(12)]
    tgt = {c: kins[i % 4] for i, c in enumerate(conds)}
    tgt["c0"] = "K7"                                   # K7: the only 'unseen' kinase, in a single group g0
    gp = pd.DataFrame({"conditionID": conds, "kinase": [tgt[c] for c in conds], "cls": "direct", "sign": -1,
                       "confidence": "high"})
    gc = pd.DataFrame({"conditionID": conds, "split": "x", "role": "gt", "n_targets": 1, "in_M1sel": True,
                       "group_id": [f"g{i // 2}" for i in range(12)], "source": "H"})
    return Evaluator(gp, gc, kins, train_pos_kinases=kins[:4], zero_policy="ignore"), conds, kins


def test_masked_delta_drops_resamples_without_subset(rt):
    import numpy as np
    import pandas as pd
    ev, conds, kins = _toy_ev()
    rng = np.random.default_rng(1)
    A = pd.DataFrame(rng.normal(size=(12, 8)), index=conds, columns=kins); A.loc["c0", "K7"] = -9.0
    Bm = pd.DataFrame(rng.normal(size=(12, 8)), index=conds, columns=kins)
    plan = rt.bs.BootstrapPlan.from_evaluator(ev, B=400, seed=0)
    uns = lambda P: P["kinase"].eq("K7").to_numpy()
    d = rt.delta(ev, A, Bm, "M1-sel", plan, uns)
    n_drawn = int((plan.W[:, list(plan.groups).index("g0")] > 0).sum())
    assert d["n_resamples"] == n_drawn < plan.B
    # single (kinase, group) unit: every drawn resample reproduces the point estimate exactly -> degenerate CI
    assert abs(d["ci_low"] - d["estimate"]) < 1e-12 and abs(d["ci_high"] - d["estimate"]) < 1e-12
    none = rt.delta(ev, A, Bm, "M1-sel", plan, lambda P: np.zeros(len(P), bool))
    assert none["n_resamples"] == 0 and np.isnan(none["p_one_sided"])


def test_secondary_block_records_failure(rt):
    def boom():
        raise RuntimeError("R died")
    r = rt.secondary_block("external", boom, lambda *a: None)
    assert "R died" in r["error"] and "Traceback" in r["traceback"]
    assert rt.secondary_block("x", lambda: {"ok": 1}, lambda *a: None) == {"ok": 1}


def test_ledger_tolerates_failed_secondary_blocks(rt, tmp_path):
    import json
    dry = ROOT / "results/test_dryrun/test_results.json"
    if not dry.exists():
        pytest.skip("no dry-run results")
    R = json.loads(dry.read_text())
    R["runtime"] = {"error": "RuntimeError('x')", "traceback": ""}
    R["external"] = {"error": "RuntimeError('y')", "traceback": ""}
    rt.write_ledger(tmp_path, R, dry=False)
    txt = (tmp_path / "LEDGER.md").read_text()
    assert txt.count("FAILED") == 2


def test_jumpsem_timeout_is_recorded(monkeypatch, tmp_path):
    import subprocess

    import numpy as np
    import pandas as pd
    import importlib
    J = importlib.import_module("kinact.methods.jumpsem")

    def hang(*a, **k):
        raise subprocess.TimeoutExpired(cmd="Rscript", timeout=7200)
    monkeypatch.setattr(J.subprocess, "run", hang)
    sites = [f"P{i:05d}_S{i}" for i in range(10)]
    Xu = pd.DataFrame(np.random.default_rng(0).normal(size=(10, 4)), index=sites, columns=list("abcd"))
    lib = pd.DataFrame({"kinase": "K1", "site": sites})
    uid, act, extra, wall, err = J._run_unit(("u", Xu, lib, 3, str(tmp_path)))
    assert act is None and "timeout" in err
