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
