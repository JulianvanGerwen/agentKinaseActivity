#!/usr/bin/env python
"""Write results/ground_truth/GT_FREEZE.json (sha256 of GT, split, rule files, scripts, sealed TEST files and the DEV
cache) - only after tests/test_gt.py passes. `sealed_gt_hash` uses the same algorithm as kinact.eval.sealed_gt_hash
(sha256 over sorted '<name>\\t<sha256>' lines of cache/sealed/*.tsv) so the PI's FREEZE.json can carry it.
Usage: python scripts/freeze_gt.py"""
import datetime
import hashlib
import json
import subprocess
import sys

from gtlib import CACHE, GT, ROOT, SEALED, sha256


def main():
    r = subprocess.run([sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider", str(ROOT / "tests" / "test_gt.py")],
                       cwd=ROOT, capture_output=True, text=True)
    print(r.stdout[-500:])
    if r.returncode != 0:
        sys.exit("tests/test_gt.py failed - GT_FREEZE.json NOT written")
    out = {"created": datetime.date.today().isoformat(), "tests": "tests/test_gt.py passed",
           "split_seed": json.loads((GT / "split_log.json").read_text())["chosen_seed"]}
    gt_files = ["gt_pairs.tsv", "gt_conditions.tsv", "kinase_universe.tsv", "harmony_bins.tsv", "harmony_bins_group.tsv",
                "condition_groups.tsv", "test_condition_hashes.txt", "split_candidates.tsv", "split_log.json",
                "condition_doses.tsv", "condition_qc.tsv", "drug_profile_sources.tsv", "paralog_pairs.tsv", "build_gt_log.json"]
    rule_files = sorted(p.name for p in GT.glob("*.tsv") if p.name not in gt_files)
    out["gt"] = {f: sha256(GT / f) for f in gt_files}
    out["rules"] = {f: sha256(GT / f) for f in rule_files}
    out["sealed"] = {p.name: sha256(p) for p in sorted(SEALED.glob("*.tsv"))}
    out["sealed_full"] = {p.name: sha256(p) for p in sorted((SEALED / "full").glob("*.tsv"))}
    lines = [f"{p.name}\t{hashlib.sha256(p.read_bytes()).hexdigest()}" for p in sorted(SEALED.glob("*.tsv"))]
    out["sealed_gt_hash"] = hashlib.sha256("\n".join(lines).encode()).hexdigest()
    out["scripts"] = {f: sha256(ROOT / "scripts" / f) for f in
                      ("gtlib.py", "build_gt.py", "make_split.py", "make_dev_cache.py", "gt_report.py", "freeze_gt.py")}
    out["tests"] = {"test_gt.py": sha256(ROOT / "tests" / "test_gt.py")}
    out["dev_cache"] = {f: sha256(CACHE / "dev" / f) for f in ("quant_dev.npy", "conditions_dev.txt", "sites_dev.txt")}
    out["inputs"] = {"cache/conditions.parquet": sha256(CACHE / "conditions.parquet"),
                     "cache/quant_wide_f32.npy": sha256(CACHE / "quant_wide_f32.npy"),
                     "cache/psp_ks_human.parquet": sha256(CACHE / "psp_ks_human.parquet")}
    (GT / "GT_FREEZE.json").write_text(json.dumps(out, indent=1) + "\n")
    print("wrote", GT / "GT_FREEZE.json")


if __name__ == "__main__":
    main()
