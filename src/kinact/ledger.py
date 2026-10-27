"""Append-only DEV ledger results/dev_ledger.tsv (DESIGN §7).

Columns: timestamp, git_commit, method, config_hash, branch, split, folds,
n_configs_tried (distinct config hashes so far in this branch, incl. this row), metrics_json.
Rows are only ever appended (O_APPEND + flock); existing bytes are never rewritten.
"""
from __future__ import annotations

import csv
import fcntl
import hashlib
import io
import json
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path

from .eval import DEFAULT_ROOT, check_split_access, to_jsonable

LEDGER_PATH = DEFAULT_ROOT / "results/dev_ledger.tsv"
COLUMNS = ["timestamp", "git_commit", "method", "config_hash", "branch", "split", "folds",
           "n_configs_tried", "metrics_json"]


def config_hash(config) -> str:
    return hashlib.sha256(json.dumps(config, sort_keys=True, default=str).encode()).hexdigest()[:16]


def git_commit(root: Path = DEFAULT_ROOT) -> str:
    try:
        out = subprocess.run(["git", "-C", str(root), "rev-parse", "--short=12", "HEAD"],
                             capture_output=True, text=True, timeout=10)
        if out.returncode != 0:
            return "nogit"
        dirty = subprocess.run(["git", "-C", str(root), "status", "--porcelain", "-uno"],
                               capture_output=True, text=True, timeout=10).stdout.strip()
        return out.stdout.strip() + ("+dirty" if dirty else "")
    except Exception:
        return "nogit"


def read_ledger(path: str | os.PathLike = LEDGER_PATH) -> list[dict]:
    path = Path(path)
    if not path.exists():
        return []
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def append(method: str, config, split, metrics: dict, branch: str = "main", folds="",
           path: str | os.PathLike = LEDGER_PATH, commit: str | None = None) -> dict:
    """Append one evaluation row; returns the row written. Refuses split 'test'."""
    splits = [split] if isinstance(split, str) else list(split)
    if "test" in splits:
        check_split_access(splits)
        raise PermissionError("the DEV ledger never records TEST evaluations")
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    ch = config_hash(config)
    fd = os.open(path, os.O_WRONLY | os.O_APPEND | os.O_CREAT, 0o644)
    with os.fdopen(fd, "a", newline="") as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        try:
            prev = read_ledger(path)
            if prev and list(prev[0].keys()) != COLUMNS:
                raise ValueError(f"{path}: unexpected header")
            seen = {r["config_hash"] for r in prev if r["branch"] == branch} | {ch}
            row = {"timestamp": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                   "git_commit": commit or git_commit(), "method": method, "config_hash": ch,
                   "branch": branch, "split": ",".join(map(str, splits)),
                   "folds": ",".join(map(str, folds)) if not isinstance(folds, str) else folds,
                   "n_configs_tried": len(seen),
                   "metrics_json": json.dumps(to_jsonable(metrics), sort_keys=True, default=str)}
            buf = io.StringIO()
            w = csv.DictWriter(buf, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
            if os.fstat(f.fileno()).st_size == 0:
                w.writeheader()
            w.writerow(row)
            f.write(buf.getvalue())
            f.flush()
        finally:
            fcntl.flock(f, fcntl.LOCK_UN)
    return row
