"""Command line interface: `kinact score input.tsv -o out.tsv [--format long|wide] [--device cpu|cuda]
[--config config/core_v3.json]`.

Input (TSV/CSV, auto-detected by extension; '-' = stdin):
  wide  first column = site key `{UniProt}_{res}{pos}` (e.g. P31749_S473), other columns = conditions (log2FC or t/z);
  long  columns uniprot, residue, position, value (+ optional condition).
Output: long (condition, kinase, z, mean, sd, n_substrates [+ z_core, z_sig with core-v3]) or wide (kinase rows).
z = core-v1 z_vif (default config) or the core-v3 fused z (--config config/core_v3.json); unscored kinases are NaN.
"""
from __future__ import annotations

import argparse
import sys
import time
import warnings


def _read(path: str):
    import pandas as pd
    sep = "," if path.lower().endswith((".csv", ".csv.gz")) else "\t"
    src = sys.stdin if path == "-" else path
    df = pd.read_csv(src, sep=sep)
    low = [str(c).lower() for c in df.columns]
    if all(c in low for c in ("uniprot", "residue", "position", "value")):
        return df
    df = df.set_index(df.columns[0])
    df.index = df.index.astype(str)
    return df


def cmd_score(a) -> int:
    from .core import KinaseActivityModel
    t0 = time.perf_counter()
    model = KinaseActivityModel.from_config(a.config)
    data = _read(a.input)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        res = model.score(data, device=a.device, n_jobs=a.threads)
    for x in w:
        print(f"warning: {x.message}", file=sys.stderr)
    out = res.long(dropna=a.drop_unscored) if a.format == "long" else res.wide()
    dst = sys.stdout if a.output in (None, "-") else a.output
    try:
        out.to_csv(dst, sep="\t", index=(a.format == "wide"), na_rep="NaN", float_format="%.6g")
    except BrokenPipeError:                       # e.g. `kinact score x.tsv | head`
        sys.stderr.close()
        return 0
    i = res.info
    if not a.quiet:
        print(f"kinact score: {i['n_conditions']} condition(s), {i['n_input_sites']} input sites "
              f"({i['n_library_sites_matched']} in library), {int(res.z.notna().to_numpy().sum())} kinase scores, "
              f"{time.perf_counter() - t0:.2f} s", file=sys.stderr)
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="kinact_mvp", description="Kinase activity inference (core-v3).")
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("score", help="score kinase activities from site-level statistics")
    s.add_argument("input", help="TSV/CSV input ('-' = stdin)")
    s.add_argument("-o", "--output", default=None, help="output TSV (default stdout)")
    s.add_argument("--format", choices=["long", "wide"], default="long")
    s.add_argument("--device", default="cpu", help="cpu | cuda | cuda:N")
    s.add_argument("--threads", type=int, default=None, help="CPU threads (default: BLAS thread limit)")
    s.add_argument("--config", default=None, help="model config (default: model/config.json = core-v3)")
    s.add_argument("--drop-unscored", action="store_true", help="long format: omit NaN (unscored) rows")
    s.add_argument("-q", "--quiet", action="store_true")
    s.set_defaults(func=cmd_score)
    a = ap.parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
