"""Shared helpers for scripts/build_gt.py and scripts/make_split.py (ground truth + split).

Everything biological lives in rule files under results/ground_truth/*.tsv; this module only parses.
"""
from __future__ import annotations

import hashlib
import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from prep_data import normalise_symbols  # noqa: E402

DATA = ROOT / "data"
CACHE = ROOT / "cache"
EXT = ROOT / "external"
GT = ROOT / "results" / "ground_truth"
SEALED = CACHE / "sealed"
HIJAZI_PMID = "31959955"
DECRYPTM_PMID = "36926954"
WILKES_PMID = "26060313"
CONF = ["low", "medium", "high"]


def norm(s) -> str:
    s = str(s).replace("Ã¶", "o").replace("ö", "o").replace("µ", "u").replace("μ", "u")
    return re.sub(r"[^a-z0-9]", "", s.lower())


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def read_rules(name) -> pd.DataFrame:
    return pd.read_csv(GT / name, sep="\t", dtype=str, comment="#", keep_default_na=False)


def write_tsv(df: pd.DataFrame, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, lineterminator="\n", float_format="%.6g")


# ----------------------------------------------------------------------------------------- kinase universe
def load_universe(hg: pd.DataFrame) -> pd.DataFrame:
    """U = kinaseSubMain kinase_list include==True (503), symbols normalised to current HGNC."""
    kl = pd.read_csv(EXT / "kinasesubmain" / "kinase_list.tsv", sep="\t", dtype=str)
    kl = kl[kl["include"] == "True"].copy()
    m = normalise_symbols(kl["gene"], hg).set_index("name")
    kl["kinase"] = [m.loc[g, "hgnc_symbol"] if m.loc[g, "match_kind"] in ("symbol", "prev_symbol", "alias_symbol")
                    else g for g in kl["gene"]]
    assert kl["kinase"].is_unique, "U symbols not unique after HGNC normalisation"
    u = pd.DataFrame({"kinase": kl["kinase"], "uniprot": kl["uniprot_acc"],
                      "manning_group": kl["manning_group"].fillna(""), "manning_family": kl["manning_family"].fillna(""),
                      "is_pseudokinase": kl["is_pseudokinase"].map({"True": True, "False": False}),
                      "ksm_gene": kl["gene"]})
    return u.sort_values("kinase").reset_index(drop=True)


class SymbolMapper:
    """Map any kinase-like token to a symbol in U (or None). Uses U symbols, KSM gene names, HGNC prev/alias."""

    def __init__(self, U: pd.DataFrame, hg: pd.DataFrame):
        self.U = set(U["kinase"])
        self.direct = {k.upper(): k for k in U["kinase"]}
        self.direct.update({g.upper(): k for g, k in zip(U["ksm_gene"], U["kinase"])})
        self.direct.update({a.upper(): k for a, k in zip(U["uniprot"], U["kinase"])})
        # HGNC lookup (same priority rule as prep_data.normalise_symbols: symbol > previous > alias; ambiguous -> None)
        pri = {"symbol": 0, "prev_symbol": 1, "alias_symbol": 2}
        h = hg.assign(p=hg["kind"].map(pri)).dropna(subset=["p"])
        h = h[h["p"] == h.groupby("query_upper")["p"].transform("min")]
        nsym = h.groupby("query_upper")["symbol"].nunique()
        h = h[h["query_upper"].map(nsym) == 1].drop_duplicates("query_upper")
        self.hgnc = dict(zip(h["query_upper"], h["symbol"]))
        self.cache = {}

    def __call__(self, tok):
        t = str(tok).strip()
        if not t:
            return None
        if t in self.cache:
            return self.cache[t]
        r = self.direct.get(t.upper())
        if r is None:
            s = self.hgnc.get(t.upper())
            r = s if s in self.U else None
        self.cache[t] = r
        return r


# ----------------------------------------------------------------------------------------- drug aliases
class DrugMatcher:
    """Find drugs (agents) in free text using results/ground_truth/drug_aliases.tsv.

    Tokens: text split on whitespace/_+,;()/: ; single tokens and adjacent pairs/triples are normalised (lower-case,
    alphanumerics only) and looked up. Class tokens (is_class=1, e.g. 'MEKi') count only when no named drug is found in
    the same text. Aliases shorter than 3 characters are never matched.
    """

    def __init__(self, aliases: pd.DataFrame | None = None):
        a = read_rules("drug_aliases.tsv") if aliases is None else aliases
        self.tab = a
        self.named, self.cls = {}, {}
        for r in a.itertuples():
            k = norm(r.alias)
            pm = str(r.pmid or "")
            if len(k) < (2 if pm else 3):
                continue
            d = self.cls if str(r.is_class) == "1" else self.named
            key = (k, pm)
            if key in d and d[key] != r.agent:
                raise ValueError(f"alias {r.alias} maps to two agents: {d[key]}, {r.agent}")
            d[key] = r.agent

    @staticmethod
    def tokens(text):
        toks = [x for x in re.split(r"[\s_+,;()/:]+", str(text)) if x]
        out = list(toks)
        out += ["".join(p) for p in zip(toks, toks[1:])]
        out += ["".join(p) for p in zip(toks, toks[1:], toks[2:])]
        return out

    def _get(self, d, k, pmid):
        return d.get((k, str(pmid))) or d.get((k, ""))

    def find(self, text, pmid="") -> list[str]:
        named, cls = [], []
        for t in self.tokens(text):
            k = norm(t)
            # also try with a trailing dose/number removed (e.g. 'caffeine0.01 uM' -> 'caffeine')
            cands = [k, re.sub(r"\d+(\d|um|nm|mm|pm|h|min)*$", "", k)]
            for kk in cands:
                if len(kk) < 2:
                    continue
                a = self._get(self.named, kk, pmid)
                if a:
                    if a not in named:
                        named.append(a)
                    break
                a = self._get(self.cls, kk, pmid)
                if a:
                    if a not in cls:
                        cls.append(a)
                    break
        return sorted(named) if named else sorted(cls)


# ----------------------------------------------------------------------------------------- time / dose parsing
_TIME = re.compile(r"(?<![A-Za-z0-9.])(\d+(?:\.\d+)?)\s*(sec|s|min|m|h|hr|d)(?![A-Za-z])")
_TUNIT = {"sec": 1 / 60, "s": 1 / 60, "min": 1, "m": 1, "h": 60, "hr": 60, "d": 1440}


def parse_minutes(condition: str, timepoint: str):
    """Treatment time in minutes from the Condition text (preferred: the Timepoint field has parse errors, e.g. 'PMA
    0.5h' -> '5h'), else from Timepoint. Returns (minutes or nan, how)."""
    for txt, how in ((condition, "condition_text"), (timepoint, "timepoint_field")):
        hits = _TIME.findall(str(txt).replace("_", " "))
        if hits:
            mins = [float(v) * _TUNIT[u] for v, u in hits]
            return max(mins), how
    return float("nan"), "unknown"


_DOSE = re.compile(r"(\d+(?:\.\d+)?(?:e[+-]?\d+)?)\s*(pM|nM|uM|µM|μM|&mu;M|mM)(?![a-zA-Z])", re.I)
_UNIT_NM = {"pm": 1e-3, "nm": 1.0, "um": 1e3, "µm": 1e3, "μm": 1e3, "&mu;m": 1e3, "mm": 1e6}


def parse_doses_nM(text: str) -> list[float]:
    return [float(v) * _UNIT_NM[u.lower()] for v, u in _DOSE.findall(str(text))]


# ----------------------------------------------------------------------------------------- union-find
class UnionFind:
    def __init__(self, items=()):
        self.p = {}
        for x in items:
            self.p[x] = x

    def find(self, x):
        self.p.setdefault(x, x)
        root = x
        while self.p[root] != root:
            root = self.p[root]
        while self.p[x] != root:
            self.p[x], x = root, self.p[x]
        return root

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            if rb < ra:
                ra, rb = rb, ra
            self.p[rb] = ra


def is_inf(x) -> bool:
    return isinstance(x, float) and math.isinf(x)
