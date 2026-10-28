"""Branch B5 — kinase self-sites (regulatory phosphosites ON the kinase protein) as an auxiliary activity signal.

Self-site table (label-free; no GT information):
  * 'psp'  : PSP Regulatory_sites rows whose (human-mapped) protein is a kinase in U and whose ON_FUNCTION
             contains '(enzymatic) activity, induced' (sign +1) or '(enzymatic) activity, inhibited' (sign -1).
             Orthologue rows (mouse/rat/...) enter via the PSP SITE_GRP_ID mapping already in the cache
             (`site_key_psp`); sites annotated both induced and inhibited are dropped.
  * 'aloop': S/T/Y residues of the activation loop (between the DFG motif and the APE motif of the kinase domain,
             kinaseSubMain phosphoformer domain table) with sign +1, unless PSP marks the site inhibitory.
Signal: s_k(c) = mean over measured self-sites i of k of sign_i * y_ic, with y = (x - med_c)/(1.4826 MAD_c) (the
B2 per-condition robust standardisation), then standardised across kinases within the condition (divided by the
condition's RMS of s over kinases with a measured self-site).

Helpers to EXCLUDE self-sites from an input vector (needed for GT-A, DESIGN §3): `self_site_keys`, `mask_self_sites`.
"""
from __future__ import annotations

import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[3]
REG_PQ = ROOT / "cache/psp_regulatory.parquet"
DOMAINS_PQ = Path(os.environ.get("KINASESUB_DIR", "kinaseSubMain"), "cache/phosformer/"
                  "kinase_domains.parquet")
ACT = re.compile(r"(?:^|;\s*)(?:enzymatic )?activity, induced")
INH = re.compile(r"(?:^|;\s*)(?:enzymatic )?activity, inhibited")


# ----------------------------------------------------------------------------- tables
def psp_self_sites(universe: pd.DataFrame, reg: pd.DataFrame | None = None) -> pd.DataFrame:
    """(kinase, site, sign, source='psp', organism) for activating / inhibitory PSP regulatory sites on U kinases."""
    reg = pd.read_parquet(REG_PQ) if reg is None else reg
    u2k = universe.dropna(subset=["uniprot"]).set_index("uniprot")["kinase"]
    r = reg[reg["site_key_psp"].notna() & reg["hs_acc"].isin(u2k.index) & reg["ON_FUNCTION"].notna()].copy()
    r["act"] = r["ON_FUNCTION"].str.contains(ACT)
    r["inh"] = r["ON_FUNCTION"].str.contains(INH)
    r = r[r["act"] | r["inh"]]
    g = r.groupby("site_key_psp").agg(act=("act", "any"), inh=("inh", "any"), hs_acc=("hs_acc", "first"),
                                      organisms=("ORGANISM", lambda s: ",".join(sorted(set(s)))))
    g = g[g["act"] != g["inh"]]                       # drop conflicting annotations
    out = pd.DataFrame({"kinase": g["hs_acc"].map(u2k).to_numpy(), "site": g.index.to_numpy(),
                        "sign": np.where(g["act"], 1, -1), "source": "psp", "organisms": g["organisms"].to_numpy()})
    return out.sort_values(["kinase", "site"]).reset_index(drop=True)


def activation_loop_sites(universe: pd.DataFrame, domains: pd.DataFrame | None = None) -> pd.DataFrame:
    """S/T/Y between DFG (D[FLYW]G) and the next APE-like motif ([ASPG][PLI]E, 12-45 aa later) in the kinase domain.
    Position = dom_start + 0-based index (checked: MAPK1 T185/Y187, AKT1 T308)."""
    if domains is None:
        if not DOMAINS_PQ.exists():
            return pd.DataFrame(columns=["kinase", "site", "sign", "source"])
        domains = pd.read_parquet(DOMAINS_PQ)
    u2k = universe.dropna(subset=["uniprot"]).set_index("uniprot")["kinase"]
    rows = []
    for acc, seq, st in zip(domains["uniprot_acc"], domains["domain_seq"], domains["dom_start"]):
        if acc not in u2k.index or not isinstance(seq, str):
            continue
        for m in re.finditer(r"D[FLYW]G", seq):
            a = re.search(r"[ASPG][PLI]E", seq[m.end() + 12: m.end() + 45])
            if a is None:
                continue
            lo, hi = m.end(), m.end() + 12 + a.start()
            for i in range(lo, hi):
                if seq[i] in "STY":
                    rows.append((u2k[acc], f"{acc}_{seq[i]}{int(st) + i}", 1, "aloop"))
            break
    return pd.DataFrame(rows, columns=["kinase", "site", "sign", "source"]).drop_duplicates(["kinase", "site"])


def self_site_table(universe: pd.DataFrame, which: str = "psp") -> pd.DataFrame:
    """which: 'psp' | 'aloop' | 'psp+aloop' (PSP sign wins where both exist)."""
    parts = []
    psp = psp_self_sites(universe)
    if "psp" in which:
        parts.append(psp)
    if "aloop" in which:
        al = activation_loop_sites(universe)
        al = al[~al["site"].isin(psp.loc[psp["sign"] < 0, "site"])]
        if "psp" in which:
            al = al[~al["site"].isin(psp["site"])]
        parts.append(al)
    return pd.concat(parts, ignore_index=True).drop_duplicates(["kinase", "site"])


# ----------------------------------------------------------------------------- exclusion helpers (GT-A)
def self_site_keys(universe: pd.DataFrame, which: str = "psp+aloop", include_all_regulatory: bool = True) -> set:
    """Site keys to remove from method input when computing GT-A. include_all_regulatory=True also removes every PSP
    regulatory site on a U kinase that mentions 'activity' (incl. 'altered'), the conservative choice."""
    keys = set(self_site_table(universe, which)["site"])
    if include_all_regulatory:
        reg = pd.read_parquet(REG_PQ)
        u = set(universe["uniprot"].dropna())
        m = reg["hs_acc"].isin(u) & reg["ON_FUNCTION"].fillna("").str.contains("activity") & reg["site_key_psp"].notna()
        keys |= set(reg.loc[m, "site_key_psp"])
    return keys


def mask_self_sites(X, keys: set):
    """Return a copy of X (DataFrame sites x conditions, or Series) with self-sites set to NaN (= not measured)."""
    X = X.copy()
    hit = X.index.isin(list(keys))
    if isinstance(X, pd.Series):
        X[hit] = np.nan
    else:
        X.loc[hit] = np.nan
    return X


# ----------------------------------------------------------------------------- signal
def self_matrix(table: pd.DataFrame, site_index, kinases) -> sp.csr_matrix:
    """sites x kinases signed incidence matrix (sign of each self-site)."""
    si = pd.Index(site_index).get_indexer(table["site"])
    ki = pd.Index(kinases).get_indexer(table["kinase"])
    ok = (si >= 0) & (ki >= 0)
    return sp.csr_matrix((table["sign"].to_numpy(float)[ok], (si[ok], ki[ok])), shape=(len(site_index), len(kinases)))


def self_signal(Y: np.ndarray, D: np.ndarray, Wself: sp.csr_matrix, standardise: bool = True):
    """Y, D: sites x conditions (standardised values, D>0 = measured). Returns s (cond x kin, NaN = no measured
    self-site) and n (cond x kin measured self-site counts)."""
    M = (D > 0).astype(float)
    num = (Wself.T @ (M * Y)).T
    n = (abs(Wself).T @ M).T
    with np.errstate(all="ignore"):
        s = np.where(n > 0, num / n, np.nan)
        if standardise:
            rms = np.sqrt(np.nanmean(s * s, axis=1, keepdims=True))
            s = s / np.where(rms > 0, rms, 1.0)
    return s, n


def per_condition_scale(S: np.ndarray) -> np.ndarray:
    """Divide each condition's scores by their RMS over scored kinases (late-combination scale matching)."""
    with np.errstate(all="ignore"):
        rms = np.sqrt(np.nanmean(S * S, axis=1, keepdims=True))
    return S / np.where(rms > 0, rms, 1.0)


def late_combine(core: np.ndarray, s: np.ndarray, w: float, mode: str) -> np.ndarray:
    """core, s: cond x kin (NaN = unscored). Modes (all use only the input data):
      add   : core + w s
      agree : core + w s 1[sign s = sign core]          (only reinforce)
      mag   : core + w |s| sign(core)                   (sign-agnostic 'this kinase changed')
      up    : core + w max(s, 0)                        (only increases of activating self-sites)
    Kinases unscored by the core but with a self-site get w*s ('add','up') so self-sites can add coverage."""
    c = per_condition_scale(core)
    s0 = np.nan_to_num(s)
    if mode == "add":
        t = s0
    elif mode == "agree":
        t = s0 * (np.sign(s0) == np.sign(np.nan_to_num(c)))
    elif mode == "mag":
        t = np.abs(s0) * np.sign(np.nan_to_num(c))
    elif mode == "up":
        t = np.maximum(s0, 0)
    else:
        raise ValueError(mode)
    out = np.nan_to_num(c) + w * t
    scored = ~np.isnan(c) | (~np.isnan(s) & (mode in ("add", "up")) & (w > 0))
    return np.where(scored, out, np.nan)


def design_with_self(design, table: pd.DataFrame, site_index, lam: float):
    """B2 Design with self-sites added as extra observations of a_k: y_i = lam * sign_i * a_k + eps_i (same noise as
    other sites). Added to the curated part (so W(alpha) keeps them at full weight)."""
    from .bayes import Design
    site_index = pd.Index(site_index)
    Ws = self_matrix(table, site_index, design.kinases)            # all sites x kin
    srows = np.flatnonzero(np.asarray(abs(Ws).sum(1)).ravel() > 0)
    rows = np.union1d(design.lib_rows, srows).astype(int)

    def remap(Wm, old_rows):
        R = sp.csr_matrix((np.ones(len(old_rows)), (np.searchsorted(rows, old_rows), np.arange(len(old_rows)))),
                          shape=(len(rows), len(old_rows)))
        return (R @ Wm).tocsr()
    Wc = remap(design.Wc, design.lib_rows) + lam * Ws[rows]
    Wp = remap(design.Wp, design.lib_rows)
    Wc.eliminate_zeros()
    return Design(design.kinases, rows, Wc.tocsr(), Wp)
