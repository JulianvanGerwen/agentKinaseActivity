"""Shared setup for branch B7 (perturbation response signatures; DESIGN §12 C10). DEV only (cache/dev +
results/ground_truth DEV files). Design / HP grid pre-declared in kinact.explore.signatures (module docstring).

Nested CV with fold-dependent signatures: for outer fold f and HP config h,
  inner_f(h) = M1-sel over folds g != f, where fold g is scored with signatures fitted on folds not in {f, g};
  the pick h*_f = argmax inner_f; fold f is then scored with signatures fitted on the 4 folds != f.
Core HPs (tau, rho_f, rho_g) per outer fold = core-v1 z_vif nested picks (results/audit_core_v1/dev_audit.json), used
for both inner and outer scoring of that outer fold, so alpha_sig = 0 reproduces core-v1.
"""
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))
from exp_b21_common import Ctx, FOLDS  # noqa: E402
from kinact.explore import bayes2 as B2  # noqa: E402
from kinact.explore import signatures as SG  # noqa: E402
from kinact.explore.bayes import report_mask  # noqa: E402

OUT = ROOT / "results/b7"; OUT.mkdir(parents=True, exist_ok=True)
CACHE = ROOT / "cache/b7"; CACHE.mkdir(parents=True, exist_ok=True)
BRANCH = "B7"
GT_PREFIX = "1e6f791c"
ALPHAS = [0.25, 0.5, 1.0, 2.0]
NTOPS = [20, 50, 100, 200]
KAPPAS = [0.0, 1.0, 4.0]
VARIANTS = {"clean+contrast": dict(nt_max=10, contrast=True),     # PRIMARY (pre-declared)
            "all+contrast": dict(nt_max=None, contrast=True),
            "clean-nocontrast": dict(nt_max=10, contrast=False),
            "all-nocontrast": dict(nt_max=None, contrast=False)}
PRIMARY = "clean+contrast"


def core_picks():
    d = json.loads((ROOT / "results/audit_core_v1/dev_audit.json").read_text())
    ch = d["results"]["core-v1 z_vif (tau,rho_f<=3,rho_g nested)"]["chosen"]
    return {f: (float(c["tau"]), float(c["rho_f"]), float(c["rho_g"])) for f, c in ch.items()}


class B7Ctx:
    def __init__(self):
        t0 = time.time()
        self.c = c = Ctx(ks={"r2np": [100]})
        assert c.gt_hash.startswith(GT_PREFIX), c.gt_hash
        self.gt_hash = c.gt_hash
        self.picks = core_picks()
        # standardised data for ALL DEV conditions (column-wise -> identical to Ctx.blocks standardisation)
        Y, M, _ = B2.standardise(c.X.to_numpy(np.float64), "wsd", None)
        Yn = np.where(M > 0, Y, np.nan).astype(np.float32)
        self.Yfit = Yn                                  # NaN = missing (signature fitting)
        self.Y = Y.astype(np.float32); self.M = M.astype(np.float32)
        self.allconds = pd.Index(c.X.columns)
        self.cond_group = c.meta["group_id"]
        # permuted GT conditions (calibration): same standardisation of c.Xp
        Yp, Mp, _ = B2.standardise(c.Xp.to_numpy(np.float64), "wsd", None)
        self.Yp, self.Mp = Yp.astype(np.float32), Mp.astype(np.float32)
        # base design on the full site index
        md = c.md
        Wl = md.W({"r2np@100": 1.0}).tocoo()
        import scipy.sparse as sp
        self.Wb = sp.csr_matrix((Wl.data, (md.lib_rows[Wl.row], Wl.col)), shape=(len(c.X.index), len(c.kin)))
        # base blocks
        self.base = {}
        for which in ("gt", "null", "perm"):
            bl = c.blocks(which, ["cur", "r2np@100"], "wsd", None)
            G, h, g0, n = B2.combine(bl, {"r2np@100": 1.0})
            self.base[which] = dict(G=G, h=h, g0=g0, n=n)
            c.drop_blocks() if which == "perm" else None
        self.K = {}
        self.fold_pos = {f: np.flatnonzero(c.cf == f) for f in FOLDS}
        nf = pd.Series(c.nulls).map(c.fold_of).to_numpy()
        self.null_fold = nf
        print(f"[b7ctx] {time.time()-t0:.0f}s sites {len(c.X.index)} conds {len(self.allconds)}", flush=True)

    def kernel(self, rf, rg):
        if (rf, rg) not in self.K:
            self.K[(rf, rg)] = B2.kernel2(self.c.Udf, self.c.kin, rf, rg)
        return self.K[(rf, rg)]

    # ------------------------------------------------------------------ signatures
    def fit(self, train_folds, eval_conds, variant, kappas=KAPPAS):
        v = VARIANTS[variant]
        c = self.c
        tr = SG.training_conditions(c.C, c.meta, train_folds, eval_conds)
        pos = SG.positives(c.P, c.C, tr, v["nt_max"])
        return SG.fit_signatures(self.Yfit, self.allconds, tr, pos, self.cond_group, eval_conds=eval_conds,
                                 kappa=list(kappas), n_max=max(NTOPS), contrast=v["contrast"])

    # ------------------------------------------------------------------ scoring
    def _sub(self, which, idx):
        b = self.base[which]
        return dict(G=b["G"][idx], h=b["h"][idx], g0=b["g0"][idx], n=b["n"][idx])

    def _YM(self, which, idx):
        if which == "gt":
            ci = self.allconds.get_indexer(np.asarray(self.c.conds)[idx]); return self.Y[:, ci], self.M[:, ci]
        if which == "null":
            ci = self.allconds.get_indexer(np.asarray(self.c.nulls)[idx]); return self.Y[:, ci], self.M[:, ci]
        return self.Yp[:, idx], self.Mp[:, idx]

    def score(self, idx_gt, Ws_list, core, alphas, which="gt", null_idx=None, return_n=False):
        """Score GT conditions idx_gt (positions in c.conds) — or permuted ('perm') — with signature design(s).
        Ws_list: list of (label, Ws csr or None). For each label and alpha: z_vif (len(idx) x p) with rho_e fitted on
        the DEV nulls under the same design. Returns dict (label, alpha) -> (Z, rho_e[, n])."""
        tau, rf, rg = core
        K = self.kernel(rf, rg)
        null_idx = np.arange(len(self.c.nulls)) if null_idx is None else null_idx
        bE, bN = self._sub(which, idx_gt), self._sub("null", null_idx)
        YE, ME = self._YM(which, idx_gt); YN, MN = self._YM("null", null_idx)
        out = {}
        for lab, Ws in Ws_list:
            if Ws is None:
                sE = sN = None
            else:
                sE = SG.sig_blocks(self.Wb, Ws, YE, ME); sN = SG.sig_blocks(self.Wb, Ws, YN, MN)
            for a in (alphas if Ws is not None else [0.0]):
                res = []
                for b, s in ((bN, sN), (bE, sE)):
                    G, h, g0, n = SG.combine_sig(b, s, a) if s is not None else (b["G"], b["h"], b["g0"], b["n"])
                    post = B2.posterior2(G, h, g0, K, tau)
                    rep = report_mask(n, K, 1, True)
                    res.append((post, rep, n))
                (pN, rN, _), (pE, rE, nE) = res
                rho, _ = B2.fit_rho_null(B2.output_stat(pN, "z_freq"), B2.m_eff(pN["q"], pN["suma"]), rN)
                Z = np.where(rE, B2.output_stat(pE, "z_vif", rho_e=rho), np.nan).astype(np.float32)
                out[(lab, a)] = (Z, rho, nE) if return_n else (Z, rho)
            if sE is not None:
                del sE, sN
        return out


def frame(ctx, parts: dict):
    """parts: fold -> array (fold rows). Full conds x kinases DataFrame (NaN where missing)."""
    A = np.full((len(ctx.c.conds), len(ctx.c.kin)), np.nan, np.float32)
    for f, arr in parts.items():
        A[ctx.fold_pos[f]] = arr
    return pd.DataFrame(A, index=ctx.c.conds, columns=ctx.c.kin)


def cfg_key(variant, alpha, n_top, kappa):
    return (variant, float(alpha), int(n_top), float(kappa)) if alpha else ("core", 0.0, 0, 0.0)


def config_dict(variant, alpha, n_top, kappa):
    v = VARIANTS.get(variant, {})
    return dict(model="bayes2+sig", base="core-v1 z_vif (per-fold core picks)", variant=variant,
                nt_max=v.get("nt_max"), contrast=v.get("contrast"), alpha_sig=alpha, n_top=n_top, kappa=kappa,
                n_min_groups=SG.N_MIN_GROUPS, clip=SG.CLIP, min_site_groups=SG.MIN_SITE_GROUPS)


def sha(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, default=str).encode()).hexdigest()[:16]
