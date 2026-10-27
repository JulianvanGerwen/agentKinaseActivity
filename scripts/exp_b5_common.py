"""B5 shared setup: DEV GT conditions, B2-standardised Y/D, self-site tables, B2 reference design."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts"))
import exp_b2_grid as E2  # noqa: E402
from kinact.explore import bayes as B  # noqa: E402
from kinact.explore import selfsites as SS  # noqa: E402
from kinact.library import load_library  # noqa: E402

FOLDS = E2.FOLDS
OUT, CACHE = ROOT / "results/b5", ROOT / "cache/b5"
CORE = dict(tau=0.125, rho_f=3.0, rho_g=0.25, alpha=1.0, intercept="median", borrow=True)
SETS = ["psp", "aloop", "psp+aloop"]


class Ctx:
    def __init__(self):
        self.P, self.C, self.Udf, self.ev, self.meta, X = E2.setup()
        self.conds = list(self.ev.cond_index)
        self.kin = list(self.Udf["kinase"])
        self.fold_of = self.meta["split"]
        self.site_index = X.index
        self.Y, self.D = E2.standardise_by_fold(X, self.conds, self.meta, {}, noise=False)
        self.X = X
        self.tables = {w: SS.self_site_table(self.Udf, w) for w in SETS}
        self.design = B.build_design(X.index, self.kin, load_library(E2.LIB_CUR),
                                     load_library(E2.LIB_PRED[0], k=E2.LIB_PRED[1]))
        self.K = B.kinase_kernel(self.Udf, self.kin, CORE["rho_f"], CORE["rho_g"])

    def signal(self, which, standardise=True):
        W = SS.self_matrix(self.tables[which], self.site_index, self.kin)
        s, n = SS.self_signal(self.Y, self.D, W, standardise)
        return pd.DataFrame(s, index=self.conds, columns=self.kin), pd.DataFrame(n, index=self.conds, columns=self.kin)

    def core(self, design=None, alpha=None):
        design = self.design if design is None else design
        G, h, n = B.suff_stats(design, CORE["alpha"] if alpha is None else alpha, self.Y, self.D)
        mu, sd = B.posterior(G, h, self.K, CORE["tau"], CORE["intercept"])
        rep = B.report_mask(n.cpu().numpy(), self.K, 1, CORE["borrow"])
        mk = lambda a: pd.DataFrame(np.where(rep, a, np.nan), index=self.conds, columns=self.kin)
        return {"z": mk(mu / sd), "mean": mk(mu)}
