"""Baseline kinase-activity methods (DESIGN §6). Common interface:

    score(X: sites x conditions DataFrame (log2FC, NaN = not measured), lib, min_substrates=3, return_n=False)
        -> conditions x kinases DataFrame (larger = more active; unscored = NaN, never 0)

See the module docstrings for formulas and reference checks (scripts/check_baseline_refs.py).
"""
from .enrich import gsea, viper
from .jumpsem import jumpsem, run_jumpsem, study_units
from .linear import ksea, mlm, norm_wmean, ulm, wmean, zscore
from .motif import kl_enrich, motif_zscore, phosx_like
from .rokai import rokai

METHODS = {"zscore": zscore, "ksea": ksea, "wmean": wmean, "norm_wmean": norm_wmean, "ulm": ulm, "mlm": mlm,
           "viper": viper, "gsea": gsea, "rokai": rokai, "motif_zscore": motif_zscore, "phosx_like": phosx_like,
           "kl_enrich": kl_enrich, "jumpsem": jumpsem}

__all__ = list(METHODS) + ["METHODS", "run_jumpsem", "study_units"]
