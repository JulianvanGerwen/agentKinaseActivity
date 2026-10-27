# Reference agreement checks (scripts/check_baseline_refs.py)

decoupler 2.3.1 run per condition on measured sites (tmin=3), library `psp_all`, 60 DEV conditions (seed 0). Spearman over kinases scored by both, per condition (min / median) and pooled over all pairs.

| method | conditions | decoupler failed | rho min | rho median | rho pooled | max abs diff | ours s | decoupler s |
|---|---|---|---|---|---|---|---|---|
| zscore | 60 | 0 | 1.0000 | 1.0000 | 1.0000 | 3.55e-15 | 0.13 | 2.2 |
| ksea | 60 | 0 | 1.0000 | 1.0000 | 1.0000 | 1.78e-15 | 0.12 | 2.1 |
| wmean | 60 | 0 | 1.0000 | 1.0000 | 1.0000 | 0 | 0.16 | 4.0 |
| norm_wmean | 60 | 0 | 0.9959 | 0.9989 | 0.9995 | 1.64 | 0.13 | 25.6 |
| ulm | 60 | 0 | 1.0000 | 1.0000 | 1.0000 | 1.47e-12 | 0.13 | 2.2 |
| mlm | 41 | 19 | 0.9371 | 1.0000 | 0.9986 | 4.94e+05 | 0.68 | 8.6 |
| viper | 60 | 0 | 1.0000 | 1.0000 | 1.0000 | 1.51e-14 | 0.64 | 2.4 |
| gsea | 60 | 0 | 0.9967 | 0.9989 | 0.9996 | 0.106 | 0.84 | 14.8 |

KSEAapp 2.0 `KSEA.Scores` (R) vs `kinact.methods.ksea(min_substrates=1)`, library `psp_human_only` (as KSData, PSP source), 12 conditions: Spearman min 1.000000, median 1.000000, max |Δz| 3.2e-14; kinase-count difference 0.

Notes: norm_wmean is computed analytically (infinite-permutation limit of decoupler's 1000-permutation z) and gsea uses its own permutation draws, so these differ from decoupler by Monte-Carlo noise only. mlm: decoupler raises on rank-deficient designs (e.g. two kinases with identical measured substrate sets); those conditions are counted as 'decoupler failed' and our pseudo-inverse solution is used there; among the full-rank conditions, rho < 0.99 occurs only on near-singular designs (1 of 41 here) where decoupler's explicit inverse returns numerically meaningless |t| (up to ~5e5). viper calls decoupler's own `_func_viper` per condition (no re-implementation).
