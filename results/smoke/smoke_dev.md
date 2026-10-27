# Smoke run on DEV (folds 0-4 pooled; no devconfirm/test)

Library `psp_all`, z-score = mean(x_s)*sqrt(n)/sd(x_c), min 3 substrates. 372 eligible GT conditions, 89 M1-sel groups. Controls: random (seed 0); prior-only and shift-prior use out-of-fold frequencies (other folds only); oracle. Runtime 2 s. Script: `scripts/smoke_dev.py`.

| method | M1-sel | M1-all | M1-sel-scored | pair coverage | harmony |
|---|---|---|---|---|---|
| zscore | 0.515 | 0.512 | 0.633 | 0.578 | 0.517 |
| random | 0.505 | 0.501 | 0.505 | 1.000 | 0.510 |
| prior_only | 0.687 | 0.543 | 0.687 | 1.000 | 0.726 |
| shift_prior | 0.546 | 0.514 | 0.546 | 1.000 | 0.561 |
| oracle | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |

## By source / sign (M1-sel)

| method | H | D | O | sign -1 | sign +1 |
|---|---|---|---|---|---|
| zscore | 0.518 | 0.560 | 0.507 | 0.513 | 0.536 |
| random | 0.542 | 0.493 | 0.493 | 0.507 | 0.521 |
| prior_only | 0.753 | 0.536 | 0.709 | 0.695 | 0.769 |
| shift_prior | 0.559 | 0.433 | 0.551 | 0.556 | 0.549 |
| oracle | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |

## Harmony per bin (M1-sel) and coverage-significant kinases

| method | 0 | 1-4 | 5-19 | >=20 | n sig (q<=0.10) |
|---|---|---|---|---|---|
| zscore | 0.501 | 0.508 | 0.531 | 0.527 | 24 |
| random | 0.544 | 0.505 | 0.489 | 0.502 | 0 |
| prior_only | 0.625 | 0.705 | 0.735 | 0.840 | 53 |
| shift_prior | 0.510 | 0.534 | 0.608 | 0.592 | 11 |
| oracle | 1.000 | 1.000 | 1.000 | 1.000 | 106 |

## Paired bootstrap (groups, B=2000, seed 0), M1-sel

- zscore_vs_prior_only: Δ=-0.172, SE=0.024, 95% basic CI [-0.179, -0.089], one-sided shifted-null p=1.0000 (C8; diagnostics: bootstrap bias -0.030, percentile CI [-0.255, -0.166], percentile p=1.0000)
- zscore_vs_random: Δ=+0.010, SE=0.013, 95% basic CI [-0.018, +0.035], one-sided shifted-null p=0.2169 (C8; diagnostics: bootstrap bias +0.002, percentile CI [-0.015, +0.039], percentile p=0.1749)

## MDE

SE(ΔM1-sel, z-score vs prior-only) = 0.0239; G_DEV = 89, G_TEST = 87; MDE = 2.5·SE·sqrt(G_DEV/G_TEST) = 0.0605 (see results/smoke/mde.json).

Raw metric dicts (incl. coverage-significant): `smoke_dev.json`. T-unseen not computed.
