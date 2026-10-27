# Smoke run on DEV (folds 0-4 pooled; no devconfirm/test)

Library `psp_all`, z-score = mean(x_s)*sqrt(n)/sd(x_c), min 3 substrates. 445 eligible GT conditions, 85 M1-sel groups. Controls: random (seed 0); prior-only and shift-prior use out-of-fold frequencies (other folds only); oracle. Runtime 3 s. Script: `scripts/smoke_dev.py`.

| method | M1-sel | M1-all | M1-sel-scored | pair coverage | harmony |
|---|---|---|---|---|---|
| zscore | 0.523 | 0.512 | 0.677 | 0.651 | 0.528 |
| random | 0.500 | 0.489 | 0.500 | 1.000 | 0.507 |
| prior_only | 0.707 | 0.580 | 0.707 | 1.000 | 0.736 |
| shift_prior | 0.525 | 0.522 | 0.525 | 1.000 | 0.555 |
| oracle | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |

## By source / sign (M1-sel)

| method | H | D | O | sign -1 | sign +1 |
|---|---|---|---|---|---|
| zscore | 0.517 | 0.569 | 0.520 | 0.515 | 0.547 |
| random | 0.483 | 0.501 | 0.522 | 0.510 | 0.471 |
| prior_only | 0.749 | 0.785 | 0.723 | 0.721 | 0.781 |
| shift_prior | 0.509 | 0.489 | 0.562 | 0.545 | 0.451 |
| oracle | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |

## Harmony per bin (M1-sel) and coverage-significant kinases

| method | 0 | 1-4 | 5-19 | >=20 | n sig (q<=0.10) |
|---|---|---|---|---|---|
| zscore | 0.501 | 0.516 | 0.544 | 0.550 | 30 |
| random | 0.504 | 0.497 | 0.486 | 0.541 | 0 |
| prior_only | 0.711 | 0.710 | 0.772 | 0.751 | 61 |
| shift_prior | 0.516 | 0.574 | 0.542 | 0.589 | 9 |
| oracle | 1.000 | 1.000 | 1.000 | 1.000 | 108 |

## Paired bootstrap (groups, B=2000, seed 0), M1-sel

- zscore_vs_prior_only: Δ=-0.184, SE=0.015, 95% CI [-0.241, -0.181], one-sided p=1.0000
- zscore_vs_random: Δ=+0.023, SE=0.013, 95% CI [-0.000, +0.052], one-sided p=0.0265

## MDE

SE(ΔM1-sel, z-score vs prior-only) = 0.0154; G_DEV = 85, G_TEST = 87; MDE = 2.5·SE·sqrt(G_DEV/G_TEST) = 0.0381 (see results/smoke/mde.json).

Raw metric dicts (incl. coverage-significant): `smoke_dev.json`. T-unseen not computed.
