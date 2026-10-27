# Smoke run on DEV (folds 0-4 pooled; no devconfirm/test)

Library `psp_all`, z-score = mean(x_s)*sqrt(n)/sd(x_c), min 3 substrates. 433 eligible GT conditions, 74 M1-sel groups. Controls: random (seed 0); prior-only and shift-prior use out-of-fold frequencies (other folds only); oracle. Runtime 2 s. Script: `scripts/smoke_dev.py`.

| method | M1-sel | M1-all | M1-sel-scored | pair coverage | harmony |
|---|---|---|---|---|---|
| zscore | 0.521 | 0.513 | 0.666 | 0.626 | 0.525 |
| random | 0.476 | 0.499 | 0.476 | 1.000 | 0.487 |
| prior_only | 0.706 | 0.544 | 0.706 | 1.000 | 0.734 |
| shift_prior | 0.547 | 0.526 | 0.547 | 1.000 | 0.564 |
| oracle | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |

## By source / sign (M1-sel)

| method | H | D | O | sign -1 | sign +1 |
|---|---|---|---|---|---|
| zscore | 0.522 | 0.566 | 0.514 | 0.512 | 0.559 |
| random | 0.480 | 0.510 | 0.484 | 0.471 | 0.531 |
| prior_only | 0.768 | 0.788 | 0.728 | 0.711 | 0.779 |
| shift_prior | 0.514 | 0.480 | 0.586 | 0.534 | 0.625 |
| oracle | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |

## Harmony per bin (M1-sel) and coverage-significant kinases

| method | 0 | 1-4 | 5-19 | >=20 | n sig (q<=0.10) |
|---|---|---|---|---|---|
| zscore | 0.501 | 0.511 | 0.545 | 0.541 | 28 |
| random | 0.491 | 0.484 | 0.443 | 0.528 | 0 |
| prior_only | 0.590 | 0.733 | 0.791 | 0.823 | 46 |
| shift_prior | 0.538 | 0.502 | 0.610 | 0.608 | 13 |
| oracle | 1.000 | 1.000 | 1.000 | 1.000 | 94 |

## Paired bootstrap (groups, B=2000, seed 0), M1-sel

- zscore_vs_prior_only: Δ=-0.184, SE=0.016, 95% CI [-0.244, -0.184], one-sided p=1.0000
- zscore_vs_random: Δ=+0.045, SE=0.011, 95% CI [+0.026, +0.070], one-sided p=0.0010

## MDE

SE(ΔM1-sel, z-score vs prior-only) = 0.0156; G_DEV = 74, G_TEST = 87; MDE = 2.5·SE·sqrt(G_DEV/G_TEST) = 0.0359 (see results/smoke/mde.json).

Raw metric dicts (incl. coverage-significant): `smoke_dev.json`. T-unseen not computed.
