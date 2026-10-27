# B1 — kinase–substrate library variants (DEV folds 0–4 pooled)

GT_FREEZE sha256 `1e6f791c6b85f15d`; 372 GT conditions, 408 configs (all in `results/dev_ledger.tsv`, branch B1). Scripts: `scripts/exp_b1_grid.py` → `scripts/exp_b1_report.py`; code `src/kinact/explore/lib_variants.py` (definitions in its docstring). Single-sample scoring, no GT information in any method. `wz` = Σwx/√Σw²/σ_c, `mean` = Σwx/Σw, `ridge` = ridge t (λ=1). Predicted libs: `meas` = top-k among sites measured in the condition, `glob` = global top-k then restricted to measured; weight (1−(r−1)/k)^γ. Mixtures: curated (w=1) + predicted top-k measured (γ=1, w·α) for kinases with < m measured curated substrates.

## 1. Reference, controls, best per library family (M1-sel; harmony bins b*; source; sign)

| config | M1-sel | M1-all | scored | cov | harm | b0 | b1-4 | b5-19 | b≥20 | H | D | O | s− | s+ | nsig |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| PSP z-score (psp_all, min3) | 0.515 | 0.512 | 0.633 | 0.578 | 0.517 | 0.501 | 0.508 | 0.531 | 0.527 | 0.518 | 0.560 | 0.507 | 0.513 | 0.536 | 24 |
| prior-only (OOF freq, oracle sign) | 0.687 | 0.543 | 0.687 | 1.000 | 0.726 | 0.625 | 0.705 | 0.735 | 0.840 | 0.753 | 0.536 | 0.709 | 0.695 | 0.769 | 53 |
| wz|r2np@200/glob/g0|min1 | 0.664 | 0.621 | 0.668 | 0.999 | 0.662 | 0.619 | 0.659 | 0.741 | 0.630 | 0.573 | 0.890 | 0.667 | 0.641 | 0.798 | 53 |
| wz|r2np@50/meas/g1|min1 | 0.663 | 0.635 | 0.663 | 1.000 | 0.661 | 0.625 | 0.654 | 0.732 | 0.631 | 0.575 | 0.908 | 0.653 | 0.644 | 0.779 | 51 |
| wz|curated_union+r2np@100/g1*a1/minf|min1 | 0.660 | 0.639 | 0.660 | 1.000 | 0.662 | 0.632 | 0.640 | 0.724 | 0.653 | 0.600 | 0.883 | 0.614 | 0.644 | 0.771 | 51 |
| wz|psp_all+r2np@100/g1*a1/minf|min1 | 0.658 | 0.639 | 0.658 | 1.000 | 0.660 | 0.635 | 0.641 | 0.725 | 0.639 | 0.596 | 0.884 | 0.610 | 0.642 | 0.772 | 49 |
| mean|r3@200/glob/g0|min1 | 0.626 | 0.614 | 0.637 | 0.970 | 0.639 | 0.599 | 0.621 | 0.707 | 0.630 | 0.533 | 0.875 | 0.643 | 0.602 | 0.783 | 29 |
| wz|curated_union+r3@50/g1*a0.1/m10|min1 | 0.625 | 0.615 | 0.625 | 1.000 | 0.628 | 0.575 | 0.649 | 0.669 | 0.619 | 0.587 | 0.878 | 0.589 | 0.604 | 0.751 | 28 |
| wz|r3@200/meas/g4|min1 | 0.625 | 0.640 | 0.625 | 1.000 | 0.639 | 0.582 | 0.612 | 0.728 | 0.632 | 0.526 | 0.870 | 0.585 | 0.606 | 0.769 | 37 |
| wz|psp_all+r3@100/g1*a1/minf|min1 | 0.622 | 0.643 | 0.622 | 1.000 | 0.635 | 0.578 | 0.616 | 0.711 | 0.636 | 0.545 | 0.873 | 0.572 | 0.604 | 0.759 | 32 |
| wz|curated_union+cantley@50/g1*a0.3/m10|min1 | 0.611 | 0.592 | 0.637 | 0.997 | 0.616 | 0.576 | 0.622 | 0.660 | 0.606 | 0.572 | 0.802 | 0.605 | 0.596 | 0.707 | 38 |
| wz|psp_all+cantley@100/g1*a1/minf|min1 | 0.597 | 0.593 | 0.622 | 0.997 | 0.604 | 0.580 | 0.584 | 0.645 | 0.608 | 0.553 | 0.767 | 0.558 | 0.580 | 0.709 | 27 |
| wz|cantley@100/meas/g0|min1 | 0.589 | 0.583 | 0.616 | 0.997 | 0.589 | 0.570 | 0.572 | 0.636 | 0.578 | 0.528 | 0.686 | 0.562 | 0.578 | 0.669 | 23 |
| wz|cantley@200/glob/g0|min1 | 0.564 | 0.540 | 0.591 | 0.929 | 0.565 | 0.554 | 0.569 | 0.576 | 0.561 | 0.531 | 0.662 | 0.538 | 0.548 | 0.654 | 16 |
| wz|curated_union|min1 | 0.563 | 0.543 | 0.650 | 0.860 | 0.565 | 0.521 | 0.573 | 0.594 | 0.572 | 0.557 | 0.693 | 0.552 | 0.554 | 0.627 | 40 |
| wz|psp_all|min1 | 0.545 | 0.534 | 0.630 | 0.822 | 0.547 | 0.514 | 0.557 | 0.561 | 0.557 | 0.546 | 0.653 | 0.526 | 0.537 | 0.602 | 37 |
| wz|invitro_sugiyama|min1 | 0.540 | 0.551 | 0.567 | 0.938 | 0.535 | 0.511 | 0.545 | 0.555 | 0.529 | 0.501 | 0.663 | 0.525 | 0.538 | 0.573 | 23 |
| wz|psp_invivo|min1 | 0.538 | 0.526 | 0.636 | 0.773 | 0.540 | 0.511 | 0.551 | 0.553 | 0.546 | 0.538 | 0.634 | 0.523 | 0.532 | 0.582 | 30 |

## 2. Best per statistic and top-10 overall

| config | M1-sel | M1-all | scored | cov | harm | b0 | b1-4 | b5-19 | b≥20 | H | D | O | s− | s+ | nsig |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| best wz: wz|r2np@200/glob/g0|min1 | 0.664 | 0.621 | 0.668 | 0.999 | 0.662 | 0.619 | 0.659 | 0.741 | 0.630 | 0.573 | 0.890 | 0.667 | 0.641 | 0.798 | 53 |
| best mean: mean|r2np@50/meas/g1|min1 | 0.660 | 0.633 | 0.660 | 1.000 | 0.655 | 0.601 | 0.655 | 0.731 | 0.631 | 0.575 | 0.908 | 0.661 | 0.646 | 0.755 | 51 |
| best ridge: ridge|psp_all+r2np@50/g1*a0.3/m3|min1 | 0.617 | 0.558 | 0.617 | 1.000 | 0.623 | 0.559 | 0.625 | 0.654 | 0.655 | 0.594 | 0.769 | 0.615 | 0.595 | 0.732 | 31 |
| wz|r2np@200/glob/g0|min1 | 0.664 | 0.621 | 0.668 | 0.999 | 0.662 | 0.619 | 0.659 | 0.741 | 0.630 | 0.573 | 0.890 | 0.667 | 0.641 | 0.798 | 53 |
| wz|r2np@50/meas/g1|min1 | 0.663 | 0.635 | 0.663 | 1.000 | 0.661 | 0.625 | 0.654 | 0.732 | 0.631 | 0.575 | 0.908 | 0.653 | 0.644 | 0.779 | 51 |
| wz|r2np@200/meas/g4|min1 | 0.663 | 0.638 | 0.663 | 1.000 | 0.665 | 0.636 | 0.643 | 0.742 | 0.637 | 0.583 | 0.892 | 0.624 | 0.644 | 0.782 | 56 |
| wz|r2np@100/meas/g4|min1 | 0.661 | 0.635 | 0.661 | 1.000 | 0.661 | 0.612 | 0.663 | 0.731 | 0.638 | 0.580 | 0.911 | 0.651 | 0.644 | 0.776 | 52 |
| wz|r2np@100/meas/g1|min1 | 0.661 | 0.635 | 0.661 | 1.000 | 0.663 | 0.633 | 0.640 | 0.742 | 0.636 | 0.578 | 0.882 | 0.622 | 0.643 | 0.776 | 52 |
| mean|r2np@50/meas/g1|min1 | 0.660 | 0.633 | 0.660 | 1.000 | 0.655 | 0.601 | 0.655 | 0.731 | 0.631 | 0.575 | 0.908 | 0.661 | 0.646 | 0.755 | 51 |
| wz|curated_union+r2np@100/g1*a1/minf|min1 | 0.660 | 0.639 | 0.660 | 1.000 | 0.662 | 0.632 | 0.640 | 0.724 | 0.653 | 0.600 | 0.883 | 0.614 | 0.644 | 0.771 | 51 |
| wz|r2np@50/meas/g0|min1 | 0.660 | 0.629 | 0.660 | 1.000 | 0.657 | 0.643 | 0.637 | 0.733 | 0.617 | 0.574 | 0.879 | 0.630 | 0.643 | 0.763 | 46 |
| mean|r2np@100/meas/g4|min1 | 0.659 | 0.634 | 0.659 | 1.000 | 0.656 | 0.593 | 0.665 | 0.729 | 0.636 | 0.580 | 0.911 | 0.662 | 0.646 | 0.756 | 52 |
| mean|r2np@200/glob/g0|min1 | 0.659 | 0.617 | 0.663 | 0.999 | 0.655 | 0.619 | 0.658 | 0.727 | 0.614 | 0.574 | 0.887 | 0.690 | 0.638 | 0.785 | 39 |

## 3. Factor effects (mean M1-sel over the grid cells sharing the factor; wz+mean, predicted only)

- **pred**: cantley: 0.555 (b0 0.541, b1-4 0.555, b5-19 0.571, b≥20 0.557); r2np: 0.633 (b0 0.588, b1-4 0.644, b5-19 0.677, b≥20 0.629); r3: 0.607 (b0 0.560, b1-4 0.591, b5-19 0.698, b≥20 0.608)
- **mode**: global: 0.583 (b0 0.553, b1-4 0.587, b5-19 0.617, b≥20 0.586); measured: 0.614 (b0 0.573, b1-4 0.606, b5-19 0.680, b≥20 0.611)
- **k**: 25: 0.585 (b0 0.553, b1-4 0.589, b5-19 0.623, b≥20 0.582); 50: 0.595 (b0 0.557, b1-4 0.596, b5-19 0.643, b≥20 0.595); 100: 0.604 (b0 0.568, b1-4 0.600, b5-19 0.659, b≥20 0.606); 200: 0.609 (b0 0.574, b1-4 0.600, b5-19 0.669, b≥20 0.610)
- **gamma**: 0: 0.600 (b0 0.564, b1-4 0.594, b5-19 0.653, b≥20 0.602); 1: 0.600 (b0 0.565, b1-4 0.598, b5-19 0.652, b≥20 0.600); 4: 0.595 (b0 0.561, b1-4 0.598, b5-19 0.640, b≥20 0.594)
- **stat**: mean: 0.597 (b0 0.558, b1-4 0.597, b5-19 0.647, b≥20 0.595); wz: 0.600 (b0 0.568, b1-4 0.596, b5-19 0.650, b≥20 0.601)
- **mixture alpha**: 0.1: 0.614; 0.3: 0.616; 1.0: 0.613
- **mixture m**: 10: 0.618; 3: 0.604; None: 0.622

## 4. Selection optimism: nested CV (choose config by M1-sel on 4 folds, score the 5th; pooled OOF)

- Apparent best (all 5 folds): **wz|r2np@200/glob/g0|min1** M1-sel 0.664.
- **Nested OOF M1-sel 0.638** (M1-all 0.631, scored 0.639, cov 1.000, harmony 0.636); optimism = +0.026.
- Chosen per held-out fold: fold0: wz|curated_union+r2np@50/g1*a1/minf|min1 (train 0.679, held-out 0.617); fold1: wz|r2np@200/meas/g4|min1 (train 0.680, held-out 0.701); fold2: wz|r2np@200/glob/g1|min1 (train 0.670, held-out 0.685); fold3: wz|r2np@200/glob/g0|min1 (train 0.670, held-out 0.665); fold4: wz|r2np@200/glob/g0|min1 (train 0.661, held-out 0.682)
- Nested OOF restricted to one family: curated_static_only 0.563; predicted_only 0.641; mixtures_only 0.646; wz_only 0.638

## 5. Shuffled-library control (best config; site-label bijection, per-kinase counts/weights kept)

| config | M1-sel | M1-all | scored | cov | harm | b0 | b1-4 | b5-19 | b≥20 | H | D | O | s− | s+ | nsig |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| best: wz|r2np@200/glob/g0|min1 | 0.664 | 0.621 | 0.668 | 0.999 | 0.662 | 0.619 | 0.659 | 0.741 | 0.630 | 0.573 | 0.890 | 0.667 | 0.641 | 0.798 | 53 |
| shuffled seed 0 | 0.480 | 0.490 | 0.480 | 1.000 | 0.480 | 0.487 | 0.468 | 0.483 | 0.480 | 0.486 | 0.452 | 0.499 | 0.507 | 0.374 | 1 |
| shuffled seed 1 | 0.494 | 0.492 | 0.494 | 0.999 | 0.495 | 0.502 | 0.497 | 0.472 | 0.511 | 0.532 | 0.480 | 0.467 | 0.488 | 0.512 | 0 |
| shuffled seed 2 | 0.521 | 0.502 | 0.521 | 1.000 | 0.522 | 0.534 | 0.495 | 0.528 | 0.529 | 0.514 | 0.566 | 0.496 | 0.533 | 0.520 | 3 |

## 6. M2 (matched universe) vs PSP z-score and paired group bootstrap (B=2000, seed 0, C8 basic CI)

- M2 best vs PSP z: 0.699 vs 0.633 (Δ +0.066; 507 pairs, 67 kinases, pair fraction 0.578)
- M2 nested_OOF vs PSP z: 0.678 vs 0.633 (Δ +0.045; 507 pairs, 67 kinases, pair fraction 0.578)
- ΔM1-sel best_vs_ref_psp_z: +0.149 [95% CI +0.106, +0.186], one-sided p=0.0005
- ΔM1-sel best_vs_prior_only: -0.023 [95% CI -0.058, +0.083], one-sided p=0.7461
- ΔM1-sel nested_OOF_vs_ref_psp_z: +0.123 [95% CI +0.072, +0.161], one-sided p=0.0005
- ΔM1-sel nested_OOF_vs_prior_only: -0.049 [95% CI -0.088, +0.057], one-sided p=0.8751
