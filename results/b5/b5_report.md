# Branch 5 — kinase self-sites as an auxiliary activity signal (DEV, folds 0-4 pooled OOF)

GT: `GT_FREEZE.json` sha256 `1e6f791c6b85f15dfe51d7eb045bf392a6b25be576b10974f7bba6c5d887baac`. Used only `cache/dev/` and
the DEV files in `results/ground_truth/`. Code: `src/kinact/explore/selfsites.py`, `scripts/exp_b5_{common,signal,combine,exclude}.py`.
Outputs: `results/b5/{signal,combine,exclude}.json`, `combine_tables.md`, `grid.tsv` (152 cfgs), `inhibitor_target_selfsignal.tsv`; ledger 172 rows (B5).
**Verdict: [-] as a core add-on**: real signal, but the B2 core already uses most self-sites as substrates (best nested
ΔM1-sel +0.003, n.s.). Keep the exclusion helper for GT-A.

## Self-site table (label-free)
- **psp**: PSP Regulatory_sites on a U kinase (human, or mapped via SITE_GRP_ID); "(enzymatic) activity, induced" → +1,
  "…inhibited" → −1, conflicting sites dropped: 1,073 sites / 276 kinases (234 inhibitory).
- **aloop**: S/T/Y between DFG and APE (kinaseSubMain domain table), +1; checked on MAPK1 T185/Y187, AKT1 T308:
  2,021 sites / 389 kinases. **psp+aloop**: 2,811 sites / 427 kinases.
- **Signal**: y = (x − med_c)/(1.4826·MAD_c) (B2 standardisation); s_k(c) = mean of sign·y over measured self-sites,
  divided by the RMS of s over kinases in c.
- **Overlap with libraries**: 44% of psp self-sites are curated edges of the *same* kinase, 50% r2np@100 edges of the
  same kinase, 83% curated substrates of *some* kinase.

## (1) Coverage
| set | U kinases in table | site in DEV compendium | measured in ≥1 DEV GT cond | median kinases/GT cond | M1-sel pairs measured | M1-sel kinases | all pairs measured |
|---|---|---|---|---|---|---|---|
| psp | 276 | 216 | 192 | 50 | 304/877 (35%) | 41/106 | 512/1906 |
| aloop | 389 | – | 183 | 28 | 178 (20%) | 35 | 277 |
| psp+aloop | 427 | – | 245 | 61 | 339 (39%) | 47 | 563 |

Measured share of M1-sel target pairs (psp): H 19%, D 45%, O 37%; s=−1 32%, s=+1 43%.

## (2) Standalone value (scored-only AUROC on pairs where the target's self-site is measured)
| set | M1-sel (½ rule) | M1-sel-scored | s=−1 | s=+1 | H | D | O | M1-all-scored |
|---|---|---|---|---|---|---|---|---|
| psp | 0.505 (cov .35) | **0.726** | 0.728 (205 pairs) | 0.644 (99) | 0.594 | 0.734 | 0.750 | 0.624 |
| aloop | 0.503 (cov .20) | 0.705 | 0.667 | 0.723 | 0.425 (12) | 0.707 | 0.729 | 0.625 |
| psp+aloop | 0.508 (cov .39) | 0.706 | 0.697 | 0.677 | 0.607 | 0.731 | 0.726 | 0.608 |

- **Sign agreement** (psp, all measured pairs): sign(s_k) = GT sign in 60% of inhibitor pairs (mean s −0.48) and in 73% of
  stimulus pairs (mean s +0.94).
- **Matched universe** (M2, psp-measured pairs): core z 0.781 vs self-sites 0.726. By sign, s−: 0.733 vs 0.728 and
  s+: 0.815 vs 0.644. By source, H 0.73 vs 0.59, D 0.94 vs 0.73, O 0.746 vs 0.750. The core is never worse except a tie on O.
- **Correlation**: Spearman(s, core z) over measured cells is only 0.18; where they differ, the core is usually right.

**Paradoxical increases** (inhibitor targets, psp; frac_up = fraction of pairs with s>0; 50% means no information):
| class | pairs | mean s | frac_up | note |
|---|---|---|---|---|
| TK (Src, Abl, EGFR, Eph, ALK, InsR) | 107 | −1.28 | 0.24 | dasatinib→LYN/EPHA2 strongly down; autophosphorylation sites are a faithful readout |
| AGC overall | 110 | −0.06 | 0.53 | uninformative |
| — PKC (Go6983 etc.) | 43 | +0.17 | 0.63 | paradoxical (known: ATP-site binders protect priming/A-loop sites) |
| — PKA (H89, AT13148) | 9 | +0.57 | 1.00 | paradoxical |
| — AKT (GDC0068, GSK690693, AT13148, Akt-1-2) | 12 | −0.79 | 0.17 | **no** paradox at the measured sites (mostly not T308/S473); only MCF7_AT13148 goes up |
| — RSK/S6K | 23 | −0.50 | 0.39 | weakly informative |
| STE7 = MEK1/2 (allosteric MEKi, refametinib, cobimetinib) | 35 | −0.17 | 0.31 | mixed; several MEKi raise S218/S222 (feedback relief) |
| ERK1/2 (GDC0994) | 6 | ≈0 | 0.67 | ERK inhibitor raises pT202/Y204 in 2/3 cell lines |
| MAP3K7/TAK1 (Hijazi off-target calls) | 15 | +0.37 | 0.80 | paradoxical or label-dubious |
| CMGC (CDK, GSK3, DYRK) | 72 | −0.26 | 0.39 | weak |
| MTOR | 17 | −0.85 | 0.18 | faithful |

- **RAF**: not assessable (GT drops paradoxical RAF positives). **By source**: frac_up H 0.52 (no information), D 0.19, O 0.28.
- **Stimuli**: MAPK family 0.94 up (36 pairs), EGFR 0.78, PIKK/CAMKL 1.0; MAP3K7 0.38 and JAK 0.0 go against the label.

## (3) Combination with the B2 reference core
Reference = fixed B2 config (τ .125, ρ_f 3, ρ_g .25, α 1, borrow); nested CV with the reference (w=0/λ=0) as a candidate.
- **late**: core/RMS_c + w·T(s); 3 HPs (set, mode, w ∈ {.1,.25,.5,1,2}); input-only modes: `add` s; `agree` s where its
  sign agrees with the core; `mag` |s|·sign(core) (sign-agnostic "this kinase changed"); `up` max(s,0) (ignores decreases).
- **bayes**: extra observations y_i = λ·sign_i·a_k + ε_i (λ ∈ {.25,.5,1,2,4}), propagated through K; 2 HPs (set, λ).

| variant (nested) | M1-sel | M1-all | harmony | H | D | O | s− | s+ | T-unseen | ΔM1-sel vs ref [95% CI], p | ΔM1-all | Δharmony | M2 A vs B |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ref_z | 0.673 | 0.619 | 0.673 | 0.647 | 0.911 | 0.681 | 0.649 | 0.807 | 0.650 | – | – | – | – |
| bayes_z (psp, λ 1–4) | 0.675 | 0.620 | 0.674 | 0.646 | 0.901 | 0.683 | 0.649 | 0.804 | 0.660 | +0.002 [−0.001,+0.007], .19 | +0.001 | +0.001 | .675/.673 |
| late_z (all modes) | 0.671 | 0.617 | 0.669 | 0.646 | 0.911 | 0.678 | 0.647 | 0.793 | 0.656 | −0.002 [−0.007,+0.002], .86 | −0.002 | −0.004 | .671/.673 |
| late_z_add | 0.674 | 0.619 | 0.673 | 0.647 | 0.912 | 0.682 | 0.649 | 0.796 | 0.654 | +0.001 [−0.003,+0.005], .35 | +0.000 | −0.000 | .674/.673 |
| late_z_mag | 0.670 | 0.618 | 0.669 | 0.641 | 0.913 | 0.681 | 0.647 | 0.793 | 0.652 | −0.003, .85 | −0.002 | −0.004 | .670/.673 |
| ref_mean | 0.674 | 0.617 | 0.674 | 0.645 | 0.911 | 0.683 | 0.650 | 0.809 | 0.653 | – | – | – | – |
| **bayes_mean (psp, λ=2 in all folds)** | **0.677** | 0.619 | 0.676 | 0.647 | 0.902 | 0.685 | 0.652 | 0.810 | 0.661 | +0.003 [−0.001,+0.008], .079 | +0.002 | +0.002 | .677/.674 |
| late_mean_add | 0.676 | 0.615 | 0.675 | 0.644 | 0.905 | 0.687 | 0.650 | 0.796 | 0.655 | +0.002 [−0.004,+0.008], .32 | −0.002 | +0.001 | .676/.674 |
| late_mean_up | 0.675 | 0.616 | 0.675 | 0.644 | 0.909 | 0.684 | 0.648 | 0.820 | 0.654 | +0.001, .32 | −0.002 | +0.001 | .675/.674 |
| any (all families) | 0.669 | 0.615 | 0.667 | 0.644 | 0.910 | 0.675 | 0.647 | 0.789 | 0.655 | −0.004, .85 | −0.004 | −0.006 | .669/.673 |

- **Harmony bins** (0|1-4|5-19|≥20): bayes_mean .600|.694|.753|.659 vs ref_mean .597|.696|.747|.656. Pooled non-nested optimum 0.680.
- **Adoption**: none passes (p<.05). bayes_mean is fold-stable (psp, λ=2) and non-inferior (worst D −0.009 < 0.02
  margin): a harmless optional switch, not a gain. Late rules are fold-unstable; `mag`/`agree` never help; `up` helps s+
  slightly (0.820 vs 0.809) at a cost on inhibitor pairs.

## Why it does not help
1. **Already in the core**: half of the self-sites are library edges of their own kinase, most others substrates of
   upstream kinases. Masking them from the input costs the core −0.010 [−0.022,+0.002] (psp) and −0.015 [−0.029,−0.003]
   (conservative set).
2. **The core is better where both exist** (matched 0.78 vs 0.73); a weaker correlated readout adds little.
3. **Paradox/feedback** (PKC/PKA, MEK/ERK, TAK1) and the uninformative Hijazi signal cancel the gain for inhibitors;
   input-only rules (`mag`, `up`, `agree`) cannot separate these cases without the GT class.

## GT-A exclusion helper (needed later)
`selfsites.self_site_keys(universe, which="psp+aloop", include_all_regulatory=True)` → 2,871 keys (1,055 in DEV cache):
psp + A-loop sites + every PSP regulatory site on a U kinase mentioning "activity". `selfsites.mask_self_sites(X, keys)`
sets them to NaN (DataFrame or Series) before scoring. Use the conservative set for GT-A: methods with curated/r2np
libraries use these sites as substrates, so unmasked GT-A would be circular.

## Caveats
- DEV only; coverage 41/106 M1-sel kinases; paradox classes small (≤43 pairs/family, often one drug): descriptive.
- The A-loop table is a DFG→APE regex on the sibling kinaseSubMain domain table (label-free); motif variants may be
  missed or mis-bounded, and A-loop sites are assumed activating.
