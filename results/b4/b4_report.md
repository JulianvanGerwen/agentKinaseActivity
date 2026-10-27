# B4 — data-driven substrate/reporter expansion from compendium co-regulation (DEV folds 0–4, pooled OOF)

GT_FREEZE sha256 `1e6f791c6b85f15dfe51d7eb045bf392a6b25be576b10974f7bba6c5d887baac`; data `cache/dev/quant_dev.npy` only
(+ `results/ground_truth/` DEV files). Code: `src/kinact/explore/coreg.py`; scripts `scripts/exp_b4_precompute.py` →
`scripts/exp_b4_grid.py` → `scripts/exp_b4_eval.py`; tests `tests/test_coreg.py` (5 pass; full suite 114 pass).
Raw: `results/b4/grid.tsv` (615 evals), `results/b4/b4_eval.json`, `results/b4/b4_tables.md`. Ledger: branch `B4`, 624 rows.
**Verdict: [-] for the B2 core (nested Δ −0.010, significantly worse); small gain only for plain weighted z, which is
already matched by the B1 library. Kill per DESIGN §7 (no DEV gain on the current best).**

## Leakage control (DESIGN §4, r1_data item 11)
- To score DEV fold f, all co-regulation statistics use ALLOWED = DEV conditions with qc ok, not in f, and not sharing
  (PMID, Control, Sample) with any fold-f condition (GT or not). This removes 114–160 extra conditions per fold;
  759–968 allowed conditions remain (fold0 is the largest fold). `coreg.allowed_mask` refuses TEST/DEV-confirm IDs.
- `tests/test_coreg.py::test_heldout_values_do_not_enter_fit` puts noise in every held-out, shared-control and
  qc-failed column and requires the fit (profile–site r, n, site–site neighbours; q = 0 and 2 factors) to stay the
  same. It also checks that perturbing one allowed column changes the fit, so the test is sensitive.
  `test_exp_b4_cache_is_fold_safe` re-checks the actual `cache/b4/allowed.json`.
- No GT labels anywhere in the fits. All hyperparameters are shared; there are no per-kinase or per-edge parameters.
- Per-condition standardisation (median/MAD) uses only that condition.
- Not covered: DEV conditions that share a control with a TEST condition. TEST control metadata is outside the files I
  may read. The DEV cache already drops TEST groups, but the final pipeline must apply this rule at G2.
- Nested-selection nuance (same as B2's v): the inner score for fold g ≠ f uses an expansion fitted on data that
  includes fold f's unlabelled values. The outer-fold score never uses fold-f data.

## Methods (all label-free; fixed B2 core τ .125, ρ_f 3, ρ_g .25, α 1 [curated_union + r2np@100], median, borrow, no noise)
1. **Kinase-anchored co-regulation expansion (ideas 1+2).**
   - Profile_k(c) = weighted z of k's seed substrates over allowed conditions. Seeds = curated_union, or
     curated_union + r2np@100.
   - Pairwise-complete Pearson r between each profile and every site with ≥30 allowed observations (GPU).
   - Rank by shrunk |r|·n/(n+30) and keep the top m ∈ {10, 25, 50, 100, 200} non-seed sites.
   - Edge weight = β·(shrunk r), β ∈ {.25, .5, 1}, sign ∈ {signed, positive-only}. Edges are added to W of the B2 core.
   - Global-factor variant q ∈ {0, 1, 3}: partial out the top-q SVD condition factors of the kinase-profile matrix
     from both profiles and sites. This is the low-rank / kinase-anchored factor idea.
2. **Data-driven propagation (idea 3, RoKAI-like).**
   - Graph: top-20 co-regulated neighbours of each of the 26,600 library sites.
   - Per condition, an unmeasured library site gets the signed r-weighted mean of its measured neighbours
     (|r| ≥ rmin ∈ {.3, .5}, ≥ min_nbr ∈ {1, 3} neighbours).
   - Imputed values enter with precision γ ∈ {.25, .5, 1}.
3. Each config is scored three ways: B2 posterior z, B2 posterior mean, and plain weighted z (wz = Σwy/√Σw²) on the
   same library. Nested CV picks, per outer fold, the best config on the other 4 folds; the no-expansion reference
   is always one of the options.

Observations on the fits:
- Correlations are high: median |r| is 0.70 at rank 10 and 0.56 at rank 200.
- Hubs dominate. The top-50 lists of 503 kinases contain only 4,993 distinct sites, and one site appears in 171
  kinases' lists. The expansions mostly capture global response axes, not kinase-specific ones.
- Factor removal (q = 3) cuts hubs to max 68 lists but grid mean drops 0.675 → 0.672.

## Results (nested OOF, 372 GT conditions; ref = same pipeline, no expansion)
| variant | M1-sel | M1-all | cov | harm | b0 | b1-4 | b5-19 | b≥20 | H | D | O | T-unseen | nsig |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ref B2 fixed core, z | 0.673 | 0.619 | 1.000 | 0.673 | 0.596 | 0.690 | 0.749 | 0.656 | 0.647 | 0.911 | 0.681 | 0.650 | 55 |
| B4 expansion, z | 0.663 | 0.615 | 1.000 | 0.668 | 0.597 | 0.685 | 0.736 | 0.652 | 0.641 | 0.884 | 0.663 | 0.642 | 53 |
| B4 imputation, z | 0.671 | 0.606 | 1.000 | 0.671 | 0.578 | 0.692 | 0.750 | 0.662 | 0.656 | 0.862 | 0.669 | 0.640 | 54 |
| ref B2 fixed core, mean | 0.674 | 0.617 | 1.000 | 0.674 | 0.597 | 0.696 | 0.747 | 0.656 | 0.645 | 0.911 | 0.683 | 0.653 | 53 |
| B4 expansion, mean | 0.665 | 0.611 | 1.000 | 0.669 | 0.596 | 0.689 | 0.737 | 0.652 | 0.633 | 0.893 | 0.668 | 0.652 | 50 |
| B4 all (exp+imp+none; z/mean) | 0.662 | 0.613 | 1.000 | 0.666 | 0.593 | 0.686 | 0.736 | 0.648 | 0.632 | 0.895 | 0.662 | 0.646 | 50 |
| ref wz (cur+r2np@100) | 0.634 | 0.612 | 0.983 | 0.634 | 0.567 | 0.666 | 0.666 | 0.639 | 0.585 | 0.912 | 0.637 | 0.639 | 40 |
| **B4 expansion, wz** | 0.660 | 0.631 | 1.000 | 0.664 | **0.645** | 0.662 | 0.691 | 0.658 | 0.631 | 0.913 | 0.648 | **0.684** | 47 |
| B4 imputation, wz | 0.644 | 0.617 | 1.000 | 0.645 | 0.590 | 0.664 | 0.683 | 0.643 | 0.616 | 0.912 | 0.608 | 0.666 | 38 |
| B2 nested (report ref) | 0.655 | 0.604 | 1.000 | 0.658 | 0.582 | 0.664 | 0.741 | 0.646 | 0.644 | 0.768 | 0.674 | 0.589 | 51 |
| B1 best wz r2np@200 glob (apparent) | 0.664 | 0.621 | 0.999 | 0.662 | 0.619 | 0.659 | 0.741 | 0.630 | 0.573 | 0.890 | 0.667 | 0.670 | 53 |
| prior-only | 0.687 | 0.543 | 1.000 | 0.726 | 0.625 | 0.705 | 0.735 | 0.840 | 0.753 | 0.536 | 0.709 | 0.324 | 53 |

T-unseen is the pooled-OOF definition from B2. The bins are the frozen harmony bins (quantified PSP substrates).

Paired group bootstrap (B = 2000, seed 0, C8 basic CI, one-sided p for Δ > 0); M2 is the matched-universe M1-sel-scored:
| comparison | ΔM1-sel | 95% CI | p | M2 A vs B (pairs) |
|---|---|---|---|---|
| B4 expansion z vs ref z | −0.010 | [−0.021, −0.003] | 0.99 | 0.663 vs 0.673 (877) |
| B4 expansion mean vs ref mean | −0.009 | [−0.020, −0.002] | 0.97 | 0.665 vs 0.674 (877) |
| B4 imputation z vs ref z | −0.002 | [−0.010, +0.003] | 0.71 | 0.670 vs 0.673 (877) |
| B4 all vs ref z | −0.011 | [−0.023, −0.002] | 0.98 | 0.662 vs 0.673 (877) |
| **B4 expansion wz vs ref wz** | **+0.026** | [+0.008, +0.044] | **0.001** | 0.656 vs 0.641 (862) |
| B4 expansion wz vs B1 best | −0.004 | [−0.038, +0.027] | 0.60 | 0.660 vs 0.668 (876) |
| B4 expansion z vs B1 best | −0.001 | [−0.036, +0.024] | 0.52 | 0.664 vs 0.668 (876) |
| B4 expansion z vs B2 nested | +0.008 | [−0.006, +0.028] | 0.20 | 0.663 vs 0.655 (877) |
| (ref z fixed core vs B2 nested) | +0.018 | [+0.010, +0.037] | 0.0005 | (core HPs chosen on all DEV) |
| B4 all vs prior-only | −0.025 | [−0.066, +0.075] | 0.76 | 0.662 vs 0.687 (877) |

The fixed B2 core was selected on all of DEV, so it carries B2's own selection optimism. The fair comparison for
B4 is therefore against the same fixed core, the "ref" rows. Against that core B4 is worse. The +0.008 against B2
nested is smaller than the fixed core's own +0.018 against B2 nested.

## Findings
1. **No gain on the Bayesian core.**
   - The apparent grid maximum is +0.006 (0.679 vs 0.673). That is below the selection noise, and nested CV turns
     it into −0.010 (CI excludes 0).
   - Folds disagree on the chosen configs (seed, q, m, sign).
   - Decryptm loses most (D 0.911 → 0.884) and O also drops (0.681 → 0.663).
   - Grid means: β .25 (0.675) > .5 (0.674) > 1 (0.671). Smaller expansions hurt less; none help reliably.
2. **Why:** the B2 core already pools signal through the family kernel and r2np@100. The co-regulated sites mostly
   add global response axes (hubs), not kinase-specific signal, so they add between-kinase correlation and reduce
   discrimination.
3. **Plain wz is the exception.**
   - Expansion helps wz significantly: +0.026, p = 0.001, with M2 also up (0.656 vs 0.641).
   - Gains are largest for understudied kinases: harmony bin 0 rises 0.567 → 0.645, and T-unseen reaches 0.684, the
     highest of any method so far. Coverage goes 0.983 → 1.0.
   - It only brings wz up to the B1 r2np@200 level (Δ −0.004 vs B1 best), so there is nothing new to adopt.
   - Possible lead: as a separate input for B6 stacking, since it is strong on T-unseen and bin 0.
4. **Imputation (idea 3) is neutral to harmful.**
   - z: −0.002 (n.s.); 4 of 5 folds chose "none". mean: −0.003.
   - D drops 0.91 → 0.86. min_nbr 3 and γ .25 are least bad.
5. **Factor residualisation (q = 1, 3) hurts.** The global axes it removes are partly real kinase biology (e.g. the
   MAPK/PI3K axes).
6. **Understudied kinases (B2 core):** bins 0 and 1–4 are unchanged (0.597 vs 0.596; 0.685 vs 0.690). No gain for
   understudied kinases on top of r2np + family kernel.

## Runtime (64-core node, RTX 4090)
- **Precompute.** About 37–49 s per fold for all variants: profiles + profile–site r for 2 seeds × 3 q, plus the
  site–site graph (7–9 s). Without factor variants it takes about 4–7 s per fold.
- **Production.** The fit runs once on all allowed conditions (≈ 8 s).
- **Per condition** (CPU, 1 thread, 200 DEV conditions, `SingleSample.score_one`): base 9.0 ms; with the chosen
  expansion 10.3 ms (+7.2k library rows, 107k vs 57k edges).
- Imputation adds one neighbour gather per condition (<1 ms). Grid 615 evals 8.6 min; eval + bootstrap 30 s.

## Caveats
- DEV only (205 configs × 3 outputs, nested CV); Pearson on clipped robust-z, Spearman not tried.
- Factor variants used SVD of the kinase-profile matrix only; NMF / masked PCA of the full site matrix not tried.
- The shared-control rule against TEST conditions is still open for G2 (see Leakage control).
