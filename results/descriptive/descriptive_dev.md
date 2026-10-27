# Descriptive analyses on DEV (folds 0–4). Not used for model selection.

DESIGN §5 "Descriptive", §12 C4 / B13 / B15. Data: `cache/dev` (1,358 conditions) plus the DEV rows of
`results/ground_truth` (GT_FREEZE `1e6f791c…`). Nothing was read from `cache/sealed/` or `external/`, and no
ledger rows were written.
Code: `scripts/desc_common.py` (score matrices, cached in `cache/descriptive/`), `desc_m3m4.py`, `desc_gta.py`,
`desc_kk.py`, `desc_resnull.py`. Raw outputs: `results/descriptive/*.json|tsv`.

**Methods.**
- **core-v1**: the `kinact.core` package with `config/core_v1.json`, output z_vif. Its HPs were fixed on all of DEV, so the core-v1 rows are slightly in-sample (package M1-sel 0.668 vs 0.663 nested).
- **core-v2**: core-v1 plus late-fusion signatures (D11, α = 1). A condition in fold f uses the B7 `clean+contrast` signature fitted on the other 4 folds, with B7's (n_top, κ) pick for fold f, so it is fold-OOF. ρ_sig is fitted on the 7 DEV nulls. The fixed α = 1 differs from B7b's nested α (fold 4 used 0.5), which explains 0.686 here vs 0.680 in B7b.
- **PSP-z min3**: `zscore`, psp_all, at least 3 substrates.
- **cur-z min1**: `zscore`, curated_union, at least 1 substrate.
- **prior-only**: s_c·freq_k, with freq computed on the other folds' GT conditions.

## 1. M3, benchmarKIN-compatible (DEV GT conditions only, C4)
Reference code: `ksp.downstream.evaluate.perturbation_eval` (kinaseSubMain, imported read-only).
- Experiments: the 372 DEV GT conditions; targets = direct high/medium positives; sign = the condition's majority sign. benchmarKIN protocol: R `scale()` per condition × sign; every other *scored* kinase is a negative (no masking); conditions with no scored target are dropped.
- AUROC: MC-subsampled, balanced negatives, 1,000× (mean ± sd). Scaled rank: rank / n_scored, mean over kinases (lower = better).

| method | M3 AUROC (own coverage) | n pos / exps / kinases | scaled rank | M3 AUROC, M1-sel conds | **M3 matched cells** | scaled rank (matched) |
|---|---|---|---|---|---|---|
| core-v1 | 0.658 ± .005 | 1906 / 372 / 254 | 0.386 | 0.726 | **0.714** | 0.353 |
| core-v2 | **0.688** ± .005 | 1906 / 372 / 254 | 0.379 | **0.773** | **0.740** | 0.342 |
| PSP-z min3 | 0.668 ± .007 | 951 / 315 / 111 | **0.371** | 0.722 | 0.668 | 0.371 |
| cur-z min1 | 0.648 ± .005 | 1547 / 350 / 194 | 0.402 | 0.724 | 0.682 | 0.377 |
| prior-only | 0.681 ± .008 | 1205 / 323 / 254 | 0.457 | 0.720 | 0.598 | 0.559 |

"Matched cells" keeps only the (condition, kinase) cells scored by all four phospho methods: 951 positives, 111 kinases.
- **Own coverage:** benchmarKIN ignores unscored kinases (rewards low coverage): PSP-z, on 111 well-annotated kinases, edges past core-v1 (0.668 vs 0.658). **Matched cells:** core-v1 +0.046 and core-v2 +0.072 over PSP-z.
- **Prior-only:** benchmarKIN's `S==0 → NA` step drops prior-only's freq = 0 cells. That inflates its own-coverage number. On matched cells it falls to 0.598.
- The per-experiment top-k hit rate is small for every method (P(rank ≤ 20) ≤ 0.14). About 500 kinases are ranked per experiment.

## 2. M4: kinase-centric, within study (scored kinases only)
Unit: a (kinase, PMID) pair where the kinase is a positive in ≥ 1 GT condition and an evaluable negative in ≥ 1 other. AUROC of s_k·score_k (positive vs negative conditions); mean over studies, then macro over kinases.

| method | M4 (own units) | kinases / units | frac kinases > .5 | **M4, common units** (106 kinases / 147 units) |
|---|---|---|---|---|
| core-v1 | 0.628 | 248 / 320 | 0.69 | **0.691** |
| core-v2 | 0.629 | 248 / 320 | 0.69 | 0.688 |
| PSP-z min3 | 0.638 | 106 / 147 | 0.75 | 0.638 |
| cur-z min1 | 0.611 | 188 / 251 | 0.70 | 0.646 |
| prior-only | 0.186 | 248 / 320 | 0.04 | 0.205 |

- On the common units, both cores are about +0.05 over the z-scores.
- core-v2 equals core-v1 within study. The signatures help ranking *across kinases* (M1, M3) but not across conditions for the same kinase.
- prior-only M4 is a fold artefact: freq_k is lower in the positive's own fold, so the number is not interpretable.

## 3. GT-A: own activating site (self-sites masked in every method's input)
- Activating sites: PSP "(enzymatic) activity, induced" sites on U kinases; raw log2FC averaged over a kinase's measured activating sites. Input masking: `selfsites.self_site_keys` (psp+aloop + every PSP "activity" regulatory site on U kinases).
- Conditions: 907 qc-ok, non-inhibitor DEV conditions (no `drug_targets_manual` agent, not Hijazi); a kinase needs ≥ 20 conditions over ≥ 2 studies. Spearman ρ; null = 2,000 permutations *within study*; BH q ≤ 0.10. Control: the same site vector vs the scores of 5 random *other* kinases.

| method | kinases | macro ρ | median ρ | ρ>0 | sig (q≤.10) | macro ρ, common 88 | sig, common 88 | control: other kinase |
|---|---|---|---|---|---|---|---|---|
| core-v1 | 145 | 0.191 | 0.201 | 82% | 80 | 0.203 | 58 | 0.044 |
| core-v2 | 145 | 0.192 | 0.206 | 83% | **81** | 0.206 | **59** | 0.040 |
| PSP-z min3 | 88 | **0.214** | 0.227 | 83% | 48 | **0.214** | 48 | 0.091 |
| cur-z min1 | 130 | 0.193 | 0.191 | 85% | 71 | 0.209 | 54 | 0.084 |
| prior-only | – | n/a (no within-kinase variation outside GT conditions) | | | | | | |

- **Specificity:** the core scores carry the most kinase-specific signal. Own-site ρ minus other-kinase ρ is about 0.15 for the cores vs about 0.12 for the z-scores; the z-scores share more of the global response.
- **PSP-z** has the highest raw ρ on its 88 kinases. **Significant kinases:** on the common 88, the cores reach 58–59 significant kinases vs 48–54 for the z-scores. Over all kinases they reach 80–81, because they also cover 57 kinases that PSP-z cannot score.
- **core-v2 vs core-v1:** per kinase, v2 − v1 = +0.002 (24 kinases up, 17 down, 104 identical).

## 4. Kinase–kinase relationships (B15)
- Profiles: 1,325 qc-ok DEV conditions; kinases with ≥ 100 scored conditions; pairwise-complete Spearman over ≥ 100 shared conditions. Jaccard on the core-v1 library (curated_union + r2np@100).
- Gold: 790 unordered U pairs from PSP KS edges A → site on kinase B (A ≠ B) where the site is a PSP regulatory "activity" site of B. Pairs with any other PSP KS link are removed from the non-gold set.

| method | kinases / pairs / gold | AUROC abs(ρ) | signed ρ | diff. family | J=0 (gold n) | 0<J≤.05 | .05<J≤.2 | J>.2 |
|---|---|---|---|---|---|---|---|---|
| core-v1 | 503 / 125,956 / 790 | 0.682 | 0.664 | 0.668 | 0.549 (110) | 0.601 (580) | 0.566 (96) | 0.38 (4) |
| core-v2 | 503 / 125,956 / 790 | 0.661 | 0.645 | 0.645 | 0.533 (110) | 0.588 (580) | 0.537 (96) | 0.31 (4) |
| PSP-z min3 | 176 / 13,145 / 537 | 0.676 | 0.689 | 0.683 | 0.655 (68) | 0.619 (428) | 0.573 (40) | 0.32 (1) |
| cur-z min1 | 369 / 64,538 / 768 | **0.741** | **0.755** | **0.739** | **0.665** (101) | **0.673** (570) | **0.710** (94) | 0.70 (3) |

- **Within strata:** the overall AUROCs (0.66–0.74) shrink by 0.08–0.13 for the cores, because part of the overall signal is library overlap. Gold pairs are enriched in small-J strata, and those strata carry a different background.
- **Within strata, ranking:** curated-z leads, the PSP-z J=0 value rests on 68 gold pairs, and core-v1 is weakest at J=0 (0.549).
- **Smoothing:** the family kernel, borrow and dense r2np library smooth the core profiles; this hurts pairwise analysis though it helps M1. **core-v2 vs core-v1:** the two agree closely (ρ-vs-ρ r = 0.975). Signatures add a shared perturbation axis that lowers gold AUROC by 0.02.

**Top 20 novel high-correlation pairs (core-v2).** Candidates: J = 0, different Manning family, no PSP KS edge in either direction; 82,800 pairs. Full list: `kk_top20_novel.tsv`.

| pair | ρ v2 (v1) | pair | ρ v2 (v1) |
|---|---|---|---|
| MAP2K2–MAPK15 | .58 (.38) | MAP2K2–NLK | .58 (.40) |
| MAP2K2–MAPK6 | .57 (.35) | MAP2K1–MAPK15 | .57 (.37) |
| MAP2K1–MAPK4 | .57 (.36) | CLK4–SRPK3 | .56 (.56) |
| MAP2K1–MAPK6 | .56 (.34) | ROR1–ZAP70 | .54 (.54) |
| CAMK2A–STK33 | .52 (.52) | ATR–CDK9 | −.51 (−.08) |
| CAMK2B–STK33 | .51 (.51) | ALK–RPS6KB1 | .49 (.24) |
| ERBB2–MAP2K2 | .49 (.33) | MAP2K1–RPS6KA2 | .48 (.27) |
| ATR–CDK3 | −.48 (−.12) | ERBB2–MAP2K1 | .48 (.34) |
| INSRR–TNK1 | .47 (.47) | PRKG1–STK32B | .47 (.47) |
| KSR1–PRKCH | −.47 (−.47) | PRKCH–RAF1 | −.47 (−.47) |

Caveats:
- **Signature axis:** half the pairs involve signature kinases (MAP2K1/2, ERBB2, ALK, ATR), and v2 inflates ρ there (ATR–CDK9 −.08 → −.51): mostly a shared MEK/ERK or DNA-damage perturbation axis, not new KS edges.
- **Biology:** MEK → atypical MAPK/NLK, CLK–SRPK and RAF/KSR–PKC anti-correlation are plausible; ROR1–ZAP70 and INSRR–TNK1 likely reflect the generic pY axis. None is validated (hypotheses only).

## 5. Residue-matched null library (r2np audit follow-up)
Every DEV-measured core-v1 edge of a kinase (56,831 edges; S 63%, T 21%, Y 16%) is replaced by a random measured site *of the same residue* (without replacement; weight and source kept), then scored with the same core-v1 model on the DEV GT conditions, 3 seeds.

| library / method | M1-sel | M1-sel-scored | cov | H | D | O |
|---|---|---|---|---|---|---|
| core-v1 (real library) | 0.668 | 0.668 | 1.00 | 0.630 | **0.911** | 0.683 |
| core-v2 | 0.686 | 0.686 | 1.00 | 0.657 | 0.911 | 0.706 |
| **residue-matched null (core-v1 model), mean of 3** | **0.529** ± .013 | 0.529 | 1.00 | 0.540 | 0.583 | 0.497 |
| PSP-z min3 | 0.515 | 0.633 | 0.58 | 0.518 | 0.560 | 0.507 |
| cur-z min1 | 0.563 | 0.650 | 0.86 | 0.557 | 0.693 | 0.552 |
| prior-only | 0.687 | 0.687 | 1.00 | 0.753 | 0.536 | 0.709 |

Per seed, the null gives M1-sel 0.544 / 0.526 / 0.518 and D 0.546 / 0.612 / 0.593.
- **Model plus residue only** (full coverage, Y kinases get pY sites): 0.529. That is above PSP-z's coverage-penalised 0.515 but far below core-v1.
- **Site identity** adds about 0.14 M1-sel overall (+0.09 H, +0.33 D, +0.19 O). **D stratum:** about 0.58 of D comes from residue type alone, which is less than the 0.70 of the audit's class-level null. The rest comes from the library.

## Summary
- **M3 (benchmarKIN-style):** ranks core-v2 > core-v1 > PSP-z on matched cells (0.740 / 0.714 / 0.668). With own coverage, benchmarKIN's coverage bias brings PSP-z and prior-only level with or above core-v1.
- **M4 (within study):** both cores are about +0.05 over the z-scores on common units; v2 = v1.
- **GT-A:** all methods track the kinase's own activating site (macro ρ about 0.2). The cores have the most significant kinases and the most kinase-specific signal.
- **Kinase–kinase:** recovering regulatory KS edges is modest once library overlap is controlled. curated-z is best there and the cores are smoothed by their prior; v2 is slightly worse than v1.
- **Residue-matched null:** core-v1's advantage is not explained by residue type or model structure. The null scores 0.53 vs 0.67.
