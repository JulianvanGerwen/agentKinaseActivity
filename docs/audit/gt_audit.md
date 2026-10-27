# Ground-truth audit (independent, pre-freeze) — 2026-10-27

Scope: DESIGN.md §3–4 (v3), GT_REPORT §9, every manual rule TSV in `results/ground_truth/`, the pipeline
(`scripts/gtlib.py`, `build_gt.py`, `make_split.py`, `make_dev_cache.py`, `gt_report.py`, `freeze_gt.py`) and
`tests/test_gt.py`. DEV rows only (`gt_pairs.tsv` / `gt_conditions.tsv`). Nothing under `cache/sealed/` was opened
and I did not read the git history. All fixes are in rule TSVs. Each changed row's note starts with `AUDIT 2026-10-27`
and gives a citation. In round 1 no code was changed. Round 2 (§3) added dose-parsing mechanics to `build_gt.py`; the biology stays in the TSVs.

## 1. Sample check (before fixes)

I drew a stratified random sample (pandas `sample`, random_state 0) of **68 DEV conditions** with role gt or
specificity: 15 decryptM dose conditions (incl. a combination ratio, antibodies, SHP099 and sub-potency doses),
12 Hijazi, 14 stimulus, 12 other inhibitor, 6 genetic or withdrawal, 4 reversed contrasts and 5 paradoxical-RAF
or specificity conditions. For each one I traced treatment text → rule → direct, downstream and masked kinases,
sign and confidence.

**Pass rate before fixes: 61/68 = 89.7 %.** The 7 failures:

| condition | problem |
|---|---|
| 21659604 `Starved` vs `Insulin` | INSR/IGF1R/PDPK1 were negatives, but the control was insulin-stimulated (reversed insulin) |
| 34086370 `MCF7_FGF7_1min`, `MCF7_FGF10_1min`, `MDAMB415_FGF7_1min` | ERK1/2 labelled direct medium at 1 min after FGF |
| 37188730 `IFN_5m` | the IFN is IFN-γ (paper abstract), so TYK2 should not be a positive and JAK2 should not be low |
| 20071362 `HB-EGF` vs `LPA` | ERK/AKT cascade labelled although LPA in the control activates ERK too |
| 34663829 `3h_TNF_CDK12IN3` | TNF 3 h labelled medium, but `DMSO_TNF_3h` (R081) is low; IKK activity is back near baseline by 3 h |

All other sampled conditions were judged correct, including signs of reversed contrasts (`Ctr` vs mitosis, DMSO vs
ALK-inhibitor, TNF vs SC-514), paradox drops (HCT116 / organoid vemurafenib and BRAFi → specificity), V600 lines
(WiDr, COLO205), decryptM dose gating and the antibody conversion. The Hijazi rows are correct *with respect to the
user kinhib matrix*. The matrix itself is a design input and I did not re-audit it.

After the fixes the failing DEV conditions are relabelled as intended. Some of them, e.g. `Starved`, `HB-EGF` and
`3h_TNF`, moved to TEST in the new split, so for those I checked the rule logic only.

## 2. Issues and fixes

Severity: **High** = a wrong positive or negative enters M1 (high or medium confidence). **Medium** = a systematic
over-confidence or a missing mask that creates wrong negatives. **Low** = confidence/consistency only, or affects
sensitivity labels only.

| # | sev | issue | fix (file / rule) | evidence |
|---|---|---|---|---|
| 1 | High | LY294002 (32034161, inferred 0.1 / 1 µM) had **MTOR as a direct, high-confidence target**; mTOR IC50 is ~5 µM. AKT downstream was medium although PI3K IC50 is 1.4 µM | R041 targets limited to class-I PI3K (outside U; MTOR is still masked by the PI3K mask); new R044b gives AKT downstream **low** for LY294002 | Brunn 1996 EMBO J (PMID 8895571); Vlahos 1994 JBC (PMID 8106507) |
| 2 | High | HB-EGF vs LPA: MAPK1/3 direct, plus AKT/MEK/RSK downstream, although LPA (the control) activates ERK and transactivates EGFR | R126 now labels EGFR only (medium); the cascade is masked through mask_class EGF | Daub 1996 Nature, doi:10.1038/379557a0 |
| 3 | High | Reversed insulin (`Starved`/`Starve` vs insulin, R247/R177): INSR, IGF1R and PDPK1 were negatives | added INSR, IGF1R, PDPK1 to the STARVATION mask (`pathway_masks.tsv`) | Manning & Toker 2017 Cell, doi:10.1016/j.cell.2017.04.001 |
| 4 | High | Dasatinib in A431/A549 (D002): **LCK** (T-cell restricted) and **KIT** were inhibited positives in epithelial lines that do not express them | D002 targets are now ABL1;SRC;LYN;YES1;FYN;EPHA2. LCK is masked as an off-target (Kd_app ≤ 3× dose). KIT is masked at ≥ 20 nM and is a negative below that dose, which is correct for an absent kinase | Human Protein Atlas cell-line RNA (LCK, KIT) |
| 5 | Medium | IFN in 37188730 was treated as "type unspecified" (TYK2/JAK2 low), but the paper uses **IFN-γ** | new stimulus class IFNG (JAK1, JAK2 direct medium) with an IFNG mask (TYK2, JAK3, ERK, AKT1 masked); R096 and R097b now use IFNG | 37188730 abstract (Commun Biol 2023); Schroder 2004 J Leukoc Biol, doi:10.1189/jlb.0603252 |
| 6 | Medium | FGF7/FGF10 at **1 min**: ERK1/2 were direct medium, but ERK activation and substrate phosphorylation lag the receptor by 2–5 min | R060 restricted to 8 min. New R060b (FGFR2 direct medium) and R060c (MAPK1/3 direct **low**) for 1 min | Kholodenko 1999 JBC, doi:10.1074/jbc.274.42.30169 |
| 7 | Medium | SP600125 at 0.1 / 1 µM (32034161) was medium, but its cellular IC50 is 5–10 µM | R031 → low | Bennett 2001 PNAS, doi:10.1073/pnas.251194298 |
| 8 | Medium | SB216763 at 0.1 µM was medium, but its cellular EC50 is 3.6 µM | new R029b gives 0.1 µM low; 1 µM stays medium (R029) | Coghlan 2000 Chem Biol, doi:10.1016/S1074-5521(00)00025-9 |
| 9 | Medium | SHP099 downstream ERK labels were gated at the biochemical SHP2 IC50 (71 nM); the cellular p-ERK IC50 is ~250 nM | `drug_potency_manual.tsv` SHP099 * 71 → 250 nM. SHP099_100nM becomes a specificity condition | Chen 2016 Nature, doi:10.1038/nature18621 |
| 10 | Medium | ATMi / DNA-PKi on the **HU** background were labelled high like on NCS; HU activates mainly ATR | R053 and R056 restricted to NCS. New R053b (HU+ATMi medium) and R056b (HU+DNA-PKi low) | Shiloh & Ziv 2013, doi:10.1038/nrm3546; Saldivar 2017, doi:10.1038/nrm.2017.67 |
| 11 | Low | TNF 3 h on a CDK-inhibitor background (R091) was medium, inconsistent with R081 (TNF 3 h low) | R091 restricted to 15 min. New R091b gives 3 h low | Hoffmann 2002 Science, doi:10.1126/science.1071914 |
| 12 | Low | UV mask missed the p38 outputs MK3, MSK1/2 and the paralogs p38β and JNK3 | added MAPK11, MAPKAPK3, RPS6KA4, RPS6KA5, MAPK10 to the UV mask | Deak 1998 EMBO J, doi:10.1093/emboj/17.15.4426 |

### Checked and confirmed (no change)

- **Potency table:** afatinib (0.5 / 14 / 1 nM), MK2206 (8 / 12 / 65), PD0325901 (0.33), refametinib (19 / 47),
  AZD8055 (0.8), dactolisib (mTOR 20.7, PI3Kα 4), AZD4547 (0.2 / 2.5 / 1.8), nintedanib (Hilberg 2008 values),
  CUDC-101 (2.4 / 15.7), lapatinib (10.8 / 9.2) and dasatinib match the primary papers.
- **Selumetinib MEK1 14 nM:** confirmed. For MEK2 the sources conflict (vendors give 53 nM or 530 nM). The Klaeger
  Kd_app is 52 nM, so I kept 50 nM.
- **Cell-line genotypes:** HL-60 NRAS Q61L, A549 KRAS G12S, HCT116 KRAS G13D, WiDr/COLO205/A375 BRAF V600E, WM239A
  BRAF-mutant and MCF7 BRAF/RAS-WT are correct.
- **§9 items 1–2, 4–8, 10:** reasonable as pre-registered. Stimulus first-tier lists are canonical. MEK is kept
  downstream, which is acceptable. The TK mask is conservative, it only removes negatives. The 10 µM default is
  conservative. Treating sorafenib as type II and LY3009120 as pan-RAF is correct. Treating the unnamed organoid
  BRAFi as paradox-inducing is conservative. The paralog rule only masks.
- **Literature off-target table:** errors in it can only remove negatives, never create a false label. Spot checks
  (Bain 2007, Davies 2000 pattern, MRT67307 → ULK1/2 and MARKs, BI-D1870 → PLK1/AURKB, Go6976 → PKD/CHK1) are
  consistent with the literature.

## 3. PI-requested changes (round 2, same day)

**(13) decryptM combination doses (High).** The labels have the form `<Drug1>+<Drug2>_<a>to<b>_<dose>`, where a:b is
the molar ratio of first to second drug. Zecha et al. 2023 (Fig. 4E legend, selumetinib:MK-2206 "concentration
ratio of 3:1") confirms this reading. Neither the paper nor the text I could access states which drug, or the
total, the label dose refers to.

The old parser assumed the label was the first drug's dose. Under that reading AZD4547 reaches 80 µM at the 1000 nM
gefitinib+AZD4547 label.

What changed:
- New rule file `combination_doses_manual.tsv` assigns a dose reading per series: `label_is_first`,
  `label_is_second`, `label_is_total` or `ambiguous`. All four series (gefitinib+AZD4547 1:80, lapatinib+AZD4547 1:1,
  selumetinib+MK2206 1:2 and 3:1) are set to **ambiguous**.
- `build_gt.py` (parsing mechanics only) handles ambiguous series as follows:
  - Each drug's dose is the maximum over the three readings, so off-target masking stays conservative.
  - Every direct and downstream label in those conditions is capped at **low**, which keeps them out of the metrics.
- 36 conditions are affected. None of them is in M1-sel any more.
- A new test (`test_combination_ambiguous_low_and_cellular_gate`) enforces this.
- If the decryptM supplement documents the reference drug, set the reading in the TSV and rebuild.

**(14) Cellular-potency override (Medium).** New rule file `drug_cellular_potency.tsv`. For a drug's primary target,
a documented cellular IC50 below the kinobead Kd_app or potency-table value becomes the decryptM dose gate. The rule
is applied to every primary target where I could document such a value:

| drug / target | cellular IC50 | replaces | source |
|---|---|---|---|
| gefitinib / EGFR | 54 nM | Kd_app 413 nM | Wakeling 2002, PMID 12384534 |
| dasatinib / ABL1 | 1 nM | Kd_app 5 nM | O'Hare 2005, PMID 15930265 (325× imatinib in Bcr-Abl cells) |

Checked but no row needed, because the cellular IC50 is at or above Kd_app / the potency value: lapatinib, AZD4547,
imatinib, MK-2206, AZD8055.

Not added, because I found no target-specific cellular IC50: nintedanib (Hilberg 2008 reports only a pooled
10–80 nM proliferation EC50), CUDC-101, dactolisib (Maira 2008 abstract has no values) and erlotinib (not
verifiable; not in decryptM).

Effect:
- Gefitinib 100 and 300 nM now count EGFR as a positive (previously specificity conditions).
- Dasatinib 1 and 3 nM now count ABL1 as a positive.

Unexpressed targets in the Hijazi set are left as given (user-provided GT); see the limitations below.

## 4. Pipeline re-run and final split

I ran `build_gt.py` → `make_split.py` → `make_dev_cache.py` → `gt_report.py` → `freeze_gt.py`. Results:

- pytest `tests/test_gt.py`: **17 passed** (one new test).
- `GT_FREEZE.json` was rewritten.
- Roles: gt 881, specificity 128, excluded 22, null 16.
- Pairs: direct 3759, downstream 705, masked 24135.

**Final split:** seed **90** (originally 79; 29 after round 1), d = 0.0275.

| | TEST | DEV-CV (folds 0–4) | DEV-confirm |
|---|---|---|---|
| conditions | 468 | 1358 | 111 |
| GT conditions | 402 | 548 | 97 |
| M1-sel pairs | 470 | 877 | 119 |

TEST has 67 M1-sel kinases; T-unseen 15, T-family-unseen 7. `cache/dev` and `GT_REPORT.md` were regenerated, and
GT_REPORT §9 items 9 and 11 describe the new rules.

## 5. Residual uncertainties (not changed)

1. **Combination series:** these stay out of the metrics until the reference drug of the label dose is documented
   (decryptM supplementary methods).
2. **Remaining Kd_app gates:** kinobead Kd_app is still the gate wherever no cellular IC50 could be documented
   (e.g. nintedanib/FGFR1 Kd_app 1.9 µM, CUDC-101/EGFR). This errs towards masking, not mislabelling.
3. **Unexpressed designated targets:** the Hijazi kinhib matrix is user-provided, all high, and left as given (PI
   decision). Other cases: dasatinib in K562 (LCK, KIT, EPHA2; R121 and D001), crizotinib/lorlatinib ROS1 in
   neuroblastoma, and BIRB796 p38γ/δ.
4. **Dasatinib ABL1 override** uses the Bcr-Abl cellular IC50. In A431/A549, ABL1 is endogenous c-Abl (medium
   confidence).
5. **32034161 doses** ("_0-1" = 0.1 µM, "_1" = 1 µM) are inferred. Fixes 1, 7 and 8 assume these doses.
6. **ERK inhibitor (Calbiochem):** the old docking-site inhibitor in Wilkes (R260) is medium; low may be more
   appropriate.
7. **CDK4/6 inhibitors:** they have no pathway mask, so an indirect CDK2 change is possible. Minor at 1 h.
8. **Inferred-but-medium rules** (R171/R172 SC-514 contrasts, R098 GLP1R HEK293, R060b FGFR2 expression per line)
   still enter M1.
