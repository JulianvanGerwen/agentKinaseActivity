# Review R1: prior work and ground-truth biology (reviewer: prior-work / kinase-signaling)

## Verdict: MAJOR revisions (GT construction). The baseline list needs only minor changes.

Checks I ran (read-only) on `results/ground_truth/*` and `scripts/build_candidate_gt.py`:
- Hijazi rows with a `kinases_downstream` mask: **0 of 183**.
- Stimulus conditions with a downstream mask: 18 of 194. Cell-cycle conditions: 0 of 9.
- 141 of 321 direct GT kinases (44%) appear *only* in conditions with ≥20 annotated targets. For Hijazi alone the figure is 189 of 306. The largest such drugs are PF3758309 (98 targets), LY2835219 (83), AZD1480 (52), Go6976 (46), Dasatinib (44) and AZD5438 (43).
- In Hijazi, Vemurafenib's targets are `MAP2K5;MAP3K20;RAF1;RIPK2`, with no BRAF. Klaeger Kd_app for BRAF is 3.4 µM. Dabrafenib's targets are BRAF/ARAF/RAF1 (−1). All three Hijazi lines are BRAF-WT, and HL60 is NRAS-Q61L.
- Lipid kinases appear as GT positives in 253 rows, across 15 genes (PIK3CA/B/D/G, PIK3C3, PIK3C2B/G, PIP4K2B/C, PIP5K1A/C, PI4KB, PIKFYVE, ADCK1, TRPM6). The non-HGNC tokens `BCR/ABL` (3 rows) and `Q6ZSR9` also occur.
- The script never builds negatives or masks. It implements no "Kd ≤ 3× dose" rule. For decryptM it uses manual target lists plus a confidence downgrade (dose < Kd/10 gives low). That is not the §3 rule ("target iff Kd ≤ d").
- Klaeger has no Kd for MTOR under Torin1 or AZD8055 (kinobeads cover PIKKs poorly). Allosteric drugs (MEKi, MK2206, rapalogs) are also poorly represented by kinobeads Kd_app.

## Proposed changes

1. **§3 Negatives: add downstream and feedback masks for every inhibitor class, Hijazi included (must).** Write a class-level mask table (`pathway_masks.tsv`), keyed by target class, and apply it to every condition:
   - MEKi/ERKi/RAFi: mask MAPK1/3, RPS6KA1–6, MKNK1/2, RPS6KA4/5 and, through feedback relief, RAF1/BRAF/ARAF.
   - PI3Ki/PDPK1i/AKTi: mask AKT1–3, SGK1–3, RPS6KB1/2, MTOR and GSK3A/B. GSK3 is expected to go **up** when AKT is inhibited.
   - mTORi: mask RPS6KB1/2, AKT (mTORC2 or IRS1 feedback) and ULK1 (expected up).
   - p38i: mask MAPKAPK2/3, RPS6KA4/5.
   - CDK7i: mask CDK1/2/9.

   Why: Hijazi currently has zero masks, so these kinases count as negatives (LITERATURE §"polypharmacology, downstream and feedback").
2. **§3: handle paradoxical RAF activation (must).** In BRAF-WT or RAS-mutant cells (Hijazi HL60/MCF7/NTERA2 for Vemurafenib/Dabrafenib, R049, R213, R100 status unknown), drop RAF-family positives and mask MAP2K1/2, MAPK1/3 and RSKs. Treating them as −1 or as negatives has the wrong sign biologically. Also note that Hijazi Vemurafenib has no BRAF target, so its GT is RAF1 only.
3. **§3 GT-H: tier positives (must).** The kinhib list uses >50% inhibition at 200 nM (DiscoverX), but the cells were dosed at 1–10 µM. benchmarKIN's TTD-curated list raised AUROC by about 0.1 over this list (LITERATURE l.365).
   - Define GT-H-primary as the designated primary targets: TTD/benchmarKIN or `drug_targets_manual`, ≤5 per drug.
   - Define GT-H-full as the kinhib list.
   - Report M1 on both. Make GT-H-primary the M1-primary set, or weight each pair by 1/(n targets in c).

   Otherwise 44% of macro-averaged GT kinases rest only on weak off-target calls from promiscuous drugs. Harmony bins would then mostly measure annotation noise.
4. **§3: off-target masking for drugs without a usable Klaeger Kd (must).** 32 of 61 Hijazi drugs and most of the literature drugs (SP600125, U0126, KN-93, Y-27632, H-1152, BIM1, SB203580, BI-D1870, Torin1, wortmannin, LY294002, 5Z-7-oxozeaenol, BI2536, VX680) have no Kd. Add a pre-registered secondary source (Davis 2011 Kd, Anastassiadis 2011, Bain 2007, ChEMBL) plus the per-drug doses from Hijazi's SI. Literature off-targets to cover at minimum: Torin1 → PRKDC/ATM/ATR, LY294002 → CSNK2/PIM1/MTOR/PRKDC, Go6976 → PRKD1–3/CHEK1, Y-27632 → PKN/PRKX, BI-D1870 → PLK1/AURKB.
5. **§3: reconcile the GT-D rule with the code (must).** Either implement "Kd ≤ d" or document the manual-list-plus-downgrade rule that the code uses. Do not apply the Kd rule to allosteric inhibitors (MEK, MK2206, rapalogs) or to PIKKs. Use literature cellular IC50 for those.
6. **§3/§5 universe: exclude lipid kinases and non-HGNC tokens from positives (must).** PI3K-family effects are captured by the AKT/SGK/S6K downstream rows. Map `BCR/ABL` to ABL1 and drop Q6ZSR9.
7. **§3: masks for stimulus and cell-cycle conditions (should).** Mask the canonical cascade:
   - EGF/insulin: MAP2K/MAPK1/3/RSK/AKT/SGK/MTOR/S6K; GSK3 expected −1.
   - TNF: MAP2K3/4/6/7, MAPK8–14, MAPKAPK2, TBK1/IKBKE/CHUK, RIPK1.
   - PMA: PRKD1–3, MAPK1/3, RSK.
   - Mitosis vs interphase: CDK2, NEK2/6/7/9, TTK, BUB1, HASPIN, MASTL, WEE1, CHEK1.

   Lower LCK (TCR) and TGFBR2 (TGF-β) to medium, because both are constitutively active. Change AKT in R177/R209/R247 to relation=downstream.
8. **§3: exclude Edelfosine and U73122 from negative-only use (should).** PLC inhibition lowers PKC activity. Also state how the `negative_control` rows with direction 0 (R140, R212, R231) are used.
9. **§3/§5: state that `downstream` pairs are a secondary positive set (should).** Report them separately (M1-downstream), keep them out of M1-primary, and do not use them as negatives.
10. **§5: paralog fairness (should).** Treat a negative paralog that shares ≥80% of its PSP substrates with a positive (MAPK1/3, PRKACA/B, RPS6KA1/3) as masked, or add a family-collapsed M1.
11. **§6 baselines (should).** Add KSTAR, which supplied benchmarKIN's tyrosine set. Add KARP and PhosR kinaseActivity. Add IKAP, the comparator in the JUMPsem paper. For M3, reproduce using benchmarKIN's own target lists and filters, not our GT. Cross-check z-score/ulm/viper against decoupler-py, and JUMPsem against the original R package, on 3 conditions.
12. **§2 (could).** Audit the PSP KS dump for edges curated from compendium PMIDs (31959955, 36926954, 26060313). Remove them, or flag them for a sensitivity analysis.

No changes needed: leakage bans on predicted libraries (§2), the reversed-condition encoding (R112/R123/R134/R232/R247 checked: signs are correct), the null controls, and the use of KSEA/z-score definitions.
