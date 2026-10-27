# Review R1 — data / leakage (reviewer: data specialist agent, 2026-10-27)

**Verdict: major revisions** (GT negatives are wrong for a large part of Hijazi, the split as specified is under-defined and very seed-sensitive, and the GT-D rule in DESIGN does not match the code). Check scripts: `/tmp/claude-588982/-workspace/4cf38946-739b-40b2-907e-e597b4a480f7/scratchpad/uf.py` and `/tmp/claude-588982/-workspace/4cf38946-739b-40b2-907e-e597b4a480f7/scratchpad/sim.py` (read-only; union-find + 20 simulated splits over `candidate_gt.tsv`).

## Concrete changes

1. **[must] §3 Negatives — add downstream masks to Hijazi.** All 2,286 Hijazi GT rows are `direct`; no `kinases_downstream`. So MAPK1/3 and RSKs are *negatives* in AZD6244/Trametinib conditions, AKT/S6K are negatives under GDC0941/CAL101/TGX221/HS173/PIK294/Torin1/PF4708671, etc. A correct method gets penalised. Apply the same downstream rule set used for GT-O (MEKi→MAPK1/3, RPS6KA1-3; ERKi→RPS6KA*; PI3Ki→AKT1-3, (mTORC1, RPS6KB1); AKTi/mTORi→RPS6KB1…) as masks for Hijazi and Wilkes, built before G1.

2. **[must] §3/§5 Universe and positives.** 17 direct GT targets are not protein kinases or not symbols: PIK3CA/B/D/G, PIK3C2B/C2G/C3, PIKFYVE, PI4KB, PIP4K2B/C, PIP5K1A/C, Q6ZSR9 (raw accession in the kinhib matrix; target of AZD6244, GDC0994, MK2206), PRKY, BCR/ABL. CAL101, TGX221, PIK294 (9 conditions) then have zero protein-kinase positives. Define U = KSM `kinase_list.tsv` include==True (503) explicitly. Drop GT positives outside U, and report how many were dropped. Keep the conditions as negative-only/specificity controls, together with Edelfosine and U73122 (6 conditions).

3. **[must] §3 GT-D: the rule in DESIGN ≠ the code.** DESIGN says "k is target iff Kd_app ≤ d". `build_candidate_gt.py` instead uses the manual target list (≤8 kinases/drug, 30 kinases total) and *downgrades confidence* by dose/Kd. 485/1,086 decryptM direct pairs are `low`. Pick one rule, write it in §3, and state which confidence tiers enter primary M1 for GT-D (GT-O states this; GT-D does not).

4. **[must] §4 Grouping must cover all 1,937 conditions, not only GT ones.**
   (a) Link drugs found in the *control* (e.g. `ctrl:Gefitinib`, `ctrl:AZD6244`).
   (b) Link the 29 non-GT drug conditions (Tofacitinib, Sorafenib, Erlotinib, Lapatinib, VE-821, PD184352… per KSM `f6b_conditions.tsv`). Otherwise they leak into co-regulation and site priors.
   (c) Use a canonical alias table (AZD6244=selumetinib, GDC0941=pictilisib, GF109203X=BIM-I …).
   In my union-find, 21 Hijazi drugs link to other studies. This gives 208 groups. Groups mix sources: the largest has 131 conditions over 8 PMIDs (gefitinib/lapatinib/AZD4547/vemurafenib), and Hijazi groups are pulled into decryptM groups. "Stratify by source H/D/O" therefore needs a rule (the group's majority source).

5. **[must] §4 Pre-register a split-selection rule; seed 0 alone is not enough.** Over 20 seeds with the stated stratified-30% rule, TEST varies a lot: 305–418 conditions, 147–276 GT kinases, 13–96 T-unseen kinases, 37–59 drugs. One PF3758309 group adds 98 kinases. Generate seeds 0–99. Pick the split closest to label-only targets (30%±3% conditions per source, ≥30 T-unseen kinases, ≥40 TEST kinases with ≥3 PSP substrates). Commit it before any method is run.

6. **[must] §5 M4 on TEST must use TEST conditions only.** The current wording is "DEV/TEST conditions".

7. **[must] §2 QC exclusion rule is too broad.** `qc_flag` removes 40 GT conditions: 28 for n_sites<2500 and 14 for |median|>0.5. Only 3–4 are scale errors (sd>3 or max>20: 25247763_1/3, 27732954_3, 30652415_1). Exclude only scale errors. Keep low-coverage and median-shift conditions, with a sensitivity analysis.

8. **[should] §3 Promiscuous-drug dominance.** 189/306 Hijazi GT kinases come only from the 13 drugs with ≥20 targets. 108 kinases come from a single drug. Of the 100 PSP-evaluable Hijazi kinases, 62 come only from drugs with ≥10 targets. Pre-register a secondary M1 on selective drugs (≤10 targets) and per-drug weighting, so the macro score is not driven by DiscoverX off-target calls.

9. **[should] §3 Off-target masking spec.** Hijazi dose is not in the metadata, so name the source (Hijazi Supp. Table). 32/61 drugs have no Klaeger profile, so state the fallback. The Klaeger parser keeps only the first gene of protein groups (14 rows, incl. PRKCE;PRKCH, PRKX;PRKY, PIP4K2A;B); expand them for masking.

10. **[should] §2 Data errors.** Two pairs of conditions hold byte-identical vectors:
    - 35688186 caffeine_k ≡ quercetin_k (6 pairs).
    - 22115753_qPTM_1 ≡ _2.

    Exclude them from every use. Add a duplicate-vector test to `tests/`. No duplicates were found across PMIDs.

11. **[should] §4 Co-regulation exclusion.** "Never TEST-group studies" leaves only 547–639/1,937 conditions. The shared-control effect is negligible: Hijazi conditions with different drugs and the same DMSO control have median r=0.013, vs 0.001 across cell lines. Instead, exclude TEST-group conditions plus any condition sharing (PMID, Control, Sample) with one (≈1,300 conditions kept). Apply the same rule inside the CV folds; it is currently asymmetric.

12. **[should] §4 Kinase-novelty tiers.** Add T-unseen-family, using the KSM paralog map (≥60% KD identity or same subfamily). For example, ROCK inhibitors differ by drug between Wilkes and Hijazi (Y-27632/H-1152 vs Ripasudil/KD025) but share targets.

13. **[should] §5 M3 harness validation.** It uses benchmarKIN Hijazi (118) and Hernandez (22 PMIDs in the compendium), which include TEST conditions. Validate the harness only on DEV conditions, or do it before the split with locked PSP+z-score, and log it.

14. **[should] §2/§7 Do not reuse KSM downstream-tuned settings.** This covers fusion 0.6/0.4, `cache/downstream/*/measured_topk`, and KL top-15; KSM chose them on benchmarKIN, which overlaps our TEST. Re-select them in DEV CV. The R2-nophospho/R3 bans are otherwise correct.

15. **[should] §2 External sets are small; report them by stratum.** The independent sets are 15 Hernandez experiments from 3 PMIDs (not "2" as the audit says) and 18 tyrosine experiments. Together they cover about 25 kinases. Their drugs overlap DEV (dasatinib, imatinib, bosutinib, EGF), so report drug-seen and drug-unseen separately. Fix the PMID count in the audit.

16. [could] §2 Join predicted libraries by 15-mer. About 10.7k compendium sites are not in the KSM site index by exact key. Disclose PSP circularity: KS rows may be curated from compendium studies and carry no PMIDs.

17. [could] GT completeness: the 19651622 SB202190 conditions lack MAPK14/11. Send this to the domain review.

No change needed: site-key and isoform mapping (99.96% window agreement), Wilkes↔Hijazi grouping (works through drug identity), and the bans on predicted libraries (L3).
