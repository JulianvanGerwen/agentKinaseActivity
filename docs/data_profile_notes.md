## Findings and caveats (hand-written; numbers come from the generated sections above)

**Site identifiers and mapping**
- `peptide_UID` = `{UniProt acc}_{position}_M{multiplicity}`. All 189,962 sites are `_M1` (singly phosphorylated site, the MaxQuant-style multiplicity) and none carry an isoform suffix. The canonical key used everywhere is `{acc}_{residue}{pos}`, with the residue taken from the centre of the 15-mer window. It is unique.
- 54.8% of compendium sites have their exact key in the PSP human site catalogue. For those, the PSP and compendium 15-mers agree in 99.96% of cases, so positions are on the same (canonical UniProt) coordinates. Another 10,815 compendium windows occur in PSP under a different key (isoform or paralog). They are not used.
- PSP KS (human kinase) has 11,255 unique human substrate sites. 5,540 map to the compendium by exact key, with 100% window agreement. 316 more map through the window on the same base accession; these are PSP isoform accessions such as `P18433-2`, i.e. position shifts. 38 match only on another accession and are not used. Total: 5,856 (52%) are quantified somewhere in the compendium. 2,560 KS rows come from mouse/rat/other substrates mapped to human through `SITE_GRP_ID`.
- Kinase names are normalised to current HGNC symbols. Note that HGNC 2026 renames STK19→WHR1, PRPF4B→PRP4K and AATK→LMTK1. The Hijazi matrix uses old symbols (ADRBK1, ZAK, GSG2, KIAA0999, MST4, PAK7, CDC2L2, SgK110→SBK3, …), which are also normalised.

**Values / normalisation**
- The values look like centred log2 ratios. The per-condition median is close to 0 (82% have |median| < 0.1) and the overall median is 0.000. Hijazi conditions have SD about 0.9 and about 8.8k sites each. decryptM conditions have SD about 0.36 (dose series, mostly small effects).
- **Data problems:** 20 conditions (all qPTM) have SD > 3, and 26 have max |log2FC| > 20 (max 197.8; HSV-1 28179408, 25247763 PIK3CA knock-ins, 26787912, 30869898 …). These look unlogged or wrongly scaled and should be excluded or winsorised. Four GT conditions are affected (25247763_1/3, 27732954_3, 30652415_1). There are also 35 conditions with |median| > 0.5 and 59 with fewer than 2,500 sites. Flags are in `conditions_profile.tsv` (`qc_flag`).
- Metadata quirks:
  - `Timepoint` is parsed wrongly in places ("PMA 0.5h" gives `5h`; 37783352_01 has Sample and Timepoint swapped).
  - There are 88 duplicated condition descriptors within studies (replicate series, e.g. decryptM Dasatinib A549 ×2, PD325901).
  - Some conditions are "reversed": the control is the perturbed state (e.g. `DMSO` vs `Ceritinib`, `Ctr` vs `mitosis`, `Starve` vs `Insulin+aa`). The GT encodes these with expected_direction +1 or −1.
- decryptM (36926954, 790 conditions = 41% of the compendium) includes 27 "0 nM vs 0 nM" conditions. These are useful as **null controls**.

**Overlap / leakage**
- Wilkes 2015 (26060313, MCF7, Cutillas lab, via Ochoa2016) uses many of the same drugs as Hijazi MCF7. The values are not identical (r ≈ 0.02–0.30; highest for mTOR/PI3K inhibitors), so these are independent experiments. Still, they share lab, cell line and drugs, so group them with Hijazi when splitting train and test.
- 32034161 (RPE1) and decryptM re-use classic inhibitors (MEKi, mTORi, AKTi), which gives extra replication across labs.

**Coverage bias**
- 430 PSP human kinases have at least one substrate. 383 have at least one substrate site quantified anywhere, 221 have at least 5, and 82 have at least 20.
- Inside Hijazi conditions, the median per-condition count is ≥3 sites for 115 kinases and ≥10 for only 45.
- The top-10 kinases (CDK1, PRKACA, CDK2, CSNK2A1, MAPK1/3, AKT1, PRKCA, SRC, GSK3B) hold 37% of all kinase–site annotations that are quantified.
- Of the 306 Hijazi target kinases, 38 have no PSP substrates at all and only 98 have a median of at least 3 quantified PSP substrates per condition. A PSP-based benchmark on Hijazi therefore evaluates about a third of the annotated targets, and those are biased toward well-studied kinases.
