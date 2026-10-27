# Review R2: prior work and ground-truth biology (reviewer: prior-work / kinase signaling)

Verdict: **minor revisions.** Most round-1 items are addressed in the design text. Two gaps would still corrupt labels: the dose source for off-target masks, and leakage in GT-A. The rest is specification detail that has to be fixed before G1.

Checks I ran (read-only):
- U has 503 genes and includes MTOR/ATM/ATR/PRKDC/SMG1/TRRAP, PDK1–4, BCKDK, TRPM6 and ADCK1. Direct GT tokens outside U: the PI3K/PI4K/PIP4K/PIP5K genes, BCR/ABL, PRKY and Q6ZSR9 (16 in total). This is consistent with A5.
- M1-sel coverage:
  - Hijazi: 114 of 183 conditions and 38 of 55 drugs. Only 81 of 290 Hijazi GT kinases remain.
  - All sources (high/medium): 702 of 753 conditions, 140 of 304 GT kinases.
- `data/01` has no dose for Hijazi (or for most non-decryptM studies). The doses are only in the paper/SI.
- `results/ground_truth/*` is still the pre-v1 candidate build: it has no masks or negatives and still contains lipid kinases. That is expected before G1 and is not a design flaw.

## 1. Round-1 must items

| r1 # | Status |
|---|---|
| 1 Downstream/feedback masks incl. Hijazi | Addressed (§3 (i)). Small gaps: CDK7i→CDK1/2/9 and p38i→MSK1/2 (RPS6KA4/5) are missing from the table. |
| 2 Paradoxical RAF | Addressed for Hijazi. Its scope beyond Hijazi is underspecified (see N3). |
| 3 Tier GT-H positives | Addressed via M1-sel (≤10 targets) + M1-all 1/n. Accepted. |
| 4 Off-target masking for drugs without Klaeger | **Partly.** The profile precedence is added, but (a) Bain 2007/Anastassiadis/ChEMBL and the named literature off-targets (Torin1→PIKKs, LY294002→CK2/PIM1, Go6976→PKD/CHEK1, Y-27632→PKN/PRKX, BI-D1870→PLK1/AURKB) are not included, and (b) the "Kd ≤ 3× dose" rule has no dose source for Hijazi or literature studies (N1). |
| 5 GT-D rule reconciled | Addressed. |
| 6 Lipid kinases / bad tokens | Addressed. |

Should items from r1:
- Not addressed: #7 (LCK/TGFBR2 downgrade, AKT relation in R177/R209/R247) and #12 (PSP compendium-PMID audit). #12 can be dropped: the PSP KS dump carries no PMIDs, so it is not auditable. State this in DATA.md.
- Addressed: #8 to #11.

## 2. New or remaining problems

**N1 (must) §3(iii) / §2 dose source.** The off-target mask and every dose-aware rule need a per-condition dose.
- Fix: add `condition_doses.tsv`, built from the Hijazi SI (per drug) and from the paper methods for GT-O, with provenance per row.
- If the dose is unknown, use a pre-registered default of 10 µM, which is conservative and gives more masking.
- Add the literature off-target table from r1 #4 as an extra mask source for drugs that have no profile.
- Test before G1: every inhibitor condition has a dose or is flagged as having the default dose.

**N2 (must) §3 GT-A leakage and definition.** Removing the self-sites from "its library" only is not enough. Method D uses them directly, and co-regulation/factorisation (§7C) can carry them into other kinases' scores.
- Fix: for GT-A, delete the kinase's PSP activating sites from the **input vector** before every method runs, including D, the co-regulation features and the baselines.
- Define the label explicitly:
  - per (condition, kinase), s = sign of the mean log2FC of its activating sites when |mean| ≥ 1;
  - otherwise the kinase is a negative only when |mean| < 0.3;
  - the band between the two is masked;
  - pairs are restricted to conditions ≤ 4 h, to limit protein-abundance confounding.
- Alternatively, freeze the per-kinase correlation that the r1 PI proposed. Pick one before G1.

**N3 (should) §3(ii) paradox scope.** "BRAF-wt cells (all Hijazi lines)" does not cover RAF inhibitors in decryptM or in other studies.
- Fix: add a cell-line genotype column (BRAF V600 / RAS / WT / unknown) and a drug-type column (paradox-inducing type I/I½: vemurafenib, dabrafenib, encorafenib, PLX4720, SB590885 vs pan-RAF/type II: LY3009120, belvarafenib, AZ628, TAK632).
- Apply the paradox rule to type I/I½ drugs in WT, RAS-mutant or unknown lines.
- In BRAF-V600 lines, keep BRAF −1 and MEK/ERK as downstream −1.

**N4 (should) §3 label precedence.** The design gives no precedence for conflicts. For example, the RAFi mask lists RAFs while RAFs are also the annotated targets, and a paralog mask can hit a co-target.
- Fix: pre-register the order paradox-drop > direct positive > downstream/off-target/paralog mask > negative, and add a test.

**N5 (should) §3 GT-O stimuli: add a time window.** A first-tier kinase gets s=+1 only at ≤60 min (ERK/AKT at ≤30 min, unless the study shows sustained activation). At later times it is masked. Feedback and receptor internalisation/degradation (EGFR after 1–2 h) flip or erase the signal. Apply the r1 #7 LCK/TGFBR2→medium change.

**N6 (should) §3 GT-D sub-Kd doses and covalent drugs.**
- At doses below Kd, the annotated target is **masked**, not negative. Kd_app is measured in lysate, so cellular potency can be higher.
- Covalent inhibitors (ibrutinib, afatinib, osimertinib, neratinib) follow the IC50 rule, as allosteric drugs do, because kinobeads Kd_app underestimates time-dependent inhibition.

**N7 (should) §6 baseline fidelity.** Run each must-baseline in two configurations:
- (a) its published default (KSEA: PSP+NetworKIN cutoff m≥5; RoKAI: PSP+SIGNOR with its network; z-score: PSP);
- (b) on the shared-library grid, which must include our final library.

Report both. This separates library gains from method gains and answers the strawman objection. The Holm comparison uses the better of (a) and (b) on DEV.

**N8 (could) §5 M1-sel coverage.** M1-sel covers only 140 of 304 GT kinases (81 for Hijazi), so the harmony bins under M1-sel will be thin. Report the n kinases per bin at G1, and show harmony for M1-all next to it.

**N9 (could) rule files.** Use PDPK1 rather than "PDK1" in the mask tables. PDK1–4 (pyruvate dehydrogenase kinases) and BCKDK are in U, so the bare token would mis-map.

**N10 (could).** State how the r1 #8 `negative_control` direction-0 rows (R140/R212/R231) are used: as specificity controls only.

No changes needed: the leakage bans, the banned-source list, the KS-source list, the null-control split, the M3 protocol, the decoupler/JUMPsem cross-checks, and the could-tier baselines (KSTAR, KARP, PhosR, IKAP).
