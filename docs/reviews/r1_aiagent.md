# Review R1 — AI-agent failure modes (reviewer: AI-agent-process expert)

## Verdict
**Major revisions to the process guards. The scientific core is sound.** The split, freeze and Holm-tested TEST protocol are good. Three things are missing. (a) The candidate GT has hygiene and semantic gaps that an agent would fill silently. (b) Nothing stops metric gaming through kinase priors learned from the labels. (c) There is no ledger or budget for DEV iterations. Enforcement is mostly "agents never touch TEST", which is a promise, not a mechanism.

Quick checks (read-only, `results/ground_truth/candidate_gt.tsv`):
- Non-kinase or unmapped GT tokens: `BCR/ABL` (Hijazi Dasatinib, direct), `Q6ZSR9` (UniProt acc, Hijazi LY2835219, 9 rows), `ALL` (decryptM 0 nM null marker, 16 rows) and 22 rows with a blank kinase (e.g. Edelfosine).
- Lipid kinases are GT targets: PIK3CA/B/D/G, PIK3C2B/G, PIK3C3, PI4KB, PIP4K2B/C, PIP5K1A/C and PIKFYVE, each with 3–53 rows. §5 says U = "human protein kinases", so their status is undefined.
- `relation` column: 4,277 direct and 513 downstream rows. MAPK1/3 have 168 rows each and AKT1 has 99, almost all downstream. §3 says downstream kinases are *masked*, but the table carries them as signed GT rows. The design does not say which reading M1 uses.
- A **data-free prior-only score** gives M1 ≈ **0.63** over 306 Hijazi target kinases. The score is −(number of *other* drugs, by name, that list kinase k as a direct inhibitor target). It is computed leave-drug-out with U = the 321 GT kinases and uses no phospho data. The same score gives 0.35 on decryptM and 0.39 on other GT because those include s=+1 stimuli. Any learned per-kinase weight or intercept can harvest this.
- `which R` is empty and `decoupler` is not installed. "Faithful JUMPsem reimplementation" and "decoupler-style" baselines therefore currently have no reference to verify against.

## Proposed changes

1. **GT validation gate before G1** (§3, §9). Priority: **must**.
   - Add `tests/test_gt.py`. It asserts that every GT kinase ∈ U (HGNC protein-kinase list), has no blank tokens, and has no non-gene tokens (BCR/ABL→ABL1, Q6ZSR9→map or drop).
   - The `ALL` null marker moves to a separate null-condition file.
   - Decide and record explicitly whether lipid kinases (PI3K/PIP4K/PIP5K/PI4K/PIKFYVE) are in U. Recommendation: exclude them from U and from the primary GT, and keep them in a sensitivity analysis.
   - The G1 freeze hash covers rule files + validator output.

2. **Define `relation=downstream` semantics** (§3, §5). Priority: **must**.
   - Primary M1/M3/M4 use `relation=direct` pairs only.
   - Downstream kinases are masked from negatives in that condition and scored in a separate secondary metric (M1-downstream).
   - Without this, MAPK1/3 and AKT1 (435 rows) would dominate whichever metric silently includes them. An agent could switch between the two readings toward the better number.

3. **Prior-only negative control and anti-intercept rule** (§5 controls, §7 adoption). Priority: **must**.
   - Add a "kinase-prior-only" baseline: training-fold target frequency, no data. Every learned method must beat it on M1 and report M4.
   - §7 adoption: no learned per-kinase term (intercept, kinase-specific weight, library weight per kinase) is adopted unless M4 also improves (paired CI).
   - Report M1 separately for s=−1 and s=+1 pairs.
   - Evidence: prior-only M1 ≈ 0.63 on Hijazi (above).

4. **DEV iteration ledger and adoption threshold** (§7, §8). Priority: **must**.
   - The harness appends every DEV evaluation to `results/dev_ledger.jsonl`. Each entry holds commit, config hash, GT/split/library hashes, metrics and wall-time. Evaluations that bypass the harness API are invalid.
   - REPORT states the number of configs tried.
   - Replace "improves M1" with: the paired fold-cluster bootstrap 95% CI for ΔM1 has its lower bound > 0, and the sign agrees on ≥2 of {GT-H, GT-D, GT-O}.
   - Carve a **DEV-confirm** set (~15% of DEV groups, seed 0, grouped like TEST). It is not used in iteration and is evaluated only at milestones (≤5 times, logged). A drop of more than the CI from DEV-CV to DEV-confirm triggers a review before TEST.

5. **Physical TEST sealing** (§4, §11). Priority: **must**.
   - Build `cache/dev/` (quant matrix, conditions, GT) with TEST conditions **and TEST-group studies** removed. Explorer and baseline agents receive only these paths.
   - Store the TEST ID list only as salted hashes, with the clear list outside the explorers' briefs.
   - A pytest check fails if any script or config references the TEST file. A second check fails if the co-regulation, site-prior or empirical-null caches were built from inputs other than the `cache/dev` hash.
   - External sets: benchmarKIN Hijazi and 22/24 Hernandez PMIDs *are in the compendium* (`docs/kinasesub_audit.md` §3). Only the tyrosine set (18), the 3 non-compendium Hernandez studies and CPTAC are external. Overlapping benchmarKIN experiments follow their compendium group's split. State this in §2.

6. **Baseline parity and fidelity gates** (§6). Priority: **must**.
   - One shared library builder (PSP mapping, orthologue rule, min-size) is used by all methods. Its library hash goes into the ledger.
   - Baselines get the same DEV-CV budget for their natural knobs (min-size, library choice), or REPORT shows both default and tuned results. A tuned method versus untuned baselines is not a fair "beats".
   - Fidelity: install R + JUMPsem/decoupleR, or the decoupler Python package (try for a time-boxed 1 h, then continue). Require Spearman ≥0.99 between per-kinase scores on a toy set and on one Hijazi cell line.
   - The z-score must reproduce the KSM harness value of 0.782 on benchmarKIN.
   - If no reference can be run, label the baseline "reimplementation, unverified vs reference" in REPORT.

7. **Calibration circularity** (§5, §7B). Priority: **should**.
   - The 27 decryptM 0-vs-0 conditions are the only real nulls. If §7B's empirical null is tuned on them, the calibration metric is circular.
   - Split the nulls by drug, half DEV and half TEST, or use site-permutation nulls only for development.

8. **Number provenance** (§9, §11). Priority: **should**.
   - REPORT and status tables are generated by script from committed result files.
   - `scripts/check_report_numbers.py` matches every numeric claim in REPORT.md and the tree's status summaries to a result file.
   - LITERATURE.md: a script resolves every PMID/DOI (NCBI esummary) and checks the title. Literature anchor numbers used in M3 (0.78–0.79, 0.721) carry table or page locators.

9. **Predicted-KS provenance** (§2). Priority: **should**.
   - Record the sha of R2-nophospho, R3 and Cantley in DATA.md.
   - Confirm from KSM records that no KSM model or variant selection used its downstream Hijazi/compendium evaluation (`KSM/docs/downstream.md` evaluates Hijazi). Disclose the result.
   - Add a per-kinase shuffled predicted-library control.

10. **Bug-driven-win guards** (§9). Priority: **should**.
    - Harness invariance tests for every method: negating the input negates the scores; permuting site IDs gives M1 ≈ 0.5 per GT source; duplicating a condition gives identical scores.
    - Reversed conditions (`expected_direction` +1 for inhibitor-labelled rows) get a unit test.
    - Any adoption with ΔM1 ≥ 0.02 triggers a non-author agent audit. The audit recomputes 3 conditions independently and lists the top-10 kinases driving the gain.

11. **Context-loss guard** (§11). Priority: **should**.
    - `scripts/state_check.py` runs at each agent start and before each ledger write. It verifies the GT, split, harness and library hashes against the FREEZE files.
    - A pytest check fails if `kinact/eval` or the GT rules change after G1 without an §12 amendment ID.

12. **Scope control** (§6). Priority: **should**.
    - Tier the baselines. Must: z-score, KSEA, JUMPsem, ulm/mlm/wmean/norm_wmean/viper, Kinase Library. Could: fgsea, PTM-SEA, RoKAI propagation, INKA.
    - Time-box each tier in the tree.

13. **GT-O rule provenance** (§3). Priority: **could**.
    - Each stimulus or genetic rule cites a PMID.
    - The rule author is different from the domain reviewer, and rules are frozen before any kinase score is computed.

Areas with no change needed: grouping by drug/study via union-find, kinase-macro M1 with M2 matched comparison, the single logged TEST run with independent reproduction, and the seed list.
