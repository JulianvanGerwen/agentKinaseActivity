# Review R2 — data / leakage + agent-process (data specialist, 2026-10-27)

**Verdict: minor-to-moderate revisions.** Most r1 items are addressed. The extra curated KS sources (§2) are not usable as specified. One source is mislabelled, and the drop rule for compendium-derived edges does not catch the main leak. Check script (read-only): `/tmp/claude-588982/-workspace/4cf38946-739b-40b2-907e-e597b4a480f7/scratchpad/ks_src.py`.

## 1. Round-1 must-items

Addressed:
- r1-data 1–7 (downstream masks, U, GT-D rule, grouping incl. controls, split-selection rule, M4 on TEST only, QC).
- r1-aiagent 1–4 (GT validation, downstream semantics, prior-only control, ledger, DEV-confirm).
- r1-aiagent 5, sealing: only partly addressed (see N6). The 6 sub-items of r1-aiagent are adequate.

Not or only partly carried over (r1 shoulds still open):
- r1-data 4c: the canonical drug alias table is not named as a committed file.
- r1-data 9: the Hijazi dose source for the "Kd ≤ 3× dose" mask is not named.
- r1-data 12: T-family-unseen is undefined.
- r1-data 13: M3 is still "validated on DEV/benchmarKIN", and benchmarKIN contains TEST conditions.
- r1-data 15: the audit still says "2 non-compendium Hernandez studies"; the correct number is 3.

## 2. New or remaining problems

**N1 [must] §2 Curated-KS leakage: the "only supporting reference" rule misses the main leak.** HPRD and KEA assigned kinases to sites from large-scale MS papers by motif. Several of those papers are compendium GT studies: 17081983 (EGF), 18669648 and 18691976 (mitosis), 20068231 (cell cycle), and 19651622 (EGF/U0126/dasatinib).
- OmniPath: 1,692 phospho edges cite a compendium PMID, but only 386 cite *only* compendium PMIDs. Most are CDK1 (≈230), CSNK2A1, PRKACA and MAPK1.
- iPTMnet: 1,209 edges cite a compendium PMID (HPRD), but only 21 cite only compendium PMIDs.
- Example: CDK1 substrates "discovered" in the mitosis datasets are then used to score CDK1 in those same GT-O cell-cycle conditions.

Fix:
- Drop every non-PSP/SIGNOR edge whose reference set contains *any* compendium PMID.
- Also drop HPRD/KEA edges whose only references are high-throughput MS papers. My operational rule: a paper cited by ≥150 edges. There are 15 such papers, covering 6,568 OmniPath edges.
- Add a pytest check that the intersection between library edge PMIDs and compendium PMIDs is empty.

**N2 [must] §2 Edges without references.** 14,653 of the 39,037 OmniPath phospho edges (37%) have no references. They come mostly from MIMP, ProtMapper, PSP duplicates and PhosphoNetworks, so the drop rule cannot be applied to them. MIMP is not a curated source. ProtMapper includes text mining.

Fix: whitelist the sources allowed without references (PhosphoSite already in PSP, PhosphoNetworks in vitro, Li2012). Drop edges supported only by MIMP or ProtMapper. Ignore OmniPath rows whose "enzyme" is not in U (for example "EGF" has 353 rows).

**N3 [must] §2 The "iKiP-DB" files are not iKiP-DB.**
- `external/ikipdb_zenodo_{20037245,21469258}` come from the 2026 "phosphoproteome atlas of human cell lines" record. They contain cell-line Cscores (kinase activity outputs, not KS edges) and `Table7_iKiPDB_Modified_Combined_PSP.gmt`. That file is a PTMsigDB-style set that mixes iKiP with PSP; the audit found it holds 61% of PSP, and its provenance cannot be separated.
- LITERATURE cites Mari 2022 (J Proteome Res 21:1575), which these files are not.

Fix: use Sugiyama 2019 MOESM tables directly, which are in vitro, HeLa-lysate and not in the compendium. Obtain the actual Mari 2022 supplement, or drop iKiP-DB. Never use Table7 or the Cscore files as a library.

**N4 [must] §2/DATA.md provenance.** OmniPath, iPTMnet, Sugiyama and iKiP have no DATA.md entry (version, checksum, licence).
- The iPTMnet readme contradicts itself: it says both CC BY-NC-SA 4.0 and CC BY 4.0.
- OmniPath licences are per resource.

Fix: add the entries before G1, and record the iPTMnet licence as NC-SA unless clarified.

**N5 [must] §3 Null controls vs. split.**
- Only 16 decryptM conditions have Condition==Control "0 nM vs 0 nM", not 27. Reconcile the count.
- Each null carries a drug name (e.g. MK2206_0nM), so union-find puts it in that drug's group. A seed-0 halving into N-fit and N-eval therefore puts TEST-group conditions into N-fit, which breaks sealing.

Fix: N-fit = nulls in DEV groups; N-eval = nulls in TEST groups. Report the per-bin calibration as underpowered (n≈8).

**N6 [should] §4/§11 Sealing is a convention, not a mechanism.**
- Agents can still read `data/03_*`.
- Any agent can create `FREEZE.json`.

Fix:
- Add a pytest/grep guard: no file under src/ or scripts/ except the cache builder references `data/03_` or the TEST ID file.
- FREEZE.json must carry the GT, split and config hashes and be written only by the PI. The evaluator verifies the hashes, not mere existence of the file.
- Keep the TEST ID list as salted hashes in explorer-visible paths.
- Exclude DEV-confirm conditions from the explorer cache and from unsupervised fits (co-regulation, noise variances, empirical null) during iteration. Currently only TEST and the held-out fold are excluded.

**N7 [should] §3/§7 GT-A is circular with component D.** GT-A uses a kinase's activating self-sites as the label, and §7-D uses the same sites as input. Fix: GT-A is reported only for methods that do not use self-sites, or D's self-sites are removed when GT-A is evaluated. Pre-register this. Also state that GT-A is evaluated on TEST groups only at TEST.

**N8 [should] §4 Split distance is under-specified.** State the formula before generating candidates (weights, which pairs). I suggest M1-sel pairs, so the targets match the primary metric. Also state the tie and infeasibility fallback, the alias table, and the family definition for T-family-unseen (KSM paralog map or Manning subfamily).

**N9 [could] §2 PSP KS high-throughput edges.** 807 of 14,724 human PSP KS edges sit on sites with LT_LIT=0. Most are in vitro HTP (CHEK1 110, AURKB 78, CDK1 71). Run a sensitivity analysis that excludes them, because PSP edges carry no PMIDs and cannot be checked against compendium studies.

**N10 [could] §2 "top-k measured sites".** Define "measured" per input condition at runtime, or on the DEV cache only. Never define it on the full compendium.
