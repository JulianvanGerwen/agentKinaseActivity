# Review R3 (final): consistency, implementability, fairness, scope

Reviewer: independent senior reviewer, round 3, 2026-10-27. Read: PROMPT.md, DESIGN.md v2 (§1–§13), docs/interfaces.md, docs/reviews/r2_*.md, docs/kinasesub_audit.md.

Read-only checks:
- 16 conditions with Condition == Control, all in decryptM 36926954. They carry 7 distinct drug-group tokens (including HDACi such as CUDC101 and Romidepsin), and their per-site SD is 0.20–0.39 with no exact zeros, so they are usable as nulls. About 5 of them will fall in TEST, so N-eval will be small.
- In candidate_gt, M1-sel (direct, high/medium, ≤10 targets) contains 1,760 pairs and 151 kinases. Signs: s = +1 for 254 pairs, of which 226 are stimuli; s = 0 for 3 pairs (negative_control rows).
- 513 pairs are "downstream". For example, EGF→MAPK1/3 is downstream, not direct.

**Verdict.** The design is sound and close to freezable. The remaining problems are (a) contradictory text that agents will implement literally, (b) four specification gaps that touch labels and sealing, and (c) some scope that should be trimmed. No change to the scientific plan is needed.

## Must

| # | Section | Change | Why |
|---|---|---|---|
| M1 | Whole doc | Before G1, fold §13 into §1–§11 as **v3**. Keep §12/§13 only as a history log. Stale text to fix:<br>• §3: "27 nulls, seed-0 halves" → B9.<br>• §3 GT-D: "lowest doses → negative-only" → B11: sub-Kd targets are masked, and the condition has 0 positives, so it is a specificity control.<br>• §5: "TEST family of 5" → 6 (B3).<br>• §5: "must also improve M4" → delete, since M4 is descriptive per B7.<br>• §2: "iKiP-DB" → dropped (B8).<br>• §4: "evaluator refuses TEST until FREEZE.json exists" → hash-verified (B10).<br>• §7 adoption rule → B1. | Agents read one section, not the override log. Each stale line above is a concrete wrong instruction. For example, an agent following the M4 clause would block adoptions on a metric that is not frozen. |
| M2 | B12 / §3(ii) | Pre-registered precedence: **paradox-drop > direct > downstream > masked (off-target / paralog / sub-Kd / time window) > negative**. Add a pytest: in BRAF-wt cells, a type I/I½ RAF inhibitor yields no RAF positives. | As written, B12 ("direct > masked") keeps RAF targets as direct positives, which contradicts the paradox rule in §3(ii). That is the rule's exact purpose (r2 priorwork N4). The "masked > downstream" order would also wipe out every downstream pair, because the downstream kinases are the mask-table entries. |
| M3 | §3 GT-O / interfaces | Define stimulus positives explicitly. direct = receptor kinase plus a named first-tier list in `condition_rules_manual.tsv`. Cascade kinases (ERK, AKT, …) = cls downstream: excluded from M1 and scored only in the downstream AUROC. | §3 says first-tier kinases are direct with s = +1, but the mask table masks "canonical cascades for stimuli", and candidate_gt marks EGF→MAPK1/3 as downstream. Nearly all s = +1 evidence (226 of 254 M1-sel pairs) depends on this choice, so it must be fixed before G1. |
| M4 | §5 M3 / §2 external | Before FREEZE, M3 and every benchmarKIN-based check run only on DEV conditions plus the tyrosine set. Do not run them on benchmarKIN's Hijazi or Hernandez experiments that map to TEST/DEV-confirm conditions; use `ldo.py` / `benchmarkin_studies.tsv` for the mapping. | "Validated on DEV/benchmarKIN" is an open path to TEST data, because benchmarKIN contains Hijazi 31959955 and 22 compendium Hernandez studies (r2 data #13, still open). |
| M5 | §7 (new rule) | Ban GT-derived kinase-level features in every method, stacker and kernel: no target frequency, no "was a positive in training", no per-kinase label counts. Learned parameters may use only label-free edge/site/kinase features (library size, motif, family). | The prior-only control (B3) is in the Holm family because a popularity prior already reaches about 0.63. If a method can feed freq_k back in through a stacking feature, then "beats prior-only" becomes trivial and understudied kinases (freq = 0) are penalised. That goes directly against the PROMPT's bias clause. |

## Should

| # | Section | Change | Why |
|---|---|---|---|
| S1 | B1 / §4 DEV-confirm | Spell out the DEV-confirm protocol:<br>• One look at the end of each B2 step that adopted something (≤5).<br>• Pass = point ΔM1-sel ≥ 0 vs the previous passed checkpoint.<br>• Fail = revert to that checkpoint.<br>• The final model is fit on the same data as when it passed (no refit on DEV-confirm), or a refit rule is pre-declared. | B1 says the bundle "is checked" but gives no timing, no criterion and no consequence. There are 6 B2 steps but only 5 looks. |
| S2 | §1 / interfaces | Unscored = NaN, by contract. A kinase with no information in a condition (no measured substrate and no kernel borrowing) must be NaN, never 0. Add a pytest: no kinase may have a constant non-NaN score across conditions. | §1 asks for "0 = no evidence", but the harness gives a NaN target ½ while a 0 target is ranked against shifted negatives. Under broad inhibitor-induced decreases that rank is below ½. The rule must be the same for baselines and our method. |
| S3 | interfaces.md | Extend the contract:<br>• `fit(X_allowed, gt_train) → model; model.score(X)` for learned methods and prior-only.<br>• OOF score storage: one file per config, each fold scored by the model trained without it.<br>• Library columns `refs` (PMIDs) and `n_orthologs`, needed for the B8 test.<br>• `harmony_bins.tsv` (kinase, group_id, n_psp_quantified, bin), hashed in GT_FREEZE.<br>• `gt_conditions.role` gains `excluded` (QC/duplicate).<br>• FREEZE.json fields per B10. | Without these, the harness builder and the method agents will each invent their own fold handling. B6 bins, B8 and B10 currently have no home in the contract. |
| S4 | B6 | Harmony bin per (kinase, group) = the **median** number of quantified PSP substrates over the group's conditions where the kinase is a positive (r2 stats N6). State it. | "Per (kinase, group) using … that condition" is self-contradictory. |
| S5 | B10 / §4 | Physically split the GT files: `gt_*_dev.tsv` is agent-visible, and the TEST and DEV-confirm rows go to a sealed file that only the PI TEST script reads. Agents see only a hash. | gt_pairs and gt_conditions currently list TEST labels in plain text. Knowing which kinases are TEST positives invites tuning toward them. |
| S6 | §6 / §2 JUMPsem | Multi-sample methods (JUMPsem; PhosR if run) are fit per study on the **same condition set as the evaluated split**: the DEV cache during DEV, and DEV ∪ TEST of that study only in the TEST run. | Running them "natively per study" on the full data reads TEST vectors during development (sealing). It would also give JUMPsem transductive information that our unsupervised parts are denied. |
| S7 | §6 / B14 | Trim the must-tier:<br>• Must: z-score, KSEA, JUMPsem, decoupler ulm + mlm + viper (use pip decoupler directly, with no re-implementation or Spearman check), RoKAI, Kinase Library motif.<br>• Move wmean, norm_wmean and fgsea to could.<br>• B14 double configs apply to the must-tier only. | There are about 10 must baselines × 2 configs × a grid, plus re-implementations and agreement checks. That is days of work and adds little; per KSM, z-score ≥ KSEA ≈ ulm > wmean already. |
| S8 | §2 top-k predicted libraries | "Measured" = sites non-NaN in the input condition, chosen at runtime (r2 data N10, not carried into §13). | The CLI must work on user data, and a compendium-wide "measured" set would leak TEST site coverage. |
| S9 | B8 | Make the HT-MS paper list script-only: compendium PMIDs plus papers cited by ≥150 edges (r2 data N1). Drop the PubMed MeSH/title heuristic. | It is deterministic and needs no external API. The PROMPT says not to stall on access. |
| S10 | B7 controls | Add a degree-preserving **shuffled-library** control to the frozen controls. | KSM lesson: a shuffled PSP library scores 0.52 through √n |z| inflation. This is the bias most likely to favour well-studied kinases, and neither site permutation nor prior-only isolates it. |
| S11 | §1 / §10 | Make the speed target a reported desideratum, not a gate. Note that the ~7 ms Bayesian estimate exceeds the 5 ms target and that batching conditions by missingness pattern is the planned fix. | The design is internally inconsistent on this point now. |
| S12 | B5 | State that a decryptM dose enters M1-sel per condition: its n_targets is counted after the Kd rule, so a drug can be in M1-sel at low dose and out at high dose. | This was left undecided (r2 stats N3). |

## Could
- C1, B4: the MDE multiplier should be (z₀.₉₉₂ + z₀.₈) ≈ 3.5, not 2.5. α/6 under Holm makes 2.5 optimistic.
- C2: CPTAC tumour benchmark and GT-A → could-tier, post-G2.
- C3: decision rule for the winner's-curse re-split: descriptive only.
- C4: key libraries internally by UniProt accession (KSM lesson 2).
- C5: report N-eval calibration as underpowered (n ≈ 5).

No other changes needed. The metric definitions, Holm family, controls, split selection, and leakage bans on predicted sources (R2-nophospho / R3 deploy / Cantley) are adequate.
