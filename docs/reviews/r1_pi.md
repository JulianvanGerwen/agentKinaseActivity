# Review R1 — external PI (strategy)

Reviewer: external PI agent, 2026-10-27. Read: PROMPT.md, DESIGN.md v0, LITERATURE.md, docs/data_profile*.md, docs/kinasesub_audit.md, results/ground_truth/*, DATA.md. Quick read-only checks noted as [chk].

## 1. Verdict
**Major revisions.** The skeleton (frozen GT, grouped splits, sealed TEST, coverage-aware kinase-macro metric, library-first priority) is right. Three strategic problems: (a) the primary metric as defined is dominated by the noisiest labels; (b) the library plan ignores the leakage-safe non-PSP KS sources already on disk, which LITERATURE ranks as the #1 lever; (c) the "explore relationships between kinases" use in PROMPT.md has no deliverable or evaluation. The baseline list is also longer than it needs to be.

## 2. Proposed changes

1. **§5 M1: stop the metric from being dominated by labels from promiscuous drugs.** (must)
   [chk] Hijazi kinhib targets per condition: median 4, mean 12.5, max 98. **214 of the 306 Hijazi GT kinases appear only in conditions with ≥10 targets.** Those labels come from DiscoverX >50% inhibition at 200 nM in vitro, or from literature. Kinase-macro averaging gives each of these kinases the same weight as AKT1 in an MK-2206 condition. So M1 mostly measures off-target labels of dasatinib/staurosporine-class drugs, which are the least reliable and the most confounded (when everything drops, any kinase "wins").
   Change: pre-register M1 over pairs from conditions with ≤10 direct targets (or Klaeger/kinhib selectivity filter), with pair weight 1/n_targets(c) inside the per-kinase average. Report "M1-promiscuous" (>10 targets) separately as a secondary metric. Freeze this before G1.
2. **§2/§7A: add the curated and in vitro KS sources already in `kinaseSubMain/external/`.** (must)
   §7A only lists PSP plus predicted edges. LITERATURE §0.1 and §5.1 found that the curated union (PSP+SIGNOR/OmniPath+GPS) was the best library in benchmarKIN, and that in vitro screens are leakage-safe.
   On disk: `omnipath`, `iptmnet`, `sugiyama2019`, `ikipdb_zenodo_*`, `kinase_library`, `phosformer`. Add them as edge-evidence tiers, and register each in DATA.md.
   Leakage rule: drop any edge whose only evidence PMID is a compendium PMID (OmniPath and iPTMnet carry references; this matters for dbPTM/HTP-derived edges). Never include PTMsigDB PERT sets or `selphi2`.
3. **§2/§3: broaden the inhibitor target profiles beyond Klaeger.** (should)
   `davis2011`, `reinecke2024` and `kinoplex_*` are on disk. Klaeger covers only 29/61 Hijazi drugs. These profiles would extend dose-aware positives and off-target masking to Hijazi, RPE1 (32034161) and organoid (38047086) drugs. Use a fixed precedence rule, frozen at G1.
4. **§3: add a secondary GT based on activating sites (GT-A), for coverage beyond inhibitor targets.** (should)
   [chk] Of the PSP "enzymatic activity, induced" sites on kinases, 184 kinases have such a site quantified in the compendium, and 146 in ≥20 conditions. 121 of these are GT kinases, and 27 have <5 quantified PSP substrates.
   Use these sites as a kinase-centric truth in **non-inhibitor** conditions only: the per-kinase correlation between the inferred score and the A-loop site log2FC across conditions. The kinase's own sites must be removed from every library for this evaluation, which makes it non-circular. Report it separately; never fold it into M1. This also tests the understudied bins that GT-H barely reaches.
5. **§1/§5/§8: make output usability and kinase relationships a deliverable.** (should)
   PROMPT.md says the output "will be used to explore relationships between kinases".
   - Ship a `kinact` CLI/Python API. Input: UniProt_site, gene_site, or a 15-mer (mapped); multi-site peptides split. Output per kinase: score, p, BH-FDR, n_eff, top contributing sites, and a family-level score.
   - Ship a precomputed compendium activity atlas (non-TEST conditions labelled).
   - Secondary evaluation: AUROC with which atlas-wide kinase–kinase score correlation recovers known kinase→kinase edges (PSP KS where the substrate is a kinase activating site; MAP2K→MAPK etc.), compared with z-score.
6. **§6: cut the baseline list to the ones that matter.** (should)
   - Keep: z-score (RoKAI form), KSEA, JUMPsem, ulm, mlm, norm_wmean, VIPER, fgsea, RoKAI+network, KL/PhosX motif-only, and the random-library control.
   - Cut: PTM-SEA/ssGSEA (≈fgsea, slow), the separate Wilcoxon and KS tests (dominated in HA2017), and INKA-lite.
   - [chk] R is not installed (`which R` returns nothing). JUMPsem: first try a ≤30 min conda install of R+lavaan and run the original package. Otherwise reimplement it in Python (semopy/lavaan-equivalent CFA) and validate on a toy matrix with a known latent factor. Document either way.
7. **§3/§5: add specificity controls from non-kinase drugs.** (should)
   Hijazi conditions with 0 kinase targets (e.g. Edelfosine) and the decryptM HDAC/curcumin/non-kinase conditions (`negative_control`, 16) should go next to the 0-vs-0 nulls in calibration: no kinase |z|>3 beyond 5%. These conditions are already in the GT table.
8. **§7 adoption rule: require the gain to hold beyond the T-seen tier.** (should)
   With only ~60 Hijazi drugs plus ~34 decryptM drugs, learned components (§7C) will mostly fit popular kinases. Adopt a learned component only if its DEV-CV gain holds on fold-level kinase-unseen pairs (the T-unseen analogue inside CV), with the paired-bootstrap lower bound >0.
9. **§3 GT-O: separate direct from downstream positives.** (could)
   513 of 4,790 GT rows are `downstream`. Report M1 for direct-only (primary) and direct+downstream, because downstream labels depend on time point and cell line.
10. **§4: check that union-find yields no giant group before freezing the split.** (could)
    Combination arms (e.g. Gefitinib+AZD4547, Selumetinib+MK2206) and stimulus co-treatments (EGF in every RPE1 arm) can chain drugs together. Log the size of the largest group. If it holds >25% of conditions, treat co-stimuli as non-linking.

Areas needing no change: §4 leakage taxonomy and the bans in §2 (consistent with kinasesub_audit.md), §8 TEST protocol, §9 verification, §10 compute.
