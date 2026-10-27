# LITERATURE.md: Kinase activity inference from phosphoproteomics

Compiled 2026-10-27 for the kinase-activity-estimator project (see PROMPT.md).
Scope: methods, kinase–substrate (KS) resources, ground-truth datasets, evaluation methodology, design implications.

**Verification tags** (please keep them when editing):
- **[V]**: checked this session against the PDF in `data/` or the publisher, PMC, PubMed or GitHub page.
- **[M]**: from reviewer memory or secondary sources only. Treat the numbers or IDs as unverified until someone checks them.

Where a number could not be confirmed, the text says so. Nothing has been filled in by guesswork.

---

## 0. Executive summary (decision-relevant)

1. **The KS library matters more than the algorithm.**
   - The largest independent benchmark, benchmarKIN (Müller-Dott et al. 2025), compared 19 algorithms × 8 libraries. Most algorithms performed similarly.
   - The **z-score (RoKAI/KSEA-style)** was the best or near-best algorithm on perturbation data: mean scaled rank 0.24, mean AUROC ≈0.79 with PhosphoSitePlus [V].
   - A **combination of curated libraries** (PSP + PTMsigDB + GPS-gold) was the best library: perturbation AUROC 0.73 averaged across methods [V].
   - Adding **NetworKIN** predictions helped on tumour benchmarks. Adding Kinase Library or Phosformer predictions (top-15 kinases per site) **hurt**: AUROC ≈0.65 on perturbation and ≈0.58 on tumour [V].
2. **The 2017 benchmark agrees.**
   - Hernandez-Armenta et al. 2017 found Z-test ≈ GSEA > KS-test ≈ Wilcoxon ≈ ridge MLR. Median AUC was 0.721 for Z and 0.734 for GSEA [V, PDF].
   - Weighting substrates by motif match gave a small but significant gain [V].
   - In-vivo-supported substrates > in vitro > in silico (NetworKIN) [V].
   - Every method was worse for kinases with ≤5 quantified substrates. MLR, KS-test and Wilcoxon degraded most [V].
3. **JUMPsem is not a single-sample method.**
   - It fits a one-factor confirmatory factor analysis per kinase across samples with `lavaan` (`kinase =~ s1 + s2 + …`) and returns `lavPredict` factor scores [V, source code].
   - So it needs a matrix of conditions. Its scores are centred across those samples and are transductive. The sign is anchored to the first indicator's loading (lavaan marker-variable default) [V code; sign consequence inferred].
   - Its paper benchmarks only against IKAP and KSEA, on 31 conditions, by "precision at rank threshold", with no AUROC [V, PDF].
4. **Leakage hot-spots in predicted KS sources.**
   - Hijazi PDTs/EBDT are derived from the *very* Hijazi inhibitor data we evaluate on [V].
   - SELPHI2, CoPheeKSA (CPTAC co-regulation) and PhosR-type co-regulation predictors learn from phosphoproteomics quantifications.
   - Some PTMsigDB perturbation signatures, and possibly PSP high-throughput (HTP) curation, may come from evaluation studies.
   - Leakage-safe predicted sources:
     - in vitro specificity: Kinase Library PSPA matrices, PhosX, Kinex
     - in vitro substrate screens: Sugiyama 2019 / iKiP-DB
     - motif-based predictors: NetPhorest/NetworKIN, GPS with care
5. **The provided compilation already contains two of the best ground-truth sources.**
   - Hijazi 2020: 61 conditions × HL60/MCF7/NTERA2, 1 h.
   - decryptM (Zecha et al. 2023 Science, PMID 36926954): 790 conditions, dose-resolved, 13 cell lines, ~34 drugs. Together with **Klaeger 2017 kinobeads** apparent Kd values (provided in `data/Klaeger/`), it allows *dose-aware* positives: a kinase counts as inhibited only when the dose exceeds its Kd_app.
6. **Biggest openings to beat baselines.** Each is detailed in §5.
   - confidence-weighted union libraries with leakage-safe predicted edges for understudied kinases
   - correcting the √n study-bias of z-scores via per-kinase calibrated nulls or shrinkage
   - multivariate or regularised joint models for overlapping substrate sets
   - learning per-edge reliability weights from perturbation data with inhibitor-grouped splits
   - family-level hierarchical priors

---

## 1. Kinase activity inference methods

Notation:
- x_i: log2FC of site i in one condition
- S_k: the measured substrates of kinase k, with m = |S_k|
- μ, σ: the mean and SD of x over all measured sites in the condition

### 1.1 KSEA: Casado et al. 2013; KSEAapp: Wiredja et al. 2017
- **Refs.**
  - Casado P et al. *Sci Signal* 2013;6(268):rs6, PMID 23532336 [M; the PMID matches a study already in our compilation].
  - Wiredja DD, Koyutürk M, Chance MR. *Bioinformatics* 2017;33(21):3489–3491, PMID 28655153, doi:10.1093/bioinformatics/btx415 [V].
- **Statistic.** One-sample z-test: z_k = (mean_{i∈S_k} x_i − μ) · √m / σ, with a two-sided normal p-value and BH FDR [M; the standard KSEAapp formula].
  - KSEAapp also reports "Enrichment" = mean(S_k)/mean(all) on the FC scale.
  - It has an `m.cutoff` minimum substrate count and an optional NetworKIN source with a score cutoff (commonly 5).
- **Inputs.** Site-level FC table (gene, residue, peptide, FC). Peptides mapping to the same site are averaged.
- **Runtime.** Closed form, milliseconds per condition.
- **Implementation.** R/CRAN `KSEAapp`, version 2.0 released March 2026 [V]; Shiny web app.
- **Performance.**
  - Hernandez-Armenta 2017: the Z-test (Casado implementation) had median AUC 0.721 and precision@recall0.5 of 0.708 [V].
  - benchmarKIN 2025: "KSEA" and "z-score (RoKAI)" are listed separately. The z-score was best overall [V].
- **Known issue.** |z| grows with √m, so kinases with many annotated substrates (well-studied ones) get systematically larger scores. This is study bias, discussed in the RoKAI paper.

### 1.2 z-score and RoKAI: Yılmaz et al. 2021
- **Ref.** Yılmaz S, Ayati M, Schlatzer D, Çiçek AE, Chance MR, Koyutürk M. *Nat Commun* 2021;12:1177, PMID 33608514, doi:10.1038/s41467-021-21211-6 [V].
- **z-score as implemented in RoKAI.** a_k = √m · mean_{S_k}(x) / σ, equivalent to Σ x_i / (σ√m). It is **not mean-centred**, unlike KSEA [V, as read from Eq. 3; minor uncertainty about the rendering].
- **RoKAI propagation.**
  - Builds a heterogeneous functional network with these edge types [V]:
    - KS links: PSP, 3,877 links, 261 kinases
    - kinase–kinase PPI: STRING, 13,031 weighted edges
    - site–site structural distance: PTMcode, 7,821 edges
    - site–site co-evolution: PTMcode, rRCS ≥0.9, 41,799 edges
  - Phosphosite values are "refined" by solving an electric-circuit linear system Cv = b with v = (CᵀC)⁻¹Cᵀb and τ = 1e-8 ridge on the diagonal. Measured sites act as current sources; unmeasured sites and kinases act as bridges [V].
  - The refined site vector is then fed to any activity method (z-score, KSEA, linear model).
- **Linear-model baseline in RoKAI.** Ridge with λ = 0.1, following Hernandez-Armenta [V].
- **Evaluation.**
  - 80 perturbations from the Ochoa 2016 atlas / Hernandez-Armenta set: 53,636 sites; 128 kinase–perturbation annotations covering 25 kinases [V].
  - Metric: P_hit(k = 10), i.e. the share of perturbed kinases ranked in the top 10.
  - Robustness was tested by **hiding a fraction of known substrates** (100 random runs) [V].
  - RoKAI improved all base methods significantly (p < 0.05). Exact numbers appear only in figures [V].
- **Important observation.** Annotated (benchmark) kinases have a median of 29 substrates versus 4 for non-annotated kinases [V]. Benchmarks therefore reward methods biased toward well-studied kinases.
- **Implementation.** MATLAB (github.com/serhan-yilmaz/RoKAI), web app rokai.io. benchmarKIN re-implements the z-score in R [V].
- **Runtime.** Network solve is a sparse linear system; seconds.

### 1.3 Hernandez-Armenta et al. 2017: benchmark and methods compared (Paper_KSEA.pdf)
- **Ref.** Hernandez-Armenta C, Ochoa D, Gonçalves E, Saez-Rodriguez J, Beltrao P. *Bioinformatics* 2017;33(12):1845–1851, doi:10.1093/bioinformatics/btx082 [V]. PMID not checked (believed to be 28200105 [M]).
- **Methods compared** [V]:
  - Z-test (Casado)
  - Kolmogorov–Smirnov test
  - Wilcoxon rank-sum (substrates vs non-substrates)
  - GSEA: weighted KS with 10,000 random same-size sets, `evocellnet/ksea` R package
  - MLR / ridge: Y = βX + ε, with X the binary KS matrix and β the activities, ridge α = 0.1 in scikit-learn
- **Scoring.** For the tests, score = −log10(p) signed by the mean substrate FC. For MLR, β is used [V].
- **Weighted variants.** Fold changes are multiplied by the motif similarity score (MSS, 0–1) from substrate-derived PWMs (±7 residues; MATCH-like). Requires ≥10 known substrates [V].
- **Library.** 2,818 curated interactions for 30 kinases from PSP, Phospho.ELM and HPRD. **150 auto-regulatory sites were excluded** [V].
- **Results** [V]:
  - Median AUC: GSEA 0.734, Z 0.721.
  - Precision@recall0.5: 0.725 and 0.708 respectively.
  - Weighted Z and weighted GSEA were significantly better than unweighted (Wilcoxon p = 8.4e-10 and 1.26e-6), but the gains were small.
  - For ≤5 quantified substrates, Wilcoxon, KS and MLR were significantly worse than GSEA/Z (p < 2.2e-16).
  - In vivo > in vitro > in silico (NetworKIN 3.0) substrates, and the ordering held after down-sampling. Stringent NetworKIN cutoffs help mainly because they re-enrich known substrates.
  - Degree-preserving randomisation of the KS network gave near-random performance, which is a good negative control to copy.

### 1.4 GSEA / fgsea / KS-type enrichment
- GSEA: Subramanian et al. *PNAS* 2005 [M].
- fgsea: Korotkevich et al., bioRxiv 060012 [M]. It uses adaptive multilevel permutations to give fast and accurate small p-values.
- Statistic: weighted running sum over sites ranked by FC, normalised ES against a permutation null.
- Strong in HA2017 (§1.3). In benchmarKIN, fgsea was among the methods with similar performance [V].
- Runtime: permutation-based; fgsea is fast (about seconds per condition for hundreds of sets).
- Ochoa 2016 used a GSEA-type approach to build its 215-kinase × 399-condition atlas [V].

### 1.5 PTM-SEA / PTMsigDB: Krug et al. 2019
- **Ref.** Krug K et al. *Mol Cell Proteomics* 2019;18(3):576–593, PMID 30563849, doi:10.1074/mcp.TIR118.000943 [V].
- **Statistic.** ssGSEA2 adapted to **site-specific, directional** signatures, with ES normalised by permutations. Signatures carry up/down direction per site [V].
- **Signatures.** Curated from more than 2,500 publications [V]. Categories:
  - KINASE: from PSP
  - PATH
  - PERT: perturbation signatures, including LINCS P100 and literature
  - DISEASE
  - iKiP-derived sets (excluded in benchmarKIN)
- **Leakage note.** PERT signatures come from phosphoproteomics experiments. Check whether any come from our evaluation studies (e.g. Hijazi, Wilkes 2015) before using them.
- **Implementation.** R/GCT, github.com/broadinstitute/ssGSEA2.0. Also in benchmarKIN.
- **Performance.** benchmarKIN tumour benchmark with curated + NetworKIN library: protein-based AUROC 0.681, activating-site AUROC 0.673 [V].

### 1.6 KinSwingR: Engholm-Keller et al. 2019
- **Ref.** Engholm-Keller K, Waardenberg AJ et al. *PLoS Biol* 2019;17(3):e3000170 [M]. Bioconductor `KinSwingR`, GPL-3 [V].
- **Algorithm** [V]:
  - Builds PWMs from known substrates (`buildPWM`) and scores all measured peptides (`scoreSequences`), which gives a probabilistic kinase × site network.
  - The **swing score** combines, per kinase, the counts of positively and negatively regulated significant sites among predicted substrates, weighted by substrate number and local network size.
  - Significance comes from permutation of the network.
- **Inputs.** Peptide sequence windows, FC and p-value (requires significance calls).
- **Performance.** No independent large benchmark found. Not in benchmarKIN.

### 1.7 KARP: Wilkes et al. 2017
- **Ref.** Wilkes EH, Casado P, Rajeeve V, Cutillas PR. *Mol Cell Proteomics* 2017 (PMC5587867) [V title; volume and pages M]. PRIDE PXD005073 [V].
- **Statistic.** Single-sample "K-score": the sum of substrate intensities divided by the sum of all site intensities, with a correction for the number of known substrates [V qualitatively; exact formula not verified].
- **Inputs.** Intensities, not FC.
- **Performance.**
  - benchmarKIN: KARP and the "number of targets" baseline were the only methods with average PHit(k) below 0.43; KARP was significantly worse on the tumour benchmarks [V].
  - KSTAR paper: KARP was second-best by rank on tyrosine data [V].

### 1.8 IKAP: Mischnik et al. 2016
- **Ref.** Mischnik M et al. *Bioinformatics* 2016;32(3):424–431, PMID 26628587, doi:10.1093/bioinformatics/btv699 [V].
- **Model.** Each phosphosite's level is the sum of kinase effects (activity × site-specific affinity). Activities and affinities are fitted jointly across conditions by a heuristic optimisation [V qualitatively; exact functional form M].
- **Behaviour.** Non-convex fitting, so results can differ between runs; this criticism is made in the JUMPsem paper [V].
- **Implementation.** MATLAB.
- **Inputs.** Needs time or condition series.

### 1.9 INKA: Beekhof et al. 2019
- **Ref.** Beekhof R et al. *Mol Syst Biol* 2019;15(4):e8250, doi:10.15252/msb.20188250 [V].
- **Approach.** Single-sample score that integrates four components [V]:
  - kinase-centric evidence: phosphopeptides on the kinase itself, especially **activation-loop** sites
  - substrate-centric evidence from PSP
  - substrate-centric evidence from NetworKIN motif predictions
  - label-free spectral counts
- **Performance.** Reported better than its own components and than KARP [V].
- **Relevance.** Single-sample tumour ranking; uses intensities rather than FC. Shows that combining kinase-centric and substrate-centric evidence helps.

### 1.10 KEA3: Kuleshov et al. 2021
- **Ref.** *Nucleic Acids Res* 2021;49(W1):W304–W316, doi:10.1093/nar/gkab359 [V]. PMID believed 34019655 [M].
- **Method.** Gene-level **over-representation** (Fisher) against several KS, kinase–protein-interaction, co-expression and co-occurrence libraries, integrated by MeanRank/TopRank [V].
- **Fit for our task.** Not site-resolved and not directional (input is a gene list). Benchmarked on recovering the kinase from L1000 expression and from kinase-inhibitor phosphoproteomics [V]. Useful only as a weak baseline.

### 1.11 PhosR: Kim et al. 2021
- **Ref.** Kim HJ et al. *Cell Rep* 2021;34:108771, doi:10.1016/j.celrep.2021.108771 [V].
- **Method.**
  - `kinaseSubstrateScore` combines a PWM motif score with the **correlation of a site's temporal profile** to the profiles of the kinase's known substrates [V qualitatively].
  - `kinaseSubstratePred` is an ensemble classifier (SVM with adaptive sampling) [M].
  - Requires multi-condition or time-course data.
- **Leakage note.** It learns from the data matrix itself, so it is transductive.

### 1.12 VIPER: Alvarez et al. 2016
- **Ref.** *Nat Genet* 2016;48:838–847 [M].
- **Method.** aREA: analytic rank-based enrichment using mode of regulation (signed) and per-edge likelihood weights, giving a three-tailed normalised ES. Optional "pleiotropy correction" for shared targets.
- **Implementation.** R `viper`, and decoupler `viper`.
- **Performance.** Comparable to others in benchmarKIN [V: "most methods perform similarly"].

### 1.13 decoupler: Badia-i-Mompel et al. 2022
- **Ref.** *Bioinformatics Advances* 2022;2(1):vbac016 [M].
- **Framework.** Python/R; one API for many enrichment methods over a weighted network (source, target, weight).
- **Methods** [M]:
  - **ulm**: per kinase, regress x on the weight vector w (0 for non-targets) and report the slope's *t*-value.
  - **mlm**: multivariate regression of x on all kinases' weight columns jointly and report the *t*-values. It handles shared substrates but becomes unstable with collinearity or small m.
  - **wsum / wmean**: Σ w_i x_i (wmean divides by Σ|w|).
  - **norm_wmean / corr_wmean**: compare against a permutation null (default 1,000 permutations) to give a z.
  - **viper**: aREA.
  - **gsva**, **aucell**, **ora**: rank-based.
  - **udt/mdt**: decision trees.
  - **consensus**: mean z across methods.
- **Performance.**
  - In decoupler's own benchmark (TFs and pathways), mlm, ulm, norm_wmean and consensus were among the best [M].
  - In benchmarKIN on kinases, ulm and mlm were comparable to the z-score but not better [V: all methods similar; z-score best].
- **Runtime.** ulm and mlm are closed-form matrix operations; very fast.

### 1.14 JUMPsem: Kong et al. 2025 (Paper_JumpSEM.pdf + source code)
- **Ref.** Kong D, Zhang A, Li L, … Peng J, Wang X. *Commun Biol* 2025;8:103, doi:10.1038/s42003-025-07548-4 [V].
- **Algorithm** (verified in `R/singleEnzymeSEM.R`, Wanglab-UTHSC/JUMPsem):
  1. Intersect PSP substrates (optionally plus MEME/MAST motif-derived extra substrates) with measured sites. Requires >1 overlapping site.
  2. Optional KMO filtering (default `kmo.off = 0`) and nearPD repair of the correlation matrix.
  3. If #substrates ≥ #samples, keep the top (n_samples − 1) substrates by correlation with PC1.
  4. Fit `lavaan::sem("KIN =~ s1 + s2 + …")`, a one-factor CFA across **samples**. Optionally refit after adding residual covariances with modification index > 4.
  5. **Activity = `lavPredict` factor scores per sample.** "Affinity" = |loadings| with z and p. Fit indices CFI/TLI/RMSEA/SRMR are reported.
- **Implications.**
  - It needs a sites × samples matrix. Factor scores are mean-zero across the supplied samples, so the activity of one condition depends on which other conditions are co-analysed. It is transductive.
  - The latent sign follows the first indicator (lavaan fixes its loading to 1). If that substrate is anti-correlated with the rest, the sign flips. This is inferred from lavaan defaults and worth testing.
  - It needs at least 2 substrates and enough samples for CFA identifiability (≥3 indicators for a just-identified one-factor model; with 2 indicators the model is under-identified without constraints).
- **Reported performance** [V, PDF]:
  - Compared only with IKAP and KSEA.
  - On 16 datasets / 31 conditions taken from two prior benchmarks, the target kinase's rank among inferred kinases was used to compute "precision" at rank cut-offs. JUMPsem was "slightly" better than IKAP (Suppl. Fig. S5c).
  - No AUROC. Speed and kinase coverage are claimed advantages; no runtime numbers are given in the main text.
- **Fair benchmarking recipe for us.** Run it per study (or per cell line) on the full condition matrix, with log2FC as input and `input.log2.norm = FALSE`. Read off each condition's factor score. Report runtime.

### 1.15 KSTAR: Crowl et al. 2022
- **Ref.** Crowl S, Jordan BT, Ahmed H, Ma CX, Naegle KM. *Nat Commun* 2022;13:4283, doi:10.1038/s41467-022-32017-5 [V].
- **Algorithm** [V]:
  - Build **50 pruned NetworKIN networks**, each with 2,000 substrates per Tyr kinase or 2,500 per S/T kinase, at most 10 (Y) or 20 (S/T) kinases per site, and study-bias-matched sampling.
  - Binarise evidence: a site counts as "observed/changed" based on a threshold.
  - Per network, run a hypergeometric test, giving 50 p-values per kinase.
  - Compare these p-values to 150 study-bias-matched random experiments with a Mann–Whitney test. Score = −log10 p. FPR is computed by treating random experiments as real.
- **Design point.** Explicitly *does not* use quantitative FC or activation-loop sites. It is single-sample and handles study bias by matched nulls.
- **Benchmark.**
  - 15 studies / 51 conditions: Tyr 20 conditions / 19 kinases; S/T 31 conditions / 38 kinases.
  - Metric: P_hit as top-10 or FDR ≤0.05, with global accuracy **averaged per kinase**.
  - Reported about 50% better than the next-best method on Tyr binary calls. No AUROC.
  - The 18 Tyr experiments were reused in benchmarKIN.
- **Implementation.** Python (`kstar`), Nextflow for large runs. Networks are on figshare ("NETWORKS_KSTAR_v1") [V].

### 1.16 Kinase Library (KL) enrichment: Johnson et al. 2023 and Yaron-Barir et al. 2024
- **Refs.**
  - Johnson JL, Yaron TM et al. *Nature* 2023;613:759–766, PMID 36631611, doi:10.1038/s41586-022-05575-3 [V]. Covers 303 S/T kinases (more than 84% of predicted-active human S/T kinases) via positional scanning peptide arrays (PSPA).
  - Yaron-Barir TM et al. *Nature* 2024;629:1174–1181, PMID 38720073, doi:10.1038/s41586-024-07407-y [V]. Covers all human Tyr kinases (≈78 canonical plus non-canonical) [M for the count].
- **Scoring** [V]:
  - Site score = product of PSSM values, normalised by a random-peptide probability.
  - The **percentile** comes from each kinase's score distribution over a reference phosphoproteome of 82,735 S/T sites.
  - Validation: 969 sites with ≥3 reports had a median percentile of 95% for the reported kinase.
- **Enrichment** [V]:
  - Each site is assigned to the **top-15** kinases by percentile rank.
  - One-sided Fisher exact test of up-regulated vs unchanged and down-regulated vs unchanged sites, with |log2FC| thresholds of 1.5 (0.5 in one panel), Haldane correction and BH adjustment.
  - Kinases significant in both directions are dropped.
- **Leakage.** None from cellular phosphoproteomics, because these are in vitro peptide libraries. The reference phosphoproteome only defines which sites exist.
- **Caveats.**
  - Captures only intrinsic specificity: close paralogs and families are nearly indistinguishable, and docking, localisation and scaffolds are ignored.
  - In benchmarKIN, adding KL top-15 targets lowered AUROC [V].
- **Implementation.** kinase-library.phosphosite.org; Python package `kinase-library` [M].

### 1.17 PhosX: Lussana & Petsalaki 2024
- **Ref.** *Bioinformatics* 2024;40(12):btae697, doi:10.1093/bioinformatics/btae697 [V]. Python, Apache-2.0, `pip install phosx` [V].
- **Algorithm** [V]:
  - Uses Kinase Library PSSMs (S/T and Y) on 10-mers. The score is the product of matrix values, quantile-normalised (10⁴ quantiles) against more than 200k PSP phosphopeptides.
  - Each peptide is assigned to its **top-n kinases (default n = 10)**.
  - GSEA-like weighted running-sum ES with permutation p-values. Score = −log10 p × sign(ES).
- **Benchmark.**
  - Hernandez-Armenta 103 experiments, 28 kinases, 173 positives (100 up, 73 down) [V].
  - Reported "marginally" better than GSEApy and KSTAR, and "clearly" better than z-score, Kinex and PTM-SEA on AUROC [V qualitatively]. Numbers are in figures only.
  - The authors note the package has since changed substantially from the paper [V].
- **Relevance.** A prior-knowledge-free baseline that is leakage-safe. Strong evidence that motif-only methods can be competitive on HA2017.

### 1.18 Kinex: Roche, bioRxiv 2023
- **Ref.** bioRxiv doi:10.1101/2023.11.23.568445; GitHub bedapub/kinex; GPL [V].
- **Method.** Scores sites with Kinase Library matrices (82,755 pre-scored S/T peptides, 303 kinases), considering acceptor favourability and phospho-priming. Up/down sets are defined by FC thresholds (±0.5 in the PhosX benchmark), followed by enrichment [V qualitatively].
- **Performance.** Below PhosX in the PhosX paper [V].

### 1.19 CoPheeKSA / CoPheeMap: Jiang et al. 2025
- **Ref.** Jiang W, Jaehnig EJ, Liao Y, Shi Z, Yaron-Barir TM, Johnson JL, Cantley LC, Zhang B. *Nat Commun* 2025;16:2766, doi:10.1038/s41467-025-57993-2 [V].
- **Method.**
  - CoPheeMap is a machine-learned co-regulation network of 26,280 sites, built from **1,195 CPTAC tumours across 11 cancer types**.
  - Network-embedding features feed a KS classifier, which predicts 24,015 KS links between 9,399 sites and 104 S/T kinases. Predictions are checked against Kinase Library specificities [V].
- **Leakage.** Trained on CPTAC phospho data, so it must not be used with any CPTAC-based evaluation. For Hijazi, risk is low but non-zero (general co-regulation).

### 1.20 SELPHI / SELPHI2: Petsalaki lab
- **SELPHI** (2015) [M]: Petsalaki E et al. *Nucleic Acids Res* 2015;43:W276. Correlation-based KS and pathway inference from the user's own dataset.
- **SELPHI2** [V]: Maier, Petursson, Lussana, Petsalaki. *Mol Cell Proteomics* 2025. ML model producing a probabilistic human KS network from "omics datasets"; trained with PSP (download 2024-09-03); server selphi2.com. Claimed better and broader than other KS predictors.
- **Leakage.** Features include phosphosite co-regulation from public phosphoproteomics compendia (which may include Hijazi and Ochoa sets). Treat as **high leakage risk** unless training data provenance is checked.

### 1.21 NetworKIN / NetPhorest / KinomeXplorer
- **Refs.** Linding R et al. *Cell* 2007;129:1415 [M]. NetworKIN 3.0 / KinomeXplorer: Horn H et al. *Nat Methods* 2014;11:603 [M].
- **Method.** Combines NetPhorest motif probabilities with STRING-based contextual proximity between kinase and substrate protein. Gives a score per (kinase, site).
- **Usage.**
  - Human predictions 3.1 are available at netphorest.science/download/networkin_human_predictions_3.1.tsv.xz [V URL, per benchmarKIN].
  - Typical cutoff: score ≥5 (benchmarKIN: 22,788 interactions) [V].
- **Leakage.** STRING includes co-expression and text mining, but not our phospho quantifications. Low risk.
- **Performance.**
  - HA2017: worst evidence class [V].
  - benchmarKIN: adding NetworKIN to curated libraries **improved** tumour AUROC (z-score protein-based 0.69, activating-site 0.677) [V]. With a kinase-matched comparison, combined 0.689 vs curated 0.668 on tumour, and 0.803 vs 0.800 on perturbation [V].
  - KSTAR builds on it.

### 1.22 KinasePA: Yang et al. 2016
- **Ref.** *Proteomics* 2016 [M]. R package `directPA` / KinasePA web.
- **Method.** "Direction analysis": tests whether a kinase's substrates move jointly in a specified direction across 2–3 perturbations (multivariate rank-based).
- **Relevance.** Designed for multi-perturbation hypotheses, not single-condition ranking. Low priority.

### 1.23 Other or related approaches
- **PHONEMeS, CARNIVAL, PHOTON (Rudolph 2016)** [M]: network or pathway contextualisation, not per-kinase activity estimators.
- **EGF signalling-pathway inference benchmark** (Garrido-Rodriguez et al.? *Nat Commun* 2026, doi:10.1038/s41467-026-69332-0; bioRxiv 2024.10.21.619348):
  - Combines OmniPath and Phosformer KS networks with EGF meta-analysis plus three new deep EGF time courses [V].
  - The literature network recovered ground-truth interactions best.
  - KS networks cover less than 5% of measured sites [V].
  - Possible source of EGF ground truth. Author order not verified.
- **Matrix factorisation / latent-factor kinase activity.** No published dedicated method was found despite searching. The closest are JUMPsem (CFA per kinase), IKAP (bilinear activity × affinity) and Ochoa-style co-regulation. This looks like an **open niche**: a low-rank or multi-task model across our ~1,900 conditions with KS-structured loadings.
- **Deep learning / learned estimators (2023–2026).** As of Oct 2026 no peer-reviewed, benchmarked DL kinase-activity estimator trained on perturbation phosphoproteomics was found. DL is used for **KS prediction** (Phosformer/Phosformer-ST, GPS 6.0 transfer learning, LinkPhinder knowledge-graph embedding, SELPHI2), not for activity estimation. This is also an open niche, with high leakage risk.
- **Bayesian / hierarchical models.** No dedicated, maintained Bayesian kinase-activity estimator was found. Uncertainty-aware occupancy models exist, but not activity models. Empirical-Bayes shrinkage across kinases or families is unexplored in the literature.

---

## 2. Kinase–substrate resources

| Resource | Type | Size / coverage | Uses phospho quantifications? | Leakage risk for our eval |
|---|---|---|---|---|
| **PhosphoSitePlus (PSP)** Kinase_Substrate_Dataset (provided) | curated LTP+HTP literature | ~390 human kinases with ≥1 substrate (benchmarKIN) [V]; heavy study bias | Some KS curated from HTP papers | Low–moderate: check whether HTP-only edges come from evaluation PMIDs |
| PSP **Regulatory_sites** (provided) | curated site function (ON_FUNCTION: "enzymatic activity, induced/inhibited") | activating sites on 280 kinases (benchmarKIN: 787 activating sites, PSP + literature) [V] | literature | Label source for tumour-style benchmark; do not also use as features in that benchmark |
| **SIGNOR 3.0** (Lo Surdo et al. *NAR* 2023) [M] | curated causal signed interactions | several thousand phosphorylation edges | no | low |
| **OmniPath** (Türei et al. *Mol Syst Biol* 2021) [M] | meta-resource (PSP, SIGNOR, HPRD, phosphoELM, dbPTM, ProtMapper/REACH text-mining…) | 467 kinases (enzyme–substrate, phosphorylation only) [V via benchmarKIN] | text-mining edges are noisy | low; filter by curation effort / ≥2 sources |
| **iPTMnet** (Huang et al. *NAR* 2018) [M] | curated + text-mined (RLIMS-P, eFIP) | — | no | low |
| **Phospho.ELM** (Dinkel et al. *NAR* 2011) [M] | curated, mostly LTP; frozen ~2011 | ~3,000 KS | no | low (largely subsumed by PSP) |
| **HPRD** [M] | curated; frozen 2010 | — | no | low |
| **PTMsigDB** (Krug 2019) [V] | site-specific signed signatures (KINASE, PERT, PATH) | KINASE sets from PSP etc. | **PERT sets are derived from phospho experiments** | moderate: exclude PERT sets that overlap evaluation studies |
| **GPS 5.0 "gold"** (Wang C et al. 2020, *GPB*) [M] / GPS 6.0 (Chen M et al. *NAR* 2023;51:W243) [V] | curated training set / ML predictor (PLR + DNN + LightGBM, transfer learning to 577 kinase predictors; pretrained on 490,762 p-sites) [V] | 352 kinases (gold) [V] | site *existence* only (not FCs) | low |
| **NetworKIN 3.1 / NetPhorest** [V/M] | motif + STRING context | 22,788 edges at score ≥5 [V] | no | low |
| **Kinase Library** S/T (Johnson 2023) + Tyr (Yaron-Barir 2024) [V] | in vitro PSPA specificity matrices → percentile per site | 303 S/T + ~all Tyr kinases | no | **none**; best source for understudied kinases, but family-level resolution |
| **Sugiyama et al. 2019** *Sci Rep* 9:10503 [V] | in vitro kinase reactions on dephosphorylated HeLa lysate peptides, 385 kinases (354 WT protein kinases); 175,574 potential direct substrates; jPOST PXD011366 | — | in vitro MS only | none for cellular eval; contains many non-physiological sites |
| **iKiP-DB** (Mari T et al. *J Proteome Res* 2022;21:1575) [V] | filtered reprocessing of Sugiyama 2019 | 26,786 interactions, >300 kinases [V via benchmarKIN] | in vitro | none; but **lowered** benchmarKIN tumour performance when added [V] |
| **Phosformer / Phosformer-ST** (Zhou et al., Kannan lab; Phosformer *Bioinformatics* 2023 [M]; Phosformer-ST venue/year not verified) | protein-LM transformer on (kinase, peptide) pairs | kinome-wide | trained on known KS (PSP-like) | low; lowered benchmarKIN AUROC when added as top-15 [V] |
| **LinkPhinder** (Nováček et al. *PLoS Comput Biol* 2020;16:e1007578) [V] | knowledge-graph embedding | "highest kinome coverage" vs 6 tools [V] | no (KG from curated DBs) | low; but trained on PSP-like edges, so overlap is circular |
| **CoPheeKSA** (Jiang 2025) [V] | ML on CPTAC co-regulation | 24,015 links, 104 S/T kinases | **yes (CPTAC)** | high for CPTAC eval |
| **SELPHI2** (Maier et al. *MCP* 2025) [V] | ML on omics co-regulation + features | broad | **likely yes** | high unless provenance checked |
| **Hijazi PDTs / EBDT** (Hijazi 2020) [V] | inhibitor-fingerprint-derived "putative downstream targets", direct + indirect; 6,206 sites, 103 kinases, 19,410 K–P relations | — | **derived from the Hijazi inhibitor phosphoproteomes we evaluate on** | **direct leakage: never use for Hijazi evaluation** |
| **Ochoa 2016 / Invergo 2020 co-regulation KS predictions** [M] | learned from phospho compendia | — | yes | high |
| **KinPred** (Xu et al. 2021, Naegle lab) [M] | harmonised multi-predictor resource | — | depends on source predictor | check |
| Dark Kinase Knowledgebase (Berginski et al. *NAR* 2021) [M] | IDG list of understudied kinases | ~160 "dark" kinases | — | use to define the understudied stratum |

**Key facts on coverage and bias.**
- Fewer than 5% of measured sites have a known kinase. The 2026 EGF benchmark [V] and JUMPsem (<10% of quantified phosphopeptides in PSP) [V] both report this.
- About 30% of UniProt kinases have no known target (SELPHI2 preprint) [V].
- Hijazi's P31/Fuj test illustrates the gap: of 8,100 phosphopeptides, PSP mapped 403 sites to 58 kinases, while PDTs mapped 1,746 sites to 82 kinases [V].

**Regulatory/activation-loop sites as activity readouts.**
- Activation-loop (A-loop) phosphorylation is a canonical activity marker. INKA uses it, and benchmarKIN's tumour benchmark uses 787 activating sites as labels [V].
- **Pitfall:** ATP-competitive inhibitors can *increase* A-loop or hydrophobic-motif phosphorylation of their own target. This is paradoxical hyperphosphorylation, documented for AKT inhibitors (Okuzumi et al. *Nat Chem Biol* 2009 [M]) and for RAF inhibitors in RAS-WT cells [M]. Feedback relief (e.g. MEK inhibitors raise upstream RAF/MEK phosphorylation [M]) causes the same problem.
- **Consequence:** A-loop sites are **not** valid inhibition labels in inhibitor experiments and may be anti-informative as features there. HA2017 and benchmarKIN both **removed autophosphorylation sites** from KS sets to avoid leakage [V].

---

## 3. Ground-truth and benchmark datasets

### 3.1 Already in our compilation (`data/01_phosphoproteomics_conditions.csv`: 1,937 conditions, 178 PMIDs)
- **Hijazi et al. 2020.**
  - Ref: *Nat Biotechnol* 38:493–502, doi:10.1038/s41587-019-0391-9; PMID 31959955; PRIDE **PXD015943** [V].
  - Compilation: 183 conditions = 61 treatments × HL60/MCF7/NTERA2, all at 60 min [V, our file]. The paper's ancillary P31/Fuj 6-inhibitor validation set is not in the compilation.
  - Protocol [V]: cells treated with **1–10 µM compound for 1 h**; label-free MS1 quantification (Pescal); phosphopeptides with the same site summed and quantile-normalised; sites with FDR <0.01 and control CV <50% kept (9,326 / 9,416 / 12,591 sites per line).
  - **Target annotation:** DiscoverX profiling of 20 inhibitors at **200 nM** against 403 kinases, plus literature/vendor data for 40 compounds. Remaining activity was min–max normalised, and >50% inhibition was taken as target [V]. This explains `data/kinhib_inhib_comb__1.csv`.
  - **Caveats:**
    - Profiling was at 200 nM, but cells were dosed at 1–10 µM, so cellular polypharmacology is broader than the annotation. Unannotated inhibited kinases will appear as "false positives".
    - Vendor-reported inhibition without numbers was set to full inhibition [V].
  - benchmarKIN's handling [V]:
    - used HL60 + MCF7 (61 experiments)
    - removed sites with FDR >0.05 and split multiply phosphorylated peptides
    - used a **TTD-based manually curated target list**, which raised average AUROC by about 0.1 versus the >50%-inhibition list (Suppl. Table 1 of benchmarKIN)
- **decryptM (Zecha et al. 2023).**
  - Ref: *Science* 380:93–101 (2023), PMID 36926954 [V citation; PRIDE accession not verified].
  - Compilation: **790 conditions**, 13 cell lines (A549, ARH-77, Ramos, SU-DHL-4, HeLa, RPMI8226, K562, PC-9, A431, BT-474, …) [V, our file].
  - Drugs include dasatinib, imatinib, gefitinib, afatinib, lapatinib, AZD4547, MK2206, pictilisib, dactolisib, AZD8055, selumetinib, refametinib, PD325901, nintedanib, staurosporine, tideglusib, SHP099, plus non-kinase drugs (HDAC/proteasome inhibitors, antibodies, chemotherapy). These non-kinase drugs serve as **negative-control conditions**.
  - Dose series run from about 0.01 nM to µM, at 30 min and other timepoints.
  - **Strength:** dose-resolved data allow graded truth. Kinobeads Kd_app values (Klaeger 2017; many of the same Kuster-lab drugs) predict at which dose each kinase becomes inhibited.
  - CurveCurator (Bayer FP et al. *Nat Commun* 2023 [M]) is the companion dose-response classifier.
- **Other perturbation studies in the compilation** (examples from a scan):
  - PMID 32034161: RPE1 cells, EGF plus MEK/AKT inhibitor titrations, 10–30 min
  - PMID 38047086: CRC organoids, BRAFi/MEKi, 2 h
  - PMID 34086370: FGF7/FGF10 stimulation time course
  - PMID 34663829: TNF ultra-fast time course in U937
  - "Ochoa2016" source: 39 conditions, likely overlapping HA2017
  - Each needs manual curation of expected kinase changes; see §4.
- **Casado 2013** (PMID 23532336) is also in the compilation.

### 3.2 Hernandez-Armenta 2017 / Ochoa 2016 perturbation gold standard
- Original: 184 positive kinase–condition pairs (30 kinases, 62 conditions, 24 publications); EGF, inhibitors, DNA damage, mitosis, etc. Table 1 of Paper_KSEA.pdf lists them [V].
- benchmarKIN's version: 103 quantile-normalised experiments, 27 kinases, 173 positives (100 up / 73 down) [V]. Zenodo **records/5645208** [V URL from benchmarKIN].
- Ochoa 2016 atlas: *Mol Syst Biol* 12:888, doi:10.15252/msb.20167295, 399 conditions, 215 kinases [V]. Data were originally on phosfate.com, which may be offline [M].
- **Overlap warning.** Several HA2017 studies (e.g. Wilkes 2015 inhibitor panel) may also be in our compilation. De-duplicate by PMID when building train/test splits.

### 3.3 benchmarKIN resources
- Ref: Müller-Dott S, Jaehnig EJ, Munchic KP, Jiang W, Yaron-Barir TM, Savage SR, Garrido-Rodriguez M, Johnson JL, Lussana A, Petsalaki E, Lei JT, Dugourd A, Krug K, Cantley LC, Mani DR, Zhang B, Saez-Rodriguez J. *Nat Commun* 2025;16:4771, PMID 40404650, doi:10.1038/s41467-025-59779-y [V for citation, PMID and DOI; author list partly M].
- Code: github.com/saezlab/benchmarKIN (R, GPL-3; docs benchmarkin.readthedocs.io) and github.com/saezlab/kinase_benchmark. Archived as Zenodo 10.5281/zenodo.15118822. Benchmark data at zenodo 12566560 [V].
- Perturbation set: 230 experiments, about 80 kinases. It combines HA2017, Hijazi HL60/MCF7, and 18 Tyr-kinase experiments from KSTAR (dasatinib, imatinib, bosutinib, nilotinib, EGF, HRG, TCR) [V].
- **CPTAC tumour benchmark** [V]:
  - 7 cohorts: BRCA, CCRCC, GBM, HNSCC, LUAD, LSCC, UCEC. COAD, OV and PDAC were excluded for low site–protein correlation.
  - Protein-based labels: top/bottom 5% of tumours by kinase protein abundance; 12,850 pairs, 388 kinases.
  - Activating-site labels: 787 activating sites on 280 kinases; 3,446 pairs, 96 kinases.
  - Scores are z-scored per cohort. AUROC on 80% subsamples ×1,000.
  - CPTAC data: Proteomic Data Commons / LinkedOmics, pan-cancer harmonised (Li Y et al. *Cancer Cell* 2023 [M]).
- **Relevance.** The tumour benchmark supplies orthogonal, non-perturbation truth for many more kinases (388), including less-studied ones. It is usable as a *secondary* test, provided CoPheeKSA and any CPTAC-trained edges are excluded.

### 3.4 Kinase inhibitor target profiles (who is truly inhibited)
- **Klaeger et al. 2017.**
  - Ref: *Science* 358:eaan4368, PMID 29191878, doi:10.1126/science.aan4368 [V].
  - Kinobeads competition profiling of **243 clinical kinase drugs** in a cell-lysate mix, giving apparent Kd (Kd_app) per drug–kinase pair.
  - Provided in `data/Klaeger/` (aan4368_Table_S1/S2, Table_S3). Also ProteomicsDB. PRIDE accession believed to be PXD005336 [M].
  - **Use:** at concentration c, a kinase is predicted to be inhibited by fraction ≈ c / (c + Kd_app), ignoring ATP competition in cells. This gives a graded, dose-aware target matrix for decryptM and for any compilation inhibitor at a known dose.
- **Davis et al. 2011** (*Nat Biotechnol* 29:1046; 72 inhibitors × 442 kinases, Kd) and **Anastassiadis et al. 2011** (*Nat Biotechnol* 29:1039; 178 inhibitors × 300 kinases, % activity) [M]. These are the literature sources Hijazi used [V that literature sources were used].
- **Therapeutic Target Database (TTD)**: the source of the benchmarKIN "curated" Hijazi targets [V].

### 3.5 Other candidate ground-truth sources (not in compilation; acquire only if cheap)
- **LINCS P100** (Abelin et al. *MCP* 2016; Litichevskiy et al. *Cell Syst* 2018) [M]: targeted ~96-site panel, ~90 drugs × 6 cell lines. Too few sites for kinome-wide inference. Low priority.
- **Mitchell et al. 2023.** Ref: *Nat Biotechnol* 41:845–857, doi:10.1038/s41587-022-01539-0 [V]. 875 compounds in HCT116, mainly **protein abundance** after 24 h. Phospho profiling is limited or absent [M], so it is of low value for kinase activity.
- **Sugiyama 2019** jPOST PXD011366: an in vitro *KS source*, not ground truth.
- **EGF time courses**:
  - Olsen et al. 2006 *Cell* (HeLa) [M]
  - Francavilla et al. 2016 *NSMB* [M]
  - the 2026 *Nat Commun* EGF benchmark with three new deep time courses [V]
  - These are expected EGFR→MAPK/AKT activation truth, partly in HA2017 already.
- **Insulin** (Humphrey et al. 2015 *Nat Biotechnol*): mouse, so not human. **Analog-sensitive kinase** and **CRISPR-KO** phosphoproteomes in human cells are scattered single-kinase studies [M]. Useful for understudied kinases, but curation-heavy. Flag as follow-up.
- **Ochoa 2020** functional landscape: *Nat Biotechnol* 38:365–373, PMID 31819260 [V]. Reference phosphoproteome of **119,809 sites** from 6,801 reanalysed experiments, with a 59-feature functional score (R package `funscoR`) [V]. Useful as site prior (functional sites may be better readouts), not as ground truth.

---

## 4. Evaluation methodology

### 4.1 How published benchmarks define positives, negatives and metrics
| Benchmark | Positives | Negatives | Metric(s) | Aggregation |
|---|---|---|---|---|
| HA2017 [V] | 184 curated kinase–condition pairs, signed by expected direction | equal number of random kinase–condition pairs from the positive set's kinases and conditions, 60 randomisations | AUROC, precision@recall 0.5 on |signed activity| | pooled across pairs |
| RoKAI 2021 [V] | 128 annotations (25 kinases) | all other kinases in a condition | P_hit(k = 10); KS statistic | pooled |
| KSTAR 2022 [V] | expected perturbed kinases | — | P_hit (top-10 or FDR ≤0.05) | **average of per-kinase accuracy** (macro) |
| benchmarKIN 2025 [V] | perturbed kinase per experiment, with activity multiplied by expected sign | other kinases scored in that experiment, **subsampled to match positives, 1,000×** | scaled rank (rank/#kinases), PHit(k = 5, 10, 20), AUROC | averaged; experiments with fewer than 5 targets in the library removed (230 → 135 for PSP) |
| PhosX 2024 [V] | HA2017 173 signed positives | (a) up vs down; (b) random kinases, 100× | AUROC, AUPR, normalised rank, top/bottom 5% recovery | pooled |
| JUMPsem 2025 [V] | target kinase | — | precision at rank cut-offs | pooled |

### 4.2 Known pitfalls
1. **Study bias / "rich get richer".**
   - Benchmark positives are dominated by well-studied kinases (RoKAI: median 29 vs 4 substrates) [V].
   - √m-scaled statistics favour these kinases, so pooled AUROC rewards study bias.
   - Mitigations:
     - macro-average per kinase (KSTAR style)
     - stratify by substrate count (≤5, 6–20, >20; HA2017 used ≤5 vs >5)
     - subsample substrates to fixed m
     - hide-substrate robustness tests (RoKAI)
     - report coverage: the number of kinases with a score per condition
2. **Coverage-dependent comparisons.**
   - Methods or libraries that score different kinase sets cannot be compared on raw ranks. benchmarKIN uses scaled rank plus a "reduced common-kinase subset" comparison; there NetworKIN 0.77, OmniPath 0.76, curated 0.75 [V].
   - Always also report results on the **intersection** of scorable kinases, and on the union with an explicit "not scored = rank of median/neutral" rule.
3. **Kinase family correlation.**
   - Paralogs (AKT1/2/3, MAPK1/3, PRKACA/B, ROCK1/2, CDK1/2, PKC isoforms, RSK1–4) share substrates and inhibitors.
   - Hijazi grouped AKT1/2 and MAPK1/3 [V].
   - Inhibitors almost always hit paralogs, so a paralog scored as "negative" is a label error. Use paralog-aware negatives: exclude same-family kinases inhibited by the drug, or evaluate at family level as a secondary metric.
4. **Inhibitor polypharmacology, downstream effects and feedback.**
   - Hijazi dosed at 1–10 µM but annotated at 200 nM [V].
   - Downstream kinases (e.g. MEKi → ERK → RSK, p90S6K; mTORi → S6K) legitimately decrease. Counting them as negatives penalises correct biology.
   - Treat known-downstream kinases as **ambiguous (excluded from negatives)** rather than negatives.
   - Use Klaeger Kd_app for dose-aware target sets.
5. **Activation-loop paradox** (§2): never use target A-loop sites as inhibition labels.
6. **Leakage.**
   - (a) KS predictors learned from the evaluation data (Hijazi PDTs, SELPHI2, CoPheeKSA vs CPTAC, PTMsigDB PERT sets, PhosR/co-regulation learned on the compilation).
   - (b) PSP edges curated from the evaluation papers themselves (filter PSP by source PMID where available, or at least test sensitivity).
   - (c) Supervised models trained on conditions sharing the **same inhibitor** (Hijazi: the same drug in 3 cell lines) or the same study/batch. Split by **inhibitor/compound** (grouped CV) and ideally also by **target kinase**, so that test kinases are unseen as positives.
   - (d) Autophosphorylation sites in the KS set (excluded in HA2017 and benchmarKIN) [V].
   - (e) Hyperparameter tuning on the test benchmark.
7. **Normalisation and site mapping.**
   - Isoform-specific positions and multiply phosphorylated peptides need consistent handling: HA2017 used monophospho-only; benchmarKIN split multi-sites.
   - Quantile normalisation was used in HA2017 and Hijazi.
   - Directional statistics (log2FC) make centring matter: KSEA centres, RoKAI does not.
8. **Sign conventions.** For inhibitor experiments, expected activity is negative. Benchmarks multiply by the expected sign so that the correct kinase should be top-ranked [V benchmarKIN]. Two-sided "regulated vs not" metrics (HA2017) and signed metrics give different rankings, so report both.

### 4.3 Recommendations for our evaluation
- **Primary metric.** Per-experiment signed rank of each true target (scaled rank; PHit@k), and AUROC with matched negatives, **macro-averaged first within kinase, then across kinases**. Also report the pooled (micro) value for comparability with benchmarKIN and HA2017.
- **Strata.** Report all of the following:
  - by kinase substrate count (PSP m: 0 / 1–4 / 5–19 / ≥20)
  - by IDG "dark" status
  - by kinase group (AGC, CAMK, CMGC, STE, TK, TKL, CK1, atypical)
  - by dataset (Hijazi, decryptM, other)
- **Coverage.** Report the number of kinases scored per condition, and the fraction of true targets that are scorable. A method that scores more kinases at equal accuracy is better (the project goal).
- **Uncertainty.** Bootstrap over **experiments grouped by inhibitor** for CIs. Use a paired comparison against the z-score baseline on identical kinase sets.
- **Controls.**
  - degree-preserving randomised KS network (should give ≈0.5) [V HA2017 design]
  - non-kinase drug conditions in decryptM (should show no strong kinase decrease)
  - permuted FC within condition
- **Runtime.** Wall-clock per condition and for the full compilation, on the same hardware, for every method.

---

## 5. Implications for our design (ranked by expected gain × evidence)

1. **Library engineering first: a confidence-weighted union library.**
   - *Evidence:*
     - benchmarKIN: the library dominates the algorithm, and the curated combination is best [V].
     - HA2017: in vivo > in vitro > in silico, and motif weighting gives a small gain [V].
     - benchmarKIN: NetworKIN helps on tumour, while naive top-15 KL/Phosformer hurt [V].
   - *Design:*
     - Start from the union of PSP + SIGNOR + OmniPath-curated + GPS-gold (+ PTMsigDB KINASE sets).
     - Give each edge a weight from evidence type (LTP/in vivo > HTP > in vitro), number of sources, and motif percentile (Kinase Library) for that kinase–site pair.
     - Add **predicted edges only as low-weight priors**, with thresholds tuned on the training split: KL percentile ≥ ~90–99 *and* rank-in-top-k, NetworKIN ≥ cutoff, iKiP/Sugiyama.
     - Prefer adding predicted edges mainly for kinases with few curated substrates.
     - Exclude PDTs, SELPHI2, CoPheeKSA and PERT sets from Hijazi/CPTAC evaluation.
2. **Fix the statistic's study bias and calibrate per kinase.**
   - *Evidence:* z ∝ √m (KSEA/RoKAI formula). RoKAI and KSTAR both identify study bias. KSTAR's matched nulls and macro accuracy address it.
   - *Design:* options include:
     - weighted mean with an analytic or permutation null per kinase (norm_wmean-like)
     - empirical-Bayes shrinkage of mean substrate FC toward 0 (or a family prior) with variance σ²/m_eff
     - correcting for correlated substrates (same protein, multi-site peptides) by collapsing per protein or estimating effective m
     - all give comparable scales across kinases with few vs many substrates
3. **Joint (multivariate) modelling of overlapping substrates.**
   - *Evidence:* mixed. Ridge MLR was worse than Z in HA2017 for small m [V]. decoupler mlm and RoKAI's linear model are competitive but not dominant.
   - *Design:* use a ridge/elastic-net or hierarchical regression x ≈ W·a with **edge weights W from step 1** and a prior a ~ N(0, τ²), possibly with family-structured covariance. Combine with the univariate score (e.g. average of standardised scores). Explicitly test on paralog-heavy cases.
4. **Network propagation / borrowing from unannotated sites.**
   - *Evidence:* RoKAI improves all base methods and robustness to hidden substrates (p < 0.05) [V].
   - *Design:* include RoKAI-style refinement (PPI, structural distance, co-evolution) as an option. It is cheap: one sparse solve per condition, or batched across conditions.
5. **Learn edge or site reliability weights from perturbation data, with strict grouped splits.**
   - *Evidence:* curated target lists change AUROC by about 0.1 (benchmarKIN), so labels and edges matter. HA2017 suggests "rational selection of target sites that best proxy kinase activity" [V].
   - *Design:*
     - On training inhibitors only, estimate per-edge "responsiveness" (how consistently a site drops when its kinase is inhibited), then shrink toward the edge prior.
     - Evaluate with leave-inhibitor-out and leave-kinase-out CV.
     - Ochoa 2020 functional scores can serve as a site prior.
   - *Risk:* overfitting to well-studied kinases. Measure separately on kinases unseen in training.
6. **Family-level hierarchical priors for understudied kinases.**
   - *Evidence:* KL motif similarity clusters the kinome [V]. Paralogs share substrates.
   - *Design:* kinase activity a_k = family effect + kinase deviation, with deviation shrunk more when m is small. Output both kinase-level and family-level activities with uncertainty. This enables "harmonious" predictions for dark kinases rather than silence.
7. **Motif-only (prior-free) channel as a fallback and ensemble member.**
   - *Evidence:* PhosX is competitive with KS-based methods on HA2017 [V].
   - *Design:* a PhosX-like KL score per kinase, used for coverage and blended with curated scores, with weights learned on training folds.
8. **Kinase-centric evidence (A-loop and regulatory sites) as an auxiliary feature only outside inhibitor contexts.**
   - *Evidence:* INKA gains from it [V]. benchmarKIN tumour activating sites track activity [V]. But ATP-competitive inhibitors cause paradoxical A-loop hyperphosphorylation [M].
   - *Design:* include optionally; evaluate separately; disable by default for inhibitor data.
9. **Dose-aware evaluation and possibly training targets from decryptM + Klaeger Kd_app.**
   - *Evidence:* graded truth exists in our compilation (790 conditions) [V]. Kinobeads Kd_app is provided.
   - *Design:* a monotone dose–activity relationship (Spearman of inferred activity vs log dose for kinases with Kd_app < max dose) is a powerful, less biased secondary metric.
10. **Speed.**
    - z-score, ulm and ridge are vectorisable as sparse matrix products over all ~1,900 conditions at once (milliseconds to seconds).
    - Avoid per-condition permutations (GSEA, PTM-SEA, VIPER-permutation, KSTAR with 50 networks × 150 nulls) in the default path. Use analytic nulls.
    - JUMPsem/lavaan is per-kinase CFA across samples, so it is slower and not single-sample. Time it but do not emulate it.

**Baselines to implement (minimum).**
- KSEA z (centred) and RoKAI z (uncentred)
- RoKAI + network
- ulm, mlm, norm_wmean, VIPER (decoupler)
- fgsea
- PTM-SEA
- Kinase Library Fisher enrichment and PhosX (motif-only)
- KSTAR (optional; slow)
- JUMPsem (per study matrix)
- a degree-preserving random-library control

Run each on PSP alone and on the final library, so that algorithm gains are separated from library gains (the benchmarKIN lesson).

---

## 6. Open uncertainties to resolve (owner: data/methods agents)
- Exact benchmarKIN AUROC tables per method × library. The PMC text gives 0.79 (z-score, perturbation) in Results and 0.76 in another passage. Pull the Source Data from Zenodo 12566560 / GitHub saezlab/kinase_benchmark.
- Whether PSP Kinase_Substrate_Dataset contains HTP edges curated from Hijazi 2020, decryptM or other compilation PMIDs. The PSP download lacks per-edge PMIDs, but Phosphorylation_site_dataset / LT_LIT / MS_LIT counts can help.
- decryptM PRIDE accession and exact timepoints per drug (check the paper's data availability).
- Phosformer-ST citation details; KARP exact K-score formula; KinSwingR paper details; Kinex final publication status.
- JUMPsem sign anchoring in practice: test on a toy dataset with known latent direction.

---

## 7. Reference list (quick lookup)
- Alvarez MJ et al. 2016 *Nat Genet* 48:838: VIPER [M]
- Badia-i-Mompel P et al. 2022 *Bioinform Adv* 2:vbac016: decoupler [M]
- Beekhof R et al. 2019 *Mol Syst Biol* 15:e8250, doi:10.15252/msb.20188250: INKA [V]
- Casado P et al. 2013 *Sci Signal* 6:rs6, PMID 23532336: KSEA [M]
- Chen M et al. 2023 *NAR* 51:W243, doi:10.1093/nar/gkad383: GPS 6.0 [V]
- Crowl S et al. 2022 *Nat Commun* 13:4283, doi:10.1038/s41467-022-32017-5: KSTAR [V]
- Engholm-Keller K et al. 2019 *PLoS Biol* 17:e3000170: KinSwingR [M]
- Hernandez-Armenta C et al. 2017 *Bioinformatics* 33:1845, doi:10.1093/bioinformatics/btx082 [V]
- Hijazi M et al. 2020 *Nat Biotechnol* 38:493, doi:10.1038/s41587-019-0391-9, PMID 31959955, PXD015943 [V]
- Horn H et al. 2014 *Nat Methods* 11:603: KinomeXplorer/NetworKIN 3.0 [M]
- Jiang W et al. 2025 *Nat Commun* 16:2766, doi:10.1038/s41467-025-57993-2: CoPheeKSA [V]
- Johnson JL et al. 2023 *Nature* 613:759, PMID 36631611: Kinase Library S/T [V]
- Kim HJ et al. 2021 *Cell Rep* 34:108771: PhosR [V]
- Klaeger S et al. 2017 *Science* 358:eaan4368, PMID 29191878: kinobeads [V]
- Kong D et al. 2025 *Commun Biol* 8:103, doi:10.1038/s42003-025-07548-4: JUMPsem [V]
- Krug K et al. 2019 *MCP* 18:576, PMID 30563849: PTM-SEA/PTMsigDB [V]
- Kuleshov MV et al. 2021 *NAR* 49:W304, doi:10.1093/nar/gkab359: KEA3 [V]
- Linding R et al. 2007 *Cell* 129:1415: NetworKIN [M]
- Lussana A, Petsalaki E 2024 *Bioinformatics* 40:btae697: PhosX [V]
- Maier et al. 2025 *MCP*: SELPHI2 [V; volume/pages M]
- Mari T et al. 2022 *J Proteome Res* 21:1575, doi:10.1021/acs.jproteome.2c00198: iKiP-DB [V]
- Mischnik M et al. 2016 *Bioinformatics* 32:424, PMID 26628587: IKAP [V]
- Mitchell DC et al. 2023 *Nat Biotechnol* 41:845, doi:10.1038/s41587-022-01539-0 [V]
- Müller-Dott S et al. 2025 *Nat Commun* 16:4771, PMID 40404650, doi:10.1038/s41467-025-59779-y: benchmarKIN [V]
- Nováček V et al. 2020 *PLoS Comput Biol* 16:e1007578: LinkPhinder [V]
- Ochoa D et al. 2016 *Mol Syst Biol* 12:888, doi:10.15252/msb.20167295 [V]
- Ochoa D et al. 2020 *Nat Biotechnol* 38:365, PMID 31819260 [V]
- Sugiyama N, Imamura H, Ishihama Y 2019 *Sci Rep* 9:10503, PXD011366 [V]
- Wiredja DD et al. 2017 *Bioinformatics* 33:3489, PMID 28655153: KSEAapp [V]
- Wilkes EH et al. 2017 *MCP* (PMC5587867): KARP [V title / M details]
- Yaron-Barir TM et al. 2024 *Nature* 629:1174, PMID 38720073: Kinase Library Tyr [V]
- Yılmaz S et al. 2021 *Nat Commun* 12:1177, PMID 33608514: RoKAI [V]
- Zecha J et al. 2023 *Science* 380:93, PMID 36926954: decryptM [V]
