# Reference tools for baseline validation

Installed 2026-10-27. Nothing committed to git.

## Versions
| Tool | Version | Where |
|---|---|---|
| R | 4.4.3 (conda-forge) | /workspace/external/envs/r (micromamba /home/tools/bin/micromamba, root /workspace/external/mamba) |
| JUMPsem | 1.1.0 (repo Wanglab-UTHSC/JUMPsem, commit 52919f560b9ac038f3641b8794ddacfd3eca40c6; doc/ also ships a JUMPsem_2.0.pdf manual, the code is 1.1.0) | source in /workspace/external/src/JUMPsem |
| lavaan 0.7.2, psych 2.6.9, EFAtools 1.1.0, tidyverse 2.0.0, devtools 2.5.2 | | R env |
| KSEAapp | 2.0 (CRAN, 2026-03-13) | R env |
| decoupler (py) | 2.3.1 (API: `dc.mt.ulm/mlm/viper/zscore/waggr/gsea/...`, `dc.op`, `dc.pp`) | /home/plm_env |
| gseapy | 1.3.1 | /home/plm_env |
(Python 3.12.3, numpy 2.5.3, pandas 3.0.5, scipy 1.18.1.) The pre-existing /home/micromamba/envs/r_env (R 4.5.3) has no lavaan/JUMPsem; not used.

## Reproduce
`bash /workspace/scripts/setup_reference_tools.sh`. Gotchas found: (1) `r-efatools` not on conda-forge, installed from CRAN, which needs the `compilers` package AND `PATH=$ENV/bin:$PATH`; (2) JUMPsem install fails (lazydata, "no lines available in input" on BOM CSVs in data/) unless `LC_ALL=C.UTF-8`; (3) the repo's `devtools::install_github` is not needed, local install from the clone works.

## Command-line use (all actually run; toy files in /workspace/external/toy/)
```
export PATH=/workspace/external/envs/r/bin:$PATH LC_ALL=C.UTF-8
cd /workspace/external/toy
Rscript run_toy.R      # JUMPsem on simulated data (inputs toy_input.csv, toy_db.csv; outputs Activity_*.txt, Affinity_*.txt, Eval_table_*.txt)
Rscript run_ksea.R     # KSEAapp (output toy_ksea_scores.csv; input toy_ksea_PX.csv); or: Rscript run_ksea.R px.csv out.csv
```
### JUMPsem toy (run_toy.R)
Input: 3 kinases (KA 6 substrates, KB 5, KC 5), 6 samples, raw intensities `2^(10 + b_i*act + noise)`, true activity KA=(-2,-1,0,1,2,3), KB=(2,1,0,0,-1,-2), KC random. Input columns must be: `Substrate_ACC, Substrate_GN, sites, Peptide, <samples...>` (col 2 and 4 are dropped; sites may be comma-separated, e.g. "S10,T12"). Call:
`JUMPsem(input=inp, datatype="psp", organism="human", database=db, input.log2.norm=TRUE)` where db is a data.frame with PSP-format columns (GENE, KINASE, KIN_ACC_ID, KIN_ORGANISM, SUB_ACC_ID, SUB_GENE, SUB_ORGANISM, SUB_MOD_RSD, ... 16 cols total; bundled `psp_data()` is the full PSP table, used when `database=NULL`).
Output (Activity, sample columns; true pattern recovered incl. sign):
```
KA  -4.061 -2.778 -0.772  0.917  2.253  4.440
KB   1.586  0.484 -0.091 -0.168 -0.240 -1.572
KC  -0.454  0.161 -0.525  1.175  0.302 -0.660   (noise kinase)
```
Affinity (|loading|, z, p): KA first substrate PA6 loading 1.000 (z=NA, fixed marker), PA5 0.819 z=30.1, ... Eval table: CFI/TLI/RMSEA/SRMR per kinase. Note: the function wraps everything in `try(..., silent=TRUE)`: on failure it returns an invisible "try-error" instead of stopping.
### KSEAapp toy (run_ksea.R)
`KSEA.Scores(KSData, PX, NetworKIN=FALSE, NetworKIN.cutoff=5)`; PX columns: `Protein, Gene, Peptide, Residue.Both, p, FC` (FC linear; "S1;T2" multi-sites are split). Toy = first 300 rows of the bundled PX (KSData bundled, 19757 PSP+NetworKIN rows). Output columns: Kinase.Gene, mS, Enrichment, m, z.score, p.value, FDR; top rows: ULK1 m=1 z=1.21; CSNK2A1 m=1 z=1.18; PRKCQ m=2 z=-0.95; MTOR m=2 z=-0.91.
KSEA formula (source): per-kinase mS = mean over (kinase,substrate gene,site) of log2|FC| (duplicates per site averaged first); z = (mS - mean(all log2FC)) * sqrt(m) / sd(all log2FC); p = pnorm(-|z|) (two-sided is NOT applied, i.e. p=one tail of |z|); FDR = BH. Kinases with m>=1 all reported (no min-m filter). Note log2|FC| discards sign of FC input if FC negative.

## JUMPsem algorithm (as implemented; R/JUMPsem.R, inputTransform.R, singleEnzymeSEM.R, pspAdjacency.R, catchResults.R)
1. Input: matrix sites x samples (multi-sample REQUIRED). Defaults: kmo.off=0, mdsite=TRUE, input.log2.norm=FALSE, relative.norm.p=TRUE, whole.proteome=NULL, relative.norm.w=TRUE, enzyList=NULL (all kinases), enzyme.organism=organism, database=NULL (bundled PSP), motif=NULL.
2. Preprocess (inputTransformPSP): keep cols 1,3,5.. ; split comma-separated sites into rows; id = `ACC_site`; drop duplicated ids (keeps first); optional log2; then per-site row-mean centering across samples (`x - rowMeans(x)`; so single-sample data would become all zeros; NA would propagate into rowMeans -> NA row). If whole.proteome given: protein-level (log2, row-mean-centred) is subtracted from the site matrix, sites without a protein row are dropped (complete.cases). Rownames "-" -> ".".
3. Adjacency (pspAdjacency): PSP rows with KIN_ORGANISM in enzyme.organism and SUB_ORGANISM in organism; site id = SUB_ACC_ID_SUB_MOD_RSD; kinase name = toupper(GENE); binary 1/0 matrix sites x kinases (unique pairs). mdsite=FALSE: gene-level using column 7 (SUB_ACC_ID) and matching first "_" token. With motif: pspMotifAdjacency adds motif-derived pairs.
4. For each kinase (singleEnzymeSEM): substrates = adjacency==1 sites; intersect with input rows; need >1 overlapping (else NULL, kinase silently dropped). Take the sample x site table.
5. KMO filter: cor (pairwise.complete); if >2 sites: if 1/kappa(cor)<1e-12 or non-finite -> Matrix::nearPD(corr=TRUE); `KMO()` (load warning: psych::KMO masks EFAtools::KMO, so psych::KMO is the one in effect) and keep sites with per-item MSA ([[2]]) > kmo.off (0 => essentially all). Need >1 sites left.
6. If #sites >= #samples: prcomp(scale=TRUE) PC1; rank sites by signed cor(PC1, site) descending; keep top (n_samples-1) (most positively correlated with PC1). (Quirk: lavaan data=all sites but model only uses selected; sign of PC1 is arbitrary from SVD.)
7. Model: one-factor CFA `KINASE =~ s1 + s2 + ...` fitted by `lavaan::sem(model, data=sites matrix)`; defaults: ML, std.lv=FALSE (first listed indicator's loading fixed to 1 = marker variable), no mean structure, orthogonal residuals. Indicator order = row order of the preprocessed input (or the PC1 ranking in step 6).
8. Modification indices: if any residual covariance (op "~~") MI>4, ALL of them are added at once ("a ~~ b") and the model refit (once, no iteration).
9. Activity = `lavPredict(fit)` (default method "EBM" = regression/Bartlett-like factor scores, one value per sample) for the latent kinase factor; reported raw (commented-out code for z-score / mean-center versions, so output is NOT standardised). Scale = units of first (marker) indicator; the model is identified by that marker, so SIGN = sign relative to the first indicator (positive when sample intensity of the marker goes up); this is why the sign convention depends on site order (typical: kinase "up" = marker site up). Loadings other than marker can be negative; Affinity table stores abs(loading), z, p (marker z/p = NA).
10. Evaluation: CFI, TLI, RMSEA, SRMR (fitMeasures) with text labels (>=.95 good, >=.90 ok; RMSEA <=.05/.08; SRMR <=.05/.10). No filtering on fit: poor-fit kinases are still reported. Kinases where lavaan errors are dropped. Output files Activity_/Affinity_/Eval_table_<timestamp>.txt in output.folder; return list(Activity (ENZYME + one column/sample), Affinity, Evaluations).
Reimplementation hints for Python: factor-analytic one-factor ML (semopy or own ML) with marker fixing, then regression factor scores; score is sample-relative (centered), so only multi-sample comparisons are meaningful; for validation compare Spearman/Pearson per kinase across samples against R output (sign-align via marker). Adding residual covariances (step 8) changes the scores, so reproduce it or note deviation.
