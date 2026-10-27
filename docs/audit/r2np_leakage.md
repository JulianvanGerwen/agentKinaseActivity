# Audit: can r2np / r3 / cantley encode compendium perturbation RESPONSES?

Independent audit, 2026-10-27. Scope: predicted KS libraries `r2np` (kinaseSubMain R2 "nophospho", DEV-only, seed 0,
sha `1edbd269…`), `r3` (R3 deploy, sha `2f62988c…`), `cantley` (`cantley_unprimed_score.npy`, sha `e5597327…`).
kinaseSubMain (`KSM/` = `/cluster/project/beltrao/JvG/Projects/hackathonKinaseSubs/kinaseSubMain`) was only read.
Empirical checks use DEV folds 0–4 only, through `kinact.explore.lib_variants.DevData` / `Evaluator`. Nothing under
`cache/sealed`, and no git history of `results/ground_truth`, was opened.
Scripts: `scripts/audit_r2np_leakage/a{1,3,4,5,6}.py` (outputs to `cache/audit_r2np/`).

**Verdict: SAFE with respect to response leakage. One residual risk, shared with every PSP baseline (§3).**
No code path gives any of the three libraries log2FC, curve or drug-target information. The high decryptM
number (D = 0.89–0.91) is explained empirically by (a) residue class and (b) PSP knowledge. Neither is leakage.

## 1. Code-level provenance

### 1.1 r2np: `r2_final_nophospho_uobs__20261026-213950__s0`

**Run identity.** Manifest `KSM/results/runs/r2_final_nophospho_uobs__…/manifest.json`, git `d7b5658`:
- argv: `scripts/models/r2_final.py --target uobs --seeds 0 --set … groups=["F1","F0","F2a","F3","F4","F5","F7"]`.
- `clean=false` and `comp_know=false`. Neither matters, because the F6 groups are absent.
- 162 features. The names are listed in the fold-0 run's `importance_fold0.tsv`.

**Assembler.** `src/ksp/models/r2/assemble.py::_build` loads exactly the stores below for these groups. F6a, F6b and
F6c are the only branches that touch `phospho/*` or `f6a*`, and they are not executed.

| group | stores / inputs | uses compendium? | label-derived? |
|---|---|---|---|
| F1 (7) | `motif/cantley_{unprimed,primed_soft,primed_hard}_{score,pct_universe}` (own matrix, else imputed) + impute distance | **existence only**. `primed_hard`: a flank S/T/Y that is a U_obs site is scored as phospho. `primed_soft`: the flank weight is the fraction of 4 MS channels that observed the flank (PSP-MS, compendium, CPTAC, Ochoa; binary per channel). `pct_universe`: reference = U_ms sites. The query site's own status is never used. | no |
| F0 (90) | AlphaFold structure (18), ESM-2 650M site PCA-64 (PCA fitted on U_ms sites) | existence only (which sites define the PCA) | no |
| F2a (6) | learned PSSM, cross-fitted (`xfit f2a_sep`, `features/motif_learned.py`) | existence only. Background = U_ms windows; `_ph` variant adds a flank-is-U_obs channel. | yes: PSP L2 (DEV) |
| F3 (12+6) | PPI screens (BioPlex, HuRI, OpenCell, Buljan, HCM), RWR; `f3nb` = cross-fitted PPI neighbourhood of training substrates | no (STRING/BioGRID were excluded as leaky) | f3nb: yes, PSP |
| F4 (10) | localisation (HPA IF, OpenCell, HCM, Hein, COMPARTMENTS-*experiments*, combined), UniProt topology | no | no |
| F5 (3) | HPA RNA tissue, DepMap Chronos co-dependency | no | no |
| F7 (35) | Manning group/family, ESM-2 kinase-domain PCA-16, Cantley profile PCA over random U_ms reference sites | existence only (choice of reference sites) | no |

**Labels.** PSP L2 only: human plus mouse/rat mapped by 15-mer, DEV proteins only (`ksp.labels.load_dev`). No
external or perturbation-derived labels (`L2+ext` was rejected, DESIGN A33). Negatives are U_ms sites, masked by
KP/Par/Hamming using training labels.

**Model selection.**
- Tree size and learning rate were chosen on KSM fold-0 PSP-label CV M1 (model_r2.md §4).
- The nophospho run is a post-hoc ablation that reuses the frozen config (model_r2.md §6.2b). It was not chosen on
  benchmarKIN or Hijazi data; downstream results for it were computed afterwards and were worse than clean.
- In B1, our own k/γ/mode choice was made on DEV GT. That is selection optimism (+0.026, nested), not leakage.

### 1.2 R3 deploy

Inputs (`model_r3.md`, README.json):
- ESM-2 650M per-residue window embeddings
- an amino-acid embedding
- a protein token
- a flank "phospho-channel" bit (flank ∈ U_ms; the centre bit is forced to 0)
- kinase descriptors: Manning group/family, ESM, and a Cantley profile over U_ms reference sites
- Cantley-atlas distillation on U_ms sites

Labels: PSP L2 on DEV+TEST KSM proteins (KSM's split, not ours). No phospho or response features. Selection used
KSM fold-0 M1 on PSP.

### 1.3 Cantley

We use `cantley_unprimed_score`: the published atlas matrices applied to the 15-mer. It uses no U_ms, no labels and
no data. It is fully clean.

### 1.4 Remaining channels

**Site existence (U_ms).** U_ms includes the compendium (50k compendium-only sites, 23 % of DEV candidates). It
enters only as binary "observed somewhere": the candidate universe, training negatives, flank priming, percentile and
PCA references. It carries no sign, magnitude, drug or condition. Empirically, r2np and r3 are *depleted* for
compendium-only sites (§2.4).

**F6 data.** None of these KSM artefacts are reachable from the nophospho config: `f6b_*` (692 inhibitor conditions,
Hijazi kinhib, decryptM CurveCurator), `f6a` co-regulation, and `f6c` CPTAC.

**Our build.** `src/kinact/library.py::verify_predicted_file` asserts the audited sha256 values and the banned-variant
list.

## 2. Empirical checks (DEV, M1-sel, wz statistic, per-kinase macro)

### 2.1 Where the decryptM number comes from

DEV D = 14 kinases. Of the 216 D pairs, 199 are Y-kinase targets: dasatinib dose series in A431, A549 and K562 (ABL1,
SRC-family, EPHA2, KIT, LCK). Dasatinib lowers pY globally. As a result, *any* library that gives Y kinases pY sites
scores Y targets above the mostly-ST negative kinases.

| library (measured top-k) | D | D, Y targets, **within-Y negatives** | H | H within-Y | O | O within-Y |
|---|---|---|---|---|---|---|
| null: random sites of the kinase's residue class @200 | 0.704 | 0.500 | 0.491 | 0.500 | 0.467 | 0.500 |
| null: site popularity (# DEV conds measured) @50 | 0.783 | 0.500 | 0.506 | 0.500 | 0.481 | 0.500 |
| null: popularity (# MS sources, existence) @50 | 0.680 | 0.500 | 0.519 | 0.500 | 0.480 | 0.499 |
| cantley @50 | 0.719 | 0.689 | 0.524 | 0.470 | 0.571 | 0.549 |
| psp_all (static) | 0.653 | 0.741 | 0.546 | 0.479 | 0.528 | 0.656 |
| r3 @50 | 0.881 | 0.734 | 0.531 | 0.405 | 0.597 | 0.664 |
| **r2np @50** | **0.908** | 0.809 | 0.575 | 0.486 | 0.655 | 0.614 |

Pooled over the dasatinib Y targets, the residue-class random null alone reaches AUROC 0.88 / 0.92 / 0.95 (A549 /
A431 / K562). r2np reaches 0.96–0.97.

The strongest evidence against response leakage comes from the second dasatinib study, which sits in O (PMID
19651622, K562, 50 nM). There, the residue-class null gives 0.06 and r2np gives 0.08–0.18 for the same 8 kinases.
cantley gives 0.21–0.34 and PSP about 0.5. r2np's Y-kinase lists therefore move with generic pY behaviour, not with
decryptM-specific responses: a library carrying decryptM dasatinib response information would not collapse on
another dasatinib experiment.

Existence and popularity nulls give exactly 0.50 within class. Site existence alone gives no discrimination between
kinases.

### 2.2 Contribution of PSP knowledge

We removed each kinase's own curated (curated_union) substrates from the r2np measured top-k. Only 3.1 % of top-50
edges are affected.
- M1-sel drops 0.663 → 0.646, and D drops 0.908 → **0.804**. Cantley without its curated edges: D 0.679.
- The largest D drops are EGFR (0.57 → 0.07) and MAP2K2 (0.98 → 0.34).

Most of r2np's D edge over cantley on the ST and receptor targets is therefore PSP substrates that r2np memorised from
training. That is legitimate: it is the same information as the PSP library.

Within-Y in D, r2np exceeds psp_all by +0.07. Most of this comes from A549 (0.80 vs 0.68); in K562 r2np is below PSP
(0.71 vs 0.73). This is based on 11 kinases and effectively 3 cell lines (one drug), so it is not a robust signal in
either direction.

### 2.3 Are r2np top-k sites enriched for responsive sites beyond PSP, specifically in decryptM?

Statistic: the fraction of a kinase's measured top-50 sites (curated edges removed) that are among the condition's
5 % most-changed sites in the GT sign direction. We compare targets with non-target kinases of the same residue class,
per source, macro over kinases (`a4.py`).

| enrichment ratio, target / non-target | D | H | O |
|---|---|---|---|
| cantley, non-curated | 1.11 | 1.00 | 1.28 |
| r2np, non-curated | 1.13 | 1.03 | 1.52 |
| r2np / cantley | **1.02** | 1.03 | 1.18 |
| r2np / cantley, curated included | 1.07 | 1.08 | 1.22 |

r2np's response enrichment over the label-free atlas is **smallest in decryptM** and largest in "other" studies. A
decryptM leak would show the opposite pattern.

### 2.4 Existence-only effects

Among each kinase's global top-200 DEV sites, the compendium-only fraction is:

| library | compendium-only fraction (top-200) | candidate base rate |
|---|---|---|
| r2np | 0.095 | 0.236 |
| r3 | 0.064 | 0.236 |
| cantley | 0.177 | 0.241 |

r2np and r3 favour well-studied sites (mean # MS sources 3.2 vs 2.7), not compendium-observed ones.

Two more checks:
- Measured-mode top-k (restricted to sites measured in the condition) behaves the same as global (B1 §3).
- The B1 shuffled-site control collapses to 0.48–0.52.

## 3. Verdict and residual risk

**SAFE.** No response, curve, drug-target or co-regulation information enters r2np, r3 or cantley. The decryptM
excess is (i) dasatinib/pY composition, which label-free nulls reproduce at 0.70–0.95, and (ii) memorised PSP
substrates (−0.10 D when removed). The response-enrichment advantage over cantley is smallest on D (×1.02).

Residual risks:
1. **Site existence.** U_ms includes compendium sites. It carries no response information; the within-class nulls
   score exactly 0.50. Estimated effect ≈ 0.
2. **PSP curation.** PSP KS rows may include relations curated from compendium papers. They carry no PMIDs, so this
   cannot be filtered. The risk is shared with psp_all and every PSP baseline, so it does not favour r2np relative to
   them. The upper bound on r2np's PSP-attributable share is about ΔM1-sel 0.017 overall (0.663 → 0.646 when curated
   edges are removed).
3. **Reporting.** D-source numbers should not be presented as evidence of method quality. Report
   residue-class-null-adjusted or within-class AUROC for D. The D stratum is about one drug (dasatinib) and has
   effectively 3 cell lines.

Not leakage, but noted:
- r3 and the r2np final model have seen all of their own training PSP labels. Their scores on PSP-labelled sites are
  in-sample.
- B1 config selection optimism is +0.026.
