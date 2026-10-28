# kinact MVP: core-v3 kinase activity inference

A self-contained, minimal package for the project's winning method, **core-v3**. It takes site-level directional
statistics from a human phosphoproteomics experiment (log2FC, t- or z-values) and returns a calibrated
activity-change z-score for each of 503 human protein kinases. See `../REPORT.md` for the method and its evaluation.

> Research code produced as part of an exercise in agentic AI for bioinformatics (see the top-level README).
> Use it for exploration, not for clinical or diagnostic decisions.

## 1. Install

Requires Python ≥ 3.9.

```bash
pip install -r mvp/requirements.txt
```

`torch` is optional: it is about 3× faster on CPU and enables `--device cuda`. Without it, a pure numpy/scipy path is used.

## 2. Build the kinase–substrate library (one time)

The library combines curated PhosphoSitePlus edges with predicted edges. PhosphoSitePlus (PSP) data are
CC BY-NC-SA 3.0 and need a free registration, so PSP-derived edges are **not shipped**. Build them from your own
download:

1. Download `Kinase_Substrate_Dataset.gz` and `Phosphorylation_site_dataset.gz` from
   <https://www.phosphosite.org/staticDownloads>.
2. Optionally (recommended) download the curated OmniPath enzyme–substrate table:

   ```bash
   curl -o omnipath_enzsub.tsv \
     "https://omnipathdb.org/enzsub?genesymbols=yes&fields=sources,references&format=tsv&license=academic"
   ```

3. Build the library. It writes `mvp/model/edges.tsv.gz` in a few seconds.

   ```bash
   python mvp/build_edges.py \
       --psp-ks Kinase_Substrate_Dataset.gz \
       --psp-sites Phosphorylation_site_dataset.gz \
       --omnipath omnipath_enzsub.tsv
   ```

With PSP 2025-05 and OmniPath 2026-10, the rebuilt library shares 97% of its edges with the one evaluated in REPORT.md.
The difference comes from iPTMnet and from the evaluation-only leakage filter, which the MVP omits.

## 3. Score

Run from the `mvp/` directory, or put `mvp/` on `PYTHONPATH`.

```bash
cd mvp
python -m kinact_mvp score example/toy_input.tsv -o example/toy_scores.tsv
python -m kinact_mvp score my_data.tsv --format wide --drop-unscored -o scores_wide.tsv
python -m kinact_mvp score my_data.tsv --device cuda -o scores.tsv        # with torch + GPU
```

From Python:

```python
import sys; sys.path.insert(0, "mvp")
from kinact_mvp import KinaseActivityModel
model = KinaseActivityModel.from_config()          # mvp/model/config.json (core-v3)
res = model.score(df)                              # df: rows = site keys, columns = conditions
res.long()                                         # condition, kinase, z, mean, sd, n_substrates, z_core, z_sig
```

### Input

Sites are keyed `{UniProt}_{residue}{position}`, for example `P31749_S473`. Two layouts are accepted:
- **wide**: the first column holds the site keys and each further column is one condition;
- **long**: columns `uniprot`, `residue`, `position`, `value`, plus an optional `condition` column.

Use the *whole* measured phosphoproteome, not only significant sites. Every value enters the per-condition
centring and scaling. Missing values are simply absent; nothing is imputed.

### Output

| column | meaning |
|---|---|
| `z` | core-v3 activity-change z-score; negative = decreased activity |
| `z_core` | z-score from the Bayesian library model alone |
| `z_sig` | perturbation-signature component, for the 55 kinases with a signature |
| `mean`, `sd` | posterior mean and SD of the core activity change |
| `n_substrates` | measured library substrates |

Unscored kinases are NaN. A kinase is scored if it has at least one measured substrate, or a related kinase
(family/group) with measured substrates.

## 4. How to read the scores (from the sealed TEST evaluation)

- **Ranking:** core-v3 ranks the truly perturbed kinase above others better than z-score, KSEA, JUMPsem and a
  PhosX-like method (TEST AUROC 0.71 vs 0.50–0.58).
- **Calibration:** |z| > 1.96 occurs in about 5% of kinases on permuted data but about 8% on real control-vs-control
  data. Treat |z| ≈ 2 as suggestive.
- **Related kinases:** calls between closely related kinases (paralogs, CDK/MAPK, AGC) partly reflect shared
  substrates.
- **Understudied kinases:** these are scored through predicted substrates and family borrowing. They are useful for
  generating hypotheses but less reliable.
- **Speed:** about 7 ms per condition on one CPU thread with torch, after a one-off model set-up of a few seconds.

## Contents

| path | content |
|---|---|
| `kinact_mvp/core.py` | model, a copy of `src/kinact/core.py` (only default paths changed) |
| `kinact_mvp/cli.py` | command line interface |
| `build_edges.py` | builds `model/edges.tsv.gz` from your PSP (+ OmniPath) download |
| `model/config.json` | frozen core-v3 hyperparameters (from `config/core_v3.json`) |
| `model/kinases.tsv` | 503-kinase universe with Manning group and family |
| `model/predicted_edges.tsv.gz` | top-100 predicted substrates per kinase (kinaseSubMain "R2 no-phospho" model) |
| `model/signatures.tsv.gz`, `e1_scales.tsv`, `sig_selection.tsv` | perturbation-response signatures and their calibration |
| `example/toy_input.tsv` | synthetic example: random values with AKT substrates lowered in `toy_AKT_down` |
