# Usage — scoring kinase activities with core-v1

core-v1 is the adopted model (RESEARCH_TREE D6/D7). It is a Bayesian linear-Gaussian model:
- library: the curated_union library plus r2np top-100 predicted edges (α = 1);
- prior: a Manning family/group kernel with τ = 0.125, ρ_f = 3, ρ_g = 0.25;
- output: z_vif, calibrated with ρ_e = 0.01. There is no winsorising.

All parameters are read from `config/core_v1.json`. The model artefacts are in `models/core_v1/`:
- `edges.tsv.gz`: the kinase–substrate edges;
- `kinases.tsv`: the kinase universe with Manning group and family;
- `manifest.json`: hashes and the ρ_e refit.

`scripts/build_core_artifacts.py` rebuilds the artefacts. It uses only label-free inputs and the 7 DEV null
conditions, and it checks that its ρ_e refit reproduces the config value. The implementation is in `src/kinact/core.py`.
It reproduces `kinact.explore.bayes2` (see "Equivalence and speed" below).

## Input
You give one or more conditions of site-level directional statistics, such as log2FC or t/z values. Sites are keyed
by `{UniProt}_{residue}{position}`, for example `P31749_S473`. Three layouts are accepted:
- **wide**: rows are sites and columns are conditions. In a file, the first column holds the site keys.
- **long**: columns `uniprot`, `residue`, `position`, `value`, plus an optional `condition` column.
- **Series**: a pandas Series indexed by site key, which counts as one condition. The Python API only.

How values are handled:
- Missing values are NaN or empty. A site that is not measured in a condition is simply absent from that condition's
  model. Nothing is imputed.
- Duplicate site keys are averaged.
- Keys are matched after stripping whitespace and upper-casing.
- Every finite value enters its condition's centre and scale: median-centring and the robust "wsd" scale.
- Only sites that are in the library enter the likelihood. Unknown and malformed keys are counted in
  `result.info["n_unknown_sites"]`.
- If no site matches the library, you get a warning and every kinase is NaN.

## Output (per condition × kinase)
| column | meaning |
|---|---|
| `z` | z_vif: a calibrated activity-change z-score. Positive means activated; 0 means no evidence. |
| `mean` | Posterior mean of the activity change, in units of the standardised input. |
| `sd` | Posterior SD. |
| `n_substrates` | Number of measured library substrates. |

The NaN contract:
- A kinase is scored if it has at least 1 measured substrate.
- It is also scored if a direct family or group relative has measured substrates. It then borrows strength from that
  relative and has `n_substrates = 0`.
- Every other kinase is **NaN** in `z`, `mean` and `sd`, never 0. `n_substrates` is always filled in.

## Command line
```bash
pip install -e .            # registers the `kinact` entry point (torch is required for the fast path / GPU)
kinact score input.tsv -o scores.tsv                    # long output: condition, kinase, z, mean, sd, n_substrates
kinact score input.tsv --format wide > scores_wide.tsv  # kinase rows; columns z:<cond>, mean:<cond>, ...
kinact score long_input.tsv --device cuda --drop-unscored -o scores.tsv
```
Options:
- `--threads N` sets the number of CPU threads. The default is the BLAS thread limit, so `OMP_NUM_THREADS` and
  `threadpoolctl` apply.
- `--config` selects another model config.
- `-q` suppresses the summary line that is printed to stderr.

Example input (`input.tsv`, wide):
```
site	ctrl_vs_EGF_5min	ctrl_vs_EGF_30min
P28482_T185	2.31	1.10
P28482_Y187	2.05	0.84
Q15418_S380	1.62	0.97
P31749_S473	0.40	1.25
...
```
Example output (long; the values are illustrative):
```
condition	kinase	z	mean	sd	n_substrates
ctrl_vs_EGF_5min	AAK1	1.04419	0.0542075	0.111697	7
ctrl_vs_EGF_5min	ABL1	1.192	0.0457979	0.11557	7
...
```

## Python API
```python
from kinact.core import KinaseActivityModel

model = KinaseActivityModel.from_config()        # config/core_v1.json (+ models/core_v1 artefacts)
res = model.score(df)                            # df: sites x conditions; or a Series; or a long table
res.z            # conditions x kinases DataFrame (z_vif, NaN = unscored)
res.mean, res.sd, res.n_substrates
res.long()       # tidy table; res.wide() -> kinase rows
res.info         # counts: input / library-matched / unknown sites, measured-site patterns, Woodbury updates

res_gpu = model.score(df, device="cuda")         # batched float64 GPU path (same numbers to ~1e-11)
state = model.prepare(df); res = model.score(state, n_jobs=16)   # split load/mapping from scoring (benchmarks)
```

## Equivalence and speed (DEV, `scripts/bench_core.py` → `results/core_v1/{equivalence,runtime}.json`)
The following is measured on the 200 fixed DEV bench conditions (seed 0) and on the full 1,358-condition DEV
compendium. Hardware: AMD EPYC 9554 and an RTX 4090.

Equivalence, on the 200 bench conditions:
- **vs the bayes2 float64 single-sample path:** max |Δz| = 8e-12, with an identical NaN pattern.
- **vs the bayes2 batched GPU scorer:** max |Δz| = 1.7e-6 (2.5e-6 over all 1,358 DEV conditions). This equals bayes2's
  own float32-vs-float64 difference.
- **Ranking:** per-condition Spearman is 1.0 once exact ties between symmetric borrowed kinases are resolved.

Runtime is reported, not gated (DESIGN §1). The numbers are written by `scripts/bench_core.py`; see `runtime.json`.
All figures are end-to-end ms per condition, with load/mapping plus scoring in brackets.

| setting | 200 bench conditions | full DEV compendium (1,358) |
|---|---|---|
| CPU, 1 thread | 6.6–6.8 (0.9 + 5.7) | 6.2 s |
| CPU, 16 threads, idle machine | 2.1 (1.0 + 1.1) | 2.2 s |
| CPU, 16 threads, load average ≈ 100 | 6.1–6.4 | 6.0 s |
| GPU | 1.9–2.6 | 2.1–2.8 s |

The targets are ≤ 5 ms per condition and ≤ 10 s for the compendium:
- The compendium target is met on every path.
- The 5 ms target is met single-threaded only on the full compendium (4.6 ms), where 395 of 958 patterns get a
  Woodbury update. On the random 200-condition set few patterns are close, so roughly one p≈500 SPD inverse per
  condition is the floor.

How the implementation gets its speed:
- Conditions that share a measured-site pattern share one factorisation.
- When a pattern is close to an already factorised one (≤ p/4 differing sites), it gets an exact Woodbury update of
  that inverse instead. This path is used single-threaded only.
- G = W′W is assembled from precomputed per-site kinase pairs.
- diag(Λ⁻¹P₀Λ⁻¹) uses the low-rank (family/group indicator) form of the prior precision.
- Patterns are spread over a thread pool, using the torch CPU float64 path.
- The GPU path is a batched float64 Cholesky.

## core-v3 (`config/core_v3.json`; RESEARCH_TREE D11/D12)
core-v3 = core-v1 + late-fused perturbation-response signatures (B7b) + the E1 per-kinase permutation-null scale (B8).
```python
model = KinaseActivityModel.from_config("config/core_v3.json")   # artefacts in models/core_v3/
res = model.score(df)
res.z        # fused, calibrated z (primary output)
res.z_core   # core-v1 z_vif (identical to config/core_v1.json)
res.z_sig    # standardised signature z (NaN: kinase without signature or < 5 measured signature sites)
```
CLI: `kinact score input.tsv --config config/core_v3.json -o scores.tsv` (long output gains `z_core`, `z_sig`).

For a kinase k with a signature (55 kinases) and condition c:
- `z_sig = [Σ_s w_ks clip(y_cs, ±4) / sqrt(Σ_s w_ks²)] / ρ_sig` over k's measured signature sites (≥ 5, else NaN).
  `y` is the same wsd-standardised input as the core; signature sites need not be library sites.
- `z = (z_core + α z_sig) / sqrt(v_k)` when both are defined (α = 1); `z = z_core` when only the core scores k;
  `z = z_sig / sqrt(v_sig,k)` when only the signature does. Every kinase without a signature keeps `z = z_core`.
- `mean`, `sd` and `n_substrates` stay the core posterior quantities.

Artefacts (`scripts/build_core_v3_artifacts.py`; DEV-CV folds 0–4 only, nothing from DEV-confirm or TEST):
| file | content | git |
|---|---|---|
| `edges.tsv.gz`, `kinases.tsv` | the core-v1 library (PSP-derived edges) | edges ignored |
| `signatures.tsv.gz` | kinase, site, signed weight, rank, n_groups (clean+contrast fit on all DEV-CV GT conditions; C10) | committed |
| `e1_scales.tsv` | v_k (fused) and v_sig,k from site-permuted DEV-CV GT conditions (seed 0) | committed |
| `sig_selection.tsv` | the (n_top, κ) freeze table (below) | committed |
| `signature_training_conditions.txt` | the 463 training conditions | committed |
| `manifest.json` | sha256 of every artefact and input, ρ_e / ρ_sig refits | committed |

Frozen values are taken from the DEV-confirm look-2 freeze (`results/g2_prep/devconfirm_look2_frozen.json`):
- α = 1;
- (n_top, κ) = (100, 0), the argmax of the pooled OOF M1-sel over the B7 primary grid;
- ρ_sig = RMS of z_sig over the 7 DEV nulls.

The build re-derives ρ_sig, the median of v_k, the 463 training conditions and the 55 signature kinases, and fails
if any of them disagrees with that record. `sig_selection.tsv` (late-fusion OOF per (n_top, κ), 0.676–0.683) is
descriptive only. As a reproduction anchor, it recovers the B7b per-fold-pick OOF of 0.681433 exactly.

Equivalence (`tests/test_core_v3.py`, `scripts/bench_core_v3.py` → `results/core_v3/`):
- vs the research path (bayes2 z_vif + `signatures_late.sig_project` + `final_box.fuse_calibrated`, float32 input as
  in B7b/B8): max |Δz| = 1.7e-6 on the 200 bench conditions and 2.5e-6 on all 1,358 DEV conditions; identical NaN
  pattern;
- the signatures, ρ_sig and E1 scales are refitted independently with the research functions in the test;
- GPU vs CPU: 6e-11.

Runtime (`results/core_v3/runtime.json`, final artefacts; numbers vary with machine load):
| setting | 200 bench conditions | DEV compendium (1,358) |
|---|---|---|
| CPU, 1 thread | 6.8 ms/cond | 6.4 s (4.7 ms/cond; 1,937 extrapolated 9.1 s) |
| CPU, 16 threads | 2.0 ms/cond | 2.1 s |
| GPU | 1.8 ms/cond | 2.0 s |
The fusion step itself costs about 1.1 ms per 200 conditions; the E1 scale is precomputed.
