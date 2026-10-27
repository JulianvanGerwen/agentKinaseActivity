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
