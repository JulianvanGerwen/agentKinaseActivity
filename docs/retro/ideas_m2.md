# Untested ideas after M1 method work: ranking and go/no-go (fresh reviewer, 2026-10-27)

Sources: RESEARCH_TREE.md, DESIGN §7/§12, LITERATURE §5, results/b*/b*_report.md, results/baselines/baselines_dev.md. Sealed data and git history were not read. No experiments were run.

Status: core-v1 has DEV nested M1-sel 0.663. The best baseline is 0.563 and the prior-only control is 0.687.
Weak spots: Hijazi (H ≈ 0.63), harmony bin ≥20 (≈ 0.645, where prior-only scores 0.84) and bin 0 (≈ 0.58).

## Already tested (do not re-propose)
- Library variants (B1) and r3/Cantley ensembles: flat. The similarity kernel: flat. The site-noise model: −0.008.
- Winsorising: +0.004. Posterior-mean output: not robust.
- Co-regulation edges (B4): −0.010. Learned edge weights from evidence/percentile/promiscuity features (B4b): nested θ = 0.
- Self-sites (B5): +0.003. Stacking (B6): A1 +0.012, p = .057. The learned stacker: −0.001.
- Signatures inside the kernel (B7): fail T-unseen. B7b late fusion is pre-declared and running; its result is not counted here.
- Lessons from these: every gain is ≤ 0.015, and grids larger than about 3 hyperparameters lose about 0.03 to selection noise.

## Ranked ideas (expected gain × P(passes hard test) / cost)

### 1. Remove the generic perturbation response from the input before scoring
- Hypothesis: many sites move in almost any perturbation (stress, cell-cycle, shared MAPK/mTOR feedback). They inflate many kinases at once, which hurts Hijazi especially (multi-target drugs in one cell background).
  B7 found that the generic-response contrast improved signatures by +0.005 to +0.009. It has never been applied to the core input.
- Experiment: y' = y − P_r y. P_r projects onto the top-r right singular vectors of the standardised site × condition matrix of ALLOWED conditions (DESIGN §4: non-TEST, not sharing a control with TEST, held-out fold excluded inside CV).
  - Hyperparameters: r ∈ {0, 1, 3, 5} and shrink factor η ∈ {0.5, 1}. Everything else is core-v1.
  - Selection is nested over 5 folds. Refit ρ_e on the DEV nulls.
  - The model stays single-sample: the projection is fixed preprocessing.
- Expected: H +0.01 to +0.02, M1-sel +0.005 to +0.012, P(hard test) ≈ 0.3.
- Risk: if the PCs are dominated by frequent-target responses (ERK/mTOR), projecting them out removes true signal for bin ≥20 kinases. Run a pre-declared variant with PCs fitted on non-GT conditions only, as a check and not as an extra hyperparameter.
  The generic vector must be refit per outer fold.
- Agent time: about 2 h.

### 2. Per-kinase robust compendium scale (label-free kinase-wise null)
- Hypothesis: calibration worsens with substrate count (|z| > 1.96 rises from 1.2% to 7.6% by bin). Hub kinases (CDKs, CK2, PKA, ERK) get large |z| in every condition, so as negatives they outrank true targets.
  Dividing z_k by a robust spread of z_k over allowed conditions, s_k = 1.4826·MAD_c(z_kc), makes kinases comparable. It also fixes desideratum 3.
- Experiment: z'_kc = z_kc / s_k^λ with λ ∈ {0, 0.5, 1}, nested. s_k is computed on allowed conditions only, excluding the held-out fold. No labels are used.
- Expected: the clearest effect is calibration flattening across bins. M1-sel ±0.01, with gains for bins 5–19/≥20 and Hijazi. P ≈ 0.25.
- Risk: s_k is larger for frequently perturbed kinases (true signal enters the MAD), which could cost bin ≥20.
  It is not a popularity prior (it is anti-popularity and uses no labels), but say so explicitly in the D5 audit.
- Agent time: about 1 h. It reuses core-v1 OOF z.

### 3. Hyperparameter bagging instead of argmax (reduce selection noise)
- Hypothesis: per-fold nested picks vary (for example ρ_f ∈ {3, 10}, τ), and B2.1 showed that selection noise costs about 0.03. Averaging z over the inner top-3 configs per fold (equal weights) lowers variance without adding hyperparameters.
- Experiment: from the existing bayes2 grid stores, z_bag = mean of inner-top-3 z_vif, then a recalibrated ρ_e. The only hyperparameter is n_bag ∈ {1, 3, 5}, nested.
- Expected: +0.003 to +0.008 overall, no particular weak spot. P ≈ 0.3. It also lowers the DEV→TEST drop.
- Risk: minimal leakage, since inner scores already exist. The gain may be below detectability.
- Agent time: about 1 h.

### 4. Sparse (ARD/Laplace) prior on activities via 3-step reweighted ridge
- Hypothesis: few kinases change per condition. The Gaussian family kernel spreads signal across paralogs and motif-sharing hubs (proline-directed CDK/MAPK, basophilic AGC), which hurts bin ≥20 and Hijazi multi-target ranking.
  Iterative reweighting τ_k² ← τ²·(|μ_k| + ε)/mean|μ| gives explaining-away with sparsity. The closed form stays at about 3× the runtime.
- Experiment: sparsity strength γ ∈ {0, 0.5, 1} as an exponent on the reweighting, 3 fixed iterations, core-v1 otherwise. Nested, with ρ_e refit.
- Expected: H and bin ≥20 +0.01. P ≈ 0.2.
- Risk: T-unseen and bin 0 may drop because weakly supported kinases are shrunk to 0 (check the T-unseen −0.01 margin). Calibration of reweighted z is approximate.
- Agent time: about 3 h.

### 5. Score-threshold predicted library for poorly annotated kinases
- Hypothesis: a fixed r2np top-100 gives low-confidence kinases (bin 0) as many edges as confident ones. Keeping edges above a global r2np score percentile (adaptive k per kinase, capped at 200) gives more signal per edge.
- Experiment: percentile threshold p ∈ {99, 99.5, 99.8} (one hyperparameter) replaces top-k. Nested.
- Expected: bin 0 / T-unseen +0.01 to +0.02, M1-sel +0.003. P ≈ 0.15. B1 showed top-k is flat, so the prior is weak.
- Risk: some kinases end up with fewer than 5 edges, so coverage drops (the ½ rule penalises this). Use a floor of 20 edges.
- Agent time: about 1 h.

## Not recommended
- A "gated A1" (co-regulation average only in bin 0): this is forking paths after seeing B6 by bin. Keep A1 as the pre-registered TEST secondary.
- Anything using target frequency or literature counts: banned by D5/C5.

## Verdict: near diminishing returns. Proceed to G2 after a bounded final round.
- Even an optimistic ΔM1-sel of +0.01 per idea is below what the DEV hard test reliably detects (bootstrap SE about 0.008) and far below the TEST MDE of about 0.06.
- None of these ideas plausibly closes the −0.024 gap to prior-only. That Holm comparison will probably fail on TEST whatever is done, so report it honestly.
- Core-v1 already beats every baseline by +0.10 and wins on T-unseen, D and calibration.
- Plan:
  - Let B7b finish.
  - Run at most ideas 1 and 3: pre-declared, ≤3 hyperparameters, a total time box of about 3 h, one combined DEV-confirm look.
  - Do the winner's-curse re-split (seeds 1–4) and freeze G2.
  - Ideas 2, 4 and 5 go to post-TEST exploration, labelled as such.
