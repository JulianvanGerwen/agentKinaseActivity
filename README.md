# Agentic AI for bioinformatics: a kinase activity inference case study

This repository is an **exercise in testing agentic AI for bioinformatics research**. A single human-written brief
(`PROMPT.md`, verbatim) asked an AI agent (Claude, acting as "PI") to run a complete computational biology
project with a team of AI sub-agents:
- literature review;
- study design with independent agentic review rounds;
- ground-truth construction and auditing;
- leakage-safe train/validation/test splits;
- baseline re-implementations and method development;
- a sealed one-shot test evaluation, independent reproduction and a final report.

Apart from the brief, the provided data and occasional prompts to continue, the code, analyses, decisions and
documents in this repository were produced by AI agents. They are shared to show what such a workflow produces,
including its records, safeguards and failure modes. They have **not been peer reviewed**. Read the results as
the output of an experiment in agentic research, not as a validated published method.

## The scientific task

The task was to estimate changes in human protein kinase activity from phosphoproteomics data (site-level log2 fold
changes). The goal was a method that is as accurate and fast as possible, covers as many kinases as possible
(including understudied ones), and beats z-score, KSEA and JUMPsem without data leakage or study bias.

## Outcome in one paragraph

The winning method, **core-v3**, has three parts:
- a Bayesian linear-Gaussian model over a curated plus predicted kinase–substrate library, with a kinase-family prior;
- late-fused perturbation-response signatures;
- per-kinase calibration.

On a sealed test set (468 conditions, 67 kinases), it reached a coverage-aware AUROC of 0.714. The comparators
scored 0.584 (z-score), 0.582 (KSEA), 0.549 (PhosX-like) and 0.505 (JUMPsem); each difference has Holm p = 0.003. It
scores 503 kinases in about 7 ms per condition. It did **not** beat a trivial "drug-target popularity" control
(0.712), which wins on well-studied kinases but fails on unseen ones. The full results, limitations and remaining
ideas are in [`REPORT.md`](REPORT.md).

## Run the method yourself

See [`mvp/README.md`](mvp/README.md). It is a minimal, self-contained package with the frozen core-v3 model and a
script that builds the kinase–substrate library from your own PhosphoSitePlus download.

```bash
pip install -r mvp/requirements.txt
python mvp/build_edges.py --psp-ks Kinase_Substrate_Dataset.gz --psp-sites Phosphorylation_site_dataset.gz \
    --omnipath omnipath_enzsub.tsv
cd mvp && python -m kinact_mvp score example/toy_input.tsv -o example/toy_scores.tsv
```

## What to look at (the agentic-research record)

| file | what it is |
|---|---|
| `PROMPT.md` | the human brief, verbatim; the only task specification the agents received |
| `RESEARCH_TREE.md` | the PI agent's research tree: status, decision log (D0–D14), every branch with outcome |
| `LITERATURE.md` | agent-written literature review |
| `DESIGN.md`, `docs/reviews/` | design and the three rounds of independent agentic review |
| `DATA.md` | every dataset: source, version, checksum, licence |
| `docs/audit/` | independent agent audits (ground truth, harness, leakage, method adoption, test runner) |
| `LESSONS.md`, `docs/retro/` | process lessons and milestone retrospectives written by fresh agents |
| `FREEZE.json`, `results/test/` | pre-registered test freeze (hashes pinned to a commit) and the one-shot test results |
| `REPORT.md` | final report, independently fact-checked by another agent |
| `src/kinact/`, `scripts/`, `tests/` | research code: evaluation harness, baselines, method branches, unit tests |
| `results/` | committed result tables per branch (`b1`…`b8`, baselines, descriptive, test) |
| `mvp/` | runnable package of the final method |

## Data and licences

No raw data are included. The phosphoproteomics compendium, PhosphoSitePlus files, publications and third-party
downloads are git-ignored. Their sources, checksums and licences are listed in `DATA.md`.
- **PhosphoSitePlus:** CC BY-NC-SA 3.0. PSP-derived kinase–substrate edges are not distributed; the MVP rebuilds
  them locally.
- **Committed tables:** small derived tables only (ground-truth labels, evaluation results, perturbation signatures,
  predicted substrate edges).
- **Research code:** reproducing the full study needs the original inputs, which are not included. The research code
  expects them under `data/`, `external/` and `cache/`, and the sibling kinaseSubMain checkout via the
  `KINASESUB_DIR` environment variable.

## Caveats about the agentic process

- **Corrected errors:** the agents made and caught errors, including a sign error in the metric, mislabelled ground
  truth, a bootstrap bias and overstated claims in the report draft. These are documented in `RESEARCH_TREE.md`,
  `docs/audit/` and `LESSONS.md`. Undetected errors may remain.
- **Ground-truth judgements:** ground-truth construction involved judgement calls by agents, for example inhibitor
  dose defaults and target sets. They are listed in `results/ground_truth/` and the audits.
