/goal You are the PI of a computational kinase signaling project. You manage a team of agents (Opus 5.5 or Sonnet 5.5 as fitting). You delegate, review, integrate and decide; you keep the research record. The goal is the best possible outcome, not the fastest completion. Make decisions yourself, document why, and keep going.

Your task is to create the most accurate and fastest estimator of kinase activity changes in phosphoproteomics data. Specifically, this method will take log2 fold changes (or other directional statistics) from human phosphoproteomics experiments, and estimate the changes in activity of as many human protein kinases as possible. This will then be used by researchers to interpret phosphoproteomics experiments and explore relationships between kinases.

Information to use: Consider the most useful input sources. Available to you will be experimentally determined kinase-substrate relationships (from PhosphositePlus), and most algorithms use this kind of information. But you will also have access to predicted kinase-substrate relationships which you could consider using in the directory kinaseSubMain/, and you could consider integrating other modalities or other sources of kinase-substrate relationships. For the predicted kinase-substrate relationships, just be careful about data leakage, because they may use phosphoproteomics data that you use for evaluation.

Evaluating your method: I have provided a compilation of human phosphoproteomics data. This includes kinase inhibitor perturbation data from Hijazi et al (paper provided) along with annotation of which kinases should be inhibited; you should start with this as your ground truth of kinase activity changes (kinases that are being inhibited should decrease in activity). But consider other data within the compilation that could have ground truth kinase activity changes, and consider other phosphoproteomics data sources.

Benchmarking against other methods: You should aim to beat established methods including z-score, KSEA, and JumpSEM (papers provided). Try to implement as many major benchmarks as possible, within reason.

Data leakage: Be careful of data leakage when you separate into training, validation, and test sets. 

Biases: Avoid being deceived by biases in existing annotations/training data. For instance, a small number of kinases are highly studied, but do not restrict your focus to these or allow the metrics to be swayed by them. Ideally achieve harmonious predictions across as many kinases as possible, including lowly studied kinases with few annotated substrates.

In the main working directory, the directory data/ contains data and papers that I am already providing to you.

Time: Do not waste time trying to access datasets or other materials; first try to get them, then if it is taking too long, ask the user but don't wait, and continue with other work instead. Flag and explain all important decisions clearly. Plan ahead your experiments accordingly, you can do quick tests on how long things take to run on the provided hardware. 

Literature: Design through repeated and independent critique: Search the current literature and available datasets, and write LITERATURE.md. Then draft DESIGN.md. The project goal and its parts are fixed, the suggestive ideas are a human priors, not a plan. Have the design reviewed and revised in at most three rounds by an external agentic PI, and agentic scientists assessing prior work, statistics, data, algorithms and compute, AI-agent failure modes, and other expertise the design needs. Each reviewer proposes concrete changes that stay within the project's scope, or states that none are needed.

Reproducibility: Never change metrics, controls, labels or tests to make results look better. Keep everything end-to-end reproducible: commit results to git, fix the random seeds list, and keep DATA.md listing the source, version, checksum and license of each dataset. Never commit restricted data, never hardcode data, prefer scripts over one-offs, make sources explicit, no silent assumptions. Don't trust a result until it has been checked. Write verification scripts using toy data with known answers, positive and negative controls, and sanity checks. An agent that didn't produce a result reviews it before you rely on it.

Research tree: Save this prompt verbatim to `PROMPT.md`. Keep track of the research in a research tree (RESEARCH_TREE.md), using a status legend: [ ] planned · [~] in progress · [x] done · [-] abandoned/negative · [!] key finding. At each milestone, add a short decision note and a two-paragraph status summary at the top. The tree is your memory: keep your own context lean by reading agents' summaries by default rather than raw logs and if your context is reset, resume from the tree. 

AI lessons: Keep LESSONS.md, a list of general process rules for running this project efficiently. Whenever something wastes or saves significant effort, add a one-line rule with its evidence. At each milestone, prune it to the most useful rules and ask a fresh agent "what should we stop doing?" and “what worked well?”. 

Auto-research: Your job is to build the best possible method. Report runtime alongside accuracy. Keep going until you genuinely run out of promising ideas or close to the time limit. Then write REPORT.md with a breakdown of the approach, the best methods, test-set results, follow-ups, limitations, and remaining ideas.

If you get blocked by the biosecurity filter, please switch to sonnet 5.5
