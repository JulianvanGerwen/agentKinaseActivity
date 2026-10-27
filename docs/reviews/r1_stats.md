# Review r1 — statistics / benchmark design (DESIGN.md v0)

Reviewer: statistics agent, 2026-10-27. Evidence script (read-only, scratch): plain PSP z-score (n≥3 quantified substrates) on the 183 Hijazi GT conditions, U = 437 kinases (PSP ∪ Hijazi targets). This is a quick check, not a validated baseline: its absolute values are only indicative.

## 1. Verdict
**Major revisions.** The harness design is sound overall (frozen GT, grouped split, kinase-macro, cluster bootstrap, Holm). Five things need fixing before G1: a sign error in the primary metric, the unscored = 0 rule, the averaging hierarchy, Coverage@0.6, and the DEV comparison procedure.

## 2. Proposed changes

**1. [must] §5 sign convention.** "t = −s·score (larger = correct)" is wrong. For an inhibitor (s = −1) with a correctly negative score it gives t < 0. Change it to **t = s·score**. Add a harness unit test: an oracle that returns −1 for every inhibitor target must give M1 = 1.

**2. [must] §5 handling of unscored kinases in M1.** Replace "unscored ⇒ score 0 with tie-aware ranking" with a neutral contribution. Let S_c be the kinases scored in c and N_c the evaluable negatives:
- If k ∈ S_c: AUROC_ck = [ Σ_{j∈N_c∩S_c} (1[t_j<t_k] + ½·1[t_j=t_k]) + ½·|N_c \ S_c| ] / |N_c|
- If k ∉ S_c: AUROC_ck = ½

Why: under 0-imputation every unscored negative sits at 0. A scored target then beats the whole unscored block whenever its sign is right. So the metric depends on where each method's 0 lies and on the score scale (KSEA z, ulm t and viper NES differ), which mixes a sign test with ranking. In the check, the same z-scores gave M1 = 0.555 with 0-imputation and 0.512 with the neutral rule. Shifting every score by a constant changed the 0-imputed M1 (0.555 → 0.545 at offset −3), while the neutral rule stayed at 0.512. The neutral rule does not depend on shift or scale, gives exactly 0.5 for random or unscored output, and still penalises low coverage. Also report two parts: **M1-scored** (pairs with k ∈ S_c, negatives in N_c∩S_c only) and the **pair coverage** |{(c,k): k∈S_c}|/|GT|.

**3. [must] §5 averaging hierarchy (condition counts).** Average in this order: pair → (kinase, GT group) mean → kinase mean → macro over kinases. The current order is condition → kinase. Why: decryptM dose series (340 conditions), Hijazi's 3 cell lines × drug, and duplicated descriptors (88) would otherwise let replicate-like conditions dominate a kinase. In Hijazi, 108/306 target kinases come from a single drug (always 3 conditions). Use the same hierarchy in M2 and in the harmony bins.

**4. [must] §5 Coverage@0.6.** Under the null it passes 13% of kinases (random scores with 0-imputation, Hijazi). Replace it with: the number of GT kinases whose per-kinase mean AUROC (from change 3) is >0.5 with one-sided p<0.05, BH q≤0.10 across kinases. Get the null by permuting the target's rank within each condition (uniform rank among N_c ∪ {k}), keeping the group structure. Report the expected number of false discoveries next to the count. Unscored kinases (AUROC exactly ½) can never count.

**5. [must] §5 Inference: definitions and the superiority rule.** Pin these down:
- (a) Resampling unit: GT groups (the union-find components from §4), stratified by GT source. The whole metric is recomputed on each resample, including which kinases are present and the bin membership.
- (b) Pairing: every method is evaluated on the identical resample.
- (c) One-sided p = (1 + #{Δ* ≤ 0}) / (B + 1), with B = 2000 and seed 0.
- (d) The family is fixed: final vs {z-score, KSEA, JUMPsem, best other baseline} on M1, plus final vs z-score on M2. That is 5 tests, Holm-adjusted. "Best other baseline" is chosen by DEV-CV M1 and written into FREEZE.json at G2, never chosen on TEST.
- (e) "Beats" = Holm p<0.05 **and** Δ ≥ δ_min. Pre-register δ_min = 0.01 on M1 (≈0.02 on M1-scored). Why: in the check, the cluster-bootstrap SE was 0.0016 for M1 and 0.0013 for a paired Δ between two z-score variants. That is because about two-thirds of kinases are fixed at ½. So practically trivial differences would come out "significant".
- (f) Rule for unseen kinases: the point estimate of Δ on T-unseen must be >0. This is descriptive, with no extra test.
- (g) External sets (benchmarKIN tyrosine/Hernandez/CPTAC) are secondary. They are reported with CIs and no superiority claim.

**6. [must] §7 comparisons on DEV and selection optimism.** "Paired fold-level bootstrap" with 5 folds has 5 resampling units, which is not valid. Replace it with a paired group-level bootstrap of out-of-fold predictions pooled across all 5 folds. An idea is adopted only if the one-sided bootstrap p<0.05 on DEV-CV M1 **and** Δ ≥ 0.005. Keep a ledger (`results/dev_ledger.tsv`) of every configuration ever evaluated on DEV, with a config hash, and report K in REPORT.md. Also report DEV-best vs TEST (the optimism gap). Any hyperparameter tuned inside a learned component uses nested (inner) grouped CV, never the outer validation fold.

**7. [should] §5/§6 tuning parity for baselines.** The final method gets many DEV experiments; baselines get none. Give each baseline in the superiority family the same small DEV grid: minimum substrates {1, 3, 5}, library {PSP, PSP+best extended}, and a sign/standardisation choice where it applies. Freeze each baseline's DEV-best configuration at G2. Report the default configurations too. Otherwise "beats" partly measures how much tuning each method got.

**8. [should] §3/§5 label noise from promiscuous drugs.** 75% of Hijazi GT pairs come from the 8 drugs with ≥20 kinhib targets (PF3758309 98, LY2835219 83, AZD1480 52, …). Even on selective drugs (≤5 targets), z-score M1-scored was only about 0.54–0.55 in the quick check. Pre-register a stratified sensitivity table now (drug target count 1–5 / 6–19 / ≥20; Kd-potency tiers where Klaeger exists). Within each kinase, weight conditions by 1/n_targets(c) as a secondary M1 variant (M1-w). Do this before any method is run, so the choice cannot follow the results. The choice of primary universe stays with the GT/domain reviewer.

**9. [should] §5 harmony bins.** Define the bin per GT pair from the number of PSP substrates of k quantified in c. It must be method-independent (PSP-only, even for methods using predicted libraries) and frozen with the GT hash. Average the bin-macro over kinase means within each bin (change 3). Replace the adoption rule "no bin degrades by >0.01" (§7): many bins hold few kinases, so their SE exceeds 0.01. The new rule: reject the idea if the lower 90% bootstrap bound of Δ in any bin is < −0.02.

**10. [should] §4 TEST size and power.** Drug union-find produces very unequal groups. A crude check over inhibitor conditions found about 160 groups, the largest having 68, 57, 53 and 46 conditions, and 40 singletons. "~30% of groups" can therefore be far from 30% of GT pairs. Choose groups to hit 30±3% of GT pairs **and** of GT kinases in each harmony bin (seed 0, with up to 1,000 tries; take the first that meets the constraint, and log the criterion). Report n_groups, n_kinases per bin and the bootstrap SE of z-score M1 on TEST before G1, as the minimum-detectable-Δ statement.

**11. [should] §5 M4 definition.** Use pairs within each study only (a stratified Mann–Whitney: AUROC within each study, weighted by number of pairs). Negatives are conditions where k is an evaluable negative (not masked). Include only kinases scored by the method in those conditions, and require ≥2 target conditions from ≥2 groups. Imputing 0 makes no sense in M4.

**12. [should] §5 M2 definition.** For each pair of methods, the universe in c is S_c(A) ∩ S_c(B), for both the targets and the negatives. GT pairs whose target is outside it are dropped. Report how many pairs remain. Use the same bootstrap as change 5.

**13. [could] §5 calibration.** The 27 null conditions all come from one study and are correlated, and a per-kinase rate from 27 draws has binomial SD ≈ 0.04. Report the pooled |z|>1.96 rate with a study-cluster CI, plus a KS/QQ check of null p-value uniformity pooled across kinases. Do not judge per-kinase rates individually.

**14. [could] §5 ties and harness checks.**
- Midrank ties are fine. Report the share of pairs where the target ties with ≥1 scored negative; paralogs sharing motif or PSP sets will do this.
- Rows with expected_direction = 0 (27, the null controls) must be excluded from M1 by an explicit assert.
- Give the tolerance for M3 validation, e.g. |M3 − 0.785| ≤ 0.02.

No change is needed for: grouped split logic (L1–L5), the single logged TEST run, B = 2000, and Holm as the adjustment method.
