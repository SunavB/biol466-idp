# 2. Methods

## 2.1 Dataset and label

From DisProt release 2025_12 (Quaglia et al., 2022) I extracted 1,279 non-redundant human proteins after collapsing UniProt isoform suffixes. At the protein level, y=1 if at least one disordered region carried an IDPO:0000011 annotation; y=0 otherwise. This gave **188 positives out of 1,279** (14.7% prevalence). At the region level (Week 9 extension), each Structural-state disorder region was a separate row, labelled y=1 if it overlapped ≥50% with any IDPO:0000011 annotation on the same protein (matching MoRF inclusion criteria; Disfani et al., 2012). This produced **3,231 regions with 391 positives** (12.1% prevalence).

## 2.2 Features

**Sequence features.** ESM-2 (`esm2_t33_650M_UR50D`, 1,280 dim; Lin et al., 2023) extracted from layer 33 in inference mode, mean-pooled at three scales: whole-protein (over the truncated sequence), disorder-pool (over DisProt-annotated residues), and per-region (over each region separately). Sequences longer than 1,022 residues were truncated to respect ESM-2's context window. Regions with start > 1,022 (5.5%) were mean-pooled over in-window residues. For representation-family robustness (Week 8), I also extracted ProstT5 embeddings (`Rostlab/ProstT5`, 1,024 dim; Heinzinger et al., 2023).

**GO features.** GOA human GAF file filtered to experimental evidence codes (EXP, IDA, IPI, IMP, IGI, IEP, HTP, HDA, HMP, HGI, HEP), excluding IEA. This produced 76,902 annotations across 1,266 of 1,279 proteins. I reduced to **GO Slim Generic** (Ashburner et al., 2000) using `goatools.mapslim` (Klopfenstein et al., 2018), producing three per-sub-ontology binary matrices: BP (64 terms), MF (36 terms), CC (25 terms).

## 2.3 Sequence-redundancy control and cross-validation

Two proteins that share >40% sequence identity are effectively the same evolutionary object for a classifier's purposes: if one is in the training set and its close homolog is in the test set, the classifier will look accurate without having learned anything generalizable. This is the standard failure mode of protein function benchmarks and can inflate apparent performance by 20-30 AUPRC points (Nair et al., 2022). I blocked it in two steps.

**CD-HIT** (Fu et al., 2012) groups sequences into clusters such that all members of a cluster share at least 40% identity with the cluster centroid. This produced 1,168 clusters from my 1,279 proteins. Two proteins in different clusters are guaranteed to share <40% identity.

**GroupKFold** (scikit-learn; Pedregosa et al., 2011) is a cross-validation scheme that treats each cluster ID as a "group" and guarantees that all proteins from a given group land in the same fold. This means no protein in the test set has a close homolog in the training set, so any measured performance is on genuinely novel sequences. Region-level rows inherit their parent protein's cluster ID, so within-protein and cross-protein homology are both held out.

I used five folds throughout. The same fold assignments were reused across every condition and every extension experiment, which enables paired statistical comparisons (each condition sees identical train/test splits, so differences are attributable only to the feature set change).

## 2.4 Factorial design and classifier

The 2×2×2 factorial over BP, MF, CC produces 8 conditions (C0 baseline through C7 full-GO). For each condition, the feature matrix concatenates the ESM-2 embedding with zero or more GO Slim matrices (1,280 dim in C0, up to 1,405 dim in C7). The classifier is held constant across conditions: **Random Forest** (500 trees, `class_weight='balanced'`, `min_samples_leaf=3`). This is the design's controlled-variable principle: only the feature set changes. I ran **XGBoost** (Chen and Guestrin, 2016) on the disorder-pool factorial as a classifier-family robustness check.

## 2.5 Metrics and statistical analysis

**Primary metric: AUPRC.** Area under the precision-recall curve. I chose AUPRC over the more conventional ROC-AUC because at 15% class prevalence, ROC-AUC can look impressive while a classifier gets almost no positives right; AUPRC penalizes exactly that failure mode and its chance baseline is the positive prevalence (0.147 at protein scale) rather than 0.5, so any lift is interpretable as "above pure guessing" (Saito and Rehmsmeier, 2015). Secondary metrics for triangulation: AUROC (prevalence-invariant, so it lets me compare across the 14.7% protein-scale and 12.1% region-scale prevalences); macro-F1; balanced accuracy; and Matthews correlation coefficient (Chicco and Jurman, 2020). For threshold-dependent metrics I binarized predictions at the positive-class prevalence rather than at 0.5, because Random Forest with class-balanced weighting at 15% prevalence often produces calibrated probabilities that never exceed 0.5, which makes the 0.5 threshold produce an all-negative-prediction artifact.

**Paired Wilcoxon signed-rank test.** Cross-condition comparisons used the paired Wilcoxon on per-fold AUPRC with alternative `greater` (each GO-augmented condition vs. C0). Wilcoxon is a non-parametric test that asks "is the median paired difference greater than zero?" and does not require the fold-level distribution to be normal — an important property because with only five folds the normality assumption of a paired t-test is untestable.

**Benjamini-Hochberg (BH) correction.** Running seven comparisons per experiment inflates the family-wise error rate: with seven independent tests at α = 0.05, the probability of at least one false positive is about 30%. BH controls the *false discovery rate* — the expected proportion of significant results that are false positives — at 5% by adjusting p-values so that only genuinely large effects survive (Benjamini and Hochberg, 1995). Both raw and BH-adjusted p-values are reported.

**Bootstrap 95% CI.** Effect sizes are reported as mean paired differences with 95% bootstrap confidence intervals (2,000 iterations, seed 42). The bootstrap CI is the range of mean lifts that could plausibly be produced by the underlying process given the observed data; the width of the CI is a direct measure of how much the point estimate could be off due to sample-size limitations, and is more informative than a p-value alone when the sample is small.

**A priori power analysis.** Before the primary factorial, a bootstrap simulation on the fold-level distribution showed that the paired Wilcoxon at α = 0.05 has approximately 80% power to detect a stable +0.02 AUPRC lift at n = 1,279 with 188 positives, but drops below 50% for a stable +0.01 lift. The study can therefore reject "GO adds a moderate-or-large effect" but cannot rule out very small effects. This bound is stated up front so the reader knows what claims the study can support before results are reported.

**Retrospective note on hypothesis formulation.** H1-H3 as pre-registered are directional ("each sub-ontology has a positive main effect on AUPRC") rather than magnitude-specific. Given the power analysis, a more informative pre-registration would have specified minimum detectable effects — for example, "H1: if PLMs have not subsumed GO Slim's BP content, we expect a BP main effect of at least +0.02 AUPRC." The directional formulation as written was near-certain to falsify given how small typical categorical-feature increments over a strong continuous baseline tend to be. Future studies in this design space should pre-register magnitude-specific hypotheses tied to a power calculation.

**Feature importance.** The mechanistic interpretation in §3.2 uses **permutation importance**: for each feature, its values are randomly shuffled and the drop in classifier performance is measured. Features the classifier truly depends on show large importance; features that are only weakly informative or redundant with other features show small importance. Aggregating importance by feature source (ESM-2 dimensions vs. GO Slim columns) gives a direct answer to the question "how much of the classifier's decision-making is coming from the sequence embedding versus the auxiliary GO annotations."

**Controls.** A label-shuffling negative control was run for each experiment; all conditions collapsed to within ±0.03 of chance-prevalence, confirming no cluster-boundary leakage. Feature importance was computed from the C7 Random Forest retrained on the full dataset.

## 2.6 Extensions

**Region-level factorial (Week 9).** The primary factorial was rerun at region scale on 3,231 rows, with the same cluster-aware CV and factorial protocol.

**Higher-resolution GO (Week 10) and strengthening pass.** To test whether the null was a coarseness artifact of GO Slim, I re-encoded GO at full-term resolution with ancestor propagation up `is_a` and `part_of` relations. This produced 1,950 BP terms, 618 MF terms, and 329 CC terms. Because the BP main effect at a single-seed configuration gave an initial +0.016 lift that looked interesting, I ran four pre-registered strengthening checks before writing it up: (1) bootstrap 95% CI on the fold-level paired diffs; (2) scrambled-BP negative control that permutes the protein-to-annotation mapping (an effect that persists under scrambling is dimensionality-driven, not annotation-specific); (3) top-20 BP features by permutation importance to check biological coherence; (4) stability across 9 (RF seed × CV shuffle seed) combinations. Pre-registered decision rules classified the result as a strong finding, a suggestive finding, or a spurious effect to be dropped.

## 2.7 Limitations

Four limitations bound the interpretation. First, **sample size is small** (n = 1,279 proteins, 188 positives). Bootstrap analysis suggests I can detect a stable +0.02 AUPRC lift but not a stable +0.01 lift. Second, **GO Slim is coarse**; I addressed this directly with the higher-resolution extension in §3.6. Third, **DisProt-annotated proteins are better-studied than average**, so my factorial evaluates GO under near-optimal availability conditions. Fourth, the negative class is d2o-negative-*or*-unannotated rather than confirmed d2o-negative; this MNAR ambiguity compresses effect sizes symmetrically and does not distort the relative comparison.
