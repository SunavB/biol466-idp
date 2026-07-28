# 2. Methods

## 2.1 Dataset and label

From DisProt release 2025_12 (Quaglia et al., 2022) I extracted 1,279 non-redundant human proteins after collapsing UniProt isoform suffixes. At the protein level, y=1 if at least one disordered region carried an IDPO:0000011 annotation; y=0 otherwise. This gave **188 positives out of 1,279** (14.7% prevalence). At the region level (Week 9 extension), each Structural-state disorder region was a separate row, labelled y=1 if it overlapped ≥50% with any IDPO:0000011 annotation on the same protein (matching MoRF inclusion criteria; Disfani et al., 2012). This produced **3,231 regions with 391 positives** (12.1% prevalence).

## 2.2 Features

**Sequence features.** ESM-2 (`esm2_t33_650M_UR50D`, 1,280 dim; Lin et al., 2023) extracted from layer 33 in inference mode, mean-pooled at three scales: whole-protein (over the truncated sequence), disorder-pool (over DisProt-annotated residues), and per-region (over each region separately). Sequences longer than 1,022 residues were truncated to respect ESM-2's context window. Regions with start > 1,022 (5.5%) were mean-pooled over in-window residues. For representation-family robustness (Week 8), I also extracted ProstT5 embeddings (`Rostlab/ProstT5`, 1,024 dim; Heinzinger et al., 2023).

**GO features.** GOA human GAF file filtered to experimental evidence codes (EXP, IDA, IPI, IMP, IGI, IEP, HTP, HDA, HMP, HGI, HEP), excluding IEA. This produced 76,902 annotations across 1,266 of 1,279 proteins. I reduced to **GO Slim Generic** (Ashburner et al., 2000) using `goatools.mapslim` (Klopfenstein et al., 2018), producing three per-sub-ontology binary matrices: BP (64 terms), MF (36 terms), CC (25 terms).

## 2.3 Sequence-redundancy control and cross-validation

I clustered sequences with **CD-HIT** (Fu et al., 2012) at 40% identity, producing 1,168 clusters. All cross-validation used scikit-learn's `GroupKFold(5)` (Pedregosa et al., 2011) on the CD-HIT cluster ID as the group, so no cluster was split across train and test. Region-level rows inherited their parent protein's cluster ID. The same fold assignments were reused across every condition and every extension experiment, which enables paired statistical comparisons.

## 2.4 Factorial design and classifier

The 2×2×2 factorial over BP, MF, CC produces 8 conditions (C0 baseline through C7 full-GO). For each condition, the feature matrix concatenates the ESM-2 embedding with zero or more GO Slim matrices (1,280 dim in C0, up to 1,405 dim in C7). The classifier is held constant across conditions: **Random Forest** (500 trees, `class_weight='balanced'`, `min_samples_leaf=3`). This is the design's controlled-variable principle: only the feature set changes. I ran **XGBoost** (Chen and Guestrin, 2016) on the disorder-pool factorial as a classifier-family robustness check.

## 2.5 Metrics and statistical analysis

Primary metric: **AUPRC**, chosen over ROC-AUC because it stays informative at 15% prevalence and calibrates against the prevalence baseline (Saito and Rehmsmeier, 2015). Secondary metrics: AUROC, macro-F1, balanced accuracy, and Matthews correlation coefficient (Chicco and Jurman, 2020). For threshold-dependent metrics I binarized at the positive-class prevalence rather than at 0.5.

Cross-condition tests used the paired Wilcoxon signed-rank test on per-fold AUPRC with alternative `greater`. Multiple-comparison correction across the seven pairwise comparisons per experiment used the Benjamini-Hochberg procedure (Benjamini and Hochberg, 1995) at α = 0.05. Effect sizes are reported as mean paired differences with 95% bootstrap CIs (2,000 iterations, seed 42).

**Controls.** A label-shuffling negative control was run for each experiment; all conditions collapsed to within ±0.03 of prevalence, confirming no leakage. Feature importance for the mechanistic interpretation used the C7 Random Forest retrained on the full dataset and its `feature_importances_` attribute aggregated by feature source.

## 2.6 Extensions

**Region-level factorial (Week 9).** The primary factorial was rerun at region scale on 3,231 rows, with the same cluster-aware CV and factorial protocol.

**Higher-resolution GO (Week 10) and strengthening pass.** To test whether the null was a coarseness artifact of GO Slim, I re-encoded GO at full-term resolution with ancestor propagation up `is_a` and `part_of` relations. This produced 1,950 BP terms, 618 MF terms, and 329 CC terms. Because the BP main effect at a single-seed configuration gave an initial +0.016 lift that looked interesting, I ran four pre-registered strengthening checks before writing it up: (1) bootstrap 95% CI on the fold-level paired diffs; (2) scrambled-BP negative control that permutes the protein-to-annotation mapping (an effect that persists under scrambling is dimensionality-driven, not annotation-specific); (3) top-20 BP features by permutation importance to check biological coherence; (4) stability across 9 (RF seed × CV shuffle seed) combinations. Pre-registered decision rules classified the result as a strong finding, a suggestive finding, or a spurious effect to be dropped.

## 2.7 Limitations

Four limitations bound the interpretation. First, **sample size is small** (n = 1,279 proteins, 188 positives). Bootstrap analysis suggests I can detect a stable +0.02 AUPRC lift but not a stable +0.01 lift. Second, **GO Slim is coarse**; I addressed this directly with the higher-resolution extension in §3.6. Third, **DisProt-annotated proteins are better-studied than average**, so my factorial evaluates GO under near-optimal availability conditions. Fourth, the negative class is d2o-negative-*or*-unannotated rather than confirmed d2o-negative; this MNAR ambiguity compresses effect sizes symmetrically and does not distort the relative comparison.
