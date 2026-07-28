# Appendices

## Appendix A. Per-fold AUPRC for the disorder-pool factorial

**Table A.1** — Per-fold AUPRC (whole-protein pooling, Random Forest, ESM-2). Input to the paired Wilcoxon tests in §3.1.

| Condition | Fold 1 | Fold 2 | Fold 3 | Fold 4 | Fold 5 | Mean | SD |
|---|---|---|---|---|---|---|---|
| C0 | 0.211 | 0.187 | 0.194 | 0.208 | 0.194 | 0.199 | 0.010 |
| C1 | 0.213 | 0.184 | 0.198 | 0.209 | 0.192 | 0.199 | 0.012 |
| C2 | 0.209 | 0.185 | 0.196 | 0.207 | 0.194 | 0.198 | 0.010 |
| C3 | 0.212 | 0.188 | 0.197 | 0.209 | 0.196 | 0.200 | 0.010 |
| C4 | 0.210 | 0.185 | 0.197 | 0.208 | 0.193 | 0.199 | 0.010 |
| C5 | 0.213 | 0.187 | 0.199 | 0.210 | 0.195 | 0.201 | 0.010 |
| C6 | 0.211 | 0.186 | 0.197 | 0.208 | 0.194 | 0.199 | 0.010 |
| C7 | 0.212 | 0.186 | 0.198 | 0.209 | 0.193 | 0.200 | 0.011 |

## Appendix B. ProstT5 replication (disorder-pool factorial)

**Table B.1** — ProstT5 disorder-pool factorial (Random Forest, 5-fold GroupKFold).

| Condition | AUPRC | 95% CI | Δ vs C0 | BH-adj p |
|---|---|---|---|---|
| C0 | 0.229 | [0.203, 0.256] | — | — |
| C1 +BP | 0.231 | [0.204, 0.259] | +0.002 | 0.87 |
| C2 +MF | 0.228 | [0.201, 0.256] | −0.001 | 0.92 |
| C3 +CC | 0.232 | [0.205, 0.260] | +0.003 | 0.83 |
| C4 | 0.230 | [0.202, 0.258] | +0.001 | 0.94 |
| C5 | 0.233 | [0.206, 0.261] | +0.004 | 0.83 |
| C6 | 0.230 | [0.203, 0.258] | +0.001 | 0.94 |
| C7 | 0.232 | [0.204, 0.260] | +0.003 | 0.83 |

No lift is significant after BH correction. The null under ProstT5 is qualitatively identical to the null under ESM-2.

## Appendix C. Region-level factorial

**Table C.1** — Region-level factorial (Random Forest, 5-fold GroupKFold on protein clusters, n = 3,231 regions).

| Condition | AUPRC | 95% CI | Δ vs C0 | BH-adj p |
|---|---|---|---|---|
| C0 | 0.273 | [0.251, 0.296] | — | — |
| C1 +BP | 0.278 | [0.256, 0.301] | +0.005 | 0.61 |
| C2 +MF | 0.275 | [0.253, 0.298] | +0.002 | 0.78 |
| C3 +CC | 0.292 | [0.269, 0.315] | **+0.019** | 0.14 |
| C4 | 0.279 | [0.257, 0.302] | +0.006 | 0.55 |
| C5 | 0.297 | [0.274, 0.320] | +0.024 | 0.11 |
| C6 | 0.294 | [0.271, 0.317] | +0.021 | 0.12 |
| C7 | 0.306 | [0.283, 0.329] | **+0.033** | 0.09 |

The largest single lift (+0.033 at C7) is the BP×CC interaction discussed in §3.6. None of the effects survives BH correction at α = 0.05.

## Appendix D. Higher-resolution GO factorial and strengthening pass

**Table D.1** — Full-term GO factorial at disorder-pool scale (Random Forest, seed 42). Feature dimensionality: 1,950 BP, 618 MF, 329 CC. No condition survives BH correction.

| Condition | AUPRC | Δ vs C0 | Wilcoxon raw p | BH-adj p |
|---|---|---|---|---|
| C0 | 0.234 | — | — | — |
| C1 +BP | 0.250 | +0.016 | 0.500 | 0.583 |
| C2 +MF | 0.233 | −0.001 | 0.594 | 0.594 |
| C3 +CC | 0.238 | +0.003 | 0.500 | 0.583 |
| C4 +BP+MF | 0.245 | +0.010 | 0.156 | 0.365 |
| C5 +BP+CC | 0.246 | +0.012 | 0.156 | 0.365 |
| C6 +MF+CC | 0.243 | +0.009 | 0.156 | 0.365 |
| C7 full | 0.243 | +0.009 | 0.219 | 0.383 |

**Table D.2** — Strengthening pass on the BP full-term effect (per-fold diffs and 9-seed stability).

| Fold | C0 baseline | C1 +BP | Diff |
|---|---|---|---|
| 1 | 0.240 | 0.330 | +0.091 |
| 2 | 0.328 | 0.317 | −0.012 |
| 3 | 0.202 | 0.210 | +0.008 |
| 4 | 0.166 | 0.170 | +0.004 |
| 5 | 0.234 | 0.224 | −0.010 |
| **Mean** | 0.234 | 0.250 | **+0.016** |
| **Bootstrap 95% CI** | | | **[−0.008, +0.054]** |
| **Scrambled-BP control** | | | **+0.009** |
| **Stability across 9 (RF, CV) seeds** | | | **mean +0.009, range [+0.003, +0.016]** |

**Table D.3** — Top-15 BP full-term features by permutation importance in the C1 model.

| Rank | GO ID | Term name |
|---|---|---|
| 1 | GO:0080090 | regulation of primary metabolic process |
| 2 | GO:0010557 | positive regulation of macromolecule biosynthetic process |
| 3 | GO:0120035 | regulation of plasma membrane bounded cell projection organization |
| 4 | GO:0002376 | immune system process |
| 5 | GO:0023052 | signaling |
| 6 | GO:0048518 | positive regulation of biological process |
| 7 | GO:0010604 | positive regulation of macromolecule metabolic process |
| 8 | GO:0023051 | regulation of signaling |
| 9 | GO:0050896 | response to stimulus |
| 10 | GO:0010605 | negative regulation of macromolecule metabolic process |
| 11 | GO:0048583 | regulation of response to stimulus |
| 12 | GO:0060255 | regulation of macromolecule metabolic process |
| 13 | GO:0140014 | mitotic nuclear division |
| 14 | GO:0009056 | catabolic process |
| 15 | GO:0009267 | cellular response to starvation |

These are ancestor-inherited generic regulatory terms rather than IDP-mechanism-specific terms (compare Figure 1b, which shows Slim top-15 as RNA binding, DNA binding, transcription regulation, etc.). Consistent with the dimensionality-hedging interpretation in §3.5.

## Appendix E. Reproducibility

**Environment.** Python 3.11.7, numpy 1.26.3, pandas 2.1.4, scikit-learn 1.4.0, xgboost 2.0.3, torch 2.2.0, fair-esm 2.0.1, transformers 4.37.2, scipy 1.12.0, statsmodels 0.14.1. Full lock file in the project repository.

**Random seed.** All experiments used random seed 42 for CV splits, classifier initialization, and bootstrap resampling.

**Data files** (in `data/`): `master_clean.csv` (1,279 proteins), `labels.csv` (d2o labels), `clusters.csv` (CD-HIT clusters), `features_esm2*.npz` (three pooling variants), `features_prostt5_disorder.npz`, `features_GO_{BP,MF,CC}_slim.csv`, `go_annotations_experimental.csv`, `go.obo`, `regions_master.csv`.

**Notebooks.** `01_data_acquisition.ipynb` through `18_generate_figures.py` in the project repository, one per week's work.
