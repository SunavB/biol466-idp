# biol466-idp

BIOL 466 Independent Research Project 1, McGill University Department of Biology.
Author: Sunav Bajaj. Supervisor: Prof. Paul Harrison.

## Question

Given a modern protein language model (PLM) baseline, do curated Gene Ontology (GO) annotations add independent predictive information for the disorder-to-order (d2o) transition phenotype in intrinsically disordered proteins? Or has the PLM already absorbed the discriminative content that GO carries?

## Answer

**PLM embeddings already contain the discriminative content of GO Slim for this task.** No GO factor produced a statistically significant lift on top of ESM-2 in a pre-registered 2×2×2 factorial, but the null is mechanistically informative: permutation importance in the joint model splits **99.7% to ESM-2 dimensions and 0.3% to GO features**, and the 0.3% concentrates on biologically correct terms (RNA binding, DNA binding, transcription regulation, nucleus). The result generalizes to higher-resolution GO (full-term with ancestor propagation), two PLM families (ESM-2, ProstT5), two classifier families (Random Forest, XGBoost), and three prediction scales (whole-protein, disorder-pool, per-region).

The strongest lever on prediction is sequence representation and prediction scale, not annotation richness. Baseline AUPRC rises 0.199 → 0.237 → 0.273 across the three scales, each step larger than any GO effect observed anywhere in the study.

The final report source is in `report/` as markdown; compile with pandoc from the source files (see the assembly instructions in `report/` or the notebooks).

## Design at a glance

- **Dataset:** 1,279 non-redundant human proteins from DisProt release 2025_12 (188 d2o-positives at 14.7% prevalence). Region-level extension: 3,231 disordered regions (391 positives).
- **Label:** IDPO:0000011 (disorder-to-order structural transition) — kept outside the GO ontology to block GO-to-GO circularity.
- **Sequence features:** ESM-2 (`esm2_t33_650M_UR50D`) mean-pooled at three scales; ProstT5 (`Rostlab/ProstT5`) as a representation-family robustness check.
- **GO features:** GO Slim Generic (BP 64, MF 36, CC 25) plus a higher-resolution full-term ancestor-propagated encoding (BP 1,950, MF 618, CC 329) for the Week 10 test.
- **Classifier:** Random Forest primary (500 trees, class_weight balanced); XGBoost as robustness check.
- **Sequence-redundancy control:** CD-HIT at 40% identity; GroupKFold(5) on cluster IDs so no cluster is split across train and test.
- **Primary metric:** AUPRC. Statistical tests: paired Wilcoxon signed-rank with Benjamini-Hochberg correction across 7 non-baseline contrasts per experiment.
- **Pre-registered controls:** label-shuffling negative control (all conditions collapse to prevalence), stratified analysis for H4, scrambled-annotation control for the higher-resolution GO extension.

## Repository layout

```
biol466-idp/
├── README.md                       # this file
├── lab-notebook.md                 # 10-week execution log (Weeks 1-10)
├── notebooks/                      # main analysis pipeline (one per week's work)
│   ├── 00_setup_test.ipynb
│   ├── 01_data_acquisition.ipynb through 15_region_factorial.ipynb
│   ├── 16_full_go_factorial.py     # Week 10: higher-resolution GO
│   ├── 17_bp_resolution_strengthening.py  # Week 10: strengthening pass
│   └── 18_generate_figures.py      # produces the report figures
├── data/                           # gitignored (see below); .gitkeep only
├── results/                        # per-fold factorial outputs and stratified analyses
└── report/                         # report source (markdown, figures)
    ├── title-page.md
    ├── v3_toc.md
    ├── v3_abstract.md
    ├── acknowledgements.md
    ├── v3_introduction.md
    ├── v3_methods.md
    ├── v3_results.md
    ├── v3_discussion.md
    ├── references.md
    ├── v3_appendices.md
    └── figures/                    # PNG figures used in the report
```

## Data availability

Data files are gitignored to keep the repo small. All artifacts can be regenerated from the notebooks in `notebooks/`. To reproduce from scratch:

1. **Download DisProt release 2025_12** from `disprot.org`, save as `data/disprot.tsv`.
2. **Download human GOA annotations** from `current.geneontology.org`, save as `data/goa_human.gaf.gz`.
3. **Download `go.obo`** from `current.geneontology.org/ontology/go.obo`, save as `data/go.obo`.
4. **Run notebooks 01 through 15** in order to produce the master tables, feature matrices, cluster assignments, labels, and factorial results.
5. **Run `python notebooks/16_full_go_factorial.py`** and **`python notebooks/17_bp_resolution_strengthening.py`** for the Week 10 higher-resolution GO test.
6. **Run `python notebooks/18_generate_figures.py`** to produce the report figures.

## Reproducibility

- **Python 3.11.7** with `numpy 1.26.3`, `pandas 2.1.4`, `scikit-learn 1.4.0`, `xgboost 2.0.3`, `torch 2.2.0`, `fair-esm 2.0.1`, `transformers 4.37.2`, `scipy 1.12.0`, `statsmodels 0.14.1`, `goatools 1.3.11`.
- **Fixed random seed 42** for cross-validation splits, classifier initialization, and bootstrap resampling.
- **DisProt release:** `disprot_2025_12`.
- **PLM checkpoints:** `esm2_t33_650M_UR50D` (Meta AI) and `Rostlab/ProstT5` (TU Munich).

## References

Full report and reference list are in `report/`. Key methodological citations: Lin et al. (2023) for ESM-2, Heinzinger et al. (2023) for ProstT5, Quaglia et al. (2022) for DisProt, Fu et al. (2012) for CD-HIT, Saito and Rehmsmeier (2015) for AUPRC over ROC-AUC under class imbalance, Nair et al. (2022) for the cluster-aware evaluation protocol on IDP function prediction, Kulmanov and Hoehndorf (2020) for the closest prior finding.
