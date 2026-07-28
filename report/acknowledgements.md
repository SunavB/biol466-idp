# Acknowledgements

I thank my BIOL 466 supervisor for their patient guidance across the thirteen weeks of this project, and in particular for their willingness to sanction the two mid-course design pivots (the switch from a GO-target to an IDPO-target framing in Week 2, and the region-level extension added in Week 9) that shaped the study into its final form. Their questioning at the post-Week-4 checkpoint on the specific mechanism by which GO Slim features could plausibly contribute information beyond ESM-2 embeddings led directly to the pre-registration of hypothesis H4 and to the stratified analysis reported in §3.7.

I thank Nancy Nelson for coordinating the BIOL 466 course structure and for the writing scaffolding that shaped the section-by-section reporting cadence used across the project's weekly lab-notebook entries.

I acknowledge the developers of the DisProt database, whose continued curation effort makes benchmarks of this kind possible; the Gene Ontology Consortium, whose sustained infrastructure work makes GO Slim tractable as a feature source; and the developers of the ESM-2 and ProstT5 protein language models, whose open release of pre-trained weights made the two-family robustness comparison in §3.5 achievable within a semester-length project timeline.

Computational work was performed on Apple Silicon (M-series) hardware with PyTorch Metal Performance Shaders acceleration. Random Forest and XGBoost baselines used the `scikit-learn` and `xgboost` open-source Python libraries respectively. All figures were produced with `matplotlib` and `seaborn`. I thank the open-source communities behind these tools.

Any errors of interpretation, statistical judgement, or methodological choice that remain in this report are my own.
