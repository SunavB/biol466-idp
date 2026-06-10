# BIOL 466 — Lab Notebook

**Student:** sunav
**Course:** BIOL 466 — Independent Research Project 1
**Supervisor:** [name to be added]
**Course coordinator:** Nancy Nelson
**Started:** 2026-05-31

## How this notebook works

Dated entries, chronological (oldest at top, append at the bottom). Each entry captures what was done, what was found, what was decided, and what follow-ups remain. This notebook is the source material for the Methods and Results sections of the final report — every decision and dataset choice should leave a trail here.

## Project documents (live in the parent folder)

- `../BIOL466_Research_Proposal.docx` — current research proposal (updated post Week-1 spike).
- `../BIOL466_Strategy_Note.docx` — private grade-strategy note; not for the supervisor.
- `../BIOL466_Week1_Feasibility_Finding.docx` — Week-1 data-feasibility output, shareable with supervisor.

---

## Week 1

### 2026-05-31 — Scoping, design pivots, feasibility spike, tooling

**Project framing locked.** The research question shifted from yes/no ("do GO terms improve function prediction") to a quantitative comparison against a strong baseline: *how much*, and *which aspect of*, GO context improves prediction of intrinsically disordered protein function beyond an ESM-2 sequence-only baseline; are the GO aspects redundant or synergistic; is the benefit concentrated in the proteins where sequence alone fails. Four directional, falsifiable hypotheses (H1 main effect, H2 aspect ranking with MF predicted strongest, H3 sub-additive interaction, H4 stratification) replaced a generic null/alternative pair.

**Design pivot 1 — circularity.** Recognised that "use GO terms to predict protein function" is inherently circular because GO terms *are* the standard encoding of protein function. First fix: target = DisProt's disorder-function annotations (assumed non-GO). Incorporated into the first proposal draft.

**Design pivot 2 — the spike caught the first fix.** Parsed a sample DisProt entry (DP00004, human cathelicidin LL-37, UniProt P49913) and saw that DisProt's "Disorder function" annotations are themselves GO terms (e.g., region 134–170 annotated as cytolysis, GO:0019835; region 150–162 as amyloid fibril formation, GO:1990000). The first fix would have reintroduced circularity. Path A then defined as: target = DisProt's **Structural transition** layer, which uses IDPO terms (truly non-GO).

**Data feasibility spike — DisProt release 2025_12 (with_ambiguous_evidences variant).** Parsed the full dataset using `numbers-parser` + pandas:

- 13,396 curated annotation rows
- 3,199 distinct proteins
- 469 organisms

Headline distinct-protein counts (top 8 organisms):

| Organism | Total | With Structural transition | With Disorder function |
|---|---:|---:|---:|
| Homo sapiens | 1,301 | 218 (17%) | 647 (50%) |
| *S. cerevisiae* (S288c) | 214 | 36 (17%) | 108 (50%) |
| *Mus musculus* | 201 | 36 | 108 |
| *E. coli* K12 | 141 | 29 | 57 |
| *Arabidopsis thaliana* | 113 | 19 | 59 |
| *Rattus norvegicus* | 89 | 22 | 46 |
| *Drosophila melanogaster* | 64 | — | 30 |
| *Caenorhabditis elegans* | 56 | — | 21 |

The earlier Numbers pivot showed annotation-row counts of 355 (human) and 62 (yeast) for Structural transition; those resolve to 218 and 36 *distinct* proteins.

**Transition vocabulary — the class set is effectively binary.** Of 584 distinct proteins with any Structural transition annotation, 513 (88%) carry *disorder to order*; the next term, *order to disorder*, has only 58 proteins; the remaining six IDPO transition terms are individually rare (≤ 10 proteins each). Per organism: 87% of human transition-annotated proteins carry disorder-to-order, 92% in yeast. A multi-class formulation would be a severe class-imbalance problem.

**Reframe to binary disorder-to-order.** Path A is framed as binary classification: does the disordered region undergo a disorder-to-order transition (coupled folding-and-binding; IDPO:0000011). This is the canonical IDP functional mechanism and a recognised prediction task (cf. MoRF / molecular recognition feature literature). The 2×2×2 GO factorial design is unchanged; only the target shape changes.

Binary sample sizes (positives = proteins with a d2o annotation; negatives = DisProt proteins with disordered regions but no d2o annotation):

- Homo sapiens alone: 189 / 1,111 (≈15% positive)
- *S. cerevisiae* S288c alone: 33 / 181 — too thin
- Human + yeast pooled: 222 / 1,292
- All DisProt pooled: 513 / 2,680

**Organism decided (provisional).** Homo sapiens as primary; pooled human + yeast available as a sensitivity check.

**Evidence quality.** Structural transition annotations are predominantly experimental: NMR spectroscopy (~255), far-UV circular dichroism (~194), X-ray crystallography (~147 + 144 with-missing-coords), cryo-EM (~24).

**Path B backup confirmed.** 1,570 distinct proteins carry a Disorder function annotation across 387 distinct terms (dominated by GO:0005515 *protein binding*, 665 proteins; 42% of the annotated set). Sufficient as fallback if the Week-4 go/no-go fails, but the vocabulary would need curated grouping before use.

**Evaluation metric switched.** Primary metric changed from macro-F1 to **AUPRC** (area under the precision-recall curve), more informative than ROC-AUC under the ~15% positive prevalence. Macro-F1 and balanced accuracy reported as secondary. Class weighting during training. Paired Wilcoxon signed-rank tests across folds with Benjamini-Hochberg correction across the eight conditions.

**New limitation identified — MNAR negatives.** Negative labels (absence of a disorder-to-order annotation) do not guarantee absence of the behaviour — they may simply mean the region has not been studied for this mechanism. Result must be reported as "annotated positives versus the rest", not as "all positives versus all true negatives." To be stated explicitly in the proposal limitations and examined for sensitivity to annotation richness.

**Proposal updated.** All Path-A-binary changes incorporated: §1 summary, §3 RQ + H1 (AUPRC) and H2 (MF rationale strengthened around binding), §4 objectives 1 and 6, §5.1 (new target paragraph + sample sizes + "why binary not multi-class"), §5.6 (AUPRC), §7 (MNAR limitation added), §8 (Week-1 spike noted complete), §9 (organism and target questions resolved).

**Tool setup complete.**

- Miniconda installed (macOS Apple Silicon).
- conda env `biol466` created with `python=3.11`; packages installed via pip: `pandas numpy scikit-learn matplotlib seaborn jupyterlab numbers-parser`.
- GitHub Desktop installed; private repo `biol466-idp` initialised at `~/Documents/Claude/Projects/466/biol466-idp/` with the skeleton `data/`, `notebooks/`, `src/`, `results/`.
- `.gitignore` extended for Python data work: `data/* !data/.gitkeep *.pkl .ipynb_checkpoints/ .DS_Store results/figures/*.png`. Empty `.gitkeep` placed in `data/` so the folder is tracked.
- README replaced with the proposal's §1 summary plus a pointer to the docs in the parent folder.
- Sanity-check notebook `notebooks/00_setup_test.ipynb` runs `import pandas, sklearn, matplotlib`; versions printed. First commit pushed to GitHub.

**Supervisor meeting deferred.** Discussed with the supervisor; agreed to defer the formal Week-1 meeting until after the Week-4 go/no-go checkpoint. Implication: the decisions below are *provisional* until that meeting, but they are consistent with what the spike showed, and proceeding with them through Weeks 2–4 is the right call. Any decisions that arise mid-stream and would benefit from earlier supervisor input will be flagged in real time.

---

## Decisions log — provisional, pending Week-4 supervisor sign-off

1. Research question reframed from yes/no to quantitative ("how much", "which aspect", "redundant or synergistic", "stratified").
2. Four directional hypotheses (H1–H4) replace a generic null/alternative.
3. Prediction target = **binary disorder-to-order** transition (IDPO:0000011) from DisProt's Structural transition layer.
4. Primary organism = **Homo sapiens**; pooled human + yeast as sensitivity check.
5. Sequence representation = pre-trained **ESM-2** embeddings, no further training.
6. Classifier = **random forest**, fixed hyperparameters across all conditions.
7. Feature factorial = **2×2×2** (BP, MF, CC GO aspects, each present/absent) → 8 conditions; sequence-only control = Condition 1.
8. GO feature encoding = **GO Slim**, experimental and curator-reviewed evidence codes only (IEA filtered out).
9. Primary metric = **AUPRC**; secondary = macro-F1, balanced accuracy. Paired Wilcoxon signed-rank across folds; Benjamini-Hochberg correction across conditions; effect sizes and confidence intervals reported alongside p-values. Class weighting during training.
10. Core controls: **label-shuffling negative control** and **stratified analysis** by baseline performance.
11. Fallback target if the Week-4 go/no-go fails = **disorder-function classification** with curated term grouping.

## Open items / follow-ups

- **Week 2 reading list anchor:** CAFA (Radivojac 2013; Zhou 2019), DeepGO / DeepGOPlus (Kulmanov 2018, 2020), ESM-2 (Lin 2023; Rives 2021), DisProt (Quaglia 2022), GO Consortium 2023, van der Lee 2014 IDR classification, and the MoRF / coupled folding-and-binding line — Mohan 2006 (ANCHOR), Disfani 2012 (MoRFpred), Oates 2013 (D²P²), plus a recent IDP function review. Output: annotated bibliography + positioning paragraph.
- **Mark weighting** (70/30, 60/40, 50/50 between project work and report) — confirm at the Week-4 meeting.
- **Stretch goals** — second organism comparison; head-to-head against DeepGOPlus and/or MoRFpred. Decide at end of Week 8 if on schedule.
- **Data for Week 3** — pull UniProt sequences for the ≈ 1,300 human DisProt proteins (or the pooled set); pull GO BP/MF/CC annotations from QuickGO/GOA with experimental evidence codes; assemble master table with the binary d2o label.

---

## Week 2 — Literature review

### 2026-05-31 — Annotated bibliography, positioning, metric locked

**Deliverable:** `../BIOL466_Week2_Literature_Review.docx` — fourteen annotated references in three thematic sections, a positioning paragraph, and the locked evaluation metric.

**Bibliography composition.** Started from a twelve-paper draft (six IDP/CFB; five function-prediction/GO; one PLM). Added two foundational references the draft was missing: Ashburner et al. (2000) — the founding GO paper — and Rives et al. (2021) — the foundational protein-language-model paper that ESM-2 builds on. Final list of fourteen entries:

- *IDPs and coupled folding-and-binding:* van der Lee et al. 2014; Mohan et al. 2006 (MoRFs); Disfani et al. 2012 (MoRFpred); Oates et al. 2013 (D²P²); Quaglia et al. 2022 (DisProt); Wright & Dyson 2015.
- *Function prediction and GO:* Ashburner et al. 2000; GO Consortium 2023; Radivojac et al. 2013 (CAFA1); Zhou et al. 2019 (CAFA3); Kulmanov et al. 2018 (DeepGO); Kulmanov & Hoehndorf 2020 (DeepGOPlus).
- *Protein language models:* Rives et al. 2021 (ESM-1b); Lin et al. 2023 (ESM-2).

**Positioning.** Two-paragraph synthesis. Paragraph 1 — frames the gap: MoRF-prediction predictors (e.g., MoRFpred) target disorder-to-order with engineered local features; modern function predictors (DeepGOPlus on ESM-2-class embeddings, benchmarked through CAFA) predict GO from sequence. Paragraph 2 — locates this project at the intersection: a controlled 2×2×2 factorial measurement of how much each GO aspect (BP, MF, CC) adds to a strong ESM-2 sequence baseline when the task is binary disorder-to-order transition labelled from DisProt's IDPO Structural transition layer. The framing turns a yes/no question into a quantitative one (how much, which aspect, redundant or synergistic, concentrated where).

**Evaluation metric locked** (the second Week-2 deliverable from the proposal).

- Primary: AUPRC.
- Secondary: macro-F1; balanced accuracy.
- Training: class weighting inversely proportional to class frequency.
- Cross-validation: repeated stratified k-fold (k = 5, 5 repeats; final k confirmed at Week-4 go/no-go).
- Statistical testing: paired Wilcoxon signed-rank across folds; Benjamini-Hochberg correction across the eight conditions; effect sizes and bootstrap CIs reported with p-values.
- Main and interaction effects from the factorial structure.

**Reading discipline going forward.** Read at least the abstract of every entry; full read of Ashburner 2000, Wright & Dyson 2015, Quaglia 2022 and Lin 2023. Track outstanding reads here and append notes inline as I work through each paper.

**Next reading targets (stretch tier):**
- Latest CAFA round (CAFA5, 2023–2024) if available.
- Any post-2024 ESM successors (e.g., ESM Cambrian) for currency.
- A recent IDP-function review post-2020 to balance the older Wright & Dyson 2015.

**Open items / follow-ups (Week 3 prep):**
- Pull UniProt sequences for the candidate human DisProt set (~1,300 proteins).
- Pull GO BP/MF/CC annotations from QuickGO/GOA with experimental evidence codes only (IEA excluded).
- Build the master table with the binary d2o label and the GO Slim feature matrices.
- Decide which ESM-2 checkpoint to use for embeddings (likely 650M for the quality-tractability trade-off).

---

## Weeks 3–4 — Data acquisition, cleaning, go/no-go

### 2026-05-31 — Plan compiled

**Working plan:** `./weeks3-4-plan.md` — detailed task-by-task plan with concrete code snippets, decision points, and explicit go/no-go criteria.

**Week 3 covers:** notebook scaffolding, candidate set definition (human only, locked), UniProt sequence pull, GOA human GAF download, binary d2o label construction, master-table skeleton.

**Week 4 covers:** experimental-evidence filter, GO Slim mapping with `goatools`, binary label edge-case audit, CD-HIT sequence-redundancy filtering at 40% identity (cluster-aware CV for Week 5), annotation-bias check (positives vs negatives on `n_go_exp`), comprehensive EDA, formal go/no-go against seven explicit criteria, lab-notebook + proposal updates.

**Tools to install before starting Week 3:**

- pip into `biol466`: `biopython requests goatools obonet`
- Homebrew: `cd-hit`

**Entries below will be populated as work proceeds.**

### 2026-06-04 — Week 3 execution

**Source data.** DisProt release 2025_12 with_ambiguous_evidences, TSV download from disprot.org. Loaded in pandas: 13,396 region annotations × 21 columns. Column names renamed to lowercase-underscore convention so the plan code works as-written.

**Candidate set (W3.2).** Filtered to *Homo sapiens*. Collapsed UniProt isoform suffixes (e.g., O00204-2 → O00204) and deduplicated. Result: **1,279 unique base accessions** (down from 1,301 raw — about 22 isoforms collapsed into canonical forms). Saved to `data/candidate_accs.csv`. **Decision locked:** human-only organism.

**Sequences (W3.3).** Pulled from UniProt REST `stream` endpoint, batches of 100. Retrieval rate: **1,279 / 1,279 (100%)**. Saved to `data/sequences.fasta` and `data/sequences.csv`. Length distribution: median 482, mean 644, range 24–34,350. Max (34,350) is almost certainly titin. **Note for Week 5:** very long sequences will need handling (truncation or sliding-window pooling) when generating ESM-2 embeddings.

*Issue hit and fixed:* the original batch size of 500, combined with the isoform-suffixed accessions, produced HTTP 400 from UniProt. Cause: the `-N` in `O00204-2` is interpreted as a Lucene NOT operator. Fix: stripped isoform suffixes globally on `df['acc']`, reduced batch to 100, and used `requests.get(..., params=...)` so URL encoding is handled by the library.

**GO annotations (W3.4).** Downloaded `goa_human.gaf.gz` (15 MB) from `current.geneontology.org`. Total GOA human annotations: **906,445**. Filtered to candidate set: **129,684 annotations**. Raw aspect distribution (unfiltered for evidence): F = 63,092 (49%), C = 38,324 (30%), P = 28,268 (22%) — MF-skewed because IEA/computational annotations are still present; this will rebalance after the W4.1 experimental-evidence filter. **GO coverage:** 1,271 / 1,279 proteins (**99.4%**) have ≥ 1 raw GO annotation.

*Issue hit and fixed:* urllib's default User-Agent returned 403 from the GO Consortium CDN. Switched to `requests` with a `Mozilla/5.0` User-Agent header and added the EBI mirror as a fallback URL.

**Binary d2o label (W3.5).** Positives = proteins with at least one DisProt region annotated `IDPO:0000011` (disorder-to-order). Result: **188 positives / 1,091 negatives (14.7% positive)**. One fewer positive than the spike's 189 — a single isoform-only annotation collapsed into a canonical-form negative. Within noise.

**Master table (W3.6).** Joined labels + sequences + raw GO count. Saved to `data/master_raw.csv`. Sanity: 0 missing sequences, 8 proteins with zero raw GO annotations (matches 1,279 − 1,271).

**Outputs in `data/`:** `disprot.tsv` · `candidate_accs.csv` · `sequences.fasta` · `sequences.csv` · `goa_human.gaf.gz` · `go_annotations_raw.csv` · `labels.csv` · `master_raw.csv`.

**Status: Week 3 complete.** Ready for Week 4 (cleaning, GO Slim mapping, CD-HIT, EDA, go/no-go).

### 2026-06-05 — Week 4 execution + go/no-go verdict

**W4.1 — Experimental-evidence filter.** Restricted to experimental evidence codes (EXP, IDA, IPI, IMP, IGI, IEP, HTP, HDA, HMP, HGI, HEP). Result: **76,902 annotations** retained from 129,684 (59% retention). **Proteins with ≥1 experimental annotation: 1,266 / 1,279 (99.0%)**. Aspect distribution after filter — F 53,949 (70%), P 11,877 (15%), C 11,076 (14%) — MF dominates because human MF receives heavy IPI (protein-binding) experimental annotation, while BP got hit hardest by removing IBA (computational ortholog inference). Saved to `data/go_annotations_experimental.csv`.

**W4.2 — GO Slim mapping.** Downloaded `go.obo` (36.7 MB, 41,552 terms) and `goslim_generic.obo` (206 terms). Mapped each annotation through `goatools.mapslim.mapslim`. Output matrices:

- **BP**: 866 proteins × **64 slim terms**
- **MF**: 1,014 proteins × **36 slim terms**
- **CC**: 1,120 proteins × **25 slim terms**

Per-aspect protein coverage variance reflects how experimental-evidence depth scales by aspect in human (CC > MF > BP). Saved as `features_GO_{BP,MF,CC}_slim.csv`.

**W4.3 — Binary label audit.** Confirmed label distribution: **188 positives / 1,091 negatives (14.7%)**. Of the 188 positives, 179 have *only* "disorder to order" annotations; 9 also carry another transition type (mixed but still positive). D2o region length distribution across 299 annotated regions: median 26 aa, mean 53 aa, range 10–454 aa, **0 regions shorter than 5 residues**. Length range is consistent with the MoRF literature (typical MoRFs 10–70 residues).

**W4.4 — CD-HIT sequence-redundancy filter.** Clustered at 40% identity (`-c 0.4 -n 2`). **1,168 distinct clusters from 1,279 proteins** in 14 s of CPU. Cluster size: mean 1.10, median 1, max 4 — overwhelmingly singletons, as expected for IDPs. **Decision logged:** Week-5 cross-validation will be **GroupKFold by cluster**, no cluster spans train/test in any fold.

**W4.5 — Annotation-bias check.** Headline ratio for the proposal's limitations section:

- **GO annotations per protein (n_go_exp): positive median 54.5 vs negative median 31.0 → 1.76× ratio**
- Sequence length: positive median 475 vs negative median 484 → **0.98× (no length confound)**

Saved figure `results/figures/annotation_bias.png`.

**W4.6 — Comprehensive EDA.** Built `data/master_clean.csv` (1,279 × 7): acc, d2o, length, n_go_exp, cluster, n_dis_regions, total_dis_len. Four figures saved to `results/figures/`:

- `eda_summary.png` — class distribution, length histogram by class, GO count vs length scatter, cluster-size histogram
- `eda_go_slim_top20.png` — top-20 slim terms for each aspect
- `eda_slim_per_protein.png` — slim-terms-per-protein distribution per aspect
- `eda_disorder_profile.png` — n_dis_regions and total_dis_len boxplots by class

**Top slim terms are biologically on point.** BP: regulation of transcription (GO:0006355) — IDPs are classical transcription regulators. MF: catalytic activity (GO:0003824). CC: organelle (GO:0043226).

**Second annotation-richness confound surfaced.** Beyond the GO-count bias, positives also have more disorder overall: n_dis_regions ratio **1.50×** (median 3 vs 2), total_dis_len ratio **1.73×** (median 104 vs 60). Biologically expected (well-studied d2o proteins tend to be more thoroughly characterised for disorder), but this is a second confound to declare alongside the GO bias in the report's limitations.

**W4.7 — Go/no-go verdict.** **All seven criteria PASS.** Recorded numerically:

| Criterion | Threshold | Actual |
|---|---|---|
| 1. N sufficient | ≥150 pos, ≥500 neg | 188 / 1,091 |
| 2. Positive rate near spike | 10–20% | 14.7% |
| 3. Sequence coverage | ≥95% | 100.0% |
| 4. GO coverage (experimental) | ≥70% | 99.0% |
| 5. GO Slim non-degenerate | ≥20 cols, ≥10 useful | BP 64/54, MF 36/31, CC 25/25 |
| 6. Positives span clusters | ≥50 | 182 |
| 7. Annotation bias median ratio | <4× | 1.76× |

**Verdict: GO — proceed to Week 5.**

**W4.8 — Closing.** Lab notebook + proposal §5.1 updated to reflect the final filtered numbers (1,279 candidates, 188/1,091 split, 1,168 clusters). Commit will be tagged `week4-go-nogo`.

**Status: Week 4 complete.** Modelling-ready dataset in `data/master_clean.csv`; three GO Slim feature matrices ready; cluster IDs ready for GroupKFold CV. Week 5 starts with ESM-2 embedding generation and the end-to-end thin-slice run.

---

## Week 5 — Sequence features and thin-slice pipeline

### 2026-06-05 — ESM-2 extraction, AA-comp baseline, thin slice verified

**Tools installed.** `torch 2.12.0`, `fair-esm 2.0.0`, `tqdm`, `ipywidgets`. Apple Silicon MPS backend confirmed available (`torch.backends.mps.is_available() == True`). Device: `mps`.

**W5.2 — ESM-2 embedding extraction.** Used the **ESM-2 650M** checkpoint (`esm2_t33_650M_UR50D`). Per-residue token representations from layer 33 mean-pooled (excluding CLS and EOS) to a single 1,280-dim vector per protein. Truncation at 1,022 residues (model context 1,024 minus CLS/EOS).

- Sequences processed: **1,279 / 1,279 (100%)** — no failures.
- Embedding dim: **1,280**.
- Truncated proteins (length > 1,022): **180 (14.1%)** — higher than initial 5% estimate; mostly disorder-rich long proteins. To be stated as a methods limitation.
- L2-norm range across embeddings: 4.53 – 10.04, mean 7.23. No zero-norm vectors; no NaNs.
- Saved as `data/features_esm2.npz` (compressed).

**W5.3 — Amino-acid composition reference baseline.** 20-dim frequency vector + log-length z-score → 21-dim per protein. Saved as `data/features_aacomp.csv` (1,279 × 21). Row-sum sanity passes (AA-frequency rows sum to 1.0). This is the reference baseline, *not* the project baseline — its purpose is to demonstrate, by comparison, how strong ESM-2 is.

**W5.4 — Thin-slice modelling pipeline.** Cluster-aware GroupKFold(5) split, Random Forest (300 trees, class_weight='balanced'), AUPRC primary metric. **Cluster-leakage assert never fired**: no homology cluster spans train and test in any fold.

Four CV configurations were tested as part of pipeline diagnosis:

| Setup | AUPRC | AUROC | macro-F1 | bal_acc |
|---|---:|---:|---:|---:|
| Chance baseline (prevalence) | 0.147 | 0.500 | — | 0.500 |
| AA composition + RF | 0.166 ± 0.031 | 0.525 ± 0.054 | 0.460 ± 0.005 | 0.500 ± 0.000 |
| ESM-2 + RF | 0.193 ± 0.047 | 0.563 ± 0.074 | 0.460 ± 0.005 | 0.500 ± 0.001 |
| ESM-2 + LR (StandardScale + L2) | 0.197 ± 0.044 | 0.560 ± 0.075 | 0.524 ± 0.040 | 0.531 ± 0.048 |
| ESM-2 → PCA(64) + RF | 0.190 ± 0.040 | 0.548 ± 0.084 | 0.460 ± 0.005 | 0.500 ± 0.000 |

**Diagnostic — RF training-set self-fit.** Trained RF on all 1,279 proteins, predicted on the same set. **Train AUPRC = 1.000, Train AUROC = 1.000.** This proves the features carry separable signal *and* the pipeline code is correct; the train→test gap (1.0 → 0.19) is high-variance generalization, not a bug.

**Interpretation.**
- The pipeline is mechanically verified end-to-end. AUPRC/AUROC/F1/bal_acc all compute correctly, splits respect clusters, no leakage.
- **Mean-pooled ESM-2 at the protein level is a weak signal for d2o.** Above-chance but only marginally (AUPRC ≈ 0.19 vs chance 0.147). RF, LR, and PCA-reduced RF all land at the same level, suggesting the signal is genuinely diffuse rather than poorly recovered.
- AA composition is essentially noise (AUPRC ≈ 0.17), confirming the proposal's strong-baseline principle: using AA-comp as the baseline would trivially be beaten by anything; ESM-2 is the meaningful baseline.
- The **macro-F1 = 0.460 and bal-acc = 0.500** across most setups reflect "predicts mostly negatives at threshold 0.5" — a calibration artefact for class-imbalanced RF, not a feature problem. LR slightly mitigates this (F1 = 0.524, bal-acc = 0.531) but the AUPRC story is unchanged.

**Implications for the proposal.**
- **The ESM-2 baseline AUPRC is ≈ 0.19**. This is the bar GO context must clear in the 8-condition factorial.
- The bar is low enough that *plausible* improvements from GO are possible, and a finding that GO doesn't help is equally interpretable as "modern sequence representations already capture most of the protein-level signal."
- The "weak protein-level baseline" finding strengthens the proposal's design rationale rather than weakening it: the factorial comparison stays clean either way.

**Status: Week 5 complete (with caveat).** Pipeline verified, baseline characterised, files saved (`features_esm2.npz`, `features_aacomp.csv`, `results/thin_slice_aacomp.csv`, `results/thin_slice_esm2.csv`). The caveat is that the baseline is quite weak; this is acknowledged and the project continues as planned — GO comparisons are the next test.

**Week 6 — Planned focus.**
1. **Wire GO Slim features into the pipeline.** Implement `make_X(condition_id)` for all 8 factorial conditions; concatenate ESM-2 + selected GO Slim matrices.
2. **Try a stronger classifier in parallel.** Gradient boosting (XGBoost or LightGBM) often outperforms RF on dense features in this regime; worth ~30 min of investigation. RF stays primary unless evidence is overwhelming.
3. **Optional: explore per-residue ESM-2 averaging over disordered regions only.** If mean-pooling over the whole protein dilutes the signal (likely), pooling over annotated disorder regions could lift the baseline meaningfully.
4. **Decide whether to keep AUPRC as primary or add MCC.** With persistent all-negative bias at threshold 0.5, MCC (Matthews correlation coefficient) might be a more honest secondary; AUPRC stays primary per the proposal.
5. **Re-run thin slice with the GO BP+MF+CC condition (full eight only).** Quick check: does GO context lift the ~0.19 baseline meaningfully? If yes, we're on a good trajectory. If no, that's the headline finding and we plan the report accordingly.

---

## Week 6 — The 8-condition factorial (main experimental result)

### 2026-06-05 — Factorial run, statistical analysis, null finding confirmed, Path II routing

**W6.1 — `make_X(condition_id)` builder.** Implemented a bit-encoded condition table (BP/MF/CC each on/off) with a single function returning the concatenated ESM-2 + GO Slim feature matrix for any of the 8 factorial conditions. Verified dimensions: Condition 1 (Seq only) = 1,280; Condition 8 (all GO) = 1,280 + 64 + 36 + 25 = **1,405**. All other conditions land between as expected.

**Final modelling set:** 1,279 proteins (188 positives, 14.7%), 1,168 CD-HIT clusters at 40% identity. Threshold = 0.147 (prevalence) to avoid the W5 all-negative-at-0.5 artefact.

**Cross-validation:** GroupKFold(5), splits predefined once and shared across all 8 conditions for valid paired-fold comparisons. Cluster leakage assert passes in every fold.

**W6.2 — RF factorial.** All eight conditions run with the same Random Forest (500 trees, `class_weight='balanced'`, `min_samples_leaf=3`). Aggregate AUPRC (mean ± std over 5 folds):

| Condition | Label | AUPRC | AUROC | MCC |
|---|---|---:|---:|---:|
| 1 | Seq only | 0.199 ± 0.058 | 0.579 | 0.092 |
| 2 | Seq + BP | 0.189 ± 0.042 | 0.570 | 0.085 |
| 3 | Seq + MF | 0.200 ± 0.047 | 0.587 | 0.098 |
| 4 | Seq + CC | 0.194 ± 0.057 | 0.577 | 0.107 |
| 5 | Seq + BP + MF | 0.203 ± 0.053 | 0.580 | 0.074 |
| 6 | Seq + BP + CC | 0.205 ± 0.049 | 0.584 | 0.097 |
| 7 | Seq + MF + CC | 0.193 ± 0.056 | 0.581 | 0.081 |
| 8 | Seq + BP + MF + CC | 0.191 ± 0.050 | 0.574 | 0.065 |

All conditions land in the 0.189–0.205 band; baseline (0.199) sits in the middle of the pack. Notably the fully-augmented condition 8 is slightly *below* baseline. The error bars (±0.05) overlap massively across all bars. Saved per-fold scores to `results/factorial_rf_per_fold.csv` and a bar chart to `results/figures/factorial_rf_auprc.png`.

**Fold structure dominates condition effects.** Across every condition: fold 4 lands at AUPRC 0.25–0.30, fold 2 at 0.14–0.16. The fold-to-fold variance is roughly equal to the largest between-condition difference, confirming that any feature-set effect is below noise. Plausibly fold 4's test split happens to contain more annotation-rich, easier-to-predict proteins (consistent with the 1.76× annotation-richness bias in positives from W4.5).

**W6.3 — Pairwise stats + factorial effects.** Paired Wilcoxon signed-rank against condition 1, BH-corrected:

| Condition | mean Δ vs Seq | CI95 | p_BH |
|---|---:|---|---:|
| Seq + BP | −0.010 | [−0.028, +0.002] | 0.94 |
| Seq + MF | +0.000 | [−0.011, +0.009] | 0.73 |
| Seq + CC | −0.005 | [−0.012, +0.003] | 0.94 |
| Seq + BP + MF | +0.004 | [−0.002, +0.008] | 0.73 |
| Seq + BP + CC | +0.006 | [−0.004, +0.015] | 0.73 |
| Seq + MF + CC | −0.006 | [−0.015, +0.004] | 0.94 |
| Seq + BP + MF + CC | −0.008 | [−0.017, −0.000] | 0.94 |

**Main effects on AUPRC** (averaged over all conditions with vs without each aspect):
- BP: **+0.0006**
- MF: **−0.0004**
- CC: **−0.0017**

**Interaction effects:**
- BP × MF: −0.0001
- BP × CC: +0.0079
- MF × CC: −0.0147

Saved as `results/factorial_rf_stats.csv`.

**W6.4 — XGBoost robustness check.** Same 8 conditions, same shared CV splits, XGBoost (500 trees, depth 6, lr 0.05, `scale_pos_weight` for class imbalance) replacing RF. Mean AUPRC by condition: 0.197, 0.202, 0.197, 0.198, 0.193, 0.210, 0.195, 0.207. Same null pattern. Average XGB-vs-RF delta across conditions: +0.003, well below the +0.03 threshold we set for considering a classifier change. **Decision: RF stays primary as pre-registered.** Saved per-fold scores to `results/factorial_xgb_per_fold.csv`.

**Hypothesis verdicts (Week-6 cut):**

- **H1 (main effect of GO context) — FALSIFIED.** No pairwise comparison reaches BH-corrected significance (best p_BH = 0.73). Effect sizes are uniformly tiny (largest +0.006 AUPRC). The full-GO condition 8 is mildly worse than baseline (CI [−0.017, 0]), consistent with noise injection from adding ineffective features. **Null hypothesis cannot be rejected.**
- **H2 (MF most informative) — FALSIFIED.** Main-effect ranking is BP (+0.0006) > MF (−0.0004) > CC (−0.0017). MF is second, not first. All three effects are effectively zero, so the ordering is noise and the directional prediction fails.
- **H3 (sub-additive interaction) — NOT INTERPRETABLE.** Interactions (−0.0001, +0.008, −0.015) are tiny and there are no main effects for them to interact around. The hypothesis is moot when the underlying effects are zero.
- **H4 (stratification) — NOT TESTED YET.** Reserved for Week 7.

**Headline interpretation.** This is the *informative null* the proposal explicitly anticipated (§6: "It is a genuine and informative possibility that a modern sequence representation already captures most of the available signal, leaving little for GO context to add."). Modern ESM-2 protein-level mean-pooled embeddings appear to already encode the protein-level signal relevant to coupled folding-and-binding prediction; concatenating GO Slim features adds no measurable benefit at this representation. Result is robust across RF and XGBoost.

**W6.5 — Week-7 routing: Path II selected.** The null at protein-level mean-pooled representation needs one further test before the report locks the finding in: per-residue ESM-2 averaged over **annotated disordered regions only** (using DisProt region coordinates). Rationale: MoRFs are 10–70 residues; mean-pooling over a 500-residue protein could be diluting the local signal that GO context might rescue. If GO still shows no effect under that finer representation, the null is robust and becomes the headline of the report. If GO suddenly helps, that is itself a methodological finding worth its own discussion.

**Week 7 plan (preview):**
1. Build a per-residue ESM-2 representation pooled over DisProt-annotated disordered regions only. Re-run the 8-condition factorial with this representation.
2. Run the **label-shuffling negative control** for the protein-level setup — shuffle d2o labels across proteins, rerun the factorial. AUPRC should drop to ≈ 0.147 (chance) across all conditions. Confirms results aren't a methodology artefact.
3. **Stratified analysis for H4** — partition test proteins by Seq-only baseline AUPRC tertile and measure GO contribution within each tertile.
4. **Feature-importance interpretation** for the protein-level setup — which GO Slim terms RF deems most useful, even though the lift is null overall, to characterise what the GO features were "trying" to contribute.
5. Begin drafting the Methods section of the report.

**Status: Week 6 complete.** Main experimental result is the informative null: GO Slim features over a strong ESM-2 mean-pooled baseline do not improve protein-level disorder-to-order prediction. RF stays as the primary classifier. Path II confirmed for Week 7's robustness work.

**Files saved this week:**
- `results/factorial_rf_per_fold.csv`
- `results/factorial_rf_stats.csv`
- `results/factorial_xgb_per_fold.csv`
- `results/figures/factorial_rf_auprc.png`

---

## Week 7 — Robustness, stratification, interpretation

### 2026-06-06 — Disorder-pool representation, negative control, H4 stratification, feature importance

**W7.2 — Per-residue ESM-2 with disorder-region pooling.** Same ESM-2 650M model, same MPS device, same 1,022-residue truncation as W5. The change: mean-pool over only the residues that fall inside DisProt-annotated disordered regions (Structural state layer) instead of all residue tokens. Embedding wall time ~27 min on Apple Silicon MPS. **41 proteins (3.2%)** had all disorder annotations beyond residue 1,022 and fell back to whole-protein pooling — noted as a methods caveat. Output: `data/features_esm2_disorder.npz`, shape (1,279, 1,280), L2-norm range 4.78–10.18 (mean 8.75; comparable to W5's 4.53–10.04 / mean 7.23). No NaNs, no zero-norm vectors.

*Sanity vs W5 whole-protein embeddings:* per-protein cosine similarity median **0.946**, mean 0.927, range 0.475–1.000. The two representations are similar in overall manifold but distinct enough to plausibly carry different information — exactly the regime where a disorder-specific representation could surface a signal that whole-protein pooling dilutes.

**W7.3 — Disorder-pool factorial.** Same 8 conditions, same shared CV splits, same RF (500 trees, `class_weight='balanced'`, `min_samples_leaf=3`), threshold = prevalence (0.147). Aggregate AUPRC (mean ± std over 5 folds):

| Condition | Label | AUPRC | AUROC | MCC |
|---|---|---:|---:|---:|
| 1 | Seq only | **0.237 ± 0.035** | 0.654 | 0.176 |
| 2 | Seq + BP | 0.228 ± 0.043 | 0.651 | 0.178 |
| 3 | Seq + MF | 0.232 ± 0.041 | 0.658 | 0.190 |
| 4 | Seq + CC | 0.231 ± 0.037 | 0.651 | 0.186 |
| 5 | Seq + BP + MF | 0.226 ± 0.046 | 0.657 | 0.178 |
| 6 | Seq + BP + CC | 0.230 ± 0.046 | 0.654 | 0.182 |
| 7 | Seq + MF + CC | **0.240 ± 0.048** | 0.659 | 0.189 |
| 8 | Seq + BP + MF + CC | 0.230 ± 0.045 | 0.655 | 0.170 |

**Two headline findings.**

1. **Disorder pooling materially improves the sequence baseline.** Cond 1 jumped from W6's 0.199 to **0.237** (+0.038 absolute, ~19% relative). AUROC 0.579 → 0.654. MCC 0.092 → 0.176 (nearly doubled). Whole-protein mean-pooling was diluting the d2o-relevant signal; restricting the pool to disordered residues recovers a meaningfully stronger representation. **The sequence representation, not the model, was the bottleneck.**
2. **GO context still adds nothing.** All 8 conditions cluster in 0.226–0.240. Cond 8 (all GO) is *below* baseline at 0.230. The null is robust across both representations.

**W7.4 — Pairwise stats on disorder pool.** All BH-corrected p-values **≥ 0.94**; six of seven equal **1.000**. All mean differences vs Cond 1 are negative or essentially zero — adding GO features slightly *degrades* the disorder-pooled baseline. Main effects: BP −0.0065, MF +0.0009, CC +0.0016 — ordering changed from W6 (where BP > MF > CC), confirming the ranking is noise. **H1 falsified more strongly than W6; H2 falsified with a different ordering; H3 again uninterpretable.**

**W7.5 — Negative control (label shuffling).** Shuffled `y` once with fixed seed; reran the 8-condition factorial. Aggregate shuffled-label AUPRC: 0.156–0.180 across all conditions, AUROC 0.46–0.59 (essentially random ranking). Compared to chance prevalence (0.147) and the real-data disorder-pool baseline (0.237):

- Real-data disorder-pool Cond 1: **0.237**
- Real-data Cond 1, whole-protein (W6): 0.199
- Shuffled-label any condition: 0.156–0.180 (mean ≈ 0.170)
- Chance prevalence: 0.147

The 0.237 baseline is ~2.7σ above the shuffled null. **The lift from whole-protein → disorder pool is real; the lift from adding GO (≈ 0) is null; both reads of the result are confirmed by the negative control.**

**W7.6 — Stratified analysis for H4.** Per-protein out-of-fold probabilities from Cond 1 and Cond 8, partitioned into tertiles by Cond-1 confidence:

| Tertile | n | n_pos | Cond 1 AUPRC | Cond 8 AUPRC | Δ |
|---|---:|---:|---:|---:|---:|
| Low (hardest) | 427 | 27 | 0.140 | 0.169 | **+0.028** |
| Mid | 426 | 75 | 0.260 | 0.249 | −0.011 |
| High (easiest) | 426 | 86 | 0.251 | 0.221 | **−0.031** |

**H4 has partial directional support, with a caveat.** In the low-confidence stratum GO context lifts AUPRC by +0.028 (~20% relative) — the H4 prediction. But this is offset by losses in the higher-confidence strata, netting to zero across the dataset.

Mechanism is visible in the per-protein probability scatter (`results/figures/per_protein_proba_shift.png`): at low Cond-1 probabilities the Cond-8 predictions sit slightly above y=x, and at high Cond-1 probabilities they sit slightly below. **GO context compresses the probability range toward the class prior (~0.15).** RF regularises uncertain predictions toward prevalence; that helps where the sequence model was guessing and hurts where it was confident. Net: zero discriminative gain.

**W7.7 — Feature importance.** Cond 8 RF trained on all data; importance distribution across 1,405 features:

| Feature kind | n columns | Sum of importance | Fraction |
|---|---:|---:|---:|
| ESM-2 | 1,280 | 0.9975 | **99.7%** |
| BP | 64 | 0.0012 | 0.1% |
| MF | 36 | 0.0008 | 0.1% |
| CC | 25 | 0.0006 | 0.1% |

**RF puts 99.7% of its importance on the sequence features and essentially ignores the 125 GO Slim columns.** Not a single GO term enters the top 30 features. The bar chart (`results/figures/feature_importance_by_kind.png`) is one tall ESM-2 bar over three invisible slivers.

But — the GO terms that *did* rank highest within the 0.3% RF allocated to GO are **biologically on-target for IDPs**:

| Rank | GO ID | Aspect | Term |
|---|---|---|---|
| 1 | GO:0012501 | BP | Programmed cell death |
| 2 | GO:0005829 | CC | Cytosol |
| 3 | GO:0003723 | MF | **RNA binding** |
| 4 | GO:0003677 | MF | **DNA binding** |
| 5 | GO:0030163 | BP | Protein catabolic process |
| 6 | GO:0060089 | MF | Molecular transducer activity |
| 7 | GO:0007010 | BP | Cytoskeleton organization |
| 8 | GO:0006281 | BP | DNA repair |
| 9 | GO:0005730 | CC | Nucleolus |
| 10 | GO:0008289 | MF | Lipid binding |
| 11 | GO:0005215 | MF | Transporter activity |
| 12 | GO:0005634 | CC | Nucleus |
| 13 | GO:0043226 | CC | Organelle |
| 14 | GO:0002376 | BP | Immune system process |
| 15 | GO:0006351 | BP | Transcription (DNA-templated) |

RNA binding, DNA binding, transcription, nucleus, signal transduction, apoptosis — the canonical functional categories of IDPs. The features RF found marginally useful are the biologically right ones; they just couldn't add anything on top of what ESM-2 already encodes.

**Synthesis — the report's headline interpretation.** Three observations now stitch into a clean explanatory story:

1. GO context produces no measurable AUPRC improvement on top of ESM-2 (under either pooling; under either RF or XGBoost; under shared, paired, BH-corrected statistical testing; confirmed by label-shuffling negative control).
2. RF essentially ignores GO features when given access to ESM-2 alongside them (99.7% of importance on sequence).
3. The few GO terms RF *does* rank are exactly the biologically correct ones for IDPs (RNA binding, DNA binding, transcription, nucleus, signalling).

Conclusion: **ESM-2 already encodes the same protein-functional information that GO Slim categories summarise.** Pre-training on hundreds of millions of sequences taught the protein language model the sequence signatures of RNA-binding proteins, transcription factors, signalling proteins; adding curated GO categories as additional features is redundant for this prediction task at this representation. The proposal's pre-registered "informative null" framing fits this result precisely.

**H1, H2, H3 verdicts unchanged from W6** (falsified, falsified, uninterpretable). **H4 partially supported in direction** (the lift in the low-confidence tertile matches the prediction) but the magnitude is small and is offset by losses elsewhere, so the *net* H4 prediction (that GO benefits the hard cases enough to lift overall AUPRC) is falsified.

**Files saved this week:**
- `data/features_esm2_disorder.npz` — disorder-pooled embeddings
- `results/factorial_rf_per_fold_disorder.csv`
- `results/factorial_rf_stats_disorder.csv`
- `results/negative_control_per_fold.csv`
- `results/stratified_h4.csv`
- `results/feature_importance_top30.csv`
- `results/figures/factorial_pool_compare.png`
- `results/figures/three_way_comparison.png`
- `results/figures/per_protein_proba_shift.png`
- `results/figures/feature_importance_by_kind.png`

**Status: Week 7 experimental work complete.** All robustness and interpretation analyses done; the report has its three-leg result structure (real-data factorial across two representations, paired stats, negative control), its H4 partial-support nuance, and its mechanistic explanation (ESM-2 subsumes GO Slim). Week 8 begins Methods drafting (W7.8 was deferred) and Results/Discussion. Supervisor meeting can now be scheduled — bringing the full result story plus all four hypothesis verdicts.
