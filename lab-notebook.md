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
