# Weeks 3–4 Work Plan — Data Acquisition, Cleaning, and Go/No-Go

**Owner:** sunav
**Compiled:** 2026-05-31
**Covers:** Weeks 3 and 4. Target end-state is a cleaned modelling-ready dataset and a formal Week-4 go/no-go decision ready to take to the supervisor.

This plan is intended to be lived in — tick boxes as you complete tasks, leave notes inline, and append observations to the lab notebook as you go.

## Top-level goals

By the end of Week 4 the project should have:

- A documented, reproducible data-acquisition pipeline (notebooks + saved intermediate files).
- A cleaned master dataset: one row per candidate human DisProt protein, carrying the binary disorder-to-order label, the protein sequence, curated GO annotations (BP/MF/CC, experimental codes only), and a cluster ID for sequence-redundancy-aware cross-validation.
- Three GO Slim feature matrices (one per aspect) ready to feed the modelling pipeline in Week 5.
- A short EDA report — figures saved to `results/figures/` plus a written summary in the lab notebook.
- A formal go/no-go decision, recorded in the lab notebook, against the explicit criteria listed in W4.7.

## Tools to install up front (~15 min)

From Terminal, add the new packages to the existing `biol466` env:

```
conda activate biol466
pip install biopython requests goatools obonet
```

CD-HIT for sequence redundancy filtering (needs Homebrew; install Homebrew first from `brew.sh` if you don't have it):

```
brew install cd-hit
```

Verify:

```
python -c "import Bio, requests, goatools, obonet; print('ok')"
cd-hit -h | head -1
```

---

## Week 3 — Data acquisition (~9 hours)

Mostly mechanical work: pull data, build a master table skeleton, sanity-check. Save every intermediate file under `data/` and commit notebooks regularly.

### W3.1 — Notebook scaffolding (15 min)

- [ ] In `notebooks/`, create three empty notebooks: `01_data_acquisition.ipynb`, `02_data_cleaning.ipynb`, `03_eda.ipynb`.
- [ ] Top cell of each: `import pandas as pd, numpy as np, matplotlib.pyplot as plt; pd.set_option('display.max_columns', 60)`.

### W3.2 — Define the candidate protein set (45 min)

In `01_data_acquisition.ipynb`:

```python
import pandas as pd

# Re-parse from the Numbers/TSV source, or load the pickle if you saved one.
# Replace path as appropriate; if you don't have the pickle, re-run the Week-1 spike parse.
df = pd.read_pickle("../data/disprot.pkl")

human = df[df['organism'] == 'Homo sapiens']
candidates = sorted(human['acc'].unique())
print(f"Candidate human proteins: {len(candidates)}")  # expect ~1,301

pd.Series(candidates, name='acc').to_csv("../data/candidate_accs.csv", index=False)
```

- [ ] Lock the organism decision in the lab notebook: **human only** (yeast has 33 d2o positives — too thin to add useful power; pooling complicates the controlled comparison). Mark as a Week-3 decision; revisit only if Week-4 shows N is borderline.

### W3.3 — Pull UniProt sequences (1–2 hours)

UniProt's REST `stream` endpoint returns FASTA for a query, batched 500 at a time to be safe.

```python
import requests
from io import StringIO
from Bio import SeqIO

accs = pd.read_csv("../data/candidate_accs.csv")['acc'].tolist()
batch = 500
records = []
for i in range(0, len(accs), batch):
    chunk = accs[i:i+batch]
    query = "+OR+".join(f"accession:{a}" for a in chunk)
    url = f"https://rest.uniprot.org/uniprotkb/stream?query={query}&format=fasta"
    r = requests.get(url, timeout=180); r.raise_for_status()
    records.extend(SeqIO.parse(StringIO(r.text), 'fasta'))
print(f"Retrieved {len(records)} sequences")

SeqIO.write(records, "../data/sequences.fasta", "fasta")
seq_df = pd.DataFrame([{
    'acc': rec.id.split('|')[1] if '|' in rec.id else rec.id,
    'sequence': str(rec.seq),
    'length': len(rec.seq),
} for rec in records])
seq_df.to_csv("../data/sequences.csv", index=False)
print(seq_df['length'].describe())
```

- [ ] Compare retrieved count vs candidate count; identify any missing accessions (likely obsolete or merged). Decision point: if > 5% missing, investigate before proceeding.
- [ ] Glance at the length distribution. Flag anything extreme (very short < 30 aa or very long > 5000 aa) for review in W4.

### W3.4 — Pull GO annotations (1 hour)

Download the GOA human GAF file (one ~30 MB file, much cleaner than the QuickGO API for bulk):

```python
import urllib.request

url = "http://current.geneontology.org/annotations/goa_human.gaf.gz"
urllib.request.urlretrieve(url, "../data/goa_human.gaf.gz")

cols = ['DB','DB_Object_ID','DB_Object_Symbol','Qualifier','GO_ID','DB_Reference',
        'Evidence_Code','With_From','Aspect','DB_Object_Name','DB_Object_Synonym',
        'DB_Object_Type','Taxon','Date','Assigned_By','Annotation_Extension',
        'Gene_Product_Form_ID']
go = pd.read_csv("../data/goa_human.gaf.gz", sep='\t', comment='!', names=cols,
                 compression='gzip', low_memory=False)
print(f"GOA human total annotations: {len(go)}")

go_candidates = go[go['DB_Object_ID'].isin(accs)].copy()
print(f"Annotations for our candidate set: {len(go_candidates)}")
go_candidates.to_csv("../data/go_annotations_raw.csv", index=False)

# Quick coverage check
print("Aspect distribution:")
print(go_candidates['Aspect'].value_counts())   # P (BP), F (MF), C (CC)
print("Proteins with at least one GO annotation:",
      go_candidates['DB_Object_ID'].nunique())
```

- [ ] Confirm aspect distribution matches expectation (BP usually largest).
- [ ] Note how many of the 1,301 candidates have ≥1 raw GO annotation. Expect close to 100% for human.

### W3.5 — Build the binary d2o label (45 min)

Positive = any DisProt region for this `acc` carries IDPO:0000011 (disorder to order).

```python
d2o = df[(df['organism']=='Homo sapiens') & (df['term']=='IDPO:0000011')]
positive_accs = set(d2o['acc'].unique())
print(f"Positive (d2o) proteins: {len(positive_accs)}")   # expect 189

labels = pd.DataFrame({'acc': candidates})
labels['d2o'] = labels['acc'].isin(positive_accs).astype(int)
labels.to_csv("../data/labels.csv", index=False)
print(labels['d2o'].value_counts())
```

- [ ] Verify the positive count matches the spike (189 for human).

### W3.6 — Build the master table skeleton (30 min)

```python
master = (labels
          .merge(seq_df, on='acc', how='left')
          .merge(go_candidates.groupby('DB_Object_ID').size()
                              .rename('n_go_raw').reset_index()
                              .rename(columns={'DB_Object_ID':'acc'}),
                 on='acc', how='left'))
master['n_go_raw'] = master['n_go_raw'].fillna(0).astype(int)
master.to_csv("../data/master_raw.csv", index=False)

print(master.head())
print(f"Missing sequence: {master['sequence'].isna().sum()}")
print(f"Zero GO annotations (raw): {(master['n_go_raw']==0).sum()}")
```

- [ ] Quick eyeball: do the columns look right? Anything obviously wrong?
- [ ] Commit and push the notebook + new files in `data/` (large `.gz` and `.fasta` should be gitignored — confirm via `git status`).

### Week 3 deliverable

- Files in `data/`: `candidate_accs.csv`, `sequences.fasta`, `sequences.csv`, `goa_human.gaf.gz`, `go_annotations_raw.csv`, `labels.csv`, `master_raw.csv`.
- Notebook `01_data_acquisition.ipynb` — clean, commented, top-to-bottom reproducible.
- Lab-notebook Week-3 entry: candidate N, retrieval rates, decisions (human-only locked).
- Commit + push.

---

## Week 4 — Cleaning, GO Slim, redundancy filtering, EDA, go/no-go (~9 hours)

More reasoning, more decisions. This week ends with a formal go/no-go.

### W4.1 — Evidence-code filter (45 min)

In `02_data_cleaning.ipynb`:

```python
go = pd.read_csv("../data/go_annotations_raw.csv", low_memory=False)

EXPERIMENTAL = {'EXP','IDA','IPI','IMP','IGI','IEP',
                'HTP','HDA','HMP','HGI','HEP'}     # standard experimental codes
go_exp = go[go['Evidence_Code'].isin(EXPERIMENTAL)].copy()
print(f"Annotations after experimental-evidence filter: {len(go_exp)} (from {len(go)})")
print("Proteins surviving with ≥1 experimental annotation:",
      go_exp['DB_Object_ID'].nunique())
go_exp.to_csv("../data/go_annotations_experimental.csv", index=False)
```

- [ ] Record: how many proteins drop to zero GO annotations after this filter? Those proteins will appear in the sequence-only condition but will be all-zeros in the GO-augmented conditions — that is *intended* and fine, but record the count.

### W4.2 — GO Slim mapping (~1.5 hours)

```python
from goatools.obo_parser import GODag
from goatools.mapslim import mapslim
import urllib.request

urllib.request.urlretrieve("http://current.geneontology.org/ontology/go.obo",
                           "../data/go.obo")
urllib.request.urlretrieve("http://current.geneontology.org/ontology/subsets/goslim_generic.obo",
                           "../data/goslim_generic.obo")

godag = GODag("../data/go.obo")
slimdag = GODag("../data/goslim_generic.obo")

def slim_terms_for(gid):
    try:
        direct, ancestors = mapslim(gid, godag, slimdag)
        return list(direct | ancestors)
    except Exception:
        return []

go_exp['slim_terms'] = go_exp['GO_ID'].apply(slim_terms_for)
exploded = (go_exp.explode('slim_terms')[['DB_Object_ID','slim_terms','Aspect']]
                  .dropna()
                  .rename(columns={'DB_Object_ID':'acc','slim_terms':'slim_term','Aspect':'aspect'}))

for letter, name in [('P','BP'), ('F','MF'), ('C','CC')]:
    sub = exploded[exploded['aspect']==letter]
    mat = (sub.assign(v=1)
              .pivot_table(index='acc', columns='slim_term', values='v', fill_value=0)
              .astype('int8'))
    mat.to_csv(f"../data/features_GO_{name}_slim.csv")
    print(f"{name}: {mat.shape[0]} proteins × {mat.shape[1]} slim terms")
```

- [ ] Inspect each matrix shape: expect tens of slim columns per aspect (~30–80).
- [ ] Check that most columns have non-trivial variance (≥ 5–10 proteins). Drop columns that are essentially constant.

### W4.3 — Binary label final checks (30 min)

The per-protein label was set in W3.5. Now do edge-case audits:

- [ ] Are there positives with **only** d2o (no other transition annotation), or are most also annotated with other transition types? (Mixed is fine.)
- [ ] Are there positives whose annotated d2o region is implausibly short (< 5 residues)? Decision: include for now; flag for review only if results look weird.
- [ ] Final label distribution. Expect ~189 / ~1,111; deviations of a few proteins are normal and OK.

### W4.4 — Sequence redundancy filtering with CD-HIT (~1 hour)

This is methodologically critical. Without it, homologous proteins can land in both training and test folds and inflate apparent accuracy across *every* condition equally — bad for honest reporting, fine for the comparison between conditions, but better to remove the confound.

From Terminal, in the `data/` folder:

```
cd ~/Documents/Claude/Projects/466/biol466-idp/data
cd-hit -i sequences.fasta -o sequences_nr40.fasta -c 0.4 -n 2 -d 0 -M 0
```

Then back in Python parse the `.clstr` file:

```python
clusters = {}
current = None
with open("../data/sequences_nr40.fasta.clstr") as f:
    for line in f:
        if line.startswith(">Cluster"):
            current = int(line.strip().split()[1])
        else:
            # Lines look like: 0  120aa, >sp|P49913|... at 92.50%
            acc = line.split(">")[1].split("|")[1] if "|" in line else line.split(">")[1].split("...")[0]
            clusters[acc] = current
cl_df = pd.Series(clusters, name='cluster').rename_axis('acc').reset_index()
cl_df.to_csv("../data/clusters.csv", index=False)
print(f"Distinct clusters: {cl_df['cluster'].nunique()} (from {len(cl_df)} proteins)")
```

- [ ] Record cluster count. For ~1,300 human IDPs at 40% identity expect a meaningful reduction (probably ~900–1,200 clusters; IDPs are usually low-homology).
- [ ] Decision logged: Week-5 cross-validation will be **GroupKFold by cluster** — no cluster spans train and test in any fold.

### W4.5 — Annotation-bias check (45 min)

The proposal flags residual GO-vs-DisProt annotation correlation as the headline limitation. Measure it now.

```python
import seaborn as sns
master = (pd.read_csv("../data/labels.csv")
            .merge(pd.read_csv("../data/sequences.csv")[['acc','length']], on='acc', how='left')
            .merge(go_exp.groupby('DB_Object_ID').size().rename('n_go_exp').reset_index()
                          .rename(columns={'DB_Object_ID':'acc'}), on='acc', how='left'))
master['n_go_exp'] = master['n_go_exp'].fillna(0).astype(int)

fig, ax = plt.subplots(1, 2, figsize=(12,4))
sns.boxplot(data=master, x='d2o', y='n_go_exp', ax=ax[0])
ax[0].set_yscale('log'); ax[0].set_title("Experimental GO annotations per protein, by d2o")
sns.boxplot(data=master, x='d2o', y='length', ax=ax[1])
ax[1].set_yscale('log'); ax[1].set_title("Sequence length per protein, by d2o")
plt.tight_layout()
plt.savefig("../results/figures/annotation_bias.png", dpi=150)
print(master.groupby('d2o')[['n_go_exp','length']].describe())
```

- [ ] Are positives more heavily annotated than negatives? Record the median ratio and the IQR overlap. The limitation paragraph already exists in the proposal; this populates it with a real number.

### W4.6 — EDA — the comprehensive look (~2 hours)

In `03_eda.ipynb`, produce and save figures (each into `results/figures/`), and write a 1–2 sentence interpretation in the notebook directly under each figure.

- [ ] Class distribution (positives vs negatives) — bar chart, with percentages.
- [ ] Sequence length distribution — histogram, log-x, overlaid by class.
- [ ] GO Slim term frequency per aspect — top-20 bars × 3 aspects.
- [ ] Slim-terms-per-protein distribution, per aspect.
- [ ] Annotation-richness scatter: `n_go_exp` vs `length`, colored by class.
- [ ] Per-protein disordered-region count and total disordered residues (from DisProt) — distribution by class.
- [ ] CD-HIT cluster-size distribution.

The output of this notebook is *the* basis for the report's Results / Data section later.

### W4.7 — Go/no-go decision (~1 hour)

Apply these explicit criteria. **All must be met for "go".**

- [ ] **N sufficient:** ≥ 150 positives and ≥ 500 negatives after cleaning.
- [ ] **Class set correct:** binary; observed positive rate within a few points of the spike's 15%.
- [ ] **Sequence coverage:** ≥ 95% of candidate accessions returned a usable sequence.
- [ ] **GO coverage usable:** at least ~70% of proteins have ≥ 1 experimental GO annotation in at least one aspect.
- [ ] **GO Slim matrices non-degenerate:** each aspect has ≥ 20 columns; most columns have variance across the dataset (≥ 5–10 annotated proteins).
- [ ] **Cluster integrity:** positives span ≥ 50 distinct CD-HIT clusters, so cluster-aware CV is non-degenerate.
- [ ] **No fatal annotation-bias confound:** positives are not extreme outliers in `n_go_exp` vs negatives (median ratio within ~4×; flag-and-control, not stop, otherwise).

**If all met:** record the go decision in the lab notebook with numbers attached. Proceed to Week 5.

**If any fail:** pause. Document the failure precisely. Pull the supervisor meeting forward (don't wait for end of project) to discuss either tightening criteria or switching to the **Path B** disorder-function fallback target (~1,570 proteins, 387 distinct GO terms needing curated grouping).

### W4.8 — Closing tasks (30 min)

- [ ] Append a Week-4 entry to `lab-notebook.md` summarising the cleaning steps, the EDA highlights, the annotation-bias number, and the go/no-go verdict with the full numbers.
- [ ] Update §5.1 of the proposal with the *final* filtered counts (positives, negatives, clusters). This becomes part of the post-Week-4 supervisor meeting agenda.
- [ ] Commit and push everything; tag the commit `week4-go-nogo`.

### Week 4 deliverable

- Cleaned master dataset (`master_clean.csv`, with cluster IDs joined in).
- Three GO Slim feature matrices.
- `02_data_cleaning.ipynb` and `03_eda.ipynb`, both clean and reproducible.
- `results/figures/*.png` for the EDA.
- Lab-notebook entries for Weeks 3 and 4 with full numbers.
- A clear go/no-go verdict.

---

## Notes on what could trip you up

- **UniProt sometimes returns fewer sequences than you ask for** if accessions are obsolete or have been merged. That's normal; record the rate and continue.
- **goatools occasionally chokes on obsolete GO terms.** The `try/except` in the `slim_terms_for` function handles this; the resulting empty list just means that one annotation can't be mapped — fine.
- **GO Slim ≠ GO.** Don't compare slim-feature predictions against literature numbers based on full GO; they're different feature spaces.
- **CD-HIT identity choice matters.** 40% is conservative; 50% is more conventional. If 40% leaves too few clusters, try 50% and document the choice.
- **Don't commit large data files to git.** The `.gitignore` should already cover `data/*`; double-check `git status` before committing.

## Where this hands off to Week 5

Week 5 will: produce ESM-2 embeddings for the cleaned set, build the eight-condition feature configurations, run the thin-slice end-to-end pipeline on the simplest baseline, and confirm the modelling machinery works before Weeks 6–8.
