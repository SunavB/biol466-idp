# Week 7 Work Plan — Robustness, Stratification, Interpretation

**Owner:** sunav
**Compiled:** 2026-06-05
**Covers:** Week 7. Path II routing from W6.5 — the protein-level mean-pooled null needs to be stress-tested with a finer-grained sequence representation (per-residue ESM-2 pooled over annotated disordered regions only), the negative control needs to be run for the pre-registered statistical protocol, the H4 stratified analysis needs to be completed, and the feature-importance interpretation needs to surface which GO Slim terms RF found useful even with null overall effect. The week also begins drafting the Methods section of the final report.

We are now a week ahead of the proposal's original timeline (Weeks 9–10 stats are mostly done; the spike compressed Weeks 3–5). Week 7 mops up the robustness and interpretation work that was originally Weeks 9–10, and starts on writing.

Tick boxes as you go. Append observations to the lab notebook at the end.

## Top-level goals

By the end of Week 7 the project should have:

- `data/features_esm2_disorder.npz` — per-residue ESM-2 embeddings mean-pooled over DisProt-annotated disordered residues only. Same 1,280 dimensions per protein, different pooling.
- `results/factorial_rf_per_fold_disorder.csv` — 8-condition RF factorial under the new representation, same shared splits.
- `results/factorial_rf_stats_disorder.csv` — paired stats / BH correction / effect sizes for the disorder-pooled setup.
- `results/negative_control_per_fold.csv` — full factorial under shuffled labels (sanity check that AUPRC ≈ 0.147 across the board).
- `results/stratified_h4.csv` — Seq-only vs all-GO AUPRC partitioned by per-protein Seq-only confidence tertile, addresses H4.
- `results/feature_importance_top20.csv` — top 20 features (most often ESM-2 dims, but watch the GO Slim terms that rank high) under the full-GO RF.
- `results/figures/` populated with: disorder-pooled factorial bar chart; negative-control histogram; stratified-by-tertile chart; feature-importance plot.
- A Methods-section first draft of the final report (`report/methods.md` or directly in a `report.docx`).
- Lab-notebook Week-7 entry recording: whether the null replicates under disorder pooling, the negative-control verdict, H4 status, the most informative GO terms (even with null overall effect), and the report draft progress.

## Tools (no new installs)

Everything needed is already in the env from Weeks 5–6. Just confirm:

```
python -c "import torch, esm, xgboost, sklearn, scipy, statsmodels; print('ok')"
```

---

## W7.1 — Disorder region extraction and notebook scaffolding (~30 min)

Create `notebooks/08_disorder_embeddings.ipynb`. Top cell:

```python
import pandas as pd, numpy as np, matplotlib.pyplot as plt
pd.set_option('display.max_columns', 60)
```

Then load DisProt and extract per-protein disordered region coordinates:

```python
disprot = pd.read_csv("../data/disprot.tsv", sep='\t', low_memory=False)
disprot = disprot.rename(columns={
    'UniProt ACC':'acc', 'Organism':'organism',
    'Term namespace':'term_namespace', 'Start':'start', 'End':'end',
    'Term name':'term_name',
})
disprot['acc'] = disprot['acc'].str.split('-').str[0]

# Structural state annotations = disordered regions
state = disprot[(disprot['organism']=='Homo sapiens') &
                (disprot['term_namespace']=='Structural state')].copy()
print(f"Disorder annotations (Structural state, human): {len(state)}")
print(f"Distinct proteins with disorder annotations: {state['acc'].nunique()}")

# Build per-acc list of (start, end) tuples (1-indexed, inclusive)
disorder_regions = (state.groupby('acc')
                         .apply(lambda x: [(int(s), int(e)) for s, e in zip(x['start'], x['end'])])
                         .to_dict())

# Sanity: distribution of disordered residues per protein
n_dis = {acc: sum(e - s + 1 for s, e in regs)
         for acc, regs in disorder_regions.items()}
print(f"\nDisordered residues per protein:")
print(pd.Series(n_dis).describe().round(1))
```

- [ ] Expect ~1,250–1,279 distinct proteins with at least one disorder annotation. Confirm.
- [ ] Median disordered residue count per protein is probably 30–80; this is what gets averaged in the new representation.

---

## W7.2 — Per-residue ESM-2 with disorder-region pooling (~2–3 hours, mostly compute)

Same ESM-2 650M model, same MPS device, same truncation at 1,022 residues. The change: instead of mean-pooling over **all** residue tokens, we mean-pool over **only the residues that fall inside DisProt-annotated disordered regions**.

```python
import torch, esm
from tqdm.auto import tqdm

device = "mps" if torch.backends.mps.is_available() else \
         ("cuda" if torch.cuda.is_available() else "cpu")
print("Using device:", device)

print("Loading ESM-2 650M (cached from W5) ...")
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
batch_converter = alphabet.get_batch_converter()
model.eval()
model = model.to(device)

seqs = pd.read_csv("../data/sequences.csv")
MAX_LEN = 1022

@torch.no_grad()
def embed_disorder(seq, acc):
    """Mean-pool ESM-2 token reps over DisProt disordered residues only.
       Returns (embedding, fallback_flag). fallback_flag=True means we
       fell back to whole-protein pooling because no disorder residues
       were in the truncation window."""
    s = seq[:MAX_LEN]
    _, _, toks = batch_converter([("p", s)])
    toks = toks.to(device)
    out = model(toks, repr_layers=[33])
    # representations[33] is (1, L+2, 1280); strip CLS at idx 0 and EOS at idx -1
    per_residue = out["representations"][33][0, 1:-1]   # shape (L, 1280)
    L = per_residue.shape[0]

    regions = disorder_regions.get(acc, [])
    if not regions:
        return per_residue.mean(dim=0).cpu().numpy().astype("float32"), True

    # Build residue indices (0-indexed) within truncation window
    indices = []
    for start, end in regions:
        for i in range(start - 1, min(end, L)):  # 1-indexed -> 0-indexed
            indices.append(i)

    if not indices:
        return per_residue.mean(dim=0).cpu().numpy().astype("float32"), True

    idx_t = torch.tensor(indices, device=device, dtype=torch.long)
    selected = per_residue[idx_t]
    return selected.mean(dim=0).cpu().numpy().astype("float32"), False

embeddings = {}
fallbacks = []
for _, row in tqdm(seqs.iterrows(), total=len(seqs)):
    try:
        emb, fb = embed_disorder(row['sequence'], row['acc'])
        embeddings[row['acc']] = emb
        if fb:
            fallbacks.append(row['acc'])
    except Exception as e:
        print(f"  Failed for {row['acc']}: {type(e).__name__}: {e}")

accs = list(embeddings.keys())
X = np.stack([embeddings[a] for a in accs])
np.savez_compressed("../data/features_esm2_disorder.npz",
                    accs=np.array(accs), X=X.astype("float32"))
print(f"\nSaved {X.shape[0]} disorder-pooled embeddings of dim {X.shape[1]}")
print(f"Fallback to whole-protein pooling (no disorder in truncation window): "
      f"{len(fallbacks)} proteins")
```

- [ ] Watch the tqdm bar; same ~20–30 min on MPS as W5.
- [ ] Expect ≤ ~10 proteins to fall back to whole-protein pooling — those with all disorder regions beyond residue 1022 (very long proteins).
- [ ] Final shape (1279, 1280) — same dimensionality as W5, different pooling semantics.

### Quick sanity check after save

```python
data = np.load("../data/features_esm2_disorder.npz", allow_pickle=True)
X, accs = data['X'], data['accs']
norms = np.linalg.norm(X, axis=1)
print("shape:", X.shape, "dtype:", X.dtype)
print(f"L2 norms: min={norms.min():.2f}, max={norms.max():.2f}, mean={norms.mean():.2f}")
print(f"NaNs? {np.isnan(X).any()}; zero-norm vectors? {(norms == 0).sum()}")

# Compare to whole-protein embeddings: cosine similarity per protein
data_old = np.load("../data/features_esm2.npz", allow_pickle=True)
X_old = data_old['X']
accs_old = list(data_old['accs'])
idx_old = {a: i for i, a in enumerate(accs_old)}
common = [a for a in accs if a in idx_old]
X_new_c = np.stack([X[list(accs).index(a)] for a in common])
X_old_c = np.stack([X_old[idx_old[a]] for a in common])

def cosine(a, b):
    return (a * b).sum(axis=1) / (np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1))

cs = cosine(X_new_c, X_old_c)
print(f"\nCosine similarity between disorder-pooled and whole-protein embeddings:")
print(f"  median={np.median(cs):.3f}, mean={cs.mean():.3f}, min={cs.min():.3f}, max={cs.max():.3f}")
```

- [ ] Expect median cosine ~0.85–0.95 (similar but not identical). If max = 1.0 for many entries, those are proteins where disordered region ≈ whole protein (high disorder content), which is fine. If median is below 0.7, the pooling is doing something very different — paste me the numbers.

---

## W7.3 — Re-run the 8-condition factorial under the disorder-pooled representation (~1 hour compute)

In `notebooks/09_factorial_disorder.ipynb`. Reuse the structure from `07_factorial.ipynb` — just swap the embedding source. Top cell loads everything, replaces `esm_X_all` with the disorder-pooled version, redefines `make_X`, then runs the 8-condition loop on the shared CV splits.

```python
import pandas as pd, numpy as np, matplotlib.pyplot as plt, os
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GroupKFold
from sklearn.metrics import (average_precision_score, f1_score,
                             balanced_accuracy_score, roc_auc_score,
                             matthews_corrcoef)

os.makedirs("../results/figures", exist_ok=True)

# Load disorder-pooled ESM-2
esm_data = np.load("../data/features_esm2_disorder.npz", allow_pickle=True)
esm_X_all = esm_data['X']
esm_accs = list(esm_data['accs'])
esm_idx = {acc: i for i, acc in enumerate(esm_accs)}

master = pd.read_csv("../data/master_clean.csv")
keep = master['acc'].isin(esm_accs) & master['cluster'].notna()
master = master[keep].reset_index(drop=True)
y = master['d2o'].values.astype(int)
groups = master['cluster'].astype(int).values

bp = pd.read_csv("../data/features_GO_BP_slim.csv", index_col=0)
mf = pd.read_csv("../data/features_GO_MF_slim.csv", index_col=0)
cc = pd.read_csv("../data/features_GO_CC_slim.csv", index_col=0)

CONDITIONS = {
    1: {'label': 'Seq only',          'bp': False, 'mf': False, 'cc': False},
    2: {'label': 'Seq + BP',          'bp': True,  'mf': False, 'cc': False},
    3: {'label': 'Seq + MF',          'bp': False, 'mf': True,  'cc': False},
    4: {'label': 'Seq + CC',          'bp': False, 'mf': False, 'cc': True },
    5: {'label': 'Seq + BP + MF',     'bp': True,  'mf': True,  'cc': False},
    6: {'label': 'Seq + BP + CC',     'bp': True,  'mf': False, 'cc': True },
    7: {'label': 'Seq + MF + CC',     'bp': False, 'mf': True,  'cc': True },
    8: {'label': 'Seq + BP + MF + CC','bp': True,  'mf': True,  'cc': True },
}

def make_X(condition_id):
    cfg = CONDITIONS[condition_id]
    X_seq = np.stack([esm_X_all[esm_idx[a]] for a in master['acc']])
    parts = [X_seq]
    for mat, present in [(bp, cfg['bp']), (mf, cfg['mf']), (cc, cfg['cc'])]:
        if present:
            parts.append(mat.reindex(master['acc']).fillna(0).values)
    return np.hstack(parts).astype(np.float32)

cv = GroupKFold(n_splits=5)
splits = list(cv.split(np.zeros(len(master)), y, groups))
THRESHOLD = float(y.mean())

results_dis = []
for cid in range(1, 9):
    cfg = CONDITIONS[cid]
    X = make_X(cid)
    print(f"\n=== Disorder-pool Condition {cid}: {cfg['label']} (X.shape {X.shape}) ===")
    for fold, (tr, te) in enumerate(splits):
        clf = RandomForestClassifier(
            n_estimators=500, class_weight='balanced',
            min_samples_leaf=3, n_jobs=-1, random_state=42)
        clf.fit(X[tr], y[tr])
        proba = clf.predict_proba(X[te])[:, 1]
        pred = (proba >= THRESHOLD).astype(int)
        results_dis.append({
            'condition': cid, 'label': cfg['label'], 'fold': fold,
            'auprc':    average_precision_score(y[te], proba),
            'auroc':    roc_auc_score(y[te], proba),
            'macro_f1': f1_score(y[te], pred, average='macro'),
            'bal_acc':  balanced_accuracy_score(y[te], pred),
            'mcc':      matthews_corrcoef(y[te], pred),
        })
        r = results_dis[-1]
        print(f"  Fold {fold}: AUPRC={r['auprc']:.3f}, F1={r['macro_f1']:.3f}, MCC={r['mcc']:.3f}")

results_dis = pd.DataFrame(results_dis)
results_dis.to_csv("../results/factorial_rf_per_fold_disorder.csv", index=False)

agg = (results_dis.groupby(['condition','label'])
                  [['auprc','auroc','macro_f1','bal_acc','mcc']]
                  .agg(['mean','std']).round(3))
print(agg.to_string())
```

- [ ] **Headline check:** does **Condition 1 AUPRC** under disorder pooling beat 0.199 (the W6 whole-protein baseline)? If yes, the disorder-specific representation is genuinely more informative — even before adding GO. If no, the protein-level mean-pooling wasn't the bottleneck.
- [ ] **Condition 8 vs Condition 1 under disorder pooling:** is the GO effect any different than the W6 null? If still null, the null is robust across representations.

Save a comparison plot:

```python
# Compare W6 (whole protein) vs W7 (disorder pool) AUPRC per condition
results_rf = pd.read_csv("../results/factorial_rf_per_fold.csv")
w6 = results_rf.groupby('condition')['auprc'].agg(['mean','std'])
w7 = results_dis.groupby('condition')['auprc'].agg(['mean','std'])

fig, ax = plt.subplots(figsize=(11, 5))
xs = np.arange(1, 9)
ax.bar(xs - 0.2, w6['mean'], 0.4, yerr=w6['std'], label='Whole-protein pool (W6)',
       color='#888', capsize=3)
ax.bar(xs + 0.2, w7['mean'], 0.4, yerr=w7['std'], label='Disorder pool (W7)',
       color='#3a7', capsize=3)
ax.axhline(y.mean(), ls='--', color='gray', label=f'chance ({y.mean():.3f})')
ax.set_xticks(xs); ax.set_xticklabels([CONDITIONS[c]['label'] for c in xs], rotation=20, ha='right')
ax.set_ylabel("AUPRC (mean ± std)")
ax.set_title("Factorial AUPRC — whole-protein vs disorder-only pooling")
ax.legend()
plt.tight_layout()
plt.savefig("../results/figures/factorial_pool_compare.png", dpi=150)
plt.show()
```

---

## W7.4 — Pairwise stats on the disorder-pooled factorial (~30 min)

Identical to W6.3 but on `results_dis`. Reuse the code; just point it at the new wide table.

```python
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests

wide = results_dis.pivot(index='fold', columns='condition', values='auprc')
control = wide[1].values
rows = []
for cid in range(2, 9):
    other = wide[cid].values
    diff = other - control
    w, p = wilcoxon(other, control, alternative='greater', zero_method='wilcox')
    rng = np.random.default_rng(42)
    boots = [diff[rng.integers(0, len(diff), len(diff))].mean() for _ in range(2000)]
    lo, hi = np.percentile(boots, [2.5, 97.5])
    rows.append({'condition': cid, 'label': CONDITIONS[cid]['label'],
                 'mean_diff_vs_seq': diff.mean(), 'ci95_lo': lo, 'ci95_hi': hi,
                 'p_raw': p})
pair = pd.DataFrame(rows)
pair['p_BH'] = multipletests(pair['p_raw'], method='fdr_bh')[1]
print(pair.round(4).to_string(index=False))
pair.to_csv("../results/factorial_rf_stats_disorder.csv", index=False)
```

- [ ] **Verdict:** if any condition reaches BH-p < 0.05 here when it didn't in W6, the disorder representation rescued some GO benefit — that's a meaningful methodological finding. Most likely outcome: still null.

---

## W7.5 — Negative control: label shuffling (~1.5 hours)

The proposal-promised negative control. Shuffle d2o labels across proteins **once** with a fixed seed, then re-run the 8-condition factorial. Expected behaviour: AUPRC ≈ prevalence (0.147) for every condition; significant departure from this would indicate leakage somewhere in the pipeline.

In a new notebook cell (or `notebooks/10_negative_control.ipynb`):

```python
rng = np.random.default_rng(2026)
y_shuffled = rng.permutation(y)
print(f"Shuffled labels: {y_shuffled.sum()} positives (preserved); "
      f"shuffled vs original agreement = {(y_shuffled == y).mean():.3f}")

results_neg = []
for cid in range(1, 9):
    cfg = CONDITIONS[cid]
    X = make_X(cid)
    for fold, (tr, te) in enumerate(splits):
        clf = RandomForestClassifier(
            n_estimators=500, class_weight='balanced',
            min_samples_leaf=3, n_jobs=-1, random_state=42)
        clf.fit(X[tr], y_shuffled[tr])
        proba = clf.predict_proba(X[te])[:, 1]
        pred = (proba >= THRESHOLD).astype(int)
        results_neg.append({
            'condition': cid, 'label': cfg['label'], 'fold': fold,
            'auprc':    average_precision_score(y_shuffled[te], proba),
            'auroc':    roc_auc_score(y_shuffled[te], proba),
        })

results_neg = pd.DataFrame(results_neg)
results_neg.to_csv("../results/negative_control_per_fold.csv", index=False)

print("\n=== Negative control: AUPRC under shuffled labels ===")
neg_agg = results_neg.groupby('condition')['auprc'].agg(['mean','std']).round(3)
neg_agg['condition_label'] = [CONDITIONS[c]['label'] for c in neg_agg.index]
print(neg_agg.to_string())
print(f"\nReference: chance AUPRC = {y.mean():.3f}")
```

- [ ] **What to look for:** every condition's mean AUPRC should sit within ±0.03 of 0.147. If any condition's shuffled AUPRC is meaningfully above this (say > 0.20), there is leakage somewhere in the pipeline — STOP and debug before continuing.

Plot the negative-control AUPRC distribution against the W6 real-data distribution to show visually that the real null overlaps the shuffled null (this is the strongest possible visual evidence that the real result is itself a null):

```python
fig, ax = plt.subplots(figsize=(9, 4))
real_means = results_rf.groupby('condition')['auprc'].mean()
neg_means  = results_neg.groupby('condition')['auprc'].mean()
xs = np.arange(1, 9)
ax.bar(xs - 0.2, real_means, 0.4, label='Real labels (W6)', color='#3a7')
ax.bar(xs + 0.2, neg_means, 0.4, label='Shuffled labels (negative control)', color='#aaa')
ax.axhline(y.mean(), ls='--', color='black', label=f'chance ({y.mean():.3f})')
ax.set_xticks(xs); ax.set_xticklabels([CONDITIONS[c]['label'] for c in xs], rotation=20, ha='right')
ax.set_ylabel('AUPRC (mean over 5 folds)')
ax.set_title('Real vs shuffled-label AUPRC by condition')
ax.legend()
plt.tight_layout()
plt.savefig("../results/figures/negative_control_vs_real.png", dpi=150)
plt.show()
```

---

## W7.6 — Stratified analysis for H4 (~1.5 hours)

H4 (per the proposal): the GO benefit is larger for proteins the sequence-only model predicts poorly.

The cleanest operationalisation:

1. For every protein, the W6 Condition-1 (Seq only) probability *when it was in the test fold* is its "Seq-only confidence."
2. Partition all proteins into **tertiles** by that confidence (low / mid / high).
3. Within each tertile, compute AUPRC for Condition 1 and Condition 8 separately.
4. H4 supported if the gap (Cond 8 − Cond 1) is largest in the low-confidence tertile.

```python
# Reload per-fold predictions; need to rebuild because we only saved metrics
# in W6. Rerun with proba storage for conditions 1 and 8 only.

protein_proba_seq, protein_proba_all, protein_y = [], [], []
protein_acc_order = master['acc'].tolist()

# We need the probability for each protein from the fold it was in test.
# Rebuild from scratch with stored probabilities:
for cond_to_track, store in [(1, protein_proba_seq), (8, protein_proba_all)]:
    X = make_X(cond_to_track)
    proba_per_protein = np.full(len(master), np.nan)
    for fold, (tr, te) in enumerate(splits):
        clf = RandomForestClassifier(
            n_estimators=500, class_weight='balanced',
            min_samples_leaf=3, n_jobs=-1, random_state=42)
        clf.fit(X[tr], y[tr])
        proba_per_protein[te] = clf.predict_proba(X[te])[:, 1]
    store.extend(proba_per_protein.tolist())

protein_proba_seq = np.array(protein_proba_seq)
protein_proba_all = np.array(protein_proba_all)

# Stratify by Seq-only confidence tertile
tertiles = pd.qcut(protein_proba_seq, q=3, labels=['low','mid','high'])

print("=== AUPRC by Seq-only confidence tertile ===")
print(f"{'tertile':<8} {'n':>5} {'pos':>5}  {'Seq AUPRC':>12} {'Seq+all-GO':>14} {'Δ':>8}")
strat_rows = []
for t in ['low', 'mid', 'high']:
    mask = (tertiles == t)
    n = mask.sum()
    n_pos = int(y[mask].sum())
    if n_pos < 5:
        print(f"{t:<8} {n:>5} {n_pos:>5}  (too few positives)")
        continue
    a_seq = average_precision_score(y[mask], protein_proba_seq[mask])
    a_all = average_precision_score(y[mask], protein_proba_all[mask])
    delta = a_all - a_seq
    print(f"{t:<8} {n:>5} {n_pos:>5}  {a_seq:>12.3f} {a_all:>14.3f} {delta:>+8.3f}")
    strat_rows.append({'tertile': t, 'n': n, 'n_pos': n_pos,
                       'auprc_seq': a_seq, 'auprc_all_go': a_all, 'delta': delta})

pd.DataFrame(strat_rows).to_csv("../results/stratified_h4.csv", index=False)
```

- [ ] **H4 verdict:** is the delta in the `low` tertile larger than the deltas in `mid` and `high`? If yes, partial support for H4. If deltas are all near zero (as expected from the W6 null), H4 is also falsified for this representation.

Also worth a scatter:

```python
fig, ax = plt.subplots(figsize=(7, 7))
colors = np.where(y == 1, '#d33', '#888')
ax.scatter(protein_proba_seq, protein_proba_all, c=colors, alpha=0.4, s=15)
mx = max(protein_proba_seq.max(), protein_proba_all.max())
ax.plot([0, mx], [0, mx], '--', color='gray', label='y=x')
ax.set_xlabel("Seq-only probability (Condition 1)")
ax.set_ylabel("Seq + all GO probability (Condition 8)")
ax.set_title("Per-protein probability shift, GO vs no-GO")
ax.legend(['y=x', 'positives (red)', 'negatives (grey)'])
plt.tight_layout()
plt.savefig("../results/figures/per_protein_proba_shift.png", dpi=150)
plt.show()
```

- [ ] If most points hug the y=x line, GO predictions are essentially identical to Seq-only — strong visual evidence for the null.

---

## W7.7 — Feature-importance interpretation (~1 hour)

Even with null overall effect, some GO Slim terms may individually be informative; documenting that is part of the proposal's promise.

```python
# Train Condition 8 RF on all data (not CV) for interpretation
X = make_X(8)
clf_full = RandomForestClassifier(
    n_estimators=500, class_weight='balanced',
    min_samples_leaf=3, n_jobs=-1, random_state=42)
clf_full.fit(X, y)

# Map column indices to feature names
n_esm = 1280
bp_names = ['BP_' + c for c in bp.columns]
mf_names = ['MF_' + c for c in mf.columns]
cc_names = ['CC_' + c for c in cc.columns]
feature_names = (['ESM_' + str(i) for i in range(n_esm)]
                 + bp_names + mf_names + cc_names)
assert len(feature_names) == X.shape[1]

imp = pd.DataFrame({
    'feature': feature_names,
    'importance': clf_full.feature_importances_,
    'kind': (['ESM-2']*n_esm + ['BP']*len(bp.columns)
             + ['MF']*len(mf.columns) + ['CC']*len(cc.columns)),
}).sort_values('importance', ascending=False)

imp.head(30).to_csv("../results/feature_importance_top30.csv", index=False)
print("=== Top 30 features ===")
print(imp.head(30).to_string(index=False))

print("\n=== Aggregate importance by feature kind ===")
agg_imp = imp.groupby('kind')['importance'].agg(['sum','mean','count']).round(4)
agg_imp['sum_normalised'] = (agg_imp['sum'] / agg_imp['sum'].sum()).round(3)
print(agg_imp.to_string())

# Top GO terms (excluding ESM)
top_go = imp[imp['kind'] != 'ESM-2'].head(15)
print("\n=== Top 15 GO Slim features ===")
print(top_go.to_string(index=False))
```

- [ ] Look for: are the top GO terms biologically sensible (binding-related MF terms, regulatory BP terms)? Record the top 5 with their biological meaning in the lab notebook.
- [ ] The aggregate-by-kind table tells you how much total importance RF distributes across the four feature groups. Even with null AUPRC effect, RF will split this importance somewhere; the relative distribution is informative.

Plot:

```python
fig, ax = plt.subplots(figsize=(8, 4))
ax.bar(agg_imp.index, agg_imp['sum_normalised'], color=['#37a','#3a7','#a37','#888'])
ax.set_ylabel('Fraction of total RF feature importance')
ax.set_title('Where Condition 8 RF puts its importance, by feature kind')
plt.tight_layout()
plt.savefig("../results/figures/feature_importance_by_kind.png", dpi=150)
plt.show()
```

---

## W7.8 — Begin drafting the Methods section (~1 hour)

Create `report/methods.md`. The draft writes itself from the lab notebook entries — pull the substance from:

- §5.1 of the proposal (dataset & target framing) — almost copy-paste with W3-W4 numbers
- W3 lab-notebook entry (data acquisition methods)
- W4 lab-notebook entry (cleaning, GO Slim, CD-HIT, EDA)
- W5 lab-notebook entry (sequence representations)
- W6 lab-notebook entry (factorial design, classifier, statistics)
- W7 lab-notebook entry (this week's robustness checks)

Recommended subsections, all 1–2 paragraphs each:

- [ ] **2.1 Dataset.** DisProt release 2025_12 → 1,279 candidate human proteins after isoform collapse → 188 / 1,091 binary d2o split.
- [ ] **2.2 Sequence features.** ESM-2 650M, mean-pooled per protein (1,280 dims). Truncation at 1,022. Plus disorder-pooled variant for the robustness check (W7).
- [ ] **2.3 Gene Ontology features.** GOA human GAF, experimental evidence codes only, GO Slim mapping via `goatools.mapslim`. Three matrices: BP (64 columns), MF (36), CC (25).
- [ ] **2.4 Sequence-redundancy handling.** CD-HIT at 40% identity, 1,168 clusters.
- [ ] **2.5 Factorial design.** Eight conditions, 2×2×2 BP/MF/CC, sequence representation common to all.
- [ ] **2.6 Classifier.** Random Forest (500 trees, `class_weight='balanced'`, `min_samples_leaf=3`). Threshold = prevalence (0.147). XGBoost as a robustness check (W6.4).
- [ ] **2.7 Cross-validation and evaluation.** Cluster-aware GroupKFold(5), splits shared across all 8 conditions for paired comparisons. Primary metric AUPRC. Secondary: macro-F1, balanced accuracy, MCC.
- [ ] **2.8 Statistical analysis.** Paired Wilcoxon vs Condition 1, Benjamini-Hochberg correction across 7 comparisons. Bootstrap 95% CIs. Factorial main and interaction effects.
- [ ] **2.9 Controls and stratification.** Label-shuffling negative control (W7.5). H4 stratified analysis by Seq-only confidence tertile (W7.6).
- [ ] **2.10 Limitations.** GO/DisProt residual correlation (1.76× n_go_exp ratio); disorder-richness confound (1.50×, 1.73×); MNAR negatives; 14% sequence truncation.

This is just a Methods draft for now — Results / Discussion / Introduction come in Weeks 8–11.

---

## W7.9 — Closing tasks (~30 min)

- [ ] Append a Week-7 entry to `lab-notebook.md` with: disorder-pooling cosine-similarity stats, disorder-pooled factorial table, pairwise stats verdict (null replicates?), negative-control aggregate, H4 stratified-tertile table, top 5 informative GO Slim terms, Methods-draft progress.
- [ ] Commit and push. Tag the commit `week7-robustness`.

---

## Notes on what could trip you up

- **Disorder coordinates are 1-indexed, inclusive.** The `state['start'] = 1, end = 50` annotation means residues 1 through 50 of the protein. The `embed_disorder` function above converts to 0-indexed Python slicing correctly — double-check by hand on `P49913` (the LL-37 entry we examined in W1: disorder region 134–170, length 170; should yield 37 residue indices).
- **Truncation interacts with disorder coords.** If a protein's only disorder annotation is at residues 1500–1800 and we truncate at 1022, no disorder residues are within the window and we fall back to whole-protein pooling. Record how many proteins this affects.
- **Stored probabilities for H4.** The W6 factorial code didn't save per-protein probabilities, so W7.6 has to rerun Condition 1 and Condition 8 specifically to populate `protein_proba_*`. The rerun uses the same `splits` and same random seed, so results are identical to W6 — no inconsistency. It just takes ~2 minutes extra.
- **MPS kernel restarts.** If the embedding loop dies silently on a long protein with the kernel-restart symptom, process in 200-protein chunks with `torch.mps.empty_cache()` between them.
- **Don't over-interpret feature importance.** RF importance is a contributions-to-tree-splits measure, not a causal one. The top GO terms tell you what RF found useful relative to ESM-2 in the all-features regime; they don't prove those terms drive d2o prediction in any deep sense.

---

## Where this hands off to Week 8

Week 8 finalises the experimental result and pushes hard on writing. Expected work:

- Run any additional robustness checks the Week-7 results suggest (e.g., per-aspect feature ablation if disorder pooling lifted things).
- Begin Results section with the figures from W6/W7 already in place.
- Begin Discussion section synthesising the H1–H4 verdicts and the null interpretation.
- First polished pass on Introduction and Abstract.
- Supervisor check-in (per the proposal's Week-1, 4, 7–8, 11 cadence) — this is the natural point for the post-Week-4 deferred meeting, given Week 7 has the robustness story and Week 8 has the writing in progress.
