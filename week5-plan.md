# Week 5 Work Plan — Sequence Features and Thin-Slice Run

**Owner:** sunav
**Compiled:** 2026-06-05
**Covers:** Week 5. Target end-state is two sets of saved sequence features (ESM-2 embeddings and a simple amino-acid-composition baseline) and a *working* end-to-end modelling pipeline that produces one cross-validated AUPRC score on the sequence-only baseline. The week's purpose is to **prove the machinery works** before any of the eight factorial conditions are run.

Tick boxes as you go. Append observations to the lab notebook at the end.

## Top-level goals

By the end of Week 5 the project should have:

- `data/features_esm2.npz` (or `.pkl`) — ESM-2 mean-pooled embeddings for all 1,279 proteins.
- `data/features_aacomp.csv` — amino-acid composition feature matrix (reference baseline; not the real baseline).
- `notebooks/04_features_esm2.ipynb`, `notebooks/05_features_aacomp.ipynb`, `notebooks/06_thin_slice.ipynb` — all clean, commented, reproducible.
- **A working end-to-end thin slice:** Random Forest trained on a cheap feature set under cluster-aware GroupKFold CV, producing one AUPRC + macro-F1 + balanced-accuracy number that is clearly above chance (0.147 baseline rate). This is the proof the machinery is correct — not a real result.
- Updated `data/master_clean.csv` if any rows had to be dropped (e.g., embedding failures).
- Lab-notebook Week-5 entry recording: ESM-2 checkpoint chosen, truncation policy, thin-slice CV scores, time taken.

## Tools to install up front (~15 min)

Activate the env and install the deep-learning stack:

```
conda activate biol466
pip install torch fair-esm
```

Verify everything loads and detect whether your Mac has the Metal Performance Shaders backend for GPU-accelerated PyTorch:

```
python - <<'PY'
import torch, esm
print("torch:", torch.__version__)
print("MPS available:", torch.backends.mps.is_available())
print("CUDA available:", torch.cuda.is_available())
print("Default device:", "mps" if torch.backends.mps.is_available()
                        else "cuda" if torch.cuda.is_available() else "cpu")
PY
```

On Apple Silicon (M1/M2/M3/M4) the MPS line should say `True`. If you're on Intel Mac, both will say `False` and we fall back to CPU — slower but still fine for 1,279 sequences.

---

## W5.1 — Notebook scaffolding (~5 min)

- [ ] In `notebooks/`, create three empty notebooks: `04_features_esm2.ipynb`, `05_features_aacomp.ipynb`, `06_thin_slice.ipynb`.
- [ ] Top cell of each:

```python
import pandas as pd, numpy as np, matplotlib.pyplot as plt
pd.set_option('display.max_columns', 60)
```

---

## W5.2 — ESM-2 embedding extraction (~2–3 hours total: ~30 min writing, ~20 min running, ~1 hour buffer)

We use the **ESM-2 650M** checkpoint (`esm2_t33_650M_UR50D`). Reasoning: the 650M model is the standard quality/tractability balance — bigger 3B and 15B checkpoints are marginal-quality improvements that require GPU memory we don't have; smaller 150M and 35M checkpoints are noticeably weaker on downstream tasks. 650M produces **1,280-dimensional per-residue embeddings**, which we mean-pool to a **single 1,280-dim vector per protein**.

**Truncation policy.** ESM-2 is trained with a 1,024-token context. Sequences longer than 1,022 residues (minus CLS + EOS) will be truncated to the first 1,022 residues. About 5% of our 1,279 proteins exceed this — including titin at 34,350 aa. We accept the truncation for now and note it as a methods detail. (Alternative: sliding-window pooling. Stretch goal if time allows.)

In `04_features_esm2.ipynb`:

```python
import torch, esm, pickle, os
import pandas as pd, numpy as np
from tqdm.auto import tqdm

device = "mps" if torch.backends.mps.is_available() else \
         ("cuda" if torch.cuda.is_available() else "cpu")
print("Using device:", device)

# Load ESM-2 650M (downloads ~2.5 GB on first run; cached after)
print("Loading ESM-2 650M ...")
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
batch_converter = alphabet.get_batch_converter()
model.eval()
model = model.to(device)
print("Loaded. Parameter count:", sum(p.numel() for p in model.parameters())/1e6, "M")

# Load sequences
seqs = pd.read_csv("../data/sequences.csv")
print(f"Embedding {len(seqs)} sequences ...")

MAX_LEN = 1022   # leaving room for CLS + EOS to total 1024

@torch.no_grad()
def embed_one(seq):
    s = seq[:MAX_LEN]
    _, _, toks = batch_converter([("p", s)])
    toks = toks.to(device)
    out = model(toks, repr_layers=[33])
    rep = out["representations"][33][0, 1:-1]   # strip CLS + EOS
    return rep.mean(dim=0).cpu().numpy().astype("float32")

embeddings = {}
for _, row in tqdm(seqs.iterrows(), total=len(seqs)):
    try:
        embeddings[row['acc']] = embed_one(row['sequence'])
    except Exception as e:
        print(f"  Failed for {row['acc']}: {type(e).__name__}: {e}")

# Save as .npz so we can mmap later if needed
accs = list(embeddings.keys())
X = np.stack([embeddings[a] for a in accs])
np.savez_compressed("../data/features_esm2.npz", accs=np.array(accs), X=X.astype("float32"))
print(f"\nSaved {X.shape[0]} embeddings of dim {X.shape[1]}")
print(f"  Truncated proteins (length > {MAX_LEN}): "
      f"{(seqs['length'] > MAX_LEN).sum()}")
```

- [ ] First run downloads the checkpoint (~2.5 GB). Subsequent runs use the cached copy.
- [ ] Watch progress with `tqdm`. Expect ~20–30 minutes on MPS, ~60–90 min on CPU.
- [ ] Record in lab notebook: device used (MPS/CUDA/CPU), wall-clock time, number truncated.
- [ ] Final shape should be (1279, 1280) (or 1277/1278 if a couple of embeddings failed).

### If you run into RAM issues

If memory pressure errors appear, restart the kernel and process in chunks of 200 with explicit `torch.mps.empty_cache()` (or `torch.cuda.empty_cache()`) between chunks. Ping me — happy to write the chunked version.

### Quick sanity check after the embeddings save

Open a new cell:

```python
data = np.load("../data/features_esm2.npz", allow_pickle=True)
X, accs = data['X'], data['accs']
print("shape:", X.shape, "dtype:", X.dtype)
print("first acc:", accs[0])
print("first 5 dims of first embedding:", X[0, :5])
print("L2 norm range:", np.linalg.norm(X, axis=1).min(),
                        np.linalg.norm(X, axis=1).max())
```

Embeddings should have non-zero L2 norm (typical range ~30–60). Any zero-norm vectors are failed embeddings — fix or drop them.

---

## W5.3 — Amino-acid composition baseline representation (~30 min)

This is **not** the real baseline (ESM-2 is the baseline per the proposal's strong-baseline principle) — but it's a useful reference point to demonstrate, in the report, how much stronger ESM-2 is. Also the right feature set for the thin slice in W5.4 because it's fast.

In `05_features_aacomp.ipynb`:

```python
import pandas as pd, numpy as np

AA20 = list("ACDEFGHIKLMNPQRSTVWY")

seqs = pd.read_csv("../data/sequences.csv")

def aa_composition(seq):
    seq = seq.upper()
    L = max(len(seq), 1)
    return [seq.count(a) / L for a in AA20]

rows = []
for _, r in seqs.iterrows():
    rows.append([r['acc']] + aa_composition(r['sequence']) + [len(r['sequence'])])

cols = ['acc'] + [f"aa_{a}" for a in AA20] + ['length']
aa_df = pd.DataFrame(rows, columns=cols).set_index('acc')

# Z-score the length column (the other columns are already 0-1 frequencies)
aa_df['length_z'] = (np.log1p(aa_df['length']) - np.log1p(aa_df['length']).mean()) / np.log1p(aa_df['length']).std()
aa_df.drop(columns=['length'], inplace=True)

aa_df.to_csv("../data/features_aacomp.csv")
print("AA composition matrix:", aa_df.shape)
print(aa_df.head())
print("\nSanity: row sums (should be ~1 for the 20 AA columns) — first 5:")
print(aa_df.iloc[:5, :20].sum(axis=1).round(3).tolist())
```

- [ ] Output: 1,279 × 21 matrix (20 amino-acid frequencies + log-length z-score).
- [ ] Row sums of the 20 AA columns should equal 1.0 (within rounding).

---

## W5.4 — The thin-slice modelling pipeline (~3 hours)

This is the load-bearing piece of Week 5. The goal: a fully working pipeline — load features, run cluster-aware CV, compute the stats — that produces **one** trustworthy AUPRC number on the sequence-only condition (using the cheap AA-composition features). When this works, swapping in ESM-2 in Week 7 is a one-line change.

In `06_thin_slice.ipynb`:

```python
import numpy as np, pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GroupKFold, RepeatedStratifiedKFold
from sklearn.metrics import (average_precision_score, f1_score,
                             balanced_accuracy_score, roc_auc_score)

# ---- Load everything aligned by acc ----
master = pd.read_csv("../data/master_clean.csv")
aa = pd.read_csv("../data/features_aacomp.csv", index_col=0)

# Align AA features to master order (drop any proteins without features, drop NaN clusters)
keep = master['acc'].isin(aa.index) & master['cluster'].notna()
master = master[keep].reset_index(drop=True)
X = aa.reindex(master['acc']).values
y = master['d2o'].values.astype(int)
groups = master['cluster'].astype(int).values

print(f"Modelling on {len(master)} proteins: "
      f"{y.sum()} positives, {(y==0).sum()} negatives; "
      f"{len(np.unique(groups))} unique clusters")
print(f"Feature shape: {X.shape}")

# ---- Cluster-aware 5-fold CV ----
cv = GroupKFold(n_splits=5)
fold_scores = []
for fold, (tr, te) in enumerate(cv.split(X, y, groups)):
    n_pos_tr = y[tr].sum(); n_pos_te = y[te].sum()
    n_clusters_tr = len(np.unique(groups[tr]))
    n_clusters_te = len(np.unique(groups[te]))
    # Verify no cluster spans train/test
    assert len(set(groups[tr]) & set(groups[te])) == 0, "Cluster leakage!"

    clf = RandomForestClassifier(
        n_estimators=300, class_weight='balanced',
        n_jobs=-1, random_state=42)
    clf.fit(X[tr], y[tr])
    proba = clf.predict_proba(X[te])[:, 1]
    pred = (proba >= 0.5).astype(int)

    fold_scores.append({
        'fold': fold,
        'n_train': len(tr), 'n_test': len(te),
        'pos_train': n_pos_tr, 'pos_test': n_pos_te,
        'clusters_train': n_clusters_tr, 'clusters_test': n_clusters_te,
        'auprc': average_precision_score(y[te], proba),
        'auroc': roc_auc_score(y[te], proba),
        'macro_f1': f1_score(y[te], pred, average='macro'),
        'bal_acc': balanced_accuracy_score(y[te], pred),
    })
    print(f"Fold {fold}: train {len(tr)} ({n_pos_tr} pos, {n_clusters_tr} clusters) → "
          f"test {len(te)} ({n_pos_te} pos, {n_clusters_te} clusters) → "
          f"AUPRC={fold_scores[-1]['auprc']:.3f}, "
          f"macro_F1={fold_scores[-1]['macro_f1']:.3f}, "
          f"bal_acc={fold_scores[-1]['bal_acc']:.3f}")

results = pd.DataFrame(fold_scores)
print("\n=== Aggregate (mean ± std across 5 folds) ===")
for col in ['auprc', 'auroc', 'macro_f1', 'bal_acc']:
    print(f"  {col}: {results[col].mean():.3f} ± {results[col].std():.3f}")

results.to_csv("../results/thin_slice_aacomp.csv", index=False)

# ---- Reference: random-guessing baseline ----
print(f"\nReference AUPRC for a random classifier at this prevalence: {y.mean():.3f}")
print(f"  → any AUPRC meaningfully above this means the pipeline learns something.")
```

What to look for in the output:

- [ ] **Fold-level sanity:** every fold has positives in the test set (no fold should have 0 positives — that would break AUPRC). Each fold should test ~250 proteins (1,279 / 5) with ~38 positives (188 / 5).
- [ ] **Cluster leakage assert:** the `assert` line should never fire. If it does, the GroupKFold setup is wrong.
- [ ] **AUPRC well above 0.147:** any score above ~0.25 means AA composition contains real signal. Below 0.20 and it's essentially noise (which is also informative — would tell us AA composition is a genuinely weak baseline, which is fine; ESM-2 is the real baseline).
- [ ] **macro-F1 well above 0.50:** for a 14.7%-prevalence problem, a model predicting "always negative" scores macro-F1 = 0.43. So macro-F1 ≥ 0.55 is meaningful.
- [ ] **Reasonable spread across folds:** if std is comparable to mean, you have a high-variance pipeline (small N per fold). Acceptable for thin slice; we'll do repeated CV in Week 8.

---

## W5.5 — Sanity checks and the verdict (~30 min)

Before declaring Week 5 done, work through this short list. Each one catches a different class of bug.

- [ ] **Embedding dim correct:** `features_esm2.npz` is shape (1,279 or close, 1,280).
- [ ] **Embedding norms non-zero:** all L2 norms in the 10–100 range. Zero-norm vectors are silent failures.
- [ ] **No NaNs anywhere:** `np.isfinite(X).all()` on both ESM-2 and AA composition matrices.
- [ ] **Cluster mapping intact:** `cluster` values from `master_clean.csv` exist for every protein.
- [ ] **CV makes biological sense:** when you print the test-fold acc lists, no two folds should share a cluster (the assert protects against this).
- [ ] **Thin-slice AUPRC > 0.25:** if much lower, debug before Week 6. If much higher (say > 0.6), look for leakage — too good is as suspicious as too bad.
- [ ] **Random-shuffle sanity check (optional, recommended):** rerun the thin slice with `y` shuffled. AUPRC should drop to ≈ 0.147 (chance). If it doesn't, there's leakage somewhere in the feature pipeline.

When all check out:

- [ ] Commit and push everything (no large model checkpoints — those are cached outside the repo).
- [ ] Append a Week-5 entry to `lab-notebook.md`: device used, total embedding time, truncation count, thin-slice CV scores, any sanity-check observations.
- [ ] Tag the commit `week5-thin-slice`.

---

## Notes on what could trip you up

- **ESM-2 download is large** (~2.5 GB) and is cached at `~/.cache/torch/hub/checkpoints/`. If the download is interrupted, delete the partial file and retry.
- **MPS quirks.** PyTorch's MPS backend is solid for inference but can be picky about dtype. If you see `MPSNDArrayDescriptor` errors, cast to float32 explicitly: `model.float()` after `model.to("mps")`.
- **The `cluster` column dtype.** Pandas may load it as float if any NaN exists. Cast with `.astype(int)` after dropping NaNs, or `GroupKFold` will treat each row as a unique group and split degenerate folds.
- **Don't pickle the model.** Save *embeddings*, not the ESM-2 weights themselves. We re-load weights from cache each kernel restart; that's intentional and reproducible.
- **Memory headroom on 8 GB Macs.** If MPS runs out of memory mid-loop, the kernel will silently restart in JupyterLab. Symptom: the cell stops mid-progress with no traceback. Fix: process in chunks, or fall back to CPU.

---

## Where this hands off to Week 6

Week 6 builds the **GO Slim feature matrices into the modelling pipeline**: the slim matrices are already built (W4.2), but the thin slice in W5.4 only uses AA composition. Week 6 wires them up as a function `make_X(condition_id)` that returns the concatenated feature matrix for any of the eight conditions:

| Condition | Sequence | BP | MF | CC |
|---|:---:|:---:|:---:|:---:|
| 1 (control) | ESM-2 | — | — | — |
| 2 | ESM-2 | ✓ | — | — |
| 3 | ESM-2 | — | ✓ | — |
| 4 | ESM-2 | — | — | ✓ |
| 5 | ESM-2 | ✓ | ✓ | — |
| 6 | ESM-2 | ✓ | — | ✓ |
| 7 | ESM-2 | — | ✓ | ✓ |
| 8 (all) | ESM-2 | ✓ | ✓ | ✓ |

Then Week 7 extends the thin-slice pipeline to run all eight conditions with one shared GroupKFold split. Week 8 runs the factorial end-to-end and collects per-fold scores. Weeks 9–10 do the statistics, negative control, stratified analysis and interpretation.
