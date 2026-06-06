# Week 6 Work Plan — Run the 8-Condition Factorial

**Owner:** sunav
**Compiled:** 2026-06-05
**Covers:** Week 6. Target end-state is the full 2×2×2 factorial result table — per-fold AUPRC/AUROC/macro-F1/balanced-accuracy for all 8 GO-feature configurations — plus the statistical analysis (paired tests, BH correction, effect sizes, CIs, main and interaction effects), and a gradient-boosting robustness check.

**This week is the main experimental result of the project.** We're a little ahead of the proposal timeline (the spike pulled forward Week 3, and Week 5's pipeline can immediately swap in GO features), so Week 6 runs what the proposal originally scheduled for Week 8.

Tick boxes as you go. Append observations to the lab notebook at the end.

## Top-level goals

By the end of Week 6 the project should have:

- A clean `make_X(condition_id)` function that builds the feature matrix for any of the 8 factorial conditions.
- Per-fold scores for all 8 conditions in `results/factorial_rf_per_fold.csv`.
- Statistical analysis output in `results/factorial_rf_stats.csv` — pairwise comparisons with BH correction, effect sizes, bootstrap CIs.
- Main and interaction effects estimated from the factorial structure, saved alongside.
- A figure (`results/figures/factorial_rf_auprc.png`) — AUPRC by condition with error bars, the visual headline.
- An optional gradient-boosting robustness check (`results/factorial_xgb_per_fold.csv`) — same 8 conditions, same splits, swapped classifier; results compared against RF.
- An informed call on whether per-residue ESM-2 features (over disordered regions only) need to be Week 7's focus.
- Lab-notebook Week-6 entry with the headline AUPRC numbers per condition, the main-effects ranking, and which hypotheses are supported / falsified so far.

## Tools to install up front (~5 min)

XGBoost for the robustness check (LightGBM is an equally good alternative if you prefer):

```
conda activate biol466
pip install xgboost
```

Verify:

```
python -c "import xgboost as xgb; print('xgb', xgb.__version__)"
```

---

## W6.1 — Notebook scaffolding and the `make_X` builder (~30 min)

In `notebooks/`, create `07_factorial.ipynb`. Top cell:

```python
import pandas as pd, numpy as np, matplotlib.pyplot as plt, seaborn as sns, os
pd.set_option('display.max_columns', 60)
os.makedirs("../results", exist_ok=True)
os.makedirs("../results/figures", exist_ok=True)
```

Then the loader + `make_X`:

```python
# ---- Load all artifacts ----
master = pd.read_csv("../data/master_clean.csv")

esm_data = np.load("../data/features_esm2.npz", allow_pickle=True)
esm_X_all = esm_data['X']
esm_accs = list(esm_data['accs'])
esm_idx = {acc: i for i, acc in enumerate(esm_accs)}

bp = pd.read_csv("../data/features_GO_BP_slim.csv", index_col=0)
mf = pd.read_csv("../data/features_GO_MF_slim.csv", index_col=0)
cc = pd.read_csv("../data/features_GO_CC_slim.csv", index_col=0)

print(f"ESM-2: {esm_X_all.shape}")
print(f"BP slim: {bp.shape}, MF slim: {mf.shape}, CC slim: {cc.shape}")

# Restrict to proteins with ESM-2 embedding AND valid cluster
keep = master['acc'].isin(esm_accs) & master['cluster'].notna()
master = master[keep].reset_index(drop=True)
y = master['d2o'].values.astype(int)
groups = master['cluster'].astype(int).values
print(f"\nFinal modelling set: {len(master)} proteins, "
      f"{y.sum()} positives, {len(np.unique(groups))} clusters")

# ---- The factorial-condition feature builder ----
# Bit-encoded condition_id: bit 0 = BP, bit 1 = MF, bit 2 = CC. Condition 1 = sequence only.
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
    for slim_mat, present in [(bp, cfg['bp']), (mf, cfg['mf']), (cc, cfg['cc'])]:
        if present:
            parts.append(slim_mat.reindex(master['acc']).fillna(0).values)
    return np.hstack(parts).astype(np.float32)

# Sanity-check dimensions
for cid in range(1, 9):
    X = make_X(cid)
    print(f"Condition {cid} ({CONDITIONS[cid]['label']}): X.shape = {X.shape}")
```

- [ ] Expect condition 1 to be (1279, 1280), condition 8 to be (1279, 1280 + 64 + 36 + 25) = (1279, 1405). Other conditions in between.

---

## W6.2 — Run all 8 conditions with shared CV splits (~1.5–2 hours of compute)

The critical design point: **all 8 conditions must use the same CV splits**, so the comparison is paired at the fold level. Predefine the splits once.

```python
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GroupKFold
from sklearn.metrics import (average_precision_score, f1_score,
                             balanced_accuracy_score, roc_auc_score,
                             matthews_corrcoef)

# Predefine CV splits ONCE so all conditions are paired by fold
cv = GroupKFold(n_splits=5)
splits = list(cv.split(np.zeros(len(master)), y, groups))
for fold, (tr, te) in enumerate(splits):
    assert len(set(groups[tr]) & set(groups[te])) == 0
    print(f"Fold {fold}: train {len(tr)} ({y[tr].sum()} pos, "
          f"{len(set(groups[tr]))} clusters) → "
          f"test {len(te)} ({y[te].sum()} pos, {len(set(groups[te]))} clusters)")

# Threshold note: use the prevalence (~0.147) instead of 0.5 to dodge the
# all-negative artifact we hit in W5. Optimal threshold per fold is also
# fine, but a fixed prevalence-based threshold is cleaner for comparability.
THRESHOLD = float(y.mean())
print(f"\nUsing classification threshold = {THRESHOLD:.3f} (prevalence)")
```

Now run all 8 conditions:

```python
results_rf = []
for cid in range(1, 9):
    cfg = CONDITIONS[cid]
    X = make_X(cid)
    print(f"\n=== Condition {cid}: {cfg['label']} (X shape {X.shape}) ===")
    for fold, (tr, te) in enumerate(splits):
        clf = RandomForestClassifier(
            n_estimators=500, class_weight='balanced',
            min_samples_leaf=3, n_jobs=-1, random_state=42)
        clf.fit(X[tr], y[tr])
        proba = clf.predict_proba(X[te])[:, 1]
        pred = (proba >= THRESHOLD).astype(int)
        results_rf.append({
            'condition': cid,
            'label': cfg['label'],
            'fold': fold,
            'n_test': len(te),
            'n_pos_test': int(y[te].sum()),
            'auprc':    average_precision_score(y[te], proba),
            'auroc':    roc_auc_score(y[te], proba),
            'macro_f1': f1_score(y[te], pred, average='macro'),
            'bal_acc':  balanced_accuracy_score(y[te], pred),
            'mcc':      matthews_corrcoef(y[te], pred),
        })
        r = results_rf[-1]
        print(f"  Fold {fold}: AUPRC={r['auprc']:.3f}, AUROC={r['auroc']:.3f}, "
              f"F1={r['macro_f1']:.3f}, MCC={r['mcc']:.3f}")

results_rf = pd.DataFrame(results_rf)
results_rf.to_csv("../results/factorial_rf_per_fold.csv", index=False)
print("\nDone. Saved per-fold scores.")
```

- [ ] **Time check:** roughly 5 folds × 8 conditions × 30–90 sec each = 20–60 minutes total wall time. Get a coffee.
- [ ] Each fold should have ≥ 30 positives in the test set; verify no fold reports `0 pos_test`.
- [ ] Expect AUPRC for condition 1 around the Week-5 number (~0.19); conditions 2–8 will tell you whether GO context lifts that.

### Aggregate view

```python
agg = (results_rf.groupby(['condition','label'])
                 [['auprc','auroc','macro_f1','bal_acc','mcc']]
                 .agg(['mean','std'])
                 .round(3))
print(agg.to_string())

# Quick-look bar chart
fig, ax = plt.subplots(figsize=(11, 5))
means = results_rf.groupby('condition')['auprc'].mean()
stds  = results_rf.groupby('condition')['auprc'].std()
labels = [CONDITIONS[c]['label'] for c in means.index]
ax.bar(range(1, 9), means.values, yerr=stds.values, capsize=4,
       color=['#888'] + ['#3a7']*7)
ax.set_xticks(range(1, 9)); ax.set_xticklabels(labels, rotation=20, ha='right')
ax.axhline(y[0:].mean(), ls='--', color='gray',
           label=f"chance ({y.mean():.3f})")
ax.set_ylabel("AUPRC (mean ± std over 5 folds)")
ax.set_title("Factorial AUPRC — Random Forest, cluster-aware GroupKFold(5)")
ax.legend()
plt.tight_layout()
plt.savefig("../results/figures/factorial_rf_auprc.png", dpi=150)
plt.show()
```

---

## W6.3 — Statistical analysis (~1.5 hours)

Three things to compute: pairwise comparisons against condition 1 with BH correction, effect sizes with bootstrap CIs, and main + interaction effects from the factorial structure.

```python
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests

# Wide view: rows = folds, columns = conditions
wide = results_rf.pivot(index='fold', columns='condition', values='auprc')

# ---- Pairwise: condition c vs condition 1 (baseline) ----
control = wide[1].values
rows = []
for cid in range(2, 9):
    other = wide[cid].values
    diff = other - control
    w, p = wilcoxon(other, control, alternative='greater', zero_method='wilcox')
    # Bootstrap CI of mean paired difference
    rng = np.random.default_rng(42)
    boots = [diff[rng.integers(0, len(diff), len(diff))].mean()
             for _ in range(2000)]
    lo, hi = np.percentile(boots, [2.5, 97.5])
    rows.append({
        'condition': cid,
        'label': CONDITIONS[cid]['label'],
        'mean_diff_vs_seq': diff.mean(),
        'ci95_lo': lo, 'ci95_hi': hi,
        'wilcoxon_W': w, 'p_raw': p,
    })
pair = pd.DataFrame(rows)
pair['p_BH'] = multipletests(pair['p_raw'], method='fdr_bh')[1]
pair = pair.round({'mean_diff_vs_seq': 3, 'ci95_lo': 3, 'ci95_hi': 3,
                   'p_raw': 4, 'p_BH': 4})
print("=== Pairwise vs condition 1 (Seq only), BH-corrected ===")
print(pair.to_string(index=False))

# ---- Factorial main and interaction effects (on AUPRC) ----
# Each condition is identified by a (BP, MF, CC) triple of 0/1.
flags = pd.DataFrame([
    {'condition': c, 'BP': int(CONDITIONS[c]['bp']),
                     'MF': int(CONDITIONS[c]['mf']),
                     'CC': int(CONDITIONS[c]['cc'])}
    for c in range(1, 9)])
wlong = wide.stack().rename('auprc').reset_index()
wlong = wlong.merge(flags, on='condition')

def main_effect(col):
    on  = wlong[wlong[col] == 1]['auprc'].mean()
    off = wlong[wlong[col] == 0]['auprc'].mean()
    return on - off

def interaction(col_a, col_b):
    # Difference of differences
    g11 = wlong[(wlong[col_a]==1) & (wlong[col_b]==1)]['auprc'].mean()
    g10 = wlong[(wlong[col_a]==1) & (wlong[col_b]==0)]['auprc'].mean()
    g01 = wlong[(wlong[col_a]==0) & (wlong[col_b]==1)]['auprc'].mean()
    g00 = wlong[(wlong[col_a]==0) & (wlong[col_b]==0)]['auprc'].mean()
    return (g11 - g10) - (g01 - g00)

print("\n=== Factorial effects on AUPRC ===")
print(f"  Main effect BP : {main_effect('BP'):+.4f}")
print(f"  Main effect MF : {main_effect('MF'):+.4f}")
print(f"  Main effect CC : {main_effect('CC'):+.4f}")
print(f"  Interaction BP×MF : {interaction('BP','MF'):+.4f}")
print(f"  Interaction BP×CC : {interaction('BP','CC'):+.4f}")
print(f"  Interaction MF×CC : {interaction('MF','CC'):+.4f}")

pair.to_csv("../results/factorial_rf_stats.csv", index=False)
```

- [ ] **What this tells you about the hypotheses** (record in lab notebook):
  - **H1 (main effect):** is any pairwise BH-corrected p < 0.05? If yes, H1 supported.
  - **H2 (MF most informative):** which main effect is largest — BP, MF or CC? If MF, H2 supported. If BP or CC, H2 falsified (note honestly).
  - **H3 (sub-additive interaction):** are the BP×MF, BP×CC, MF×CC interactions negative? If yes, supports H3 (redundancy). If positive, falsifies H3 (synergy).
  - **H4 (stratification):** not testable here — that's Week 7 / Week 9, after we have per-protein predictions.

---

## W6.4 — Gradient-boosting robustness check (~1.5 hours)

Run the same 8 conditions with XGBoost using the same CV splits. This is a robustness check, not the primary analysis: the proposal specifies RF, and we keep RF as the headline classifier unless XGBoost is dramatically better and we have a principled reason to switch.

```python
import xgboost as xgb

results_xgb = []
for cid in range(1, 9):
    cfg = CONDITIONS[cid]
    X = make_X(cid)
    print(f"\n=== XGB Condition {cid}: {cfg['label']} ===")
    for fold, (tr, te) in enumerate(splits):
        # scale_pos_weight handles class imbalance equivalently to class_weight
        spw = (y[tr] == 0).sum() / max(y[tr].sum(), 1)
        clf = xgb.XGBClassifier(
            n_estimators=500, max_depth=6, learning_rate=0.05,
            subsample=0.8, colsample_bytree=0.8,
            scale_pos_weight=spw,
            n_jobs=-1, random_state=42,
            eval_metric='aucpr', tree_method='hist')
        clf.fit(X[tr], y[tr])
        proba = clf.predict_proba(X[te])[:, 1]
        pred = (proba >= THRESHOLD).astype(int)
        results_xgb.append({
            'condition': cid, 'label': cfg['label'], 'fold': fold,
            'auprc': average_precision_score(y[te], proba),
            'auroc': roc_auc_score(y[te], proba),
            'macro_f1': f1_score(y[te], pred, average='macro'),
            'bal_acc': balanced_accuracy_score(y[te], pred),
            'mcc': matthews_corrcoef(y[te], pred),
        })
        r = results_xgb[-1]
        print(f"  Fold {fold}: AUPRC={r['auprc']:.3f}, MCC={r['mcc']:.3f}")

results_xgb = pd.DataFrame(results_xgb)
results_xgb.to_csv("../results/factorial_xgb_per_fold.csv", index=False)

# Compare RF vs XGB AUPRC per condition
print("\n=== RF vs XGB mean AUPRC by condition ===")
cmp = pd.DataFrame({
    'condition': range(1, 9),
    'label':     [CONDITIONS[c]['label'] for c in range(1, 9)],
    'RF_auprc':  results_rf.groupby('condition')['auprc'].mean().round(3).values,
    'XGB_auprc': results_xgb.groupby('condition')['auprc'].mean().round(3).values,
})
print(cmp.to_string(index=False))
```

- [ ] **Decision rule:** if XGB's mean AUPRC across conditions is more than +0.03 above RF, consider adding XGB as a secondary classifier in the report (with the standard "RF was pre-registered, XGB as robustness" framing). Otherwise stick with RF.

---

## W6.5 — Interim review and Week-7 routing (~30 min)

Based on the Week 6 results, decide which of three Week-7 paths to take:

- [ ] **Path I — happy path: GO lifts the baseline meaningfully (≥ +0.03 AUPRC over Seq only, at least one BH-p < 0.05).** Week 7 focuses on the **negative control** (shuffle GO labels, rerun) and the **stratified analysis** (H4 — does GO help where sequence struggles?). The report writes itself.

- [ ] **Path II — null finding: no condition is significantly above Seq only.** Equally interpretable result; Week 7 stress-tests this by trying **per-residue ESM-2 embeddings averaged over annotated disordered regions only** (rather than whole-protein mean-pooled). If the per-residue version *does* show GO benefit, that's a methodological finding. If not, the null is the headline and the report's discussion section makes the case ("modern PLM embeddings already capture protein-level signal").

- [ ] **Path III — mixed: BP lifts but MF and CC don't, or vice versa.** Week 7 does feature-importance analysis to characterise *which* GO terms drive the lift, and stratified analysis on the helped subset.

Record which path you're on in the lab notebook.

---

## W6.6 — Closing tasks (~30 min)

- [ ] Append a Week-6 entry to `lab-notebook.md` with: the 8 condition × mean-AUPRC table; the H1/H2/H3 verdicts; the main and interaction effect sizes; the RF-vs-XGB comparison; the Week-7 path chosen.
- [ ] Commit and push everything. Tag the commit `week6-factorial`.

---

## Notes on what could trip you up

- **Threshold matters for F1/bal_acc/MCC.** We're using the prevalence (0.147) as the threshold. This is honest and consistent but it will make the absolute numbers look different from any literature comparison that uses 0.5. Note in the report.
- **Shared splits, paired tests.** The whole comparison rests on every condition seeing the same train/test split per fold. If you regenerate splits inside the loop by mistake, the pairing breaks and the Wilcoxon test loses power. The `splits = list(cv.split(...))` line outside the loop is what prevents this.
- **Don't run XGBoost with `tree_method='gpu_hist'` on Apple Silicon.** CUDA-only; use `tree_method='hist'` (CPU). The code above already does this.
- **GO matrix coverage.** Not every protein has BP/MF/CC annotations. The `reindex(...).fillna(0)` keeps shapes aligned by treating "no annotation" as the all-zero feature vector. This is the correct semantics — absent ≠ "negative" but it's how the factorial structure naturally encodes "no GO info available."
- **n_jobs=-1 on Apple Silicon.** Sometimes causes weird stalls in scikit-learn with very small problems. If you see a freeze, change to `n_jobs=4` and re-run.

---

## Where this hands off to Week 7

Week 7 then runs the **negative control** (label shuffling), the **stratified analysis** for H4, and feature-importance interpretation. With Week 6's results in hand, Week 7 is mostly diagnostics and storytelling, not new modelling — *unless* W6.5 puts us on Path II, in which case Week 7 also includes the per-residue ESM-2 reframe.
