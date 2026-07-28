"""
16_full_go_factorial.py

Week 10 fast experiment: higher-resolution GO representation.

We showed in W5-W9 that GO Slim (65-100 dim per sub-ontology) does not help
on top of ESM-2. This tests whether full-term GO WITH ancestor propagation
(hundreds of dim per sub-ontology, information-preserving up the DAG) changes
the picture.

If the null holds here too: it holds across THREE levels of GO resolution
(Slim, full-term binary, and higher-resolution). Stronger null.

If the null breaks here: we have a positive result. Story becomes "GO
resolution mattered; Slim was too coarse to help but hierarchical full-term
GO does add value."

Both outcomes serve the rewritten narrative.

Run this from the notebooks/ directory:
    python 16_full_go_factorial.py

Expected runtime on M-series Mac: ~4-6 minutes.
"""

import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

import numpy as np
import pandas as pd
from pathlib import Path
from collections import defaultdict
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GroupKFold
from sklearn.metrics import average_precision_score
from scipy.stats import wilcoxon
import warnings
warnings.filterwarnings("ignore")

# -----------------------------------------------------------------------------
# 1. Load master data and existing pipeline artifacts
# -----------------------------------------------------------------------------
DATA = Path("../data")
print("Loading master tables...")
master = pd.read_csv(DATA / "master_clean.csv")
print("  master_clean columns:", list(master.columns))

# master_clean.csv already has acc, d2o, and cluster.
# Only merge in what's missing.
df = master.copy()
if "d2o" not in df.columns:
    df = df.merge(pd.read_csv(DATA / "labels.csv"), on="acc")
if "cluster" not in df.columns:
    df = df.merge(pd.read_csv(DATA / "clusters.csv"), on="acc", how="left")
print(f"  Proteins: {len(df):,}")
print(f"  Positives (d2o): {int(df['d2o'].sum()):,}")

# Load disorder-pool ESM-2 features (this is the strongest baseline scale)
print("Loading disorder-pool ESM-2 embeddings...")
esm_data = np.load(DATA / "features_esm2_disorder.npz", allow_pickle=True)
print("  npz keys:", list(esm_data.keys()))

# Auto-detect which key is the accession list and which is the feature matrix.
# The feature matrix is 2-D; the accession vector is 1-D.
acc_key = None
feat_key = None
for k in esm_data.keys():
    arr = esm_data[k]
    if arr.ndim == 2:
        feat_key = k
    elif arr.ndim == 1:
        acc_key = k
if acc_key is None or feat_key is None:
    # Fallback: try common names
    for cand in ["accs", "ids", "acc", "accessions"]:
        if cand in esm_data.keys():
            acc_key = cand; break
    for cand in ["features", "X", "embeddings", "feats"]:
        if cand in esm_data.keys():
            feat_key = cand; break

esm_accs = esm_data[acc_key]
esm_feats = esm_data[feat_key]
# Handle bytes vs str accessions
esm_accs = [a.decode() if isinstance(a, bytes) else str(a) for a in esm_accs]
esm_lookup = {a: esm_feats[i] for i, a in enumerate(esm_accs)}
print(f"  Using acc_key='{acc_key}', feat_key='{feat_key}'")
print(f"  ESM-2 shape: {esm_feats.shape}")

# Align to master
df = df[df["acc"].isin(esm_lookup)].reset_index(drop=True)
X_esm = np.stack([esm_lookup[a] for a in df["acc"]])
y = df["d2o"].values.astype(int)
groups = df["cluster"].fillna(-1).astype(int).values
print(f"  Aligned: {len(df):,} proteins, {y.sum()} positives\n")

# -----------------------------------------------------------------------------
# 2. Load raw GO annotations (experimental evidence only) + parent map from obo
# -----------------------------------------------------------------------------
print("Loading raw GO annotations (experimental evidence)...")
go_ann = pd.read_csv(DATA / "go_annotations_experimental.csv")
go_ann = go_ann[go_ann["DB_Object_ID"].isin(df["acc"])].copy()
go_ann["aspect"] = go_ann["Aspect"].map({"P": "BP", "F": "MF", "C": "CC"})
print(f"  Annotations after align: {len(go_ann):,}")
print(f"  By aspect: {dict(go_ann['aspect'].value_counts())}\n")

print("Parsing go.obo for ancestor relationships...")

def parse_obo(path):
    """Minimal .obo parser. Returns {go_id: {'name':..., 'namespace':..., 'is_a':[...], 'part_of':[...]}}."""
    terms = {}
    current = None
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line == "[Term]":
                current = {"is_a": [], "part_of": []}
            elif line == "" and current is not None:
                if "id" in current:
                    terms[current["id"]] = current
                current = None
            elif current is not None:
                if line.startswith("id: "):
                    current["id"] = line[4:]
                elif line.startswith("name: "):
                    current["name"] = line[6:]
                elif line.startswith("namespace: "):
                    current["namespace"] = line[11:]
                elif line.startswith("is_a: "):
                    parent = line[6:].split(" ! ")[0].strip()
                    current["is_a"].append(parent)
                elif line.startswith("relationship: part_of "):
                    parent = line[len("relationship: part_of "):].split(" ! ")[0].strip()
                    current["part_of"].append(parent)
                elif line.startswith("is_obsolete: true") and current is not None:
                    current["obsolete"] = True
    return terms

terms = parse_obo(DATA / "go.obo")
terms = {k: v for k, v in terms.items() if not v.get("obsolete")}
print(f"  Parsed {len(terms):,} GO terms")

def ancestors(go_id, memo={}):
    """All ancestors of go_id including itself, walking is_a + part_of."""
    if go_id in memo:
        return memo[go_id]
    if go_id not in terms:
        return {go_id}
    seen = {go_id}
    stack = list(terms[go_id]["is_a"]) + list(terms[go_id]["part_of"])
    while stack:
        p = stack.pop()
        if p in seen or p not in terms:
            continue
        seen.add(p)
        stack.extend(terms[p]["is_a"])
        stack.extend(terms[p]["part_of"])
    memo[go_id] = seen
    return seen

print("Computing ancestor closures (this is the propagation step)...")
# Test on one term
test = "GO:0006355"  # regulation of transcription
print(f"  Test: ancestors of GO:0006355 -> {len(ancestors(test))} terms")

# -----------------------------------------------------------------------------
# 3. Build propagated full-term GO feature matrices per sub-ontology
# -----------------------------------------------------------------------------
print("\nBuilding propagated GO feature matrices...")

def build_go_matrix(annotations_df, aspect_label, protein_order):
    """Return (matrix, term_list) for the given aspect with ancestor propagation."""
    sub = annotations_df[annotations_df["aspect"] == aspect_label]
    # For each protein, collect all annotated terms + all their ancestors
    protein_terms = defaultdict(set)
    for acc, gid in zip(sub["DB_Object_ID"], sub["GO_ID"]):
        protein_terms[acc].update(ancestors(gid))
    # Feature vocabulary: all terms that appear (after propagation) in >= 5 proteins
    counts = defaultdict(int)
    for terms_set in protein_terms.values():
        for t in terms_set:
            counts[t] += 1
    vocab = sorted([t for t, c in counts.items() if 5 <= c <= len(protein_order) - 5])
    term_to_col = {t: i for i, t in enumerate(vocab)}
    # Build matrix
    mat = np.zeros((len(protein_order), len(vocab)), dtype=np.float32)
    for i, acc in enumerate(protein_order):
        for t in protein_terms.get(acc, set()):
            if t in term_to_col:
                mat[i, term_to_col[t]] = 1.0
    return mat, vocab

order = df["acc"].tolist()
X_bp, vocab_bp = build_go_matrix(go_ann, "BP", order)
X_mf, vocab_mf = build_go_matrix(go_ann, "MF", order)
X_cc, vocab_cc = build_go_matrix(go_ann, "CC", order)
print(f"  BP full: {X_bp.shape[1]} terms (vs {65} in Slim)")
print(f"  MF full: {X_mf.shape[1]} terms")
print(f"  CC full: {X_cc.shape[1]} terms")
print(f"  Non-zero density BP: {X_bp.mean():.3f}")
print()

# -----------------------------------------------------------------------------
# 4. Factorial: baseline vs +BP vs +MF vs +CC vs full
# -----------------------------------------------------------------------------
print("Running 8-condition factorial with GroupKFold(5)...")

CONDITIONS = {
    "C0_baseline":      (0, 0, 0),
    "C1_BP":            (1, 0, 0),
    "C2_MF":            (0, 1, 0),
    "C3_CC":            (0, 0, 1),
    "C4_BP_MF":         (1, 1, 0),
    "C5_BP_CC":         (1, 0, 1),
    "C6_MF_CC":         (0, 1, 1),
    "C7_full":          (1, 1, 1),
}

def build_features(bp, mf, cc):
    parts = [X_esm]
    if bp: parts.append(X_bp)
    if mf: parts.append(X_mf)
    if cc: parts.append(X_cc)
    return np.hstack(parts)

gkf = GroupKFold(n_splits=5)
results = {cond: [] for cond in CONDITIONS}

for fold, (train_idx, test_idx) in enumerate(gkf.split(X_esm, y, groups)):
    y_tr, y_te = y[train_idx], y[test_idx]
    for cond, (bp, mf, cc) in CONDITIONS.items():
        X_full = build_features(bp, mf, cc)
        X_tr, X_te = X_full[train_idx], X_full[test_idx]
        clf = RandomForestClassifier(
            n_estimators=500, class_weight="balanced",
            min_samples_leaf=3, n_jobs=-1, random_state=42
        )
        clf.fit(X_tr, y_tr)
        proba = clf.predict_proba(X_te)[:, 1]
        auprc = average_precision_score(y_te, proba)
        results[cond].append(auprc)
    print(f"  Fold {fold+1}/5 done")

# -----------------------------------------------------------------------------
# 5. Report
# -----------------------------------------------------------------------------
print("\n" + "="*72)
print("RESULTS: Higher-resolution GO factorial (full-term, ancestor-propagated)")
print("="*72)
print(f"{'Condition':<18} {'AUPRC':>8} {'Δ vs C0':>10} {'p (Wilcoxon)':>14}")
print("-"*72)

baseline = np.array(results["C0_baseline"])
for cond in CONDITIONS:
    arr = np.array(results[cond])
    mean = arr.mean()
    delta = mean - baseline.mean()
    if cond == "C0_baseline":
        pstr = "-"
    else:
        try:
            _, p = wilcoxon(arr, baseline, alternative="greater")
            pstr = f"{p:.3f}"
        except ValueError:
            pstr = "n/a"
    print(f"{cond:<18} {mean:>8.4f} {delta:>+10.4f} {pstr:>14}")

# BH correction across 7 non-baseline comparisons
print("\nBenjamini-Hochberg-adjusted p-values:")
pvals = []
conds_ordered = [c for c in CONDITIONS if c != "C0_baseline"]
for cond in conds_ordered:
    try:
        _, p = wilcoxon(np.array(results[cond]), baseline, alternative="greater")
    except ValueError:
        p = 1.0
    pvals.append(p)
# Simple BH
pvals_arr = np.array(pvals)
order_p = np.argsort(pvals_arr)
ranks = np.empty_like(order_p)
ranks[order_p] = np.arange(1, len(pvals_arr) + 1)
adj = pvals_arr * len(pvals_arr) / ranks
# Enforce monotonicity
adj_sorted = adj[order_p]
for i in range(len(adj_sorted) - 2, -1, -1):
    adj_sorted[i] = min(adj_sorted[i], adj_sorted[i+1])
adj_final = np.empty_like(adj_sorted)
adj_final[order_p] = adj_sorted
for cond, p, pa in zip(conds_ordered, pvals, adj_final):
    print(f"  {cond:<18} raw p={p:.3f}  BH-adj={pa:.3f}")

# Save results
out = pd.DataFrame(results)
out.to_csv("../data/results_16_full_go.csv", index=False)
print(f"\nSaved per-fold results to ../data/results_16_full_go.csv")

# -----------------------------------------------------------------------------
# 6. Headline
# -----------------------------------------------------------------------------
print("\n" + "="*72)
best_lift = max(np.array(results[c]).mean() - baseline.mean() for c in conds_ordered)
best_cond = max(conds_ordered, key=lambda c: np.array(results[c]).mean() - baseline.mean())
if adj_final.min() < 0.05:
    print("HEADLINE: Higher-resolution GO IS significant after BH correction.")
    print(f"  Best condition: {best_cond} at +{best_lift:.4f} AUPRC")
    print("  Narrative: GO Slim was too coarse; hierarchical full-term GO does add signal.")
else:
    print("HEADLINE: Higher-resolution GO also fails to help.")
    print(f"  Best lift observed: +{best_lift:.4f} (BH-adj p={adj_final.min():.3f})")
    print(f"  BP dim: {X_bp.shape[1]}, MF dim: {X_mf.shape[1]}, CC dim: {X_cc.shape[1]}")
    print("  Narrative: null generalizes across GO resolution (Slim, full-term binary).")
print("="*72)
