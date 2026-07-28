"""
17_bp_resolution_strengthening.py

Four checks to strengthen or reject the Week 10 finding that BP main effect
at full-term GO resolution is roughly 8x larger than at GO Slim resolution.

Checks:
  1. Bootstrap 95% CI on the BP main effect (fold-level paired diffs).
  2. Scrambled-BP negative control: shuffle protein-to-BP mapping.
     If effect persists, the +0.016 is spurious.
  3. Top-20 BP full-term features by permutation importance.
     If biologically coherent, effect has mechanistic support.
  4. Stability across three RandomForest seeds and CV shuffle seeds.

Reuses the setup from 16_full_go_factorial.py.

Runtime: ~15-20 minutes.
"""

import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

import numpy as np
import pandas as pd
from pathlib import Path
from collections import defaultdict
from sklearn.ensemble import RandomForestClassifier
from sklearn.inspection import permutation_importance
from sklearn.model_selection import GroupKFold
from sklearn.metrics import average_precision_score
from scipy.stats import wilcoxon
import warnings
warnings.filterwarnings("ignore")

RNG_MAIN = 42
BOOTSTRAP_ITERS = 2000
DATA = Path("../data")

# -----------------------------------------------------------------------------
# Setup (same as 16_full_go_factorial.py)
# -----------------------------------------------------------------------------
print("Loading data...")
master = pd.read_csv(DATA / "master_clean.csv")
df = master.copy()
esm_data = np.load(DATA / "features_esm2_disorder.npz", allow_pickle=True)
acc_key = next(k for k in esm_data.keys() if esm_data[k].ndim == 1)
feat_key = next(k for k in esm_data.keys() if esm_data[k].ndim == 2)
esm_accs = [a.decode() if isinstance(a, bytes) else str(a) for a in esm_data[acc_key]]
esm_lookup = {a: esm_data[feat_key][i] for i, a in enumerate(esm_accs)}
df = df[df["acc"].isin(esm_lookup)].reset_index(drop=True)
X_esm = np.stack([esm_lookup[a] for a in df["acc"]])
y = df["d2o"].values.astype(int)
groups = df["cluster"].fillna(-1).astype(int).values

# GO annotations + ancestor closure (BP only for these checks)
go_ann = pd.read_csv(DATA / "go_annotations_experimental.csv")
go_ann = go_ann[go_ann["DB_Object_ID"].isin(df["acc"])].copy()
go_ann["aspect"] = go_ann["Aspect"].map({"P": "BP", "F": "MF", "C": "CC"})

def parse_obo(path):
    terms = {}; current = None
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line == "[Term]": current = {"is_a": [], "part_of": []}
            elif line == "" and current is not None:
                if "id" in current: terms[current["id"]] = current
                current = None
            elif current is not None:
                if line.startswith("id: "): current["id"] = line[4:]
                elif line.startswith("is_a: "):
                    current["is_a"].append(line[6:].split(" ! ")[0].strip())
                elif line.startswith("relationship: part_of "):
                    current["part_of"].append(line[len("relationship: part_of "):].split(" ! ")[0].strip())
                elif line.startswith("is_obsolete: true"): current["obsolete"] = True
    return {k: v for k, v in terms.items() if not v.get("obsolete")}

terms = parse_obo(DATA / "go.obo")

def ancestors(gid, memo={}):
    if gid in memo: return memo[gid]
    if gid not in terms: return {gid}
    seen = {gid}; stack = list(terms[gid]["is_a"]) + list(terms[gid]["part_of"])
    while stack:
        p = stack.pop()
        if p in seen or p not in terms: continue
        seen.add(p); stack.extend(terms[p]["is_a"]); stack.extend(terms[p]["part_of"])
    memo[gid] = seen
    return seen

def build_bp_matrix(order, ann_df=None):
    if ann_df is None: ann_df = go_ann
    sub = ann_df[ann_df["aspect"] == "BP"]
    protein_terms = defaultdict(set)
    for acc, gid in zip(sub["DB_Object_ID"], sub["GO_ID"]):
        protein_terms[acc].update(ancestors(gid))
    counts = defaultdict(int)
    for ts in protein_terms.values():
        for t in ts: counts[t] += 1
    vocab = sorted([t for t, c in counts.items() if 5 <= c <= len(order) - 5])
    idx = {t: i for i, t in enumerate(vocab)}
    mat = np.zeros((len(order), len(vocab)), dtype=np.float32)
    for i, acc in enumerate(order):
        for t in protein_terms.get(acc, set()):
            if t in idx: mat[i, idx[t]] = 1.0
    return mat, vocab

order = df["acc"].tolist()
print("Building BP full-term matrix...")
X_bp, bp_vocab = build_bp_matrix(order)
print(f"  BP full-term: {X_bp.shape[1]} terms, density={X_bp.mean():.3f}")

# -----------------------------------------------------------------------------
# Helper: run a single-condition factorial (baseline vs +BP) and return per-fold
# -----------------------------------------------------------------------------
def factorial_bp(X_bp_arg, rf_seed=42, cv_seed=None):
    gkf = GroupKFold(n_splits=5)
    if cv_seed is not None:
        rng = np.random.default_rng(cv_seed)
        perm = rng.permutation(len(y))
        y_p, groups_p, X_esm_p, X_bp_p = y[perm], groups[perm], X_esm[perm], X_bp_arg[perm]
    else:
        y_p, groups_p, X_esm_p, X_bp_p = y, groups, X_esm, X_bp_arg
    baseline, plus_bp = [], []
    for tr, te in gkf.split(X_esm_p, y_p, groups_p):
        for X_full, container in [(X_esm_p, baseline), (np.hstack([X_esm_p, X_bp_p]), plus_bp)]:
            clf = RandomForestClassifier(
                n_estimators=500, class_weight="balanced",
                min_samples_leaf=3, n_jobs=-1, random_state=rf_seed,
            )
            clf.fit(X_full[tr], y_p[tr])
            p = clf.predict_proba(X_full[te])[:, 1]
            container.append(average_precision_score(y_p[te], p))
    return np.array(baseline), np.array(plus_bp)

# -----------------------------------------------------------------------------
# CHECK 1: Bootstrap 95% CI on BP main effect
# -----------------------------------------------------------------------------
print("\n" + "="*72)
print("CHECK 1: Bootstrap 95% CI on BP main effect")
print("="*72)
baseline_folds, plus_bp_folds = factorial_bp(X_bp, rf_seed=42)
diffs = plus_bp_folds - baseline_folds
print(f"  Baseline AUPRC (per fold): {baseline_folds}")
print(f"  +BP AUPRC (per fold):      {plus_bp_folds}")
print(f"  Paired diffs (per fold):   {diffs}")
print(f"  Mean lift: {diffs.mean():+.4f}")

# Bootstrap
rng = np.random.default_rng(RNG_MAIN)
boot_means = np.array([rng.choice(diffs, size=len(diffs), replace=True).mean()
                       for _ in range(BOOTSTRAP_ITERS)])
lo, hi = np.quantile(boot_means, [0.025, 0.975])
print(f"  Bootstrap 95% CI on mean lift: [{lo:+.4f}, {hi:+.4f}]")
print(f"  Fraction of bootstraps > 0: {(boot_means > 0).mean():.3f}")

# -----------------------------------------------------------------------------
# CHECK 2: Scrambled-BP negative control
# -----------------------------------------------------------------------------
print("\n" + "="*72)
print("CHECK 2: Scrambled-BP negative control")
print("="*72)
scramble_rng = np.random.default_rng(RNG_MAIN)
perm = scramble_rng.permutation(len(y))
X_bp_scrambled = X_bp[perm]
_, plus_bp_scrambled = factorial_bp(X_bp_scrambled, rf_seed=42)
scrambled_diffs = plus_bp_scrambled - baseline_folds
print(f"  Real BP lift:      {diffs.mean():+.4f}")
print(f"  Scrambled BP lift: {scrambled_diffs.mean():+.4f}")
print(f"  Difference (real - scrambled): {(diffs.mean() - scrambled_diffs.mean()):+.4f}")
if scrambled_diffs.mean() < diffs.mean() - 0.005:
    print("  RESULT: Scrambled control gives smaller lift. Real BP signal is doing work.")
else:
    print("  RESULT: Scrambled control matches real BP. Effect is likely spurious.")

# -----------------------------------------------------------------------------
# CHECK 3: Top-20 BP full-term features by permutation importance
# -----------------------------------------------------------------------------
print("\n" + "="*72)
print("CHECK 3: Top-20 BP features by permutation importance")
print("="*72)
X_full = np.hstack([X_esm, X_bp])
gkf = GroupKFold(n_splits=5)
# Fit on 4 folds, compute permutation importance on the 5th
splits = list(gkf.split(X_full, y, groups))
train_idx, test_idx = splits[0]
clf = RandomForestClassifier(n_estimators=500, class_weight="balanced",
                              min_samples_leaf=3, n_jobs=-1, random_state=42)
clf.fit(X_full[train_idx], y[train_idx])
print("  Computing permutation importance (this takes a few minutes)...")
result = permutation_importance(
    clf, X_full[test_idx], y[test_idx], n_repeats=5,
    random_state=42, scoring="average_precision", n_jobs=-1,
)
imps = result.importances_mean

esm_dim = X_esm.shape[1]
bp_imps = imps[esm_dim:]
top20_idx = np.argsort(bp_imps)[::-1][:20]

# Try to get term names from the obo for readability
term_names = {t: terms[t].get("name", "") if t in terms else "" for t in bp_vocab}

print(f"\n  Rank  GO ID        Importance   Term name")
print(f"  ----  -----------  -----------  " + "-"*40)
for rank, i in enumerate(top20_idx, 1):
    gid = bp_vocab[i]
    print(f"  {rank:>4}  {gid:<11}  {bp_imps[i]:>+11.5f}  {term_names.get(gid, '')[:40]}")

# -----------------------------------------------------------------------------
# CHECK 4: Stability across 3 seeds
# -----------------------------------------------------------------------------
print("\n" + "="*72)
print("CHECK 4: Stability across 3 RF seeds x 3 CV shuffle seeds")
print("="*72)
print(f"  RF seed | CV seed | Baseline | +BP    | Lift    ")
print(f"  --------+---------+----------+--------+--------")
lifts = []
for rf_seed in [0, 42, 123]:
    for cv_seed in [None, 7, 999]:
        b, p = factorial_bp(X_bp, rf_seed=rf_seed, cv_seed=cv_seed)
        lift = (p - b).mean()
        lifts.append(lift)
        cv_disp = "-" if cv_seed is None else str(cv_seed)
        print(f"  {rf_seed:>7} | {cv_disp:>7} | {b.mean():.4f}   | {p.mean():.4f} | {lift:+.4f}")
lifts = np.array(lifts)
print(f"\n  Lifts across 9 (RF, CV) combinations:")
print(f"    Mean: {lifts.mean():+.4f}")
print(f"    Range: [{lifts.min():+.4f}, {lifts.max():+.4f}]")
print(f"    SD: {lifts.std():.4f}")
print(f"    Fraction positive: {(lifts > 0).mean():.2%}")

# -----------------------------------------------------------------------------
# HEADLINE
# -----------------------------------------------------------------------------
print("\n" + "="*72)
print("STRENGTHENING SUMMARY")
print("="*72)
print(f"BP main effect at full-term GO resolution:")
print(f"  Point estimate lift:  {diffs.mean():+.4f}")
print(f"  Bootstrap 95% CI:     [{lo:+.4f}, {hi:+.4f}]")
print(f"  Scrambled control:    {scrambled_diffs.mean():+.4f}")
print(f"  Stability (9 seed combos): mean={lifts.mean():+.4f}, range=[{lifts.min():+.4f}, {lifts.max():+.4f}], {(lifts > 0).mean():.0%} positive")
print("\nDecision rules for how to write this up:")
print("  * If CI lower bound > 0 AND scrambled < real - 0.005 AND >=80% seeds positive:")
print("      -> STRONG directional finding. Write it up as a key result.")
print("  * If CI lower bound near 0 OR scrambled ~ real OR seeds mixed:")
print("      -> WEAK finding. Report honestly but as suggestive, not a headline.")
print("  * If any check fails badly:")
print("      -> Drop the claim. The +0.016 was noise.")
print("="*72)
