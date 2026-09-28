"""Refit-aware, leakage-free crossed family x gene interval (queue item 5 remainder).
Protocol: docs/PREREG_REFIT_INTERVAL_20260928.md (committed pre-outcome).
Per draw: resample family/gene cluster counts -> weighted restandardization ->
refit both logistic arms on the SAME fixed GroupKFold(5) gene folds with
sample weights -> weighted OOF AUROC delta. CNN feature columns fixed
(declared in the prereg; TargetScan-trained, labels disjoint from CLIP)."""
import csv, json, os
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold
from expression_confound import attach_expression, A, E, C, CSV, EXPR

ROOT = Path(__file__).resolve().parents[1]
B = 300; SEED = 20260928
raw = pd.read_csv(ROOT / CSV); expr = pd.read_csv(ROOT / EXPR, sep='\t')
d, _ = attach_expression(raw, expr)
y = d.label.to_numpy(); genes = d.gene.to_numpy()
Xa = d[A + E].to_numpy(float); Xc = d[A + E + C].to_numpy(float)

fam = {r['MiRBase ID']: r['miR family'] for r in csv.DictReader((ROOT / 'data/miR_Family_Info.txt').open(), delimiter='\t') if r['Species ID'] == '9606'}
missing = sorted(set(d.mirna) - set(fam)); assert not missing, missing
families = d.mirna.map(fam).to_numpy()
uf, ifam = np.unique(families, return_inverse=True)
ug, igene = np.unique(genes, return_inverse=True)
folds = list(GroupKFold(5).split(Xa, y, genes))  # deterministic, identical across draws

def wstd(X, w):
    W = w.sum()
    mu = (X * w[:, None]).sum(0) / W
    var = (w[:, None] * (X - mu) ** 2).sum(0) / W
    sd = np.sqrt(var); sd[sd == 0] = 1
    return (X - mu) / sd

def arm_scores(X, w):
    p = np.zeros(len(y))
    for tr, te in folds:
        m = LogisticRegression(max_iter=2000).fit(X[tr], y[tr], sample_weight=w[tr])
        p[te] = m.predict_proba(X[te])[:, 1]
    return p

w1 = np.ones(len(y))
pa0 = arm_scores(wstd(Xa, w1), w1); pc0 = arm_scores(wstd(Xc, w1), w1)
base = float(roc_auc_score(y, pa0)); full = float(roc_auc_score(y, pc0))
prior = json.loads((ROOT / 'results/expression_confound.json').read_text())['pooled_auroc_gene_grouped_cv']
assert abs(base - prior['AE']) < 1e-8 and abs(full - prior['AEC']) < 1e-8, (base, full, prior)
print('integrity gate passed: committed AUROCs reproduced at uniform weights', flush=True)

ck = ROOT / 'results/refit_interval_checkpoint.json'
recs = json.loads(ck.read_text())['draws'] if ck.exists() else []
rng = np.random.default_rng(SEED)
for b in range(B):
    fc = np.bincount(rng.integers(0, len(uf), len(uf)), minlength=len(uf))
    gc = np.bincount(rng.integers(0, len(ug), len(ug)), minlength=len(ug))
    if b < len(recs):
        continue
    w = (fc[ifam] * gc[igene]).astype(float)
    degenerate = any(w[tr][y[tr] == cls].sum() <= 0 for tr, _ in folds for cls in (0, 1))
    if degenerate:
        recs.append({'b': b, 'skipped': True})
    else:
        pa = arm_scores(wstd(Xa, w), w); pc = arm_scores(wstd(Xc, w), w)
        recs.append({'b': b, 'skipped': False,
                     'delta': float(roc_auc_score(y, pc, sample_weight=w) - roc_auc_score(y, pa, sample_weight=w))})
    if (b + 1) % 25 == 0:
        tmp = ck.with_suffix('.json.tmp'); tmp.write_text(json.dumps({'partial': True, 'draws': recs}, indent=1) + '\n'); os.replace(tmp, ck)
        print('checkpoint', b + 1, flush=True)
assert len(recs) == B
vals = [r['delta'] for r in recs if not r['skipped']]
out = {'type': 'refit-aware crossed family-and-gene weighted bootstrap (evaluator + standardization refit per draw)',
       'protocol': 'docs/PREREG_REFIT_INTERVAL_20260928.md', 'n_rows': len(y), 'n_genes': len(ug),
       'n_families': len(uf), 'seed': SEED, 'bootstrap_replicates': B, 'n_skipped_degenerate': int(B - len(vals)),
       'base_auroc': base, 'full_auroc': full, 'delta': full - base,
       'ci95': [float(x) for x in np.percentile(vals, [2.5, 97.5])],
       'n_nonpositive': int(sum(v <= 0 for v in vals)),
       'fixed_prediction_comparison': {'ci95': [-0.0015680281768226001, 0.0027334023043778972], 'commit': '719db2d'},
       'limits': 'CNN feature extractor not retrained (declared cost bound); features TargetScan-trained, labels disjoint from CLIP; not an independent-library replication'}
(ROOT / 'results/refit_interval.json').write_text(json.dumps(out, indent=1) + '\n')
print(json.dumps(out), flush=True)
