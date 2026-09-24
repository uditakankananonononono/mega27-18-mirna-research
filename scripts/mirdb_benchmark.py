"""Head-to-head vs a public leader: miRDB v6.0 (MirTarget) against our DuplexCNN on the
identical miRTarBase 10.0 strong-evidence benchmark rows (results/mirtarbase_rows.csv).

miRDB lists only predictions with score >= 50; an absent (miRNA, gene) pair is scored 0
(below miRDB's own reporting threshold). Gene score = max over the gene's RefSeq transcripts
(RefSeq->symbol via MyGene.info, data/mirdb/refseq2symbol.json). Our score = -min CNN
context++ prediction (more negative = stronger repression). Evaluation restricted to
miRNAs present in miRDB. Caveat: miRDB's MirTarget was trained on external CLIP/array
data and may overlap literature-validated targets; this is a leakage risk in miRDB's favour.
Output: results/mirdb_benchmark.json
"""
import csv, json
from collections import defaultdict
from pathlib import Path
import numpy as np
from sklearn.metrics import roc_auc_score
from scipy.stats import wilcoxon

R = Path(__file__).resolve().parents[1]
m2s = json.load(open(R / 'data/mirdb/refseq2symbol.json'))
mirdb = defaultdict(dict)
for line in open(R / 'data/mirdb/rows_mirs.tsv'):
    mi, nm, sc = line.rstrip('\n').split('\t')
    g = m2s.get(nm)
    if g:
        mirdb[mi][g] = max(mirdb[mi].get(g, 0.0), float(sc))
rows = [r for r in csv.DictReader(open(R / 'results/mirtarbase_rows.csv')) if r['mirna'] in mirdb]
y = np.array([int(r['label']) for r in rows])
cnn = np.array([-float(r['min']) for r in rows])
md = np.array([mirdb[r['mirna']].get(r['gene'], 0.0) for r in rows])
per = {}
by = defaultdict(list)
for i, r in enumerate(rows):
    by[r['mirna']].append(i)
for mi, idx in by.items():
    yy = y[idx]
    if yy.min() == yy.max() or yy.sum() < 5:
        continue
    per[mi] = {'n': len(idx), 'n_pos': int(yy.sum()),
               'cnn_auroc': float(roc_auc_score(yy, cnn[idx])),
               'mirdb_auroc': float(roc_auc_score(yy, md[idx]))}
c = np.array([v['cnn_auroc'] for v in per.values()]); d = np.array([v['mirdb_auroc'] for v in per.values()])
out = {'n_rows': len(rows), 'n_pos': int(y.sum()), 'n_mirnas_in_mirdb': len(by),
       'n_mirnas_evaluated': len(per),
       'pooled_cnn_auroc': float(roc_auc_score(y, cnn)),
       'pooled_mirdb_auroc': float(roc_auc_score(y, md)),
       'frac_rows_mirdb_predicted': float((md > 0).mean()),
       'frac_pos_mirdb_predicted': float((md[y == 1] > 0).mean()),
       'frac_neg_mirdb_predicted': float((md[y == 0] > 0).mean()),
       'median_per_mirna_cnn': float(np.median(c)), 'median_per_mirna_mirdb': float(np.median(d)),
       'mirdb_minus_cnn_median': float(np.median(d - c)),
       'mirdb_wins': int((d > c).sum()),
       'wilcoxon_p': float(wilcoxon(d, c).pvalue), 'per_mirna': per}
json.dump(out, open(R / 'results/mirdb_benchmark.json', 'w'), indent=1)
print(json.dumps({k: v for k, v in out.items() if k != 'per_mirna'}, indent=1))
