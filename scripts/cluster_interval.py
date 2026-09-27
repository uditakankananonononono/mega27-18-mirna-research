"""Exploratory miRNA-cluster bootstrap of per-miRNA conditional CLIP AUROC gains.

This is a cluster-level uncertainty check on the median effect across miRNAs,
not a pooled pair-level AUROC CI or an independent replication.
"""
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent.parent
src=ROOT/'results/expression_confound.json'
d=json.loads(src.read_text())
rows=d['rows']
assert len(rows)==130 and len({r['mirna'] for r in rows})==len(rows)
deltas=np.asarray([float(r['delta']) for r in rows])
rng=np.random.default_rng(20260927)
B=10000
boots=np.empty(B)
for i in range(B):
    boots[i]=np.median(rng.choice(deltas,size=len(deltas),replace=True))
out={
    'type':'exploratory percentile bootstrap over miRNA clusters (the rows in expression_confound.json)',
    'estimand':'median within-miRNA delta AUROC (AEC minus AE); not pooled AUROC delta',
    'source':'results/expression_confound.json',
    'seed':20260927,'bootstrap_replicates':B,'n_mirna':len(rows),
    'median_delta':float(np.median(deltas)),
    'mean_delta':float(np.mean(deltas)),
    'median_ci95':[float(v) for v in np.percentile(boots,[2.5,97.5])],
    'n_positive':int(np.sum(deltas>0)),
    'limits':'miRNAs may share seed families and genes; this is not a family- or gene-clustered interval and does not identify causal model value'
}
p=(ROOT/'results/expression_cluster_interval.json')
p.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
