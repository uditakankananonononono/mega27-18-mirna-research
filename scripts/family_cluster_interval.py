"""Exploratory family-cluster bootstrap of the within-miRNA median CLIP gain.

Not a pooled AUROC CI or gene-cluster CI. TargetScan human miR-family mappings
join all 130 per-miRNA rows; resample entire families, preserving within-family
miRNA deltas and family size, then recompute the median over sampled miRNAs.
"""
import csv, json
from collections import defaultdict
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
rows=json.loads((ROOT/'results/expression_confound.json').read_text())['rows']
fam={r['MiRBase ID']:r['miR family'] for r in csv.DictReader((ROOT/'data/miR_Family_Info.txt').open(),delimiter='\t') if r['Species ID']=='9606'}
assert len(rows)==130 and len({r['mirna'] for r in rows})==130
assert all(r['mirna'] in fam for r in rows)
by=defaultdict(list)
for r in rows:by[fam[r['mirna']]].append(float(r['delta']))
keys=sorted(by)
assert len(keys)==103
rng=np.random.default_rng(20260928)
B=10000
boot=[]
for _ in range(B):
    sampled=rng.choice(len(keys),size=len(keys),replace=True)
    boot.append(float(np.median([v for j in sampled for v in by[keys[j]]])))
vals=np.array([float(r['delta']) for r in rows]); boot=np.array(boot)
out={'type':'exploratory family-cluster percentile bootstrap',
     'estimand':'median within-miRNA conditional AUROC increment; NOT pooled delta',
     'source':'results/expression_confound.json',
     'mapping':'TargetScan human miR family in data/miR_Family_Info.txt',
     'seed':20260928,'bootstrap_replicates':B,'n_mirna':len(rows),'n_families':len(keys),
     'median_delta':float(np.median(vals)),
     'median_ci95':[float(x) for x in np.percentile(boot,[2.5,97.5])],
     'family_sizes':{k:len(by[k]) for k in keys},
     'limitations':'family resampling preserves within-family miRNA ties, not shared genes across families; not a pooled AUROC interval and not independent-library replication'}
(ROOT/'results/family_cluster_interval.json').write_text(json.dumps(out,indent=1)+'\n')
print(json.dumps({k:v for k,v in out.items() if k!='family_sizes'},indent=1))
