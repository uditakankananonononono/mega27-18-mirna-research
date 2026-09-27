"""Exploratory fixed-prediction gene-cluster CI; locked design in docs/PREREG_POOLED_GENE_CI_20260928.md."""
import json, os
from pathlib import Path
import numpy as np,pandas as pd
from sklearn.metrics import roc_auc_score
from expression_confound import attach_expression,cv_scores,standardize,A,E,C,CSV,EXPR
ROOT=Path(__file__).resolve().parents[1]
B=300; SEED=20260928
raw=pd.read_csv(ROOT/CSV);expr=pd.read_csv(ROOT/EXPR,sep='\t');d,_=attach_expression(raw,expr)
y=d.label.to_numpy(); groups=d.gene.to_numpy()
pa=cv_scores(standardize(d[A+E].to_numpy(float)),y,groups)
pc=cv_scores(standardize(d[A+E+C].to_numpy(float)),y,groups)
base=float(roc_auc_score(y,pa)); full=float(roc_auc_score(y,pc))
prior=json.loads((ROOT/'results/expression_confound.json').read_text())['pooled_auroc_gene_grouped_cv']
if abs(base-prior['AE'])>1e-8 or abs(full-prior['AEC'])>1e-8:raise ValueError('recomputed AUROC differs from source')
unique,inv=np.unique(groups,return_inverse=True)
ck=ROOT/'results/pooled_gene_interval_checkpoint.json'
vals=json.loads(ck.read_text())['bootstrap_deltas'] if ck.exists() else []
rng=np.random.default_rng(SEED)
for b in range(B):
    c=np.bincount(rng.integers(0,len(unique),len(unique)),minlength=len(unique))
    if b<len(vals):continue
    w=c[inv]
    vals.append(float(roc_auc_score(y,pc,sample_weight=w)-roc_auc_score(y,pa,sample_weight=w)))
    if (b+1)%25==0:
        temp=ck.with_suffix('.json.tmp');temp.write_text(json.dumps({'partial':True,'bootstrap_deltas':vals},indent=1)+'\n');os.replace(temp,ck)
        print('checkpoint',b+1,flush=True)
if len(vals)!=B:raise ValueError('incomplete bootstrap')
out={'type':'exploratory fixed-prediction gene-cluster pooled AUROC delta',
     'protocol':'docs/PREREG_POOLED_GENE_CI_20260928.md','n_rows':len(y),'n_genes':len(unique),
     'seed':SEED,'bootstrap_replicates':B,'base_auroc':base,'full_auroc':full,
     'delta':full-base,'gene_ci95':[float(x) for x in np.percentile(vals,[2.5,97.5])],
     'n_nonpositive':int(sum(v<=0 for v in vals)),
     'limits':'predictions fixed; gene resampling alone; family dependence and global-scaling transduction remain; not external replication'}
(ROOT/'results/pooled_gene_interval.json').write_text(json.dumps(out,indent=1)+'\n')
print(json.dumps(out),flush=True)
