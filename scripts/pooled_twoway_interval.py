"""Exploratory two-way fixed-prediction gene x miR-family CLIP interval.
Protocol: docs/PREREG_POOLED_TWOWAY_CI_20260928.md (committed pre-outcome).
"""
import csv,json,os
from pathlib import Path
import numpy as np,pandas as pd
from sklearn.metrics import roc_auc_score
from expression_confound import attach_expression,cv_scores,standardize,A,E,C,CSV,EXPR
ROOT=Path(__file__).resolve().parents[1];B=300;SEED=20260928
raw=pd.read_csv(ROOT/CSV);expr=pd.read_csv(ROOT/EXPR,sep='\t');d,_=attach_expression(raw,expr)
y=d.label.to_numpy();genes=d.gene.to_numpy()
pa=cv_scores(standardize(d[A+E].to_numpy(float)),y,genes)
pc=cv_scores(standardize(d[A+E+C].to_numpy(float)),y,genes)
base=float(roc_auc_score(y,pa));full=float(roc_auc_score(y,pc))
prior=json.loads((ROOT/'results/expression_confound.json').read_text())['pooled_auroc_gene_grouped_cv']
assert abs(base-prior['AE'])<1e-8 and abs(full-prior['AEC'])<1e-8
fam={r['MiRBase ID']:r['miR family'] for r in csv.DictReader((ROOT/'data/miR_Family_Info.txt').open(),delimiter='\t') if r['Species ID']=='9606'}
missing=sorted(set(d.mirna)-set(fam));assert not missing,missing
families=d.mirna.map(fam).to_numpy();uf,ifam=np.unique(families,return_inverse=True);ug,igene=np.unique(genes,return_inverse=True)
ck=ROOT/'results/pooled_twoway_interval_checkpoint.json';vals=json.loads(ck.read_text())['bootstrap_deltas'] if ck.exists() else []
rng=np.random.default_rng(SEED)
for b in range(B):
    fc=np.bincount(rng.integers(0,len(uf),len(uf)),minlength=len(uf))
    gc=np.bincount(rng.integers(0,len(ug),len(ug)),minlength=len(ug))
    if b<len(vals):continue
    w=fc[ifam]*gc[igene]
    vals.append(float(roc_auc_score(y,pc,sample_weight=w)-roc_auc_score(y,pa,sample_weight=w)))
    if (b+1)%25==0:
        tmp=ck.with_suffix('.json.tmp');tmp.write_text(json.dumps({'partial':True,'bootstrap_deltas':vals},indent=1)+'\n');os.replace(tmp,ck)
        print('checkpoint',b+1,flush=True)
assert len(vals)==B
out={'type':'exploratory crossed-family-and-gene fixed-prediction weighted bootstrap',
     'protocol':'docs/PREREG_POOLED_TWOWAY_CI_20260928.md','n_rows':len(y),'n_genes':len(ug),'n_families':len(uf),
     'seed':SEED,'bootstrap_replicates':B,'base_auroc':base,'full_auroc':full,'delta':full-base,
     'ci95':[float(x) for x in np.percentile(vals,[2.5,97.5])],
     'n_nonpositive':int(sum(v<=0 for v in vals)),
     'limits':'fits fixed and global scaling from source retained; not refit uncertainty or an independent-library replication'}
(ROOT/'results/pooled_twoway_interval.json').write_text(json.dumps(out,indent=1)+'\n')
print(json.dumps(out),flush=True)
