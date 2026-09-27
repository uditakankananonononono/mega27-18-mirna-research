"""Rebuild missing ignored DuplexCNN weights with the original fixed training path.

Uses the same data.build, gene_split, feature normalization and fit_cnn call as
src.mirna.train.main(max_rows=80000, epochs=6); avoids the separate ridge arms.
The existing archived benchmark JSON is not overwritten.
"""
import hashlib, json
from pathlib import Path
import numpy as np
import torch
from src.mirna.data import load_mirnas, load_utrs, build, gene_split
from src.mirna.train import fit_cnn, predict, metrics
ROOT = Path(__file__).resolve().parents[1]
D = ROOT / 'data'
mirs = load_mirnas(D / 'miR_Family_Info.txt')
utrs = load_utrs(D / 'human_utrs.tsv')
U, M, ST, y, genes, meta, F = build(D / 'human_sites.tsv', utrs, mirs, max_rows=80000)
ST = ST.clip(0,3)
tr, te = gene_split(genes)
mu, sd = F[tr].mean(0), F[tr].std(0) + 1e-6
F = ((F-mu)/sd).astype(np.float32)
print('n/train/test',len(y),len(tr),len(te),flush=True)
net=fit_cnn(U,M,ST,y,tr,epochs=6,seed=0,bs=256,lr=2e-3,log=lambda s:print(s,flush=True),F=F)
p=predict(net,U,M,ST,te,F=F)
result=metrics(y[te],p)
print('rebuild metrics',result,flush=True)
original=json.loads((ROOT/'results/context_pp_benchmark.json').read_text())
assert (len(y),len(tr),len(te))==(original['n_sites'],original['n_train'],original['n_test'])
assert abs(result['pearson']-original['duplex_cnn_plus_context']['pearson'])<1e-5, 'Rebuilt weights disagree with archived metric'
out=ROOT/'results/duplex_cnn.pt'
torch.save(net.state_dict(),out)
hashval=hashlib.sha256(out.read_bytes()).hexdigest()
print('weights SHA256',hashval,flush=True)
