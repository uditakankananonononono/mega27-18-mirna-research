"""Ablation: DuplexCNN with vs without the explicit complementarity pair-map
branch, same data/split/epochs/seed as context_pp_benchmark.json.
Output: results/ablation_pairmap.json"""
import json, sys, time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import numpy as np
import torch
import torch.nn as nn
from mirna.data import load_mirnas, load_utrs, build, gene_split
from mirna.model import DuplexCNN
from mirna.train import metrics, predict

ROOT = Path(__file__).resolve().parents[1]

class NoPairCNN(nn.Module):
    def __init__(self, n_site_types=4, width=32, n_feat=0):
        super().__init__()
        self.utr = nn.Sequential(nn.Conv1d(4, width, 7, padding=3), nn.ReLU(),
                                 nn.Conv1d(width, width, 5, padding=2), nn.ReLU(),
                                 nn.AdaptiveMaxPool1d(1))
        self.mir = nn.Sequential(nn.Conv1d(4, width, 5, padding=2), nn.ReLU(),
                                 nn.AdaptiveMaxPool1d(1))
        self.st = nn.Embedding(n_site_types, 8)
        self.head = nn.Sequential(nn.Linear(2 * width + 8 + n_feat, 64), nn.ReLU(), nn.Linear(64, 1))

    def forward(self, u, m, st, f=None):
        parts = [self.utr(u).flatten(1), self.mir(m).flatten(1), self.st(st)]
        if f is not None:
            parts.append(f)
        return self.head(torch.cat(parts, 1)).squeeze(1)

def fit(net, U, M, ST, y, tr, F, epochs=6, seed=0, bs=256, lr=2e-3):
    torch.manual_seed(seed)
    opt = torch.optim.Adam(net.parameters(), lr=lr)
    Ut, Mt, St, yt, Ft = map(torch.from_numpy, (U, M, ST, y, F))
    for ep in range(epochs):
        perm = np.random.default_rng(seed + ep).permutation(tr)
        net.train(); tot = 0.0
        for i in range(0, len(perm), bs):
            b = torch.from_numpy(perm[i:i + bs])
            loss = ((net(Ut[b], Mt[b], St[b], Ft[b]) - yt[b]) ** 2).mean()
            opt.zero_grad(); loss.backward(); opt.step(); tot += loss.item() * len(b)
        print(f"epoch {ep} mse {tot/len(perm):.4f}", flush=True)
    return net

d = ROOT / "data"
mirs = load_mirnas(d / "miR_Family_Info.txt")
utrs = load_utrs(d / "human_utrs.tsv")
U, M, ST, y, genes, meta, F = build(d / "human_sites.tsv", utrs, mirs, max_rows=80000)
ST = ST.clip(0, 3)
tr, te = gene_split(genes)
mu, sd = F[tr].mean(0), F[tr].std(0) + 1e-6
F = ((F - mu) / sd).astype(np.float32)

res = {"split": "held-out genes (20%), same as context_pp_benchmark.json", "n_test": int(len(te))}
for name, net in [("full_duplexcnn", DuplexCNN(n_feat=F.shape[1])),
                  ("no_pairmap", NoPairCNN(n_feat=F.shape[1]))]:
    t0 = time.time()
    net = fit(net, U, M, ST, y, tr, F)
    p = predict(net, U, M, ST, te, F=F)
    res[name] = {**metrics(y[te], p), "train_seconds": round(time.time() - t0, 1)}
    print(name, res[name], flush=True)
res["delta_full_minus_nopair"] = {k: round(res["full_duplexcnn"][k] - res["no_pairmap"][k], 5)
                                  for k in ("pearson", "spearman", "rmse")}
json.dump(res, open(ROOT / "results" / "ablation_pairmap.json", "w"), indent=1)
print(json.dumps(res["delta_full_minus_nopair"], indent=1))
