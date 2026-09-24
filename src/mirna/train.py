"""Train DuplexCNN to predict TargetScan context++ on held-out genes; compare
to site-type-mean baseline and ridge regression on window one-hots."""
from __future__ import annotations
import json, sys, time
from pathlib import Path
import numpy as np
import torch
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import Ridge

from .data import load_mirnas, load_utrs, build, gene_split
from .model import DuplexCNN

ROOT = Path(__file__).resolve().parents[2]


def metrics(y, p):
    return {"pearson": float(pearsonr(y, p)[0]), "spearman": float(spearmanr(y, p)[0]),
            "rmse": float(np.sqrt(np.mean((y - p) ** 2)))}


def fit_cnn(U, M, ST, y, tr, epochs=6, seed=0, bs=256, lr=2e-3, log=print, F=None):
    torch.manual_seed(seed)
    net = DuplexCNN(n_feat=0 if F is None else F.shape[1])
    Ft = None if F is None else torch.from_numpy(F)
    opt = torch.optim.Adam(net.parameters(), lr=lr)
    Ut, Mt, St, yt = map(torch.from_numpy, (U, M, ST, y))
    for ep in range(epochs):
        perm = np.random.default_rng(seed + ep).permutation(tr)
        net.train(); tot = 0.0
        for i in range(0, len(perm), bs):
            b = torch.from_numpy(perm[i:i + bs])
            loss = ((net(Ut[b], Mt[b], St[b], None if Ft is None else Ft[b]) - yt[b]) ** 2).mean()
            opt.zero_grad(); loss.backward(); opt.step(); tot += loss.item() * len(b)
        log(f"epoch {ep} mse {tot/len(perm):.4f}")
    return net


def predict(net, U, M, ST, idx, bs=1024, F=None):
    net.eval(); out = []
    with torch.no_grad():
        for i in range(0, len(idx), bs):
            b = torch.from_numpy(idx[i:i + bs])
            out.append(net(torch.from_numpy(U)[b], torch.from_numpy(M)[b], torch.from_numpy(ST)[b], None if F is None else torch.from_numpy(F)[b]).numpy())
    return np.concatenate(out)


def main(max_rows=80000, epochs=6):
    d = ROOT / "data"
    mirs = load_mirnas(d / "miR_Family_Info.txt")
    utrs = load_utrs(d / "human_utrs.tsv")
    U, M, ST, y, genes, meta, F = build(d / "human_sites.tsv", utrs, mirs, max_rows=max_rows)
    ST = ST.clip(0, 3)
    tr, te = gene_split(genes)
    print("n", len(y), "train", len(tr), "test", len(te), flush=True)
    res = {"n_sites": int(len(y)), "n_train": int(len(tr)), "n_test": int(len(te)),
           "split": "held-out genes (20%)"}
    # baseline 1: site-type mean
    means = {s: y[tr][ST[tr] == s].mean() for s in np.unique(ST[tr])}
    res["site_type_mean"] = metrics(y[te], np.array([means.get(s, y[tr].mean()) for s in ST[te]]))
    # baseline 2: ridge on flattened one-hots + site type
    X = np.concatenate([U.reshape(len(y), -1), M.reshape(len(y), -1), np.eye(4)[ST]], 1)
    rg = Ridge(alpha=10.0).fit(X[tr], y[tr])
    res["ridge_onehot"] = metrics(y[te], rg.predict(X[te]))
    mu, sd = F[tr].mean(0), F[tr].std(0) + 1e-6
    F = ((F - mu) / sd).astype(np.float32)
    Xf = np.concatenate([X, F], 1)
    rg2 = Ridge(alpha=10.0).fit(Xf[tr], y[tr])
    res["ridge_onehot_plus_context"] = metrics(y[te], rg2.predict(Xf[te]))
    print("ridge", res["ridge_onehot"], res["ridge_onehot_plus_context"], flush=True)
    t0 = time.time()
    net = fit_cnn(U, M, ST, y, tr, epochs=epochs, log=lambda s: print(s, flush=True), F=F)
    p = predict(net, U, M, ST, te, F=F)
    res["duplex_cnn_plus_context"] = metrics(y[te], p)
    res["cnn_train_seconds"] = round(time.time() - t0, 1)
    print(json.dumps(res, indent=1), flush=True)
    (ROOT / "results").mkdir(exist_ok=True)
    json.dump(res, open(ROOT / "results" / "context_pp_benchmark.json", "w"), indent=1)
    torch.save(net.state_dict(), ROOT / "results" / "duplex_cnn.pt")


if __name__ == "__main__":
    main(*(int(a) for a in sys.argv[1:]))
