import sys
from pathlib import Path
import numpy as np, torch
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from mirna.data import one_hot, utr_window, gene_split
from mirna.model import DuplexCNN, pair_map
from mirna.train import fit_cnn, predict


def test_window_coords_pad():
    w = utr_window("ACGU", 1, 2, win=6)
    assert len(w) == 6 and w.startswith("N") and "AC" in w


def test_window_centres_site():
    utr = "A" * 30 + "GGGGGGG" + "A" * 30
    w = utr_window(utr, 31, 37, win=11)
    assert w == "AA" + "GGGGGGG" + "AA"


def test_pair_map_complement():
    u = torch.from_numpy(one_hot("AC", 2))[None]
    m = torch.from_numpy(one_hot("GU", 2))[None]
    pm = pair_map(u, m)[0, 0]
    assert pm[1, 0] == 1  # U pairs A
    assert pm[0, 1] == 1  # G pairs C
    assert pm[0, 0] == 0


def test_gene_split_disjoint():
    g = np.array(list("aabbccddeeffgghhiijj"))
    tr, te = gene_split(g, 0.3)
    assert not set(g[tr]) & set(g[te]) and len(te) > 0


def test_cnn_learns_seed_signal():
    rng = np.random.default_rng(0)
    n = 600
    seqs = ["".join(rng.choice(list("ACGU"), 50)) for _ in range(n)]
    mir = "UAGCUUAUCAGACUGAUGUUGA"
    site = "AUAAGCUA"  # reverse complement of seed region, 8mer
    y = np.zeros(n, dtype=np.float32)
    for i in range(0, n, 2):
        seqs[i] = seqs[i][:21] + site + seqs[i][29:]; y[i] = -1.0
    U = np.stack([one_hot(s, 50) for s in seqs]); M = np.stack([one_hot(mir, 22)] * n)
    ST = np.zeros(n, dtype=np.int64)
    idx = np.arange(n); tr, te = idx[:500], idx[500:]
    net = fit_cnn(U, M, ST, y, tr, epochs=8, log=lambda s: None)
    p = predict(net, U, M, ST, te)
    assert np.corrcoef(p, y[te])[0, 1] > 0.8
