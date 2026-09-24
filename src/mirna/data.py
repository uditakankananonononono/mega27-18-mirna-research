"""Build miRNA:UTR site windows from TargetScan vert_80 human data."""
from __future__ import annotations
import csv
from pathlib import Path
import numpy as np

BASES = "ACGU"
IDX = {b: i for i, b in enumerate(BASES)}
WIN = 50          # UTR window (nt) centred on the site
MIR_LEN = 22


def one_hot(seq: str, length: int) -> np.ndarray:
    x = np.zeros((4, length), dtype=np.float32)
    for i, c in enumerate(seq[:length]):
        j = IDX.get(c)
        if j is not None:
            x[j, i] = 1.0
    return x


def load_mirnas(path: Path) -> dict[str, str]:
    out = {}
    with open(path) as fh:
        r = csv.reader(fh, delimiter="\t"); next(r)
        for row in r:
            if row[2] == "9606":
                out[row[3]] = row[4].upper().replace("T", "U")
    return out


def load_utrs(path: Path) -> dict[str, str]:
    out = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            t, g, s = line.rstrip("\n").split("\t")
            out[t.split(".")[0]] = s.upper().replace("T", "U")
    return out


def utr_window(utr: str, start: int, end: int, win: int = WIN) -> str:
    """start/end are 1-based inclusive UTR coordinates (TargetScan)."""
    centre = (start - 1 + end) // 2
    lo = centre - win // 2
    seg = "".join(utr[i] if 0 <= i < len(utr) else "N" for i in range(lo, lo + win))
    return seg


def build(sites_path: Path, utrs: dict, mirs: dict, max_rows: int | None = None, seed: int = 0):
    rows = []
    with open(sites_path) as fh:
        r = csv.reader(fh, delimiter="\t"); next(r)
        for t, g, m, st, s, e, c in r:
            if t not in utrs or m not in mirs or int(st) < 1:
                continue
            rows.append((t, g, m, int(st), int(s), int(e), float(c)))
    rng = np.random.default_rng(seed)
    if max_rows and len(rows) > max_rows:
        rows = [rows[i] for i in rng.choice(len(rows), max_rows, replace=False)]
    U = np.stack([one_hot(utr_window(utrs[t], s, e), WIN) for t, g, m, st, s, e, c in rows])
    M = np.stack([one_hot(mirs[m], MIR_LEN) for t, g, m, st, s, e, c in rows])
    ST = np.array([st for *_, st, s, e, c in [(r[0], r[1], r[2], r[3], r[4], r[5], r[6]) for r in rows]], dtype=np.int64)
    y = np.array([r[6] for r in rows], dtype=np.float32)
    genes = np.array([r[1] for r in rows])
    meta = rows
    return U, M, ST, y, genes, meta


def gene_split(genes: np.ndarray, frac_test: float = 0.2, seed: int = 0):
    uniq = np.unique(genes)
    rng = np.random.default_rng(seed)
    test = set(rng.choice(uniq, int(len(uniq) * frac_test), replace=False))
    mask = np.array([g in test for g in genes])
    return np.where(~mask)[0], np.where(mask)[0]
