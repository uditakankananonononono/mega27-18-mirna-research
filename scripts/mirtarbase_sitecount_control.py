"""Site-count control for the miRTarBase benchmark (item 18).

The CLIP per-miRNA feature audit found the CNN min-score's univariate signal is
explained by seed-site count. This script reruns the miRTarBase evaluation rows
(results/mirtarbase_rows.csv) with the pure seed-site COUNT as the predictor, on
exactly the same (miRNA, gene, label) rows, to test whether the literature
benchmark AUROC (0.6829 pooled) survives the same confound check.
Output: results/mirtarbase_sitecount_control.json
"""
import csv, json
from collections import defaultdict
from pathlib import Path
import numpy as np
from sklearn.metrics import roc_auc_score
import sys
sys.path.insert(0, "src")
from mirna.data import load_utrs, load_mirnas
from mirna.cli import rc

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / "data"

utrs = load_utrs(D / "human_utrs.tsv")
tid2gene, gene2tids = {}, {}
with open(D / "human_utrs.tsv") as fh:
    next(fh)
    for line in fh:
        t, g, _ = line.rstrip("\n").split("\t", 2)
        tid2gene[t.split(".")[0]] = g
        gene2tids.setdefault(g, []).append(t.split(".")[0])
mirs = load_mirnas(D / "miR_Family_Info.txt")

# rows for this benchmark
rows = defaultdict(dict)  # mi -> gene -> label
with open(ROOT / "results" / "mirtarbase_rows.csv") as fh:
    for r in csv.DictReader(fh):
        rows[r["mirna"]][r["gene"]] = int(r["label"])

def seed_count(seq, g):
    seed8, seed7 = rc(seq[1:8]), rc(seq[1:7])
    n = 0
    for t in gene2tids.get(g, []):
        u = utrs[t]
        i = u.find(seed8)
        while i >= 0:
            n += 1; i = u.find(seed8, i + 1)
        i = u.find(seed7)
        while i >= 0:
            if i + 6 < len(u) and u[i + 6] == "A" and u[i:i + 7] != seed8:
                n += 1
            i = u.find(seed7, i + 1)
    return n

pooled_y, pooled_c, per_mir = [], [], {}
for mi, genes in sorted(rows.items()):
    seq = mirs.get(mi)
    if seq is None:
        continue
    y, c = [], []
    for g, lab in genes.items():
        y.append(lab); c.append(seed_count(seq, g))
    y = np.array(y); c = np.array(c)
    pooled_y += y.tolist(); pooled_c += c.tolist()
    auroc = round(float(roc_auc_score(y, c)), 4) if 0 < y.sum() < len(y) else None
    per_mir[mi] = {"auroc_site_count": auroc, "n": int(len(y))}

pooled_y = np.array(pooled_y); pooled_c = np.array(pooled_c)
out = {
    "predictor": "seed-site count per gene (8mer + 7mer-A1, same enumeration as benchmark)",
    "n_rows": int(len(pooled_y)),
    "pooled_auroc_site_count": round(float(roc_auc_score(pooled_y, pooled_c)), 4),
    "median_per_mirna_auroc_site_count": round(float(np.median([v["auroc_site_count"] for v in per_mir.values() if v["auroc_site_count"] is not None])), 4),
    "reference_cnn_min_pooled_auroc": 0.6829,
    "reference_cnn_min_median_per_mirna": 0.6926,
    "per_mirna": per_mir,
}
json.dump(out, open(ROOT / "results" / "mirtarbase_sitecount_control.json", "w"), indent=1)
print(json.dumps({k: v for k, v in out.items() if k != "per_mirna"}, indent=1))
