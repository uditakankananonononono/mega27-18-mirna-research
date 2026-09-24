"""Per-miRNA univariate AUROC: CNN scores vs seed-count/length features.

clip_falsification_136.json reports a pooled delta AUROC; here we ask WHERE
the CNN signal lives: for each miRNA, univariate AUROC of each feature
(min/sum CNN score, n8, n7, log len) against CLIP labels, then a Wilcoxon
paired comparison across miRNAs of CNN-min vs the best seed feature.
Output: results/per_mirna_feature_auroc.json + paper/figs/fig_permirna_scatter.pdf
"""
import csv, json, math
import numpy as np
from scipy.stats import rankdata, wilcoxon
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rows = {}
with open("results/clip_rows_136.csv") as fh:
    for r in csv.DictReader(fh):
        rows.setdefault(r["mirna"], []).append(r)

def auroc(vals, labels):
    v = np.asarray(vals, float); y = np.asarray(labels, int)
    npos = int(y.sum()); nneg = len(y) - npos
    if npos < 10 or nneg < 10:
        return None
    rk = rankdata(v)
    return float((rk[y == 1].sum() - npos * (npos + 1) / 2) / (npos * nneg))

FEATS = ["min", "sum", "n8", "n7", "len"]
per = {}
for mir, rs in rows.items():
    labels = [int(r["label"]) for r in rs]
    per[mir] = {}
    for f in FEATS:
        vals = [math.log1p(float(r[f])) if f in ("n8", "n7", "len") else float(r[f]) for r in rs]
        per[mir][f] = auroc(vals, labels)
    per[mir]["n"] = len(rs)
    per[mir]["npos"] = int(sum(labels))

ok = {m: d for m, d in per.items() if all(d[f] is not None for f in FEATS)}
print(f"miRNAs with usable labels: {len(ok)} / {len(per)}", flush=True)

cnn = np.array([ok[m]["min"] for m in ok])
bestseed = np.array([max(ok[m]["n8"], ok[m]["n7"], ok[m]["len"]) for m in ok])
w, p = wilcoxon(cnn - bestseed)
win = int((cnn > bestseed).sum())
top = sorted(ok.items(), key=lambda kv: -kv[1]["min"])[:10]
print(f"CNN-min AUROC median {np.median(cnn):.4f} vs best-seed median {np.median(bestseed):.4f}")
print(f"CNN-min wins in {win}/{len(ok)} miRNAs; Wilcoxon p={p:.3g}")
print("top CNN-min miRNAs:", [(m, round(d["min"], 3)) for m, d in top])

fig, ax = plt.subplots(figsize=(4.2, 4.0))
ax.scatter(bestseed, cnn, s=10, alpha=0.6)
lo, hi = 0.35, 0.75
ax.plot([lo, hi], [lo, hi], ls="--", color="gray", lw=0.8)
ax.set_xlabel("best seed-feature AUROC (per miRNA)")
ax.set_ylabel("CNN min-score AUROC (per miRNA)")
ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
fig.tight_layout(); fig.savefig("paper/figs/fig_permirna_scatter.pdf")
print("saved fig_permirna_scatter.pdf")

json.dump({"n_mirnas": len(ok), "features": FEATS,
           "cnn_min_median": round(float(np.median(cnn)), 4),
           "best_seed_median": round(float(np.median(bestseed)), 4),
           "cnn_wins": win, "wilcoxon_p": float(p),
           "top10_cnn_min": [(m, round(d["min"], 4)) for m, d in top],
           "per_mirna": {m: {f: round(d[f], 4) for f in FEATS} for m, d in ok.items()}},
          open("results/per_mirna_feature_auroc.json", "w"), indent=1)
print("saved results/per_mirna_feature_auroc.json")
