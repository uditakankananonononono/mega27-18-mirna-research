"""Per-miRNA conditional delta-AUROC: is the pooled CNN gain broad or driven by a few miRNAs?

For each miRNA in results/clip_rows_136.csv (committed; unit = (miRNA,gene),
genes with a conserved site for the same seed already excluded upstream), fit
5-fold CV logistic models on covariates only [n8, n7, len] vs covariates+CNN
[n8, n7, len, min, sum] and record within-miRNA delta AUROC. Wilcoxon
signed-rank across miRNAs tests breadth. Hermetic: reads only the committed CSV.
"""
import json
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from scipy.stats import wilcoxon

CSV = "results/clip_rows_136.csv"
OUT = "results/clip_per_mirna_delta.json"
BASE = ["n8", "n7", "len"]
CNN = ["min", "sum"]


def per_mirna_delta(df, seed=0, min_pos=50, min_neg=50):
    """Return per-miRNA (base_auroc, cnn_auroc, delta) via 5-fold CV, or None if
    the miRNA lacks class balance. Extracted for testability."""
    y = df["label"].to_numpy()
    if y.sum() < min_pos or (1 - y).sum() < min_neg:
        return None
    Xb = df[BASE].to_numpy(float)
    Xc = df[BASE + CNN].to_numpy(float)
    skf = StratifiedKFold(5, shuffle=True, random_state=seed)
    pb = np.zeros(len(df)); pc = np.zeros(len(df))
    for tr, te in skf.split(Xb, y):
        mb = LogisticRegression(max_iter=1000).fit(Xb[tr], y[tr])
        mc = LogisticRegression(max_iter=1000).fit(Xc[tr], y[tr])
        pb[te] = mb.predict_proba(Xb[te])[:, 1]
        pc[te] = mc.predict_proba(Xc[te])[:, 1]
    ab, ac = roc_auc_score(y, pb), roc_auc_score(y, pc)
    return float(ab), float(ac), float(ac - ab)


def main():
    df = pd.read_csv(CSV)
    rows = []
    for mirna, g in df.groupby("mirna"):
        r = per_mirna_delta(g)
        if r is None:
            continue
        rows.append({"mirna": mirna, "n_pairs": int(len(g)),
                     "n_pos": int(g["label"].sum()),
                     "base_auroc": r[0], "cnn_auroc": r[1], "delta": r[2]})
    d = np.array([r["delta"] for r in rows])
    stat = wilcoxon(d)
    out = {
        "design": "Per-miRNA 5-fold CV logistic AUROC, covariates [n8,n7,len] vs +CNN [min,sum]; "
                  "min_pos/min_neg >= 50. Same committed rows as results/clip_falsification_136.json.",
        "n_mirnas_tested": len(rows),
        "delta_median": float(np.median(d)),
        "delta_iqr": [float(np.percentile(d, 25)), float(np.percentile(d, 75))],
        "frac_positive": float((d > 0).mean()),
        "wilcoxon_stat": float(stat.statistic), "wilcoxon_p": float(stat.pvalue),
        "top5": sorted(rows, key=lambda r: -r["delta"])[:5],
        "bottom5": sorted(rows, key=lambda r: r["delta"])[:5],
        "rows": rows,
    }
    with open(OUT, "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != "rows"}, indent=1)[:1200])


if __name__ == "__main__":
    main()
