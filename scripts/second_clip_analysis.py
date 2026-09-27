"""Locked replication analysis on TarBase v9 labels (lane-18 item 8), per
docs/PREREG_SECOND_CLIP_20260927.md. Mirrors the primary endpoint
(scripts/expression_confound.py) exactly: same features, same models
(logistic regression, gene-grouped 5-fold CV), same expression covariate
(HPA HEK293 nTPM, log1p). No re-tuning, no feature additions.

Primary replication metric: pooled delta AUROC (AEC - AE).
Decision rule (locked): replication supports the primary claim only if the
pooled delta is positive AND its within-resource bootstrap CI (B=1000, rows
resampled by gene, seed 0) excludes 0. Alpha 0.01 one-sided in the direction
of the committed point estimate is evaluated downstream in the paper/queue;
this script reports the numbers.

Secondary: per-miRNA deltas on the 120 eligible miRNAs
(results/second_clip_label_tarbase_audit.json), then miRNA-cluster
percentile bootstrap (10,000 replicates, seed 20260927) of the median
per-miRNA delta, matching scripts/cluster_interval.py.
"""
import json
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold, StratifiedKFold
from sklearn.metrics import roc_auc_score

CSV = "results/second_clip_rows_tarbase.csv"
EXPR = "results/hpa_hek293_ntpm.tsv"
AUDIT = "results/second_clip_label_tarbase_audit.json"
OUT = "results/second_clip_analysis.json"
A = ["n8", "n7", "len"]
E = ["lexpr"]
C = ["min", "sum"]

def attach_expression(df, expr):
    e = expr.groupby("Gene name")["nTPM"].max()
    out = df.copy()
    out["nTPM"] = out["gene"].map(e)
    unmatched = int(out.loc[out["nTPM"].isna(), "gene"].nunique())
    out = out.dropna(subset=["nTPM"])
    out["lexpr"] = np.log1p(out["nTPM"].astype(float))
    return out, unmatched

def cv_scores(X, y, groups):
    p = np.zeros(len(y))
    for tr, te in GroupKFold(5).split(X, y, groups):
        m = LogisticRegression(max_iter=2000).fit(X[tr], y[tr])
        p[te] = m.predict_proba(X[te])[:, 1]
    return p

def standardize(X):
    mu, sd = X.mean(0), X.std(0)
    sd[sd == 0] = 1
    return (X - mu) / sd

def per_mirna(g, seed=0):
    y = g["label"].to_numpy()
    if y.sum() < 50 or (1 - y).sum() < 50:
        return None
    Xa = standardize(g[A + E].to_numpy(float))
    Xc = standardize(g[A + E + C].to_numpy(float))
    skf = StratifiedKFold(5, shuffle=True, random_state=seed)
    pa = np.zeros(len(y)); pc = np.zeros(len(y))
    for tr, te in skf.split(Xa, y):
        pa[te] = LogisticRegression(max_iter=1000).fit(Xa[tr], y[tr]).predict_proba(Xa[te])[:, 1]
        pc[te] = LogisticRegression(max_iter=1000).fit(Xc[tr], y[tr]).predict_proba(Xc[te])[:, 1]
    return roc_auc_score(y, pa), roc_auc_score(y, pc)

def main():
    df = pd.read_csv(CSV)
    expr = pd.read_csv(EXPR, sep="\t")
    eligible = set(json.load(open(AUDIT))["summary"]["eligible"])
    n_rows0, n_genes0 = len(df), df["gene"].nunique()
    df, unmatched = attach_expression(df, expr)
    y = df["label"].to_numpy()
    groups = df["gene"].to_numpy()
    res = {}
    for name, cols in [("A", A), ("AE", A + E), ("AC", A + C), ("AEC", A + E + C)]:
        p = cv_scores(standardize(df[cols].to_numpy(float)), y, groups)
        res[name] = p
        print("cv done", name, flush=True)
    aurocs = {k: float(roc_auc_score(y, p)) for k, p in res.items()}
    delta = aurocs["AEC"] - aurocs["AE"]
    # within-resource bootstrap CI: rows resampled by gene, B=1000, seed 0
    rng = np.random.default_rng(0)
    ug = np.unique(groups)
    gene_rows = {g: np.where(groups == g)[0] for g in ug}
    boots = np.empty(1000)
    for b in range(1000):
        gs = rng.choice(ug, size=len(ug), replace=True)
        idx = np.concatenate([gene_rows[g] for g in gs])
        boots[b] = roc_auc_score(y[idx], res["AEC"][idx]) - roc_auc_score(y[idx], res["AE"][idx])
        if b % 100 == 0: print("boot", b, flush=True)
    ci = [float(v) for v in np.percentile(boots, [2.5, 97.5])]
    ci99 = [float(v) for v in np.percentile(boots, [0.5, 99.5])]
    # secondary: per-miRNA deltas on eligible miRNAs
    rows = []
    for mirna, g in df.groupby("mirna"):
        if mirna not in eligible:
            continue
        r = per_mirna(g)
        if r is None:
            continue
        rows.append({"mirna": mirna, "auroc_AE": float(r[0]), "auroc_AEC": float(r[1]),
                     "delta": float(r[1] - r[0])})
        print("per-mirna", mirna, flush=True)
    deltas = np.array([r["delta"] for r in rows])
    rng2 = np.random.default_rng(20260927)
    B = 10000
    meds = np.empty(B)
    for i in range(B):
        meds[i] = np.median(rng2.choice(deltas, size=len(deltas), replace=True))
    out = {
        "design": __doc__.strip().splitlines()[0],
        "resource": "DIANA TarBase v9.0, AGO-CLIP whitelist (HITS-CLIP/PAR-CLIP/CLASH/qCLASH)",
        "rows_source": CSV, "n_rows_input": int(n_rows0), "n_genes_input": int(n_genes0),
        "n_rows_matched": int(len(df)), "n_genes_unmatched": unmatched,
        "n_pos": int(y.sum()), "n_neg": int((1 - y).sum()),
        "pooled_auroc_gene_grouped_cv": aurocs,
        "pooled_delta_AEC_minus_AE": float(delta),
        "bootstrap_gene_cluster_B1000_seed0": {"ci95": ci, "ci99": ci99,
            "frac_positive": float((boots > 0).mean())},
        "eligible_mirnas": len(eligible),
        "per_mirna_n": len(rows),
        "per_mirna_delta_median": float(np.median(deltas)),
        "per_mirna_delta_iqr": [float(np.percentile(deltas, 25)), float(np.percentile(deltas, 75))],
        "per_mirna_frac_positive": float((deltas > 0).mean()),
        "cluster_bootstrap": {"seed": 20260927, "B": B,
            "median_delta": float(np.median(deltas)),
            "median_ci95": [float(v) for v in np.percentile(meds, [2.5, 97.5])]},
        "rows": rows,
    }
    with open(OUT, "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != "rows"}, indent=1), flush=True)

if __name__ == "__main__":
    main()
