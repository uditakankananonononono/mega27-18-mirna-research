"""Expression confound control for the CLIP falsification (HPA HEK293 RNA).

ENCORI AGO-CLIP positives come largely from HEK293 libraries, so a gene that is
highly expressed in HEK293 is more likely to yield CLIP reads regardless of
miRNA targeting. Question: does the CNN's conditional CLIP gain survive adding
HEK293 expression (Human Protein Atlas rna_celline, nTPM) as a covariate?

Inputs (committed): results/clip_rows_136.csv, results/hpa_hek293_ntpm.tsv
(Gene name, nTPM; derived from HPA rna_celline.tsv.zip, cell line HEK293).
Models (5-fold CV logistic, grouped by gene so no gene spans train and test):
  A  = [n8,n7,len]                      site-count covariates
  AE = A + [log1p nTPM]                 + expression
  AEC= AE + [min,sum]                   + CNN
Pooled AUROCs and per-miRNA delta(AEC-AE) with Wilcoxon signed-rank.
"""
import json
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold, StratifiedKFold
from sklearn.metrics import roc_auc_score
from scipy.stats import wilcoxon, mannwhitneyu

CSV = "results/clip_rows_136.csv"
EXPR = "results/hpa_hek293_ntpm.tsv"
OUT = "results/expression_confound.json"
A = ["n8", "n7", "len"]
E = ["lexpr"]
C = ["min", "sum"]


def attach_expression(df, expr):
    """Join log1p(nTPM) by gene symbol; returns (joined df, n_unmatched_genes)."""
    e = expr.groupby("Gene name")["nTPM"].max()
    out = df.copy()
    out["nTPM"] = out["gene"].map(e)
    unmatched = int(out.loc[out["nTPM"].isna(), "gene"].nunique())
    out = out.dropna(subset=["nTPM"])
    out["lexpr"] = np.log1p(out["nTPM"].astype(float))
    return out, unmatched


def cv_scores(X, y, groups, seed=0):
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
    n_rows0, n_genes0 = len(df), df["gene"].nunique()
    df, unmatched = attach_expression(df, expr)
    y = df["label"].to_numpy()
    groups = df["gene"].to_numpy()
    res = {}
    for name, cols in [("expr_only", E), ("A", A), ("AE", A + E), ("AC", A + C), ("AEC", A + E + C)]:
        p = cv_scores(standardize(df[cols].to_numpy(float)), y, groups)
        res[name] = float(roc_auc_score(y, p))
    # effect of expression on labels
    mw = mannwhitneyu(df.loc[df.label == 1, "lexpr"], df.loc[df.label == 0, "lexpr"])
    rows = []
    for mirna, g in df.groupby("mirna"):
        r = per_mirna(g)
        if r is None:
            continue
        rows.append({"mirna": mirna, "auroc_AE": float(r[0]), "auroc_AEC": float(r[1]),
                     "delta": float(r[1] - r[0])})
    d = np.array([r["delta"] for r in rows])
    w = wilcoxon(d)
    out = {
        "design": __doc__.strip().splitlines()[0],
        "source": "Human Protein Atlas rna_celline.tsv.zip (last-modified 2025-11-05), cell line HEK293, nTPM",
        "n_rows_input": int(n_rows0), "n_genes_input": int(n_genes0),
        "n_rows_matched": int(len(df)), "n_genes_unmatched": unmatched,
        "median_log1p_ntpm_pos": float(df.loc[df.label == 1, "lexpr"].median()),
        "median_log1p_ntpm_neg": float(df.loc[df.label == 0, "lexpr"].median()),
        "mannwhitney_p_expr_pos_vs_neg": float(mw.pvalue),
        "pooled_auroc_gene_grouped_cv": res,
        "pooled_delta_cnn_given_A": res["AC"] - res["A"],
        "pooled_delta_cnn_given_AE": res["AEC"] - res["AE"],
        "per_mirna_n": len(rows),
        "per_mirna_delta_median": float(np.median(d)),
        "per_mirna_delta_iqr": [float(np.percentile(d, 25)), float(np.percentile(d, 75))],
        "per_mirna_frac_positive": float((d > 0).mean()),
        "per_mirna_wilcoxon_p": float(w.pvalue),
        "rows": rows,
    }
    with open(OUT, "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != "rows"}, indent=1))


if __name__ == "__main__":
    main()
