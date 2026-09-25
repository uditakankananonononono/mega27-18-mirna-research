"""OmniPath non-miRTarBase curated interactions as independent labels (tool 28).

PRE-REGISTERED (written after counting positives per miRNA, before any score was compared):
Labels: OmniPath /interactions?datasets=mirnatarget (human), keeping only rows whose
sources do NOT include miRTarBase (ncRDeathDB, miRecords, miR2Disease, miRDeathDB,
SIGNOR). MIMAT -> name via TargetScan miR_Family_Info. Universe per miRNA: the
CLIP-scored genes in results/clip_rows_136.csv (non-conserved seed-match genes).
Include miRNAs with >= 5 positives in their universe.
Scores: CNN (-min), site count (n8*10 + n7, ties broken by len; close to, but not identical with,
the lexicographic Enrichr site-count arm), HEK293 expression (HPA nTPM).
Gates:
 G1 >= 10 eligible miRNAs.
 G2 label-permutation control (200 seeded shuffles per miRNA): median of the
    per-miRNA permutation-mean CNN AUROC in [0.45, 0.55].
Hypotheses (one-sided sign tests across eligible miRNAs, alpha 0.05):
 H1 CNN AUROC > 0.5.  H2 CNN AUROC > site-count AUROC.  H3 CNN AUROC > expression AUROC.
"""
import csv, json, sys
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
sys.path.insert(0, "scripts")
from gprofiler_enrichment import signtest


def main(out="results/omnipath_independent.json"):
    o = pd.read_csv("data/omnipath/mirnatarget_human.tsv", sep="\t")
    n_all = len(o); o = o[~o.sources.str.contains("miRTarBase")]
    acc2name = {r["MiRBase Accession"]: r["MiRBase ID"] for r in csv.DictReader(open("data/miR_Family_Info.txt"), delimiter="\t") if r["Species ID"] == "9606"}
    o = o.assign(mirna=o.source.map(acc2name)).dropna(subset=["mirna"])
    pos = o.groupby("mirna").target_genesymbol.apply(set)
    d = pd.read_csv("results/clip_rows_136.csv", usecols=["mirna", "gene", "min", "n8", "n7", "len"]).drop_duplicates(["mirna", "gene"])
    expr = pd.read_csv("results/hpa_hek293_ntpm.tsv", sep="\t").groupby("Gene name")["nTPM"].max()
    per = {}; perm_means = []
    for mi, s in d.groupby("mirna"):
        if mi not in pos.index: continue
        y = s.gene.isin(pos[mi]).values.astype(int)
        if y.sum() < 5: continue
        cnn = -s["min"].values; sc = s.n8.values * 10 + s.n7.values + s["len"].values / 1e6
        ex = s.gene.map(expr).fillna(0).values
        rng = np.random.default_rng(sum(map(ord, mi)))
        pm = float(np.mean([roc_auc_score(rng.permutation(y), cnn) for _ in range(200)])); perm_means.append(pm)
        per[mi] = {"n": len(y), "n_pos": int(y.sum()), "auroc_cnn": roc_auc_score(y, cnn), "auroc_sitecount": roc_auc_score(y, sc),
                   "auroc_expr": roc_auc_score(y, ex), "perm_mean_cnn": pm}
        print(mi, per[mi]["n_pos"], round(per[mi]["auroc_cnn"], 3), round(per[mi]["auroc_sitecount"], 3), round(per[mi]["auroc_expr"], 3), round(pm, 3), flush=True)
    st = lambda f: signtest(sum(f(v) > 0 for v in per.values()), sum(f(v) < 0 for v in per.values()))
    vd = lambda h: "CONFIRMED" if h["p_one_sided"] < 0.05 else "FALSIFIED"
    h1 = st(lambda v: v["auroc_cnn"] - 0.5); h2 = st(lambda v: v["auroc_cnn"] - v["auroc_sitecount"]); h3 = st(lambda v: v["auroc_cnn"] - v["auroc_expr"])
    med = lambda k: float(np.median([v[k] for v in per.values()]))
    J = {"tool": "OmniPath web service (omnipathdb.org /interactions, datasets=mirnatarget)", "n_rows_all": n_all, "n_rows_non_mirtarbase": int(len(o)),
         "gates": {"G1_n_eligible": len(per), "G1_pass": len(per) >= 10, "G2_median_perm_mean": float(np.median(perm_means)),
                   "G2_pass": 0.45 <= float(np.median(perm_means)) <= 0.55},
         "median_auroc": {"cnn": med("auroc_cnn"), "sitecount": med("auroc_sitecount"), "expr": med("auroc_expr")},
         "H1_cnn_above_chance": dict(h1, verdict=vd(h1)), "H2_cnn_vs_sitecount": dict(h2, verdict=vd(h2)), "H3_cnn_vs_expr": dict(h3, verdict=vd(h3)),
         "per_mirna": per}
    json.dump(J, open(out, "w"), indent=1)
    print(json.dumps({k: J[k] for k in ["n_rows_non_mirtarbase", "gates", "median_auroc", "H1_cnn_above_chance", "H2_cnn_vs_sitecount", "H3_cnn_vs_expr"]}))


if __name__ == "__main__":
    main(*sys.argv[1:])
