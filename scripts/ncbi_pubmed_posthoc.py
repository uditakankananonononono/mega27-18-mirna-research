"""NOT PRE-REGISTERED sensitivity analysis for ncbi_pubmed_bias.py (written after gate G2 failed).

G2 (>= 95% symbol mapping) failed at 89.8%; the pre-registered script scores unmapped symbols
as 0 publications, which would make sets rich in unmapped symbols look less studied. Here:
(a) unmapped fraction in CNN top-200 vs uniform draws; (b) P1-P3 and P4/P5 recomputed after
dropping unmapped genes from the universe before any set is drawn. Seeds f"{mi}-pmph-<arm>".
"""
import json, random, sys
import numpy as np, pandas as pd
sys.path.insert(0, "scripts")
from gprofiler_enrichment import K, signtest
from string_exprmatched import matched


def main(out="results/ncbi_pubmed_posthoc.json"):
    gi = pd.read_csv("data/ncbi_gene/Homo_sapiens.gene_info.gz", sep="\t", compression="gzip", usecols=["GeneID", "Symbol"])
    pm = pd.read_csv("data/ncbi_gene/human_gene2pubmed_counts.tsv", sep="\t")
    lp = np.log10(1 + gi.merge(pm, on="GeneID", how="left").fillna({"n_pubmed": 0}).groupby("Symbol").n_pubmed.max())
    lo = pd.read_csv("data/gnomad/lof_v211.txt.bgz", sep="\t", compression="gzip", usecols=["gene", "oe_lof_upper"]).dropna().groupby("gene").oe_lof_upper.min()
    h = pd.read_csv("data/hpo/genes_to_disease.txt", sep="\t"); dis = set(h[h.association_type == "MENDELIAN"].gene_symbol)
    panel = json.load(open("results/ncbi_pubmed_bias.json"))["panel"]
    d0 = pd.read_csv("results/clip_rows_136.csv", usecols=["mirna", "gene", "min", "len"]).drop_duplicates(["mirna", "gene"])
    d0 = d0[d0.mirna.isin(panel)]
    expr = pd.read_csv("results/hpa_hek293_ntpm.tsv", sep="\t").groupby("Gene name")["nTPM"].max()
    med = lambda gs: float(np.median(lp.reindex(gs).values)); mlo = lambda gs: float(np.nanmedian(lo.reindex(gs).values)); nd = lambda gs: int(sum(g in dis for g in gs))
    per = {}
    for mi in panel:
        s0 = d0[d0.mirna == mi]; cnn0 = s0.nsmallest(K, "min").gene
        un_cnn = float((~cnn0.isin(lp.index)).mean()); un_all = float((~s0.gene.isin(lp.index)).mean())
        s = s0[s0.gene.isin(lp.index)].reset_index(drop=True); U = s.gene.tolist(); cnn = s.nsmallest(K, "min").gene.tolist()
        se, sl, sp = s.assign(x=s.gene.map(expr)), s.assign(x=s["len"].astype(float)), s.assign(x=s.gene.map(lp))
        ru, re_, rl, rp = (random.Random(f"{mi}-pmph-{a}") for a in ("unif", "expr", "len", "pub"))
        P = [matched(sp, cnn, rp) for _ in range(3)]
        per[mi] = {"unmapped_cnn": un_cnn, "unmapped_universe": un_all, "cnn": med(cnn), "uniform_mean": float(np.mean([med(ru.sample(U, K)) for _ in range(3)])),
                   "expr_mean": float(np.mean([med(matched(se, cnn, re_)) for _ in range(3)])), "len_mean": float(np.mean([med(matched(sl, cnn, rl)) for _ in range(3)])),
                   "cnn_loeuf": mlo(cnn), "pub_loeuf_mean": float(np.mean([mlo(x) for x in P])), "cnn_hpo": nd(cnn), "pub_hpo_mean": float(np.mean([nd(x) for x in P]))}
    V = list(per.values()); J = {"note": "NOT pre-registered; unmapped symbols dropped", "panel": panel,
         "unmapped_cnn_higher": dict(signtest(sum(v["unmapped_cnn"] > v["unmapped_universe"] for v in V), sum(v["unmapped_cnn"] < v["unmapped_universe"] for v in V)),
                                      median_cnn=float(np.median([v["unmapped_cnn"] for v in V])), median_universe=float(np.median([v["unmapped_universe"] for v in V])))}
    for name, a, b, sign in (("P1_less_studied_vs_uniform", "cnn", "uniform_mean", -1), ("P2_vs_expr", "cnn", "expr_mean", -1), ("P3_vs_utrlen", "cnn", "len_mean", -1),
                             ("P4_loeuf_vs_pubmatched", "cnn_loeuf", "pub_loeuf_mean", 1), ("P5_hpo_vs_pubmatched", "cnn_hpo", "pub_hpo_mean", -1)):
        w = sum((v[a] - v[b]) * sign > 0 for v in V); l = sum((v[a] - v[b]) * sign < 0 for v in V)
        J[name] = dict(signtest(w, l), median_diff=float(np.median([v[a] - v[b] for v in V])))
    J["per_mirna"] = per
    json.dump(J, open(out, "w"), indent=1); print(json.dumps({k: J[k] for k in J if k not in ("per_mirna", "panel")}))


if __name__ == "__main__":
    main(*sys.argv[1:])
