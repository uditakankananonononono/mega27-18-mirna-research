"""PRE-REGISTERED study-bias test with NCBI Gene (tool 37).

Data: NCBI Gene FTP Homo_sapiens.gene_info.gz (GeneID -> Symbol) and gene2pubmed.gz streamed
and reduced to human per-GeneID PubMed counts (data/ncbi_gene/human_gene2pubmed_counts.tsv;
each file counted once). Question: are CNN candidates less studied, and does study bias
explain the dosage/disease depletion (tools 30-33)? Committed BEFORE running.
Panel: 13 discovery + 30 held-out miRNAs; CLIP universe; K=200. Genes missing from
gene2pubmed get count 0. Statistic: median log10(1 + PubMed count).
Comparators (3 draws, seeds f"{mi}-pm-<arm>"): uniform, expression-matched, UTR-length-matched,
publication-matched (decile matching on log10(1+count), string_exprmatched.matched).
Positive control: TargetScan conserved top-200 more studied than 3 uniform random UTR genes.
Gates: G1 positive control sign test p<0.05 (one-sided, ties excluded);
       G2 >= 95% of pooled universe symbols map to a GeneID.
Hypotheses, one-sided sign tests over 43 miRNAs, ties excluded, alpha 0.05:
 P1 CNN less studied than uniform; P2 than expression-matched; P3 than UTR-length-matched.
 P4 CNN median LOEUF (gnomAD v2.1.1) > publication-matched (the LOEUF depletion survives study bias).
 P5 CNN HPO Mendelian count < publication-matched (the HPO depletion survives study bias).
"""
import json, random, sys
import numpy as np, pandas as pd
sys.path.insert(0, "scripts")
from gprofiler_enrichment import PANEL, K, signtest
from string_exprmatched import matched


def main(out="results/ncbi_pubmed_bias.json"):
    gi = pd.read_csv("data/ncbi_gene/Homo_sapiens.gene_info.gz", sep="\t", compression="gzip", usecols=["GeneID", "Symbol"])
    pm = pd.read_csv("data/ncbi_gene/human_gene2pubmed_counts.tsv", sep="\t")
    cnt = gi.merge(pm, on="GeneID", how="left").fillna({"n_pubmed": 0}).groupby("Symbol").n_pubmed.max()
    lp = np.log10(1 + cnt)
    lo = pd.read_csv("data/gnomad/lof_v211.txt.bgz", sep="\t", compression="gzip", usecols=["gene", "oe_lof_upper"]).dropna().groupby("gene").oe_lof_upper.min()
    h = pd.read_csv("data/hpo/genes_to_disease.txt", sep="\t"); dis = set(h[h.association_type == "MENDELIAN"].gene_symbol)
    panel = PANEL + json.load(open("results/gnomad_replication.json"))["panel"]
    d = pd.read_csv("results/clip_rows_136.csv", usecols=["mirna", "gene", "min", "len"]).drop_duplicates(["mirna", "gene"])
    d = d[d.mirna.isin(panel)]
    map_rate = float(pd.Series(d.gene.unique()).isin(cnt.index).mean())
    expr = pd.read_csv("results/hpa_hek293_ntpm.tsv", sep="\t").groupby("Gene name")["nTPM"].max()
    utrg = sorted(set(pd.read_csv("data/human_utrs.tsv", sep="\t", usecols=[1]).iloc[:, 0]))
    ts = pd.read_csv("data/human_sites.tsv", sep="\t", usecols=["gene", "mirna", "context_pp"])
    med = lambda gs: float(np.median(lp.reindex(gs).fillna(0).values))
    mlo = lambda gs: float(np.nanmedian(lo.reindex(gs).values))
    nd = lambda gs: int(sum(g in dis for g in gs))
    per = {}
    for mi in panel:
        s = d[d.mirna == mi].reset_index(drop=True); U = s.gene.tolist(); cnn = s.nsmallest(K, "min").gene.tolist()
        se, sl, sp = s.assign(x=s.gene.map(expr)), s.assign(x=s["len"].astype(float)), s.assign(x=s.gene.map(lp).fillna(0))
        ru, re_, rl, rp, rb = (random.Random(f"{mi}-pm-{a}") for a in ("unif", "expr", "len", "pub", "utr"))
        P = [matched(sp, cnn, rp) for _ in range(3)]
        row = {"universe": len(U), "cnn": med(cnn), "uniform_mean": float(np.mean([med(ru.sample(U, K)) for _ in range(3)])),
               "expr_mean": float(np.mean([med(matched(se, cnn, re_)) for _ in range(3)])), "len_mean": float(np.mean([med(matched(sl, cnn, rl)) for _ in range(3)])),
               "pub_matched_logpub_mean": float(np.mean([med(x) for x in P])),
               "cnn_loeuf": mlo(cnn), "pub_loeuf_mean": float(np.mean([mlo(x) for x in P])), "cnn_hpo": nd(cnn), "pub_hpo_mean": float(np.mean([nd(x) for x in P]))}
        t = ts[(ts.mirna == mi) & ts.gene.isin(utrg)].groupby("gene").context_pp.min().nsmallest(K); row["ts_n"] = len(t)
        if len(t) >= 20:
            row["ts"] = med(list(t.index)); row["ts_rand_mean"] = float(np.mean([med(rb.sample(utrg, len(t))) for _ in range(3)]))
        per[mi] = row
    tsr = [v for v in per.values() if "ts" in v]
    g1 = signtest(sum(v["ts"] > v["ts_rand_mean"] for v in tsr), sum(v["ts"] < v["ts_rand_mean"] for v in tsr))
    J = {"tool": "NCBI Gene FTP (gene_info + gene2pubmed)", "n_symbols": int(len(cnt)), "panel": panel, "K": K,
         "gates": {"G1_ts_vs_random": g1, "G1_n": len(tsr), "G1_pass": bool(g1["p_one_sided"] < 0.05 and g1["wins"] > g1["losses"]), "G2_map_rate": map_rate, "G2_pass": bool(map_rate >= 0.95)}}
    V = list(per.values())
    for name, a, b, sign in (("P1_less_studied_vs_uniform", "cnn", "uniform_mean", -1), ("P2_vs_expr", "cnn", "expr_mean", -1), ("P3_vs_utrlen", "cnn", "len_mean", -1),
                             ("P4_loeuf_vs_pubmatched", "cnn_loeuf", "pub_loeuf_mean", 1), ("P5_hpo_vs_pubmatched", "cnn_hpo", "pub_hpo_mean", -1)):
        w = sum((v[a] - v[b]) * sign > 0 for v in V); l = sum((v[a] - v[b]) * sign < 0 for v in V); r = signtest(w, l)
        J[name] = dict(r, median_diff=float(np.median([v[a] - v[b] for v in V])), verdict="CONFIRMED" if r["p_one_sided"] < 0.05 else "FALSIFIED")
    J["per_mirna"] = per
    json.dump(J, open(out, "w"), indent=1)
    print(json.dumps({k: J[k] for k in J if k not in ("per_mirna", "panel")}))


if __name__ == "__main__":
    main(*sys.argv[1:])
