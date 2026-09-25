"""PRE-REGISTERED DECIPHER haploinsufficiency-prediction test (tool 38).

Label: DECIPHER HI_Predictions_Version3.bed.gz (Ni et al.; deciphergenomics.org downloads; counted
once); name field "SYMBOL|score|pct", column 5 = predicted probability of haploinsufficiency.
This model is built from genomic, evolutionary and network features (not gnomAD LOEUF), so it is a
further, partly independent dosage label. Committed BEFORE running.
Panel: 13 discovery + 30 held-out miRNAs; CLIP universe restricted to genes with a score; K=200.
Comparators (3 draws, seeds f"{mi}-dc-<arm>"): uniform, expression-matched, UTR-length-matched,
publication-matched (NCBI PubMed deciles; unmapped symbols dropped from the universe first).
Positive control: TargetScan conserved top-200 higher median score than 3 uniform random UTR genes.
Statistic: median HI score.
Gates: G1 positive control one-sided sign test p<0.05 (ties excluded); G2 coverage of pooled
CLIP universe >= 85%.
Hypotheses (CNN lower), one-sided sign tests over 43 miRNAs, ties excluded, alpha 0.05:
 H1 vs uniform; H2 vs expression-matched; H3 vs UTR-length-matched; H4 vs publication-matched.
"""
import gzip, json, random, sys
import numpy as np, pandas as pd
sys.path.insert(0, "scripts")
from gprofiler_enrichment import PANEL, K, signtest
from string_exprmatched import matched


def load_hi(path="data/decipher/HI_Predictions_Version3.bed.gz"):
    out = {}
    for ln in gzip.open(path, "rt"):
        if ln.startswith("track"):
            continue
        f = ln.rstrip("\n").split("\t"); g = f[3].split("|")[0]; out[g] = max(out.get(g, 0.0), float(f[4]))
    return pd.Series(out)


def main(out="results/decipher_hi.json"):
    hi = load_hi()
    gi = pd.read_csv("data/ncbi_gene/Homo_sapiens.gene_info.gz", sep="\t", compression="gzip", usecols=["GeneID", "Symbol"])
    pm = pd.read_csv("data/ncbi_gene/human_gene2pubmed_counts.tsv", sep="\t")
    lp = np.log10(1 + gi.merge(pm, on="GeneID", how="left").fillna({"n_pubmed": 0}).groupby("Symbol").n_pubmed.max())
    panel = PANEL + json.load(open("results/gnomad_replication.json"))["panel"]
    d0 = pd.read_csv("results/clip_rows_136.csv", usecols=["mirna", "gene", "min", "len"]).drop_duplicates(["mirna", "gene"])
    d0 = d0[d0.mirna.isin(panel)]; cov = float(pd.Series(d0.gene.unique()).isin(hi.index).mean())
    d = d0[d0.gene.isin(hi.index)]
    expr = pd.read_csv("results/hpa_hek293_ntpm.tsv", sep="\t").groupby("Gene name")["nTPM"].max()
    utrg = sorted(set(pd.read_csv("data/human_utrs.tsv", sep="\t", usecols=[1]).iloc[:, 0]) & set(hi.index))
    ts = pd.read_csv("data/human_sites.tsv", sep="\t", usecols=["gene", "mirna", "context_pp"])
    med = lambda gs: float(np.median(hi.reindex(gs).values))
    per = {}
    for mi in panel:
        s = d[d.mirna == mi].reset_index(drop=True); U = s.gene.tolist(); cnn = s.nsmallest(K, "min").gene.tolist()
        se, sl = s.assign(x=s.gene.map(expr)), s.assign(x=s["len"].astype(float))
        sp = s[s.gene.isin(lp.index)].reset_index(drop=True); sp = sp.assign(x=sp.gene.map(lp)); cnnp = sp.nsmallest(K, "min").gene.tolist()
        ru, re_, rl, rp, rb = (random.Random(f"{mi}-dc-{a}") for a in ("unif", "expr", "len", "pub", "utr"))
        row = {"universe": len(U), "cnn": med(cnn), "uniform_mean": float(np.mean([med(ru.sample(U, K)) for _ in range(3)])),
               "expr_mean": float(np.mean([med(matched(se, cnn, re_)) for _ in range(3)])), "len_mean": float(np.mean([med(matched(sl, cnn, rl)) for _ in range(3)])),
               "cnn_pubuniv": med(cnnp), "pub_mean": float(np.mean([med(matched(sp, cnnp, rp)) for _ in range(3)]))}
        t = ts[(ts.mirna == mi) & ts.gene.isin(utrg)].groupby("gene").context_pp.min().nsmallest(K); row["ts_n"] = len(t)
        if len(t) >= 20:
            row["ts"] = med(list(t.index)); row["ts_rand_mean"] = float(np.mean([med(rb.sample(utrg, len(t))) for _ in range(3)]))
        per[mi] = row
    V = list(per.values()); tsr = [v for v in V if "ts" in v]
    g1 = signtest(sum(v["ts"] > v["ts_rand_mean"] for v in tsr), sum(v["ts"] < v["ts_rand_mean"] for v in tsr))
    J = {"tool": "DECIPHER HI predictions v3 (deciphergenomics.org)", "n_scored_genes": int(len(hi)), "panel": panel, "K": K,
         "gates": {"G1_ts_vs_random": g1, "G1_n": len(tsr), "G1_pass": bool(g1["p_one_sided"] < 0.05 and g1["wins"] > g1["losses"]), "G2_coverage": cov, "G2_pass": bool(cov >= 0.85)}}
    for name, a, b in (("H1_vs_uniform", "cnn", "uniform_mean"), ("H2_vs_expr", "cnn", "expr_mean"), ("H3_vs_utrlen", "cnn", "len_mean"), ("H4_vs_pubmatched", "cnn_pubuniv", "pub_mean")):
        r = signtest(sum(v[a] < v[b] for v in V), sum(v[a] > v[b] for v in V))
        J[name] = dict(r, median_diff=float(np.median([v[a] - v[b] for v in V])), verdict="CONFIRMED" if r["p_one_sided"] < 0.05 else "FALSIFIED")
    J["per_mirna"] = per
    json.dump(J, open(out, "w"), indent=1)
    print(json.dumps({k: J[k] for k in J if k not in ("per_mirna", "panel")}))


if __name__ == "__main__":
    main(*sys.argv[1:])
