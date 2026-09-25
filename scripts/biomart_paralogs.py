"""PRE-REGISTERED paralog-buffering test with Ensembl BioMart (tool 40).

Data: Ensembl BioMart martservice query (hsapiens_gene_ensembl, protein_coding; attributes
ensembl_gene_id, external_gene_name, hsapiens_paralog_ensembl_gene) saved as
data/biomart/paralogs.tsv (counted once). Paralog count per symbol = number of distinct paralog
gene IDs. Question: genes with paralogs are buffered against dosage change; do CNN candidates
carry more paralogs, and does the LOEUF depletion survive paralog matching? Committed BEFORE running.
Panel: 13 discovery + 30 held-out miRNAs; CLIP universe restricted to genes with a BioMart
entry and a LOEUF; K=200. Comparators (3 draws, seeds f"{mi}-bm-<arm>"): uniform,
expression-matched, paralog-count-matched (decile matcher on log2(1+paralogs)).
Positive control (direction not assumed): TargetScan conserved top-200 vs 3 uniform random UTR
genes, two-sided sign test on median paralog count.
Gates: G1 two-sided p<0.05; G2 >= 90% of pooled universe symbols found in BioMart.
Hypotheses (one-sided sign tests over 43 miRNAs, ties excluded, alpha 0.05):
 B1 CNN median log2(1+paralogs) > uniform; B2 > expression-matched.
 B3 CNN median LOEUF > paralog-matched mean (the dosage depletion is not paralog buffering).
"""
import json, random, sys
import numpy as np, pandas as pd
sys.path.insert(0, "scripts")
from gprofiler_enrichment import PANEL, K, signtest
from string_exprmatched import matched


def paralog_counts(path="data/biomart/paralogs.tsv"):
    b = pd.read_csv(path, sep="\t", dtype=str).dropna(subset=["Gene name"])
    return b.groupby("Gene name")["Human paralogue gene stable ID"].nunique()


def main(out="results/biomart_paralogs.json"):
    pc = paralog_counts(); lpc = np.log2(1 + pc)
    lo = pd.read_csv("data/gnomad/lof_v211.txt.bgz", sep="\t", compression="gzip", usecols=["gene", "oe_lof_upper"]).dropna().groupby("gene").oe_lof_upper.min()
    panel = PANEL + json.load(open("results/gnomad_replication.json"))["panel"]
    d0 = pd.read_csv("results/clip_rows_136.csv", usecols=["mirna", "gene", "min"]).drop_duplicates(["mirna", "gene"])
    d0 = d0[d0.mirna.isin(panel)]; cov = float(pd.Series(d0.gene.unique()).isin(pc.index).mean())
    d = d0[d0.gene.isin(pc.index) & d0.gene.isin(lo.index)]
    expr = pd.read_csv("results/hpa_hek293_ntpm.tsv", sep="\t").groupby("Gene name")["nTPM"].max()
    utrg = sorted(set(pd.read_csv("data/human_utrs.tsv", sep="\t", usecols=[1]).iloc[:, 0]) & set(pc.index))
    ts = pd.read_csv("data/human_sites.tsv", sep="\t", usecols=["gene", "mirna", "context_pp"])
    mp = lambda gs: float(np.median(lpc.reindex(gs).values)); ml = lambda gs: float(np.median(lo.reindex(gs).values))
    per = {}
    for mi in panel:
        s = d[d.mirna == mi].reset_index(drop=True); U = s.gene.tolist(); cnn = s.nsmallest(K, "min").gene.tolist()
        se, sp = s.assign(x=s.gene.map(expr)), s.assign(x=s.gene.map(lpc))
        ru, re_, rp, rb = (random.Random(f"{mi}-bm-{a}") for a in ("unif", "expr", "par", "utr"))
        P = [matched(sp, cnn, rp) for _ in range(3)]
        row = {"universe": len(U), "cnn": mp(cnn), "uniform_mean": float(np.mean([mp(ru.sample(U, K)) for _ in range(3)])),
               "expr_mean": float(np.mean([mp(matched(se, cnn, re_)) for _ in range(3)])), "par_matched_par_mean": float(np.mean([mp(x) for x in P])),
               "cnn_loeuf": ml(cnn), "par_loeuf_mean": float(np.mean([ml(x) for x in P]))}
        t = ts[(ts.mirna == mi) & ts.gene.isin(utrg)].groupby("gene").context_pp.min().nsmallest(K); row["ts_n"] = len(t)
        if len(t) >= 20:
            row["ts"] = mp(list(t.index)); row["ts_rand_mean"] = float(np.mean([mp(rb.sample(utrg, len(t))) for _ in range(3)]))
        per[mi] = row
    V = list(per.values()); tsr = [v for v in V if "ts" in v]
    w, l = sum(v["ts"] > v["ts_rand_mean"] for v in tsr), sum(v["ts"] < v["ts_rand_mean"] for v in tsr)
    g1 = dict(signtest(w, l), p_two_sided=min(1.0, 2 * signtest(max(w, l), min(w, l))["p_one_sided"]))
    J = {"tool": "Ensembl BioMart martservice (hsapiens paralogues)", "n_symbols": int(len(pc)), "panel": panel, "K": K,
         "gates": {"G1_ts_vs_random": g1, "G1_n": len(tsr), "G1_pass": bool(g1["p_two_sided"] < 0.05), "G2_coverage": cov, "G2_pass": bool(cov >= 0.90)}}
    vd = lambda h: "CONFIRMED" if h["p_one_sided"] < 0.05 else "FALSIFIED"
    for name, a, b in (("B1_more_paralogs_vs_uniform", "cnn", "uniform_mean"), ("B2_vs_expr", "cnn", "expr_mean"), ("B3_loeuf_vs_paralog_matched", "cnn_loeuf", "par_loeuf_mean")):
        r = signtest(sum(v[a] > v[b] for v in V), sum(v[a] < v[b] for v in V))
        J[name] = dict(r, median_diff=float(np.median([v[a] - v[b] for v in V])), verdict=vd(r))
    J["per_mirna"] = per
    json.dump(J, open(out, "w"), indent=1)
    print(json.dumps({k: J[k] for k in J if k not in ("per_mirna", "panel")}))


if __name__ == "__main__":
    main(*sys.argv[1:])
