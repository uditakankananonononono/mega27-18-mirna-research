"""PRE-REGISTERED Mendelian-disease-gene test (HPO annotations, tool 33).

Label: genes with a MENDELIAN association in the HPO release file genes_to_disease.txt
(github.com/obophenotype/human-phenotype-ontology, latest release; counted once).
Question: do CNN top-200 non-conserved candidates carry fewer Mendelian disease genes,
in line with the dosage finding (tools 30-32)? Committed BEFORE running.
Panel: 13 discovery + 30 held-out miRNAs; CLIP universe; K=200. Comparators (3 draws,
seeds f"{mi}-hpo-<arm>"): uniform, expression-matched, UTR-length-matched.
Positive control: TargetScan conserved top-200 vs 3 uniform random UTR genes
(miRNAs with >= 20 conserved-site genes).
Statistic: count of Mendelian genes per set.
Gates: G1 positive control more disease genes than random, sign test p<0.05 (ties excluded);
       G2 disease-gene prevalence in pooled universe between 10% and 40%.
Hypotheses (CNN fewer), one-sided sign tests over 43 miRNAs, ties excluded, alpha 0.05:
 M1 vs uniform; M2 vs expression-matched; M3 vs UTR-length-matched.
"""
import json, random, sys
import numpy as np, pandas as pd
sys.path.insert(0, "scripts")
from gprofiler_enrichment import PANEL, K, signtest
from string_exprmatched import matched


def main(out="results/hpo_disease_genes.json"):
    h = pd.read_csv("data/hpo/genes_to_disease.txt", sep="\t")
    dis = set(h[h.association_type == "MENDELIAN"].gene_symbol)
    panel = PANEL + json.load(open("results/gnomad_replication.json"))["panel"]
    d = pd.read_csv("results/clip_rows_136.csv", usecols=["mirna", "gene", "min", "len"]).drop_duplicates(["mirna", "gene"])
    d = d[d.mirna.isin(panel)]
    expr = pd.read_csv("results/hpa_hek293_ntpm.tsv", sep="\t").groupby("Gene name")["nTPM"].max()
    utrg = sorted(set(pd.read_csv("data/human_utrs.tsv", sep="\t", usecols=[1]).iloc[:, 0]))
    ts = pd.read_csv("data/human_sites.tsv", sep="\t", usecols=["gene", "mirna", "context_pp"])
    c = lambda gs: int(sum(g in dis for g in gs))
    per = {}
    for mi in panel:
        s = d[d.mirna == mi].reset_index(drop=True); U = s.gene.tolist(); cnn = s.nsmallest(K, "min").gene.tolist()
        se, sl = s.assign(x=s.gene.map(expr)), s.assign(x=s["len"].astype(float))
        ru, re_, rl, rb = (random.Random(f"{mi}-hpo-{a}") for a in ("unif", "expr", "len", "utr"))
        row = {"universe": len(U), "universe_dis": c(U), "cnn": c(cnn), "uniform_mean": float(np.mean([c(ru.sample(U, K)) for _ in range(3)])),
               "expr_mean": float(np.mean([c(matched(se, cnn, re_)) for _ in range(3)])), "len_mean": float(np.mean([c(matched(sl, cnn, rl)) for _ in range(3)]))}
        t = ts[(ts.mirna == mi) & ts.gene.isin(utrg)].groupby("gene").context_pp.min().nsmallest(K); row["ts_n"] = len(t)
        if len(t) >= 20:
            row["ts"] = c(list(t.index)); row["ts_rand_mean"] = float(np.mean([c(rb.sample(utrg, len(t))) for _ in range(3)]))
        per[mi] = row
    tsr = [v for v in per.values() if "ts" in v]
    g1 = signtest(sum(v["ts"] > v["ts_rand_mean"] for v in tsr), sum(v["ts"] < v["ts_rand_mean"] for v in tsr))
    prev = len(set(d.gene) & dis) / d.gene.nunique()
    J = {"tool": "Human Phenotype Ontology gene-disease annotations (genes_to_disease.txt)", "n_mendelian_genes": len(dis), "panel": panel, "K": K,
         "gates": {"G1_ts_vs_random": g1, "G1_n": len(tsr), "G1_pass": bool(g1["p_one_sided"] < 0.05 and g1["wins"] > g1["losses"]),
                   "G2_prevalence": prev, "G2_pass": bool(0.10 <= prev <= 0.40)}}
    for name, k in (("M1_vs_uniform", "uniform_mean"), ("M2_vs_expr", "expr_mean"), ("M3_vs_utrlen", "len_mean")):
        r = signtest(sum(v["cnn"] < v[k] for v in per.values()), sum(v["cnn"] > v[k] for v in per.values()))
        pc, pk = sum(v["cnn"] for v in per.values()), sum(v[k] for v in per.values())
        J[name] = dict(r, pooled_cnn=pc, pooled_comparator=pk, pooled_ratio=pc / pk, verdict="CONFIRMED" if r["p_one_sided"] < 0.05 else "FALSIFIED")
    J["per_mirna"] = per
    json.dump(J, open(out, "w"), indent=1)
    print(json.dumps({k: J[k] for k in ["n_mendelian_genes", "gates", "M1_vs_uniform", "M2_vs_expr", "M3_vs_utrlen"]}))


if __name__ == "__main__":
    main(*sys.argv[1:])
