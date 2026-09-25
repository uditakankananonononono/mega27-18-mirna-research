"""PRE-REGISTERED independent test of the gnomAD finding with ClinGen dosage curation (tool 31).

Named finding under test (results/gnomad_constraint.json + gnomad_replication.json):
CNN top-200 non-conserved candidates are depleted of dosage-sensitive genes beyond
expression and UTR length. ClinGen (expert curation, independent of gnomAD) gives an
orthogonal label: haploinsufficiency score 3 ("sufficient evidence") = HI gene.
Committed BEFORE running.
Panel: the 13 discovery miRNAs + the 30 held-out miRNAs of gnomad_replication.json (43).
Universe: CLIP non-conserved seed-match genes (results/clip_rows_136.csv). K=200.
Comparators (3 draws, seeds f"{mi}-cg-<arm>"): uniform; expression-matched; UTR-length-matched.
Positive control: TargetScan conserved top-200 (context++) vs 3 uniform random UTR genes,
  for miRNAs with >= 20 conserved-site genes.
Statistics: per-miRNA count of HI genes in the set; pooled counts over miRNAs.
Gates: G1 positive control pooled HI count > pooled random mean AND sign test p<0.05
  (ties excluded); G2 >= 150 HI genes present in the pooled CLIP universe.
Hypotheses (CNN depleted), one-sided sign tests over miRNAs (ties excluded), alpha 0.05:
 C1 CNN HI count < uniform mean; C2 < expression-matched mean; C3 < UTR-length-matched mean.
Also reported: pooled CNN/comparator count ratios. Low per-miRNA counts mean many ties;
if fewer than 15 non-tied miRNAs remain for a hypothesis it is reported INCONCLUSIVE.
"""
import json, random, sys
import numpy as np, pandas as pd
sys.path.insert(0, "scripts")
from gprofiler_enrichment import PANEL, K, signtest
from string_exprmatched import matched


def main(out="results/clingen_dosage.json"):
    cg = pd.read_csv("data/clingen/clingen_gene_curation_GRCh38.tsv", sep="\t", skiprows=5)
    cg.columns = [c.lstrip("#") for c in cg.columns]
    hi = set(cg[cg["Haploinsufficiency Score"].astype(str) == "3"]["Gene Symbol"])
    panel = PANEL + json.load(open("results/gnomad_replication.json"))["panel"]
    d = pd.read_csv("results/clip_rows_136.csv", usecols=["mirna", "gene", "min", "len"]).drop_duplicates(["mirna", "gene"])
    d = d[d.mirna.isin(panel)]
    expr = pd.read_csv("results/hpa_hek293_ntpm.tsv", sep="\t").groupby("Gene name")["nTPM"].max()
    utrg = sorted(set(pd.read_csv("data/human_utrs.tsv", sep="\t", usecols=[1]).iloc[:, 0]))
    ts = pd.read_csv("data/human_sites.tsv", sep="\t", usecols=["gene", "mirna", "context_pp"])
    c = lambda gs: int(sum(g in hi for g in gs))
    per = {}
    for mi in panel:
        s = d[d.mirna == mi].reset_index(drop=True); U = s.gene.tolist()
        cnn = s.nsmallest(K, "min").gene.tolist()
        se, sl = s.assign(x=s.gene.map(expr)), s.assign(x=s["len"].astype(float))
        ru, re_, rl, rb = (random.Random(f"{mi}-cg-{a}") for a in ("unif", "expr", "len", "utr"))
        row = {"universe": len(U), "universe_hi": c(U), "cnn": c(cnn), "uniform": [c(ru.sample(U, K)) for _ in range(3)],
               "expr": [c(matched(se, cnn, re_)) for _ in range(3)], "len": [c(matched(sl, cnn, rl)) for _ in range(3)]}
        t = ts[(ts.mirna == mi) & ts.gene.isin(utrg)].groupby("gene").context_pp.min().nsmallest(K)
        row["ts_n"] = len(t)
        if len(t) >= 20:
            row["ts"] = c(list(t.index)); row["ts_rand"] = [c(rb.sample(utrg, len(t))) for _ in range(3)]
        for k in ("uniform", "expr", "len", "ts_rand"):
            if k in row: row[k + "_mean"] = float(np.mean(row[k]))
        per[mi] = row
        print(mi, len(U), row["universe_hi"], row["cnn"], row["uniform_mean"], row["expr_mean"], row["len_mean"], row.get("ts"), row.get("ts_rand_mean"), flush=True)
    def st(a, b, rows):
        return signtest(sum(v[a] < v[b] for v in rows), sum(v[a] > v[b] for v in rows))
    tsr = [v for v in per.values() if "ts" in v]
    g1 = signtest(sum(v["ts"] > v["ts_rand_mean"] for v in tsr), sum(v["ts"] < v["ts_rand_mean"] for v in tsr))
    g1_pool = (sum(v["ts"] for v in tsr), sum(v["ts_rand_mean"] for v in tsr))
    pooled_u = len(set(d.gene) & hi)
    J = {"tool": "ClinGen Dosage Sensitivity curation (ftp.clinicalgenome.org)", "n_hi_genes": len(hi), "panel": panel, "K": K,
         "gates": {"G1_ts_vs_random": g1, "G1_n": len(tsr), "G1_pooled_ts_vs_random": g1_pool,
                   "G1_pass": bool(g1_pool[0] > g1_pool[1] and g1["p_one_sided"] < 0.05), "G2_hi_in_universe": pooled_u, "G2_pass": bool(pooled_u >= 150)}}
    for name, k in (("C1_cnn_vs_uniform", "uniform_mean"), ("C2_cnn_vs_expr", "expr_mean"), ("C3_cnn_vs_utrlen", "len_mean")):
        h = st("cnn", k, per.values()); nt = h["wins"] + h["losses"]
        pool = (sum(v["cnn"] for v in per.values()), sum(v[k] for v in per.values()))
        J[name] = dict(h, non_tied=nt, pooled_cnn=pool[0], pooled_comparator=pool[1], pooled_ratio=pool[0] / pool[1] if pool[1] else None,
                       verdict="INCONCLUSIVE" if nt < 15 else ("CONFIRMED" if h["p_one_sided"] < 0.05 else "FALSIFIED"))
    J["per_mirna"] = per
    json.dump(J, open(out, "w"), indent=1)
    print(json.dumps({k: J[k] for k in ["n_hi_genes", "gates", "C1_cnn_vs_uniform", "C2_cnn_vs_expr", "C3_cnn_vs_utrlen"]}))


if __name__ == "__main__":
    main(*sys.argv[1:])
