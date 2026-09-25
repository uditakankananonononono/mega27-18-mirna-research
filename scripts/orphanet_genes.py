"""PRE-REGISTERED Orphanet rare-disease gene test (tool 42).

Label: Orphadata en_product6.xml (Orphanet disorder-gene associations; bulk file, counted once).
A gene is a disease gene if it has an Assessed (status id 17991) association whose type name starts
with "Disease-causing germline mutation(s)"; it is a LoF disease gene if that type is
"Disease-causing germline mutation(s) (loss of function) in". Orphanet curation is manual and
independent of gnomAD/HPO pipelines; the explicit LoF mechanism is a direct haploinsufficiency-type label.
Committed BEFORE running. Panel: 13 discovery + 30 held-out miRNAs (43); K=200; statistic =
fraction of the top-200 that are (LoF) disease genes. Genes absent from Orphanet count as 0.
Comparators (3 draws, seeds f"{mi}-or-<arm>"): expression-matched, UTR-length-matched,
publication-matched (NCBI PubMed deciles; unmapped symbols dropped first).
Gates: G1 positive control: TargetScan conserved top-200 disease-gene fraction > 3 uniform random
UTR genes, one-sided sign test p<0.05. G2 parse sanity: >= 3000 disease genes and >= 800 LoF genes.
Hypotheses (CNN lower; one-sided sign tests over 43 miRNAs, ties excluded, alpha 0.05):
 O1 disease genes vs expression-matched; O2 disease genes vs publication-matched;
 O3 LoF genes vs publication-matched; O4 LoF genes vs UTR-length-matched.
"""
import json, random, sys
import xml.etree.ElementTree as ET
import numpy as np, pandas as pd
sys.path.insert(0, "scripts")
from gprofiler_enrichment import PANEL, K, signtest
from string_exprmatched import matched

LOF = "Disease-causing germline mutation(s) (loss of function) in"


def load_orphanet(path="data/orphanet/en_product6.xml"):
    dis, lof = set(), set()
    for a in ET.parse(path).getroot().iter("DisorderGeneAssociation"):
        sym = a.findtext("Gene/Symbol"); st = a.find("DisorderGeneAssociationStatus"); ty = a.findtext("DisorderGeneAssociationType/Name") or ""
        if not sym or st is None or st.get("id") != "17991" or not ty.startswith("Disease-causing germline mutation(s)"):
            continue
        dis.add(sym)
        if ty == LOF:
            lof.add(sym)
    return dis, lof


def main(out="results/orphanet_genes.json"):
    dis, lof = load_orphanet()
    gi = pd.read_csv("data/ncbi_gene/Homo_sapiens.gene_info.gz", sep="\t", compression="gzip", usecols=["GeneID", "Symbol"])
    pm = pd.read_csv("data/ncbi_gene/human_gene2pubmed_counts.tsv", sep="\t")
    lp = np.log10(1 + gi.merge(pm, on="GeneID", how="left").fillna({"n_pubmed": 0}).groupby("Symbol").n_pubmed.max())
    panel = PANEL + json.load(open("results/gnomad_replication.json"))["panel"]
    d = pd.read_csv("results/clip_rows_136.csv", usecols=["mirna", "gene", "min", "len"]).drop_duplicates(["mirna", "gene"])
    d = d[d.mirna.isin(panel)]
    expr = pd.read_csv("results/hpa_hek293_ntpm.tsv", sep="\t").groupby("Gene name")["nTPM"].max()
    utrg = sorted(set(pd.read_csv("data/human_utrs.tsv", sep="\t", usecols=[1]).iloc[:, 0]))
    ts = pd.read_csv("data/human_sites.tsv", sep="\t", usecols=["gene", "mirna", "context_pp"])
    fd = lambda gs: float(np.mean([g in dis for g in gs])); fl = lambda gs: float(np.mean([g in lof for g in gs]))
    per = {}
    for mi in panel:
        s = d[d.mirna == mi].reset_index(drop=True); cnn = s.nsmallest(K, "min").gene.tolist()
        se, sl = s.assign(x=s.gene.map(expr)), s.assign(x=s["len"].astype(float))
        sp = s[s.gene.isin(lp.index)].reset_index(drop=True); sp = sp.assign(x=sp.gene.map(lp)); cnnp = sp.nsmallest(K, "min").gene.tolist()
        re_, rl, rp, rb = (random.Random(f"{mi}-or-{a}") for a in ("expr", "len", "pub", "utr"))
        E = [matched(se, cnn, re_) for _ in range(3)]; L = [matched(sl, cnn, rl) for _ in range(3)]; P = [matched(sp, cnnp, rp) for _ in range(3)]
        row = {"universe": len(s), "cnn_dis": fd(cnn), "expr_dis": float(np.mean([fd(x) for x in E])),
               "cnnp_dis": fd(cnnp), "pub_dis": float(np.mean([fd(x) for x in P])), "cnnp_lof": fl(cnnp), "pub_lof": float(np.mean([fl(x) for x in P])),
               "cnn_lof": fl(cnn), "len_lof": float(np.mean([fl(x) for x in L]))}
        t = ts[(ts.mirna == mi) & ts.gene.isin(utrg)].groupby("gene").context_pp.min().nsmallest(K); row["ts_n"] = len(t)
        if len(t) >= 20:
            row["ts"] = fd(list(t.index)); row["ts_rand_mean"] = float(np.mean([fd(rb.sample(utrg, len(t))) for _ in range(3)]))
        per[mi] = row
    V = list(per.values()); tsr = [v for v in V if "ts" in v]
    g1 = signtest(sum(v["ts"] > v["ts_rand_mean"] for v in tsr), sum(v["ts"] < v["ts_rand_mean"] for v in tsr))
    J = {"tool": "Orphanet / Orphadata en_product6 (disorder-gene associations)", "n_disease_genes": len(dis), "n_lof_genes": len(lof), "panel": panel, "K": K,
         "gates": {"G1_ts_vs_random": g1, "G1_n": len(tsr), "G1_pass": bool(g1["p_one_sided"] < 0.05 and g1["wins"] > g1["losses"]),
                   "G2_pass": bool(len(dis) >= 3000 and len(lof) >= 800)}}
    for name, a, b in (("O1_dis_vs_expr", "cnn_dis", "expr_dis"), ("O2_dis_vs_pub", "cnnp_dis", "pub_dis"),
                       ("O3_lof_vs_pub", "cnnp_lof", "pub_lof"), ("O4_lof_vs_utrlen", "cnn_lof", "len_lof")):
        r = signtest(sum(v[a] < v[b] for v in V), sum(v[a] > v[b] for v in V))
        J[name] = dict(r, median_diff=float(np.median([v[a] - v[b] for v in V])), verdict="CONFIRMED" if r["p_one_sided"] < 0.05 else "FALSIFIED")
    J["per_mirna"] = per
    json.dump(J, open(out, "w"), indent=1)
    print(json.dumps({k: J[k] for k in J if k not in ("per_mirna", "panel")}))


if __name__ == "__main__":
    main(*sys.argv[1:])
