"""PRE-REGISTERED MGI mouse-mutant phenotype test (tool 43).

Label: MGI HMD_HumanPhenotype.rpt (informatics.jax.org reports; human symbol -> mouse orthologue ->
top-level Mammalian Phenotype systems observed in mouse mutants; bulk file, counted once).
An experimental, non-human, non-population label: does the orthologue's mutant show
mortality/aging (MP:0010768) or embryo (MP:0005380) phenotypes? Committed BEFORE running.
Universe: CLIP genes with an HMD row carrying >= 1 MP term (genes never phenotyped are dropped, so
"no annotation" is not read as "no phenotype"). Panel: 13 discovery + 30 held-out miRNAs (43); K=200.
Statistic: fraction of the top-200 with the phenotype. Comparators (3 draws, seeds f"{mi}-mg-<arm>"):
expression-matched, UTR-length-matched, publication-matched (NCBI PubMed deciles).
Gates: G1 positive control: TargetScan conserved top-200 (phenotyped UTR genes) mortality fraction >
3 uniform random phenotyped UTR genes, one-sided sign test p<0.05. G2 >= 80% of pooled CLIP-universe
symbols have an HMD row.
Hypotheses (CNN lower; one-sided sign tests over 43 miRNAs, ties excluded, alpha 0.05):
 M1 mortality vs expression-matched; M2 mortality vs publication-matched;
 M3 mortality vs UTR-length-matched; M4 embryo phenotype vs publication-matched.
"""
import json, random, sys
import numpy as np, pandas as pd
sys.path.insert(0, "scripts")
from gprofiler_enrichment import PANEL, K, signtest
from string_exprmatched import matched

MORT, EMB = "MP:0010768", "MP:0005380"


def load_hmd(path="data/mgi/HMD_HumanPhenotype.rpt"):
    """Return (all symbols with a row, {symbol: set(MP)} for phenotyped symbols)."""
    rows, ph = set(), {}
    for ln in open(path):
        f = ln.rstrip("\n").split("\t")
        if not f[0]:
            continue
        rows.add(f[0]); mp = {x.strip() for x in (f[4] if len(f) > 4 else "").split(",") if x.strip()}
        if mp:
            ph.setdefault(f[0], set()).update(mp)
    return rows, ph


def main(out="results/mgi_mouse_ko.json"):
    rows, ph = load_hmd(); mort = {g for g, s in ph.items() if MORT in s}; emb = {g for g, s in ph.items() if EMB in s}
    gi = pd.read_csv("data/ncbi_gene/Homo_sapiens.gene_info.gz", sep="\t", compression="gzip", usecols=["GeneID", "Symbol"])
    pm = pd.read_csv("data/ncbi_gene/human_gene2pubmed_counts.tsv", sep="\t")
    lp = np.log10(1 + gi.merge(pm, on="GeneID", how="left").fillna({"n_pubmed": 0}).groupby("Symbol").n_pubmed.max())
    panel = PANEL + json.load(open("results/gnomad_replication.json"))["panel"]
    d0 = pd.read_csv("results/clip_rows_136.csv", usecols=["mirna", "gene", "min", "len"]).drop_duplicates(["mirna", "gene"])
    d0 = d0[d0.mirna.isin(panel)]; cov = float(pd.Series(d0.gene.unique()).isin(rows).mean())
    d = d0[d0.gene.isin(ph.keys())]
    expr = pd.read_csv("results/hpa_hek293_ntpm.tsv", sep="\t").groupby("Gene name")["nTPM"].max()
    utrg = sorted(set(pd.read_csv("data/human_utrs.tsv", sep="\t", usecols=[1]).iloc[:, 0]) & set(ph))
    ts = pd.read_csv("data/human_sites.tsv", sep="\t", usecols=["gene", "mirna", "context_pp"])
    fm = lambda gs: float(np.mean([g in mort for g in gs])); fe = lambda gs: float(np.mean([g in emb for g in gs]))
    per = {}
    for mi in panel:
        s = d[d.mirna == mi].reset_index(drop=True); cnn = s.nsmallest(K, "min").gene.tolist()
        se, sl = s.assign(x=s.gene.map(expr)), s.assign(x=s["len"].astype(float))
        sp = s[s.gene.isin(lp.index)].reset_index(drop=True); sp = sp.assign(x=sp.gene.map(lp)); cnnp = sp.nsmallest(K, "min").gene.tolist()
        re_, rl, rp, rb = (random.Random(f"{mi}-mg-{a}") for a in ("expr", "len", "pub", "utr"))
        P = [matched(sp, cnnp, rp) for _ in range(3)]
        row = {"universe": len(s), "cnn_mort": fm(cnn), "expr_mort": float(np.mean([fm(matched(se, cnn, re_)) for _ in range(3)])),
               "len_mort": float(np.mean([fm(matched(sl, cnn, rl)) for _ in range(3)])),
               "cnnp_mort": fm(cnnp), "pub_mort": float(np.mean([fm(x) for x in P])), "cnnp_emb": fe(cnnp), "pub_emb": float(np.mean([fe(x) for x in P]))}
        t = ts[(ts.mirna == mi) & ts.gene.isin(utrg)].groupby("gene").context_pp.min().nsmallest(K); row["ts_n"] = len(t)
        if len(t) >= 20:
            row["ts"] = fm(list(t.index)); row["ts_rand_mean"] = float(np.mean([fm(rb.sample(utrg, len(t))) for _ in range(3)]))
        per[mi] = row
    V = list(per.values()); tsr = [v for v in V if "ts" in v]
    g1 = signtest(sum(v["ts"] > v["ts_rand_mean"] for v in tsr), sum(v["ts"] < v["ts_rand_mean"] for v in tsr))
    J = {"tool": "MGI HMD_HumanPhenotype (mouse mutant phenotypes)", "n_rows": len(rows), "n_phenotyped": len(ph), "n_mortality": len(mort), "n_embryo": len(emb),
         "panel": panel, "K": K, "gates": {"G1_ts_vs_random": g1, "G1_n": len(tsr), "G1_pass": bool(g1["p_one_sided"] < 0.05 and g1["wins"] > g1["losses"]),
                                            "G2_coverage": cov, "G2_pass": bool(cov >= 0.80)}}
    for name, a, b in (("M1_mort_vs_expr", "cnn_mort", "expr_mort"), ("M2_mort_vs_pub", "cnnp_mort", "pub_mort"),
                       ("M3_mort_vs_utrlen", "cnn_mort", "len_mort"), ("M4_emb_vs_pub", "cnnp_emb", "pub_emb")):
        r = signtest(sum(v[a] < v[b] for v in V), sum(v[a] > v[b] for v in V))
        J[name] = dict(r, median_diff=float(np.median([v[a] - v[b] for v in V])), verdict="CONFIRMED" if r["p_one_sided"] < 0.05 else "FALSIFIED")
    J["per_mirna"] = per
    json.dump(J, open(out, "w"), indent=1)
    print(json.dumps({k: J[k] for k in J if k not in ("per_mirna", "panel")}))


if __name__ == "__main__":
    main(*sys.argv[1:])
