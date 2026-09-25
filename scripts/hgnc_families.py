"""PRE-REGISTERED gene-family confound test with HGNC (tool 41).

Data: HGNC complete set (approved symbols, gene groups, locus_group, location), saved as
data/hgnc/hgnc_complete_set.txt (bulk file, counted once). Question: is the CNN dosage
depletion a gene-family artefact, i.e. are CNN candidates enriched in large, dosage-tolerant
families (C2H2 zinc fingers, KRAB-ZNFs), and does the LOEUF depletion survive their removal?
Committed BEFORE running. Panel: 13 discovery + 30 held-out miRNAs (43); K=200; CLIP universe
restricted to genes with an HGNC protein-coding approved symbol and a gnomAD LOEUF.
Comparators: expression-matched (3 draws, seeds f"{mi}-hg-<arm>").
Gates: G1 >= 90% of pooled CLIP-universe symbols are HGNC protein-coding symbols.
       G2 (parse sanity, known biology: ZNF clusters on 19q13) fraction of C2H2-ZNF genes on
       chr19 / fraction elsewhere >= 3.
Hypotheses (one-sided sign tests over 43 miRNAs, ties excluded, alpha 0.05):
 F1 CNN fraction C2H2-ZNF > expression-matched mean.
 F2 CNN median log2(1+largest HGNC group size) > expression-matched mean.
 F3 after removing all C2H2-ZNF genes from the universe, CNN median LOEUF > expression-matched mean.
"""
import json, random, sys
import numpy as np, pandas as pd
sys.path.insert(0, "scripts")
from gprofiler_enrichment import PANEL, K, signtest
from string_exprmatched import matched


def load_hgnc(path="data/hgnc/hgnc_complete_set.txt"):
    h = pd.read_csv(path, sep="\t", dtype=str, usecols=["symbol", "locus_group", "status", "location", "gene_group", "gene_group_id"])
    h = h[(h.locus_group == "protein-coding gene") & (h.status == "Approved")].copy()
    h["znf"] = h.gene_group.fillna("").str.contains("Zinc fingers C2H2-type", regex=False)
    ex = h.dropna(subset=["gene_group_id"]).assign(g=lambda x: x.gene_group_id.str.split("|")).explode("g")
    size = ex.g.value_counts()
    ex["n"] = ex.g.map(size)
    h["fam"] = np.log2(1 + h.symbol.map(ex.groupby("symbol").n.max()).fillna(0))
    h["chr19"] = h.location.fillna("").str.match(r"^19[pq]")
    return h.set_index("symbol")


def g2_ratio(h):
    a = h.loc[h.chr19, "znf"].mean(); b = h.loc[~h.chr19, "znf"].mean()
    return float(a / b) if b > 0 else float("inf")


def main(out="results/hgnc_families.json"):
    h = load_hgnc()
    lo = pd.read_csv("data/gnomad/lof_v211.txt.bgz", sep="\t", compression="gzip", usecols=["gene", "oe_lof_upper"]).dropna().groupby("gene").oe_lof_upper.min()
    panel = PANEL + json.load(open("results/gnomad_replication.json"))["panel"]
    d0 = pd.read_csv("results/clip_rows_136.csv", usecols=["mirna", "gene", "min"]).drop_duplicates(["mirna", "gene"])
    d0 = d0[d0.mirna.isin(panel)]; cov = float(pd.Series(d0.gene.unique()).isin(h.index).mean())
    d = d0[d0.gene.isin(h.index) & d0.gene.isin(lo.index)]
    expr = pd.read_csv("results/hpa_hek293_ntpm.tsv", sep="\t").groupby("Gene name")["nTPM"].max()
    fz = lambda gs: float(h.znf.reindex(gs).mean()); ff = lambda gs: float(np.median(h.fam.reindex(gs).values))
    ml = lambda gs: float(np.median(lo.reindex(gs).values))
    per = {}
    for mi in panel:
        s = d[d.mirna == mi].reset_index(drop=True); cnn = s.nsmallest(K, "min").gene.tolist()
        se = s.assign(x=s.gene.map(expr)); r1, r3 = random.Random(f"{mi}-hg-expr"), random.Random(f"{mi}-hg-noznf")
        E = [matched(se, cnn, r1) for _ in range(3)]
        s3 = s[~s.gene.map(h.znf)].reset_index(drop=True); c3 = s3.nsmallest(K, "min").gene.tolist()
        s3e = s3.assign(x=s3.gene.map(expr))
        per[mi] = {"universe": len(s), "cnn_znf": fz(cnn), "expr_znf": float(np.mean([fz(x) for x in E])),
                   "cnn_fam": ff(cnn), "expr_fam": float(np.mean([ff(x) for x in E])),
                   "noznf_universe": len(s3), "noznf_cnn_loeuf": ml(c3),
                   "noznf_expr_loeuf": float(np.mean([ml(matched(s3e, c3, r3)) for _ in range(3)]))}
    V = list(per.values()); g2 = g2_ratio(h)
    J = {"tool": "HGNC complete set (gene groups)", "n_protein_coding": int(len(h)), "n_znf": int(h.znf.sum()), "panel": panel, "K": K,
         "gates": {"G1_coverage": cov, "G1_pass": bool(cov >= 0.90), "G2_chr19_znf_ratio": g2, "G2_pass": bool(g2 >= 3)}}
    vd = lambda r: "CONFIRMED" if r["p_one_sided"] < 0.05 else "FALSIFIED"
    for name, a, b in (("F1_znf_vs_expr", "cnn_znf", "expr_znf"), ("F2_family_size_vs_expr", "cnn_fam", "expr_fam"),
                       ("F3_loeuf_noznf_vs_expr", "noznf_cnn_loeuf", "noznf_expr_loeuf")):
        r = signtest(sum(v[a] > v[b] for v in V), sum(v[a] < v[b] for v in V))
        J[name] = dict(r, median_diff=float(np.median([v[a] - v[b] for v in V])), verdict=vd(r))
    J["per_mirna"] = per
    json.dump(J, open(out, "w"), indent=1)
    print(json.dumps({k: J[k] for k in J if k not in ("per_mirna", "panel")}))


if __name__ == "__main__":
    main(*sys.argv[1:])
