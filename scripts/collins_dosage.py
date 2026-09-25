"""PRE-REGISTERED test of the gnomAD finding with Collins et al. 2022 rCNV dosage scores (tool 32).

Label: pHaplo and pTriplo (Zenodo record 6347673, Collins_rCNV_2022.dosage_sensitivity_scores),
probabilities of haploinsufficiency / triplosensitivity estimated from ~950k rare CNVs.
Caveat fixed in advance: the Collins model uses gene features that include gnomAD
constraint, so this label is only partly independent of LOEUF.
Committed BEFORE running.
Panel: 13 discovery + 30 held-out miRNAs (as clingen_dosage.py); CLIP universe restricted
to genes with a score; K=200. Comparators (3 draws, seeds f"{mi}-co-<arm>"): uniform,
expression-matched, UTR-length-matched. Positive control: TargetScan conserved top-200
vs 3 uniform random UTR genes with a score (miRNAs with >= 20 conserved-site genes).
Statistic: median pHaplo of the set (primary); median pTriplo (secondary).
Gates: G1 positive control higher median pHaplo than random, sign test p<0.05;
       G2 pooled universe coverage >= 80%.
Hypotheses (CNN lower), one-sided sign tests over 43 miRNAs, ties excluded, alpha 0.05:
 D1 pHaplo CNN < uniform; D2 < expression-matched; D3 < UTR-length-matched.
 T3 (secondary) pTriplo CNN < UTR-length-matched.
"""
import json, random, sys
import numpy as np, pandas as pd
sys.path.insert(0, "scripts")
from gprofiler_enrichment import PANEL, K, signtest
from string_exprmatched import matched


def main(out="results/collins_dosage.json"):
    sc = pd.read_csv("data/collins/Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz", sep="\t", compression="gzip")
    sc.columns = [c.lstrip("#") for c in sc.columns]; sc = sc.groupby("gene")[["pHaplo", "pTriplo"]].max()
    med = lambda gs, col="pHaplo": float(np.median(sc[col].reindex(gs).values))
    panel = PANEL + json.load(open("results/gnomad_replication.json"))["panel"]
    d0 = pd.read_csv("results/clip_rows_136.csv", usecols=["mirna", "gene", "min", "len"]).drop_duplicates(["mirna", "gene"])
    d0 = d0[d0.mirna.isin(panel)]; cov = float(d0.gene.isin(sc.index).mean()); d = d0[d0.gene.isin(sc.index)]
    expr = pd.read_csv("results/hpa_hek293_ntpm.tsv", sep="\t").groupby("Gene name")["nTPM"].max()
    utrg = sorted(set(pd.read_csv("data/human_utrs.tsv", sep="\t", usecols=[1]).iloc[:, 0]) & set(sc.index))
    ts = pd.read_csv("data/human_sites.tsv", sep="\t", usecols=["gene", "mirna", "context_pp"])
    per = {}
    for mi in panel:
        s = d[d.mirna == mi].reset_index(drop=True); U = s.gene.tolist(); cnn = s.nsmallest(K, "min").gene.tolist()
        se, sl = s.assign(x=s.gene.map(expr)), s.assign(x=s["len"].astype(float))
        ru, re_, rl, rb = (random.Random(f"{mi}-co-{a}") for a in ("unif", "expr", "len", "utr"))
        sets = {"uniform": [ru.sample(U, K) for _ in range(3)], "expr": [matched(se, cnn, re_) for _ in range(3)], "len": [matched(sl, cnn, rl) for _ in range(3)]}
        row = {"universe": len(U), "cnn": med(cnn), "cnn_tri": med(cnn, "pTriplo")}
        for k, v in sets.items():
            row[k + "_mean"] = float(np.mean([med(x) for x in v])); row[k + "_tri_mean"] = float(np.mean([med(x, "pTriplo") for x in v]))
        t = ts[(ts.mirna == mi) & ts.gene.isin(utrg)].groupby("gene").context_pp.min().nsmallest(K)
        row["ts_n"] = len(t)
        if len(t) >= 20:
            row["ts"] = med(list(t.index)); row["ts_rand_mean"] = float(np.mean([med(rb.sample(utrg, len(t))) for _ in range(3)]))
        per[mi] = row
    lo = lambda a, b: signtest(sum(v[a] < v[b] for v in per.values()), sum(v[a] > v[b] for v in per.values()))
    vd = lambda h: "CONFIRMED" if h["p_one_sided"] < 0.05 else "FALSIFIED"
    tsr = [v for v in per.values() if "ts" in v]
    g1 = signtest(sum(v["ts"] > v["ts_rand_mean"] for v in tsr), sum(v["ts"] < v["ts_rand_mean"] for v in tsr))
    J = {"tool": "Collins et al. 2022 rCNV dosage sensitivity scores (Zenodo 6347673)", "n_scored_genes": int(len(sc)), "panel": panel, "K": K,
         "gates": {"G1_ts_vs_random": g1, "G1_n": len(tsr), "G1_pass": bool(g1["p_one_sided"] < 0.05 and g1["wins"] > g1["losses"]), "G2_coverage": cov, "G2_pass": bool(cov >= 0.80)}}
    for name, a, b in (("D1_phaplo_vs_uniform", "cnn", "uniform_mean"), ("D2_phaplo_vs_expr", "cnn", "expr_mean"), ("D3_phaplo_vs_utrlen", "cnn", "len_mean"),
                       ("T3_ptriplo_vs_utrlen", "cnn_tri", "len_tri_mean")):
        h = lo(a, b); J[name] = dict(h, median_diff=float(np.median([v[a] - v[b] for v in per.values()])), verdict=vd(h))
    J["per_mirna"] = per
    json.dump(J, open(out, "w"), indent=1)
    print(json.dumps({k: J[k] for k in ["gates", "D1_phaplo_vs_uniform", "D2_phaplo_vs_expr", "D3_phaplo_vs_utrlen", "T3_ptriplo_vs_utrlen"]}))


if __name__ == "__main__":
    main(*sys.argv[1:])
