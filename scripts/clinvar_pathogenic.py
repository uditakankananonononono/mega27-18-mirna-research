"""PRE-REGISTERED ClinVar pathogenic-gene test (tool 34).

Label: ClinVar gene_specific_summary.txt (NCBI FTP, dated in its header; counted once).
Primary: gene has >= 5 alleles reported Pathogenic/Likely pathogenic ("PLP gene").
Secondary: >= 1 PLP allele. Known confound fixed in advance: PLP counts grow with coding
length and clinical testing effort; UTR-length matching does not address coding length.
Question: do CNN top-200 non-conserved candidates carry fewer PLP genes (dosage/disease
finding, tools 30-33)? Committed BEFORE running.
Panel: 13 discovery + 30 held-out miRNAs; CLIP universe; K=200. Comparators (3 draws, seeds
f"{mi}-cv-<arm>"): uniform, expression-matched, UTR-length-matched. Positive control:
TargetScan conserved top-200 vs 3 uniform random UTR genes (>= 20 conserved-site genes).
Gates: G1 positive control more PLP genes than random, sign test p<0.05 (ties excluded);
       G2 primary PLP prevalence in pooled universe between 5% and 40%.
Hypotheses (CNN fewer PLP genes, primary label), one-sided sign tests over 43 miRNAs,
ties excluded, alpha 0.05: V1 vs uniform; V2 vs expression-matched; V3 vs UTR-length-matched.
Secondary S3: >= 1 PLP label, CNN fewer than UTR-length-matched.
"""
import json, random, sys
import numpy as np, pandas as pd
sys.path.insert(0, "scripts")
from gprofiler_enrichment import PANEL, K, signtest
from string_exprmatched import matched


def main(out="results/clinvar_pathogenic.json"):
    cv = pd.read_csv("data/clinvar/gene_specific_summary.txt", sep="\t", skiprows=1)
    cv.columns = [c.lstrip("#") for c in cv.columns]
    n = pd.to_numeric(cv["Alleles_reported_Pathogenic_Likely_pathogenic"], errors="coerce").fillna(0)
    plp5, plp1 = set(cv.Symbol[n >= 5]), set(cv.Symbol[n >= 1])
    date = open("data/clinvar/gene_specific_summary.txt").readline().strip("#\n")
    panel = PANEL + json.load(open("results/gnomad_replication.json"))["panel"]
    d = pd.read_csv("results/clip_rows_136.csv", usecols=["mirna", "gene", "min", "len"]).drop_duplicates(["mirna", "gene"])
    d = d[d.mirna.isin(panel)]
    expr = pd.read_csv("results/hpa_hek293_ntpm.tsv", sep="\t").groupby("Gene name")["nTPM"].max()
    utrg = sorted(set(pd.read_csv("data/human_utrs.tsv", sep="\t", usecols=[1]).iloc[:, 0]))
    ts = pd.read_csv("data/human_sites.tsv", sep="\t", usecols=["gene", "mirna", "context_pp"])
    c = lambda gs, S=plp5: int(sum(g in S for g in gs))
    per = {}
    for mi in panel:
        s = d[d.mirna == mi].reset_index(drop=True); U = s.gene.tolist(); cnn = s.nsmallest(K, "min").gene.tolist()
        se, sl = s.assign(x=s.gene.map(expr)), s.assign(x=s["len"].astype(float))
        ru, re_, rl, rb = (random.Random(f"{mi}-cv-{a}") for a in ("unif", "expr", "len", "utr"))
        L = [matched(sl, cnn, rl) for _ in range(3)]
        row = {"universe": len(U), "cnn": c(cnn), "cnn_any": c(cnn, plp1), "uniform_mean": float(np.mean([c(ru.sample(U, K)) for _ in range(3)])),
               "expr_mean": float(np.mean([c(matched(se, cnn, re_)) for _ in range(3)])), "len_mean": float(np.mean([c(x) for x in L])),
               "len_any_mean": float(np.mean([c(x, plp1) for x in L]))}
        t = ts[(ts.mirna == mi) & ts.gene.isin(utrg)].groupby("gene").context_pp.min().nsmallest(K); row["ts_n"] = len(t)
        if len(t) >= 20:
            row["ts"] = c(list(t.index)); row["ts_rand_mean"] = float(np.mean([c(rb.sample(utrg, len(t))) for _ in range(3)]))
        per[mi] = row
    tsr = [v for v in per.values() if "ts" in v]
    g1 = signtest(sum(v["ts"] > v["ts_rand_mean"] for v in tsr), sum(v["ts"] < v["ts_rand_mean"] for v in tsr))
    prev = len(set(d.gene) & plp5) / d.gene.nunique()
    J = {"tool": "ClinVar gene_specific_summary (NCBI FTP)", "clinvar_date": date, "n_plp5": len(plp5), "n_plp1": len(plp1), "panel": panel, "K": K,
         "gates": {"G1_ts_vs_random": g1, "G1_n": len(tsr), "G1_pass": bool(g1["p_one_sided"] < 0.05 and g1["wins"] > g1["losses"]),
                   "G2_prevalence": prev, "G2_pass": bool(0.05 <= prev <= 0.40)}}
    for name, a, k in (("V1_vs_uniform", "cnn", "uniform_mean"), ("V2_vs_expr", "cnn", "expr_mean"), ("V3_vs_utrlen", "cnn", "len_mean"), ("S3_any_vs_utrlen", "cnn_any", "len_any_mean")):
        r = signtest(sum(v[a] < v[k] for v in per.values()), sum(v[a] > v[k] for v in per.values()))
        pc, pk = sum(v[a] for v in per.values()), sum(v[k] for v in per.values())
        J[name] = dict(r, pooled_cnn=pc, pooled_comparator=pk, pooled_ratio=pc / pk, verdict="CONFIRMED" if r["p_one_sided"] < 0.05 else "FALSIFIED")
    J["per_mirna"] = per
    json.dump(J, open(out, "w"), indent=1)
    print(json.dumps({k: J[k] for k in ["clinvar_date", "n_plp5", "gates", "V1_vs_uniform", "V2_vs_expr", "V3_vs_utrlen", "S3_any_vs_utrlen"]}))


if __name__ == "__main__":
    main(*sys.argv[1:])
