"""PRE-REGISTERED held-out replication of the gnomAD constraint reversal (tool 30 follow-up).

Discovery (post hoc, results/gnomad_constraint.json): on the 13-miRNA panel the CNN
top-200 had a higher median LOEUF (less constrained) than uniform, site-count and
expression-matched sets in 13/13 miRNAs. Candidate explanation: 3' UTR length.
This script is committed BEFORE it is run.
Panel: miRNAs in results/clip_rows_136.csv with >= 1000 LOEUF-mapped universe genes,
excluding the 13 discovery miRNAs; sorted, shuffled with random.Random("gnomad-rep-2026"),
first 30. K=200. Universe = CLIP non-conserved seed-match genes with a LOEUF.
Comparators (3 draws each, seeds f"{mi}-rep-<arm>"): uniform random; expression-matched
(HEK293 nTPM deciles, string_exprmatched.matched); UTR-length-matched (same decile
matcher on the gene's longest 3' UTR length, column len).
Statistic: median LOEUF; for each miRNA compare CNN with the mean of the 3 draws.
Gates: G1 >= 25 eligible miRNAs; G2 pooled LOEUF coverage >= 85%.
Hypotheses (one-sided sign tests, alpha 0.05, ties excluded):
 R1 CNN median LOEUF > uniform (replicates the reversal).
 R2 CNN > expression-matched.
 R3 CNN > UTR-length-matched. If R1 holds and R3 fails, the reversal is attributed to
    UTR length (a length artefact), not to the CNN's sequence preference.
"""
import json, random, sys
import numpy as np, pandas as pd
sys.path.insert(0, "scripts")
from gprofiler_enrichment import PANEL, K, signtest
from string_exprmatched import matched


def main(out="results/gnomad_replication.json"):
    g = pd.read_csv("data/gnomad/lof_v211.txt.bgz", sep="\t", compression="gzip", usecols=["gene", "oe_lof_upper"]).dropna()
    lo = g.groupby("gene").oe_lof_upper.min()
    med = lambda gs: float(np.median(lo.reindex(gs).values))
    d0 = pd.read_csv("results/clip_rows_136.csv", usecols=["mirna", "gene", "min", "len"]).drop_duplicates(["mirna", "gene"])
    d0 = d0[~d0.mirna.isin(PANEL)]
    cov = float(d0.gene.isin(lo.index).mean())
    d = d0[d0.gene.isin(lo.index)]
    n = d.mirna.value_counts()
    cand = sorted(n[n >= 1000].index); random.Random("gnomad-rep-2026").shuffle(cand); panel = cand[:30]
    expr = pd.read_csv("results/hpa_hek293_ntpm.tsv", sep="\t").groupby("Gene name")["nTPM"].max()
    per = {}
    for mi in panel:
        s = d[d.mirna == mi].reset_index(drop=True); U = s.gene.tolist()
        cnn = s.nsmallest(K, "min").gene.tolist()
        se, sl = s.assign(x=s.gene.map(expr)), s.assign(x=s["len"].astype(float))
        ru, re_, rl = (random.Random(f"{mi}-rep-{a}") for a in ("unif", "expr", "len"))
        row = {"universe": len(U), "cnn": med(cnn), "cnn_median_len": float(s[s.gene.isin(cnn)]["len"].median()), "universe_median_len": float(s["len"].median()),
               "uniform": [med(ru.sample(U, K)) for _ in range(3)], "expr": [med(matched(se, cnn, re_)) for _ in range(3)],
               "len": [med(matched(sl, cnn, rl)) for _ in range(3)]}
        for k in ("uniform", "expr", "len"): row[k + "_mean"] = float(np.mean(row[k]))
        per[mi] = row
        print(mi, len(U), round(row["cnn"], 3), round(row["uniform_mean"], 3), round(row["expr_mean"], 3), round(row["len_mean"], 3), row["cnn_median_len"], row["universe_median_len"], flush=True)
    st = lambda k: signtest(sum(v["cnn"] > v[k] for v in per.values()), sum(v["cnn"] < v[k] for v in per.values()))
    vd = lambda h: "CONFIRMED" if h["p_one_sided"] < 0.05 else "FALSIFIED"
    r1, r2, r3 = st("uniform_mean"), st("expr_mean"), st("len_mean")
    J = {"panel": panel, "K": K, "gates": {"G1_n": len(panel), "G1_pass": bool(len(panel) >= 25), "G2_coverage": cov, "G2_pass": bool(cov >= 0.85)},
         "R1_cnn_vs_uniform": dict(r1, verdict=vd(r1)), "R2_cnn_vs_expr": dict(r2, verdict=vd(r2)), "R3_cnn_vs_utrlen": dict(r3, verdict=vd(r3)),
         "cnn_shorter_utr_than_universe": int(sum(v["cnn_median_len"] < v["universe_median_len"] for v in per.values())),
         "per_mirna": per}
    json.dump(J, open(out, "w"), indent=1)
    print(json.dumps({k: J[k] for k in ["gates", "R1_cnn_vs_uniform", "R2_cnn_vs_expr", "R3_cnn_vs_utrlen", "cnn_shorter_utr_than_universe"]}))


if __name__ == "__main__":
    main(*sys.argv[1:])
