"""PRE-REGISTERED held-out replication of the Enrichr reverse-lookup positive
(written after results/enrichr_reverse_lookup.json and before running this file).

Panel: miRNAs in results/clip_rows_136.csv that are terms in miRTarBase_2017,
excluding the 13 discovery miRNAs; shuffled with random.Random("replication-2026")
and the first 30 taken. Sets (K=200): CNN top-K; site-count top-K (n8, n7, len desc);
3 random sets matched to the CNN set's HEK293 expression deciles; 3 uniform random.
Statistic: own-term rank (Enrichr p-value order).
Hypotheses (one-sided sign tests across the 30, alpha 0.05, ties excluded):
 R1 CNN rank < mean uniform-random rank.
 R2 CNN rank < site-count rank.
 R3 CNN rank < mean expression-matched rank.
Discovery claim holds only if R2 AND R3 are confirmed.
"""
import json, random, sys
import numpy as np, pandas as pd
sys.path.insert(0, "scripts")
from enrichr_reverse_lookup import enrich, own_rank, PANEL, K, signtest, LIB
from string_exprmatched import matched


def main(out="results/enrichr_replication.json"):
    terms = {l.split("\t")[0] for l in open(f"data/enrichr/{LIB}.gmt") if l.strip()}
    d = pd.read_csv("results/clip_rows_136.csv", usecols=["mirna", "gene", "min", "n8", "n7", "len"])
    cand = sorted(m for m in d.mirna.unique() if m in terms and m not in PANEL)
    random.Random("replication-2026").shuffle(cand); panel = cand[:30]
    expr = pd.read_csv("results/hpa_hek293_ntpm.tsv", sep="\t").groupby("Gene name")["nTPM"].max()
    per = {}
    for mi in panel:
        s = d[d.mirna == mi].drop_duplicates("gene").reset_index(drop=True); s = s.assign(x=s.gene.map(expr))
        U = s.gene.tolist(); cnn = s.nsmallest(K, "min").gene.tolist()
        sc = s.sort_values(["n8", "n7", "len"], ascending=False).head(K).gene.tolist()
        rm = random.Random(f"{mi}-rep-match"); ru = random.Random(f"{mi}-rep-unif")
        row = {"universe": len(U), "cnn": own_rank(enrich(cnn), mi)["rank"], "sitecount": own_rank(enrich(sc), mi)["rank"],
               "matched": [own_rank(enrich(matched(s, cnn, rm)), mi)["rank"] for _ in range(3)],
               "uniform": [own_rank(enrich(ru.sample(U, K)), mi)["rank"] for _ in range(3)]}
        row["matched_mean"] = float(np.mean(row["matched"])); row["uniform_mean"] = float(np.mean(row["uniform"]))
        per[mi] = row
        print(mi, row["cnn"], row["sitecount"], row["matched"], row["uniform"], flush=True)
        json.dump({"partial": True, "per_mirna": per}, open(out, "w"), indent=1)
    st = lambda f: signtest(sum(f(v) > 0 for v in per.values()), sum(f(v) < 0 for v in per.values()))
    vd = lambda h: "CONFIRMED" if h["p_one_sided"] < 0.05 else "FALSIFIED"
    r1 = st(lambda v: v["uniform_mean"] - v["cnn"]); r2 = st(lambda v: v["sitecount"] - v["cnn"]); r3 = st(lambda v: v["matched_mean"] - v["cnn"])
    J = {"design": "pre-registered held-out replication of enrichr_reverse_lookup H1 + post-hoc controls", "library": LIB, "panel": panel, "K": K,
         "median_rank": {"cnn": float(np.median([v["cnn"] for v in per.values()])), "sitecount": float(np.median([v["sitecount"] for v in per.values()])),
                         "matched": float(np.median([v["matched_mean"] for v in per.values()])), "uniform": float(np.median([v["uniform_mean"] for v in per.values()]))},
         "R1_cnn_vs_uniform": dict(r1, verdict=vd(r1)), "R2_cnn_vs_sitecount": dict(r2, verdict=vd(r2)), "R3_cnn_vs_exprmatched": dict(r3, verdict=vd(r3)),
         "discovery_replicates": vd(r2) == "CONFIRMED" and vd(r3) == "CONFIRMED", "per_mirna": per}
    json.dump(J, open(out, "w"), indent=1)
    print(json.dumps({k: J[k] for k in ["median_rank", "R1_cnn_vs_uniform", "R2_cnn_vs_sitecount", "R3_cnn_vs_exprmatched", "discovery_replicates"]}))


if __name__ == "__main__":
    main(*sys.argv[1:])
