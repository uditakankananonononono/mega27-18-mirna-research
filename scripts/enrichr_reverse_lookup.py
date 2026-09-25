"""Enrichr reverse-lookup audit (tool 25): can a target set name its own miRNA?

PRE-REGISTERED (written before any group comparison was run):
Enrichr API (maayanlab.cloud/Enrichr addList + enrich), library miRTarBase_2017
(3,240 miRNA terms). Same 13-miRNA panel, K=200 sets and seeds as
scripts/gprofiler_enrichment.py: CNN top-K, CLIP+ top-K, HEK293 expression top-K,
R=3 random draws from the per-miRNA CLIP universe. Statistic: rank of the query
miRNA's own term among all returned terms by Enrichr p-value (unreturned = worst+1).
Gates:
 G1 positive control: first 200 genes of the library's own hsa-miR-21-5p set
    return hsa-miR-21-5p at rank 1.
 G2 calibration: median own-term rank of random sets > 100.
 G3 all 13 panel miRNAs present as terms in the downloaded library GMT.
Hypotheses (one-sided sign tests across miRNAs, alpha 0.05):
 H1 CNN own-term rank < mean random own-term rank.
 H2 CLIP+ own-term rank < mean random rank (caveat: miRTarBase 2017 includes
    CLIP-derived evidence, so H2 may be partly circular).
 H3 CNN own-term rank < expression top-K own-term rank.
"""
import json, sys, time, random
import numpy as np, pandas as pd, requests
sys.path.insert(0, "scripts")
from gprofiler_enrichment import PANEL, K, R, signtest

BASE = "https://maayanlab.cloud/Enrichr/"
LIB = "miRTarBase_2017"


def enrich(genes, tries=4):
    for i in range(tries):
        try:
            r = requests.post(BASE + "addList", files={"list": (None, "\n".join(genes)), "description": (None, "mega27")}, timeout=60)
            r.raise_for_status(); uid = r.json()["userListId"]
            e = requests.get(BASE + "enrich", params={"userListId": uid, "backgroundType": LIB}, timeout=120)
            e.raise_for_status(); time.sleep(0.3)
            return e.json()[LIB]
        except Exception as ex:
            err = ex; time.sleep(3 * (i + 1))
    raise err


def own_rank(res, mi):
    res = sorted(res, key=lambda x: x[2])
    for i, x in enumerate(res):
        if x[1] == mi:
            return {"rank": i + 1, "p": x[2], "adj_p": x[6], "n_terms": len(res)}
    return {"rank": len(res) + 1, "p": None, "adj_p": None, "n_terms": len(res)}


def main(out="results/enrichr_reverse_lookup.json"):
    gmt = {l.split("\t")[0]: [g for g in l.rstrip("\n").split("\t")[2:] if g] for l in open(f"data/enrichr/{LIB}.gmt") if l.strip()}
    g3 = all(m in gmt for m in PANEL)
    pc = own_rank(enrich(gmt["hsa-miR-21-5p"][:200]), "hsa-miR-21-5p")
    d = pd.read_csv("results/clip_rows_136.csv", usecols=["mirna", "gene", "min", "label"])
    expr = pd.read_csv("results/hpa_hek293_ntpm.tsv", sep="\t").groupby("Gene name")["nTPM"].max()
    per = {}
    for mi in PANEL:
        s = d[d.mirna == mi].drop_duplicates("gene"); U = s.gene.tolist()
        cnn = s.nsmallest(K, "min").gene.tolist(); clip = s[s.label == 1].nsmallest(K, "min").gene.tolist()
        e = s.assign(x=s.gene.map(expr)).dropna(subset=["x"]); ex = e.nlargest(K, "x").gene.tolist()
        rng = random.Random(f"{mi}-seed"); rnd = [rng.sample(U, K) for _ in range(R)]
        row = {"library_set_size": len(gmt.get(mi, [])), "cnn": own_rank(enrich(cnn), mi), "clip": own_rank(enrich(clip), mi),
               "expr": own_rank(enrich(ex), mi), "random": [own_rank(enrich(x), mi) for x in rnd]}
        row["rand_mean_rank"] = float(np.mean([x["rank"] for x in row["random"]]))
        per[mi] = row
        print(mi, row["cnn"]["rank"], row["clip"]["rank"], row["expr"]["rank"], [x["rank"] for x in row["random"]], flush=True)
    rr = [x["rank"] for v in per.values() for x in v["random"]]
    st = lambda f: signtest(sum(f(v) > 0 for v in per.values()), sum(f(v) < 0 for v in per.values()))
    h1 = st(lambda v: v["rand_mean_rank"] - v["cnn"]["rank"])
    h2 = st(lambda v: v["rand_mean_rank"] - v["clip"]["rank"])
    h3 = st(lambda v: v["expr"]["rank"] - v["cnn"]["rank"])
    vd = lambda h: "CONFIRMED" if h["p_one_sided"] < 0.05 else "FALSIFIED"
    J = {"tool": "Enrichr API (maayanlab.cloud/Enrichr)", "library": LIB, "library_terms": len(gmt), "K": K, "R": R, "panel": PANEL,
         "gates": {"G1_miR21_rank": pc["rank"], "G1_pass": pc["rank"] == 1, "G2_median_random_rank": float(np.median(rr)),
                   "G2_pass": float(np.median(rr)) > 100, "G3_pass": g3},
         "H1_cnn_vs_random": dict(h1, verdict=vd(h1)), "H2_clip_vs_random": dict(h2, verdict=vd(h2)),
         "H3_cnn_vs_expression": dict(h3, verdict=vd(h3)), "per_mirna": per}
    json.dump(J, open(out, "w"), indent=1)
    print(json.dumps({k: J[k] for k in ["gates", "H1_cnn_vs_random", "H2_clip_vs_random", "H3_cnn_vs_expression"]}))


if __name__ == "__main__":
    main(*sys.argv[1:])
