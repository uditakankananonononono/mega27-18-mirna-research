"""g:Profiler functional-enrichment audit of miRNA target sets (tool 20).

PRE-REGISTERED (written before any group comparison was run):
Panel: 13 well-studied miRNAs present in results/clip_rows_136.csv.
Per miRNA, universe U = all genes scored for that miRNA (custom background).
Sets (size K=200): CNN top-K (lowest 'min' context score), CLIP+ top-K
(label=1 ranked by CNN score, or all if fewer), expression top-K
(HPA HEK293 nTPM, confound control), and R=3 random size-K draws (seeded).
Query: g:GOSt, sources GO:BP/KEGG/REAC, g:SCS threshold 0.05, domain custom.
Gates:
 G1 positive control: 20 canonical cell-cycle genes vs whole-genome background
    must return GO:0007049 (cell cycle) significant.
 G2 calibration: median significant-term count over all random sets <= 2.
 G3 mapping: median unmapped fraction of queried genes < 0.10.
Hypotheses:
 H1 CNN top-K sets have more significant terms than the mean of their
    matched random sets (one-sided sign test across miRNAs, alpha 0.05).
 H2 term-set Jaccard(CNN, CLIP+) exceeds mean Jaccard(random, CLIP+)
    (one-sided sign test across miRNAs).
 H3 confound: CNN top-K term count exceeds expression top-K term count
    (sign test; if not, enrichment is not attributable to targeting).
"""
import json, sys, time, random
import numpy as np, pandas as pd, requests
from scipy.stats import binomtest

API = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
PANEL = ["hsa-let-7a-5p", "hsa-miR-1-3p", "hsa-miR-122-5p", "hsa-miR-125b-5p",
         "hsa-miR-145-5p", "hsa-miR-155-5p", "hsa-miR-16-5p", "hsa-miR-17-5p",
         "hsa-miR-19a-3p", "hsa-miR-200c-3p", "hsa-miR-21-5p", "hsa-miR-34a-5p",
         "hsa-miR-92a-3p"]
CELLCYCLE = ["CDK1", "CCNB1", "CCNA2", "CDC20", "PLK1", "BUB1", "AURKA", "CCNE1",
             "CDK2", "CDK4", "CCND1", "E2F1", "MCM2", "PCNA", "CDC25A", "CHEK1",
             "WEE1", "RB1", "CDKN1A", "MAD2L1"]
K, R = 200, 3


def gost(genes, background=None, tries=4):
    body = {"organism": "hsapiens", "query": list(genes),
            "sources": ["GO:BP", "KEGG", "REAC"], "user_threshold": 0.05,
            "significance_threshold_method": "g_SCS", "no_evidences": True}
    if background is not None:
        body["background"] = list(background)
        body["domain_scope"] = "custom"
    for i in range(tries):
        try:
            r = requests.post(API, json=body, timeout=120)
            r.raise_for_status()
            return r.json()
        except Exception as e:  # network retry
            err = e; time.sleep(3 * (i + 1))
    raise err


def summarize(js, nq):
    res = [x for x in js.get("result", []) if x.get("significant")]
    meta = js.get("meta", {}).get("genes_metadata", {})
    failed = meta.get("failed", [])
    amb = meta.get("ambiguous", {}) or {}
    return {"n_sig": len(res), "terms": sorted(x["native"] for x in res),
            "top": [(x["native"], x["name"], x["p_value"]) for x in
                    sorted(res, key=lambda x: x["p_value"])[:5]],
            "unmapped_frac": len(failed) / max(nq, 1), "n_ambiguous": len(amb)}


def jaccard(a, b):
    a, b = set(a), set(b)
    return len(a & b) / len(a | b) if (a | b) else 0.0


def signtest(wins, losses):
    n = wins + losses
    p = binomtest(wins, n, 0.5, alternative="greater").pvalue if n else 1.0
    return {"wins": wins, "losses": losses, "ties_excluded": True, "p_one_sided": p}


def main(out="results/gprofiler_enrichment.json"):
    d = pd.read_csv("results/clip_rows_136.csv", usecols=["mirna", "gene", "min", "label"])
    expr = pd.read_csv("results/hpa_hek293_ntpm.tsv", sep="\t")
    expr = expr.groupby("Gene name")["nTPM"].max()
    pc = summarize(gost(CELLCYCLE), len(CELLCYCLE))
    g1 = "GO:0007049" in pc["terms"]
    per = {}
    for mi in PANEL:
        s = d[d.mirna == mi].drop_duplicates("gene")
        U = s.gene.tolist()
        cnn = s.nsmallest(K, "min").gene.tolist()
        clip = s[s.label == 1].nsmallest(K, "min").gene.tolist()
        e = s.assign(x=s.gene.map(expr)).dropna(subset=["x"])
        ex = e.nlargest(K, "x").gene.tolist()
        rng = random.Random(f"{mi}-seed")
        rnd = [rng.sample(U, K) for _ in range(R)]
        row = {"universe": len(U), "cnn": summarize(gost(cnn, U), K),
               "clip": summarize(gost(clip, U), len(clip)),
               "expr": summarize(gost(ex, U), len(ex)),
               "random": [summarize(gost(r, U), K) for r in rnd]}
        row["cnn_clip_overlap_genes"] = len(set(cnn) & set(clip))
        row["J_cnn_clip"] = jaccard(row["cnn"]["terms"], row["clip"]["terms"])
        row["J_rand_clip"] = float(np.mean([jaccard(r["terms"], row["clip"]["terms"]) for r in row["random"]]))
        row["rand_mean_nsig"] = float(np.mean([r["n_sig"] for r in row["random"]]))
        per[mi] = row
        print(mi, row["cnn"]["n_sig"], row["clip"]["n_sig"], row["expr"]["n_sig"],
              [r["n_sig"] for r in row["random"]], round(row["J_cnn_clip"], 3), flush=True)
    rand_all = [r["n_sig"] for v in per.values() for r in v["random"]]
    unm = [x["unmapped_frac"] for v in per.values() for x in [v["cnn"], v["clip"], v["expr"]] + v["random"]]
    g2 = float(np.median(rand_all)) <= 2
    g3 = float(np.median(unm)) < 0.10

    def st(f):
        w = sum(1 for v in per.values() if f(v) > 0); l = sum(1 for v in per.values() if f(v) < 0)
        return signtest(w, l)
    h1 = st(lambda v: v["cnn"]["n_sig"] - v["rand_mean_nsig"])
    h2 = st(lambda v: v["J_cnn_clip"] - v["J_rand_clip"])
    h3 = st(lambda v: v["cnn"]["n_sig"] - v["expr"]["n_sig"])
    verdict = lambda h: "CONFIRMED" if h["p_one_sided"] < 0.05 else "FALSIFIED"
    outj = {"tool": "g:Profiler g:GOSt REST API (biit.cs.ut.ee/gprofiler)",
            "gprofiler_version": pc_meta_version(), "K": K, "R": R, "panel": PANEL,
            "gates": {"G1_cellcycle_GO0007049": g1, "G1_top": pc["top"],
                      "G2_median_random_nsig": float(np.median(rand_all)), "G2_pass": g2,
                      "G3_median_unmapped": float(np.median(unm)), "G3_pass": g3},
            "H1_cnn_vs_random": dict(h1, verdict=verdict(h1)),
            "H2_termJaccard_cnn_vs_random": dict(h2, verdict=verdict(h2)),
            "H3_cnn_vs_expression": dict(h3, verdict=verdict(h3)),
            "per_mirna": per}
    json.dump(outj, open(out, "w"), indent=1)
    print(json.dumps({k: outj[k] for k in ["gates", "H1_cnn_vs_random", "H2_termJaccard_cnn_vs_random", "H3_cnn_vs_expression"]}, indent=1))


_VER = {}


def pc_meta_version():
    if "v" not in _VER:
        js = gost(["TP53"])
        _VER["v"] = js.get("meta", {}).get("version", "unknown")
    return _VER["v"]


if __name__ == "__main__":
    main(*sys.argv[1:])
