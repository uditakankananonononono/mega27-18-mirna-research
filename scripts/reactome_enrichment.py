"""Reactome AnalysisService pathway audit with controls fixed up front (tool 26).

PRE-REGISTERED (written before any group comparison was run):
Reactome AnalysisService /identifiers/projection (human, whole-Reactome background).
13-miRNA panel and K=200 as scripts/gprofiler_enrichment.py. Arms per miRNA:
CNN top-K; site-count top-K (n8, n7, len desc); 3 random sets matched to the CNN
set's HEK293 expression deciles (string_exprmatched.matched); 3 uniform random.
Statistic: number of pathways with entities FDR < 0.05 (all pages).
Gates:
 G1 positive control: 20 cell-cycle genes give >= 20 pathways at FDR<0.05 incl.
    one whose name contains "Cell Cycle".
 G2 calibration: median uniform-random pathway count <= 5.
 G3 median identifiers-not-found fraction < 0.10.
Hypotheses (one-sided sign tests across miRNAs, alpha 0.05, ties excluded):
 H1 CNN count > mean uniform count.
 H2 CNN count > site-count count.
 H3 CNN count > mean expression-matched count.
"""
import json, random, sys, time
import numpy as np, pandas as pd, requests
sys.path.insert(0, "scripts")
from gprofiler_enrichment import PANEL, CELLCYCLE, K, signtest
from string_exprmatched import matched

URL = "https://reactome.org/AnalysisService/identifiers/projection"


def analyse(genes, tries=4):
    for i in range(tries):
        try:
            r = requests.post(URL, data="\n".join(genes), headers={"Content-Type": "text/plain"},
                              params={"pageSize": 3000, "page": 1, "species": 9606, "includeDisease": "true"}, timeout=120)
            r.raise_for_status(); js = r.json(); time.sleep(0.3)
            sig = [p for p in js.get("pathways", []) if p["entities"]["fdr"] < 0.05]
            return {"n_sig": len(sig), "top": [(p["stId"], p["name"], p["entities"]["fdr"]) for p in sorted(sig, key=lambda p: p["entities"]["fdr"])[:3]],
                    "not_found_frac": js.get("identifiersNotFound", 0) / len(genes)}
        except Exception as e:
            err = e; time.sleep(3 * (i + 1))
    raise err


def main(out="results/reactome_enrichment.json"):
    pc = analyse(CELLCYCLE)
    d = pd.read_csv("results/clip_rows_136.csv", usecols=["mirna", "gene", "min", "n8", "n7", "len"])
    expr = pd.read_csv("results/hpa_hek293_ntpm.tsv", sep="\t").groupby("Gene name")["nTPM"].max()
    per = {}
    for mi in PANEL:
        s = d[d.mirna == mi].drop_duplicates("gene").reset_index(drop=True); s = s.assign(x=s.gene.map(expr))
        U = s.gene.tolist(); cnn = s.nsmallest(K, "min").gene.tolist()
        sc = s.sort_values(["n8", "n7", "len"], ascending=False).head(K).gene.tolist()
        rm, ru = random.Random(f"{mi}-reactome-match"), random.Random(f"{mi}-reactome-unif")
        row = {"cnn": analyse(cnn), "sitecount": analyse(sc), "matched": [analyse(matched(s, cnn, rm)) for _ in range(3)],
               "uniform": [analyse(ru.sample(U, K)) for _ in range(3)]}
        row["matched_mean"] = float(np.mean([x["n_sig"] for x in row["matched"]])); row["uniform_mean"] = float(np.mean([x["n_sig"] for x in row["uniform"]]))
        per[mi] = row
        print(mi, row["cnn"]["n_sig"], row["sitecount"]["n_sig"], [x["n_sig"] for x in row["matched"]], [x["n_sig"] for x in row["uniform"]], flush=True)
    un = [x["n_sig"] for v in per.values() for x in v["uniform"]]
    nf = [x["not_found_frac"] for v in per.values() for x in [v["cnn"], v["sitecount"]] + v["matched"] + v["uniform"]]
    st = lambda f: signtest(sum(f(v) > 0 for v in per.values()), sum(f(v) < 0 for v in per.values()))
    vd = lambda h: "CONFIRMED" if h["p_one_sided"] < 0.05 else "FALSIFIED"
    h1 = st(lambda v: v["cnn"]["n_sig"] - v["uniform_mean"]); h2 = st(lambda v: v["cnn"]["n_sig"] - v["sitecount"]["n_sig"])
    h3 = st(lambda v: v["cnn"]["n_sig"] - v["matched_mean"])
    ver = requests.get("https://reactome.org/ContentService/data/database/version", timeout=30).text.strip()
    J = {"tool": "Reactome AnalysisService (reactome.org)", "reactome_version": ver, "K": K, "panel": PANEL,
         "gates": {"G1_cellcycle_nsig": pc["n_sig"], "G1_top": pc["top"], "G1_pass": pc["n_sig"] >= 20 and any("Cell Cycle" in t[1] for t in pc["top"]),
                   "G2_median_uniform_nsig": float(np.median(un)), "G2_pass": float(np.median(un)) <= 5,
                   "G3_median_not_found": float(np.median(nf)), "G3_pass": float(np.median(nf)) < 0.10},
         "H1_cnn_vs_uniform": dict(h1, verdict=vd(h1)), "H2_cnn_vs_sitecount": dict(h2, verdict=vd(h2)),
         "H3_cnn_vs_exprmatched": dict(h3, verdict=vd(h3)), "per_mirna": per}
    json.dump(J, open(out, "w"), indent=1)
    print(json.dumps({k: J[k] for k in ["reactome_version", "gates", "H1_cnn_vs_uniform", "H2_cnn_vs_sitecount", "H3_cnn_vs_exprmatched"]}))


if __name__ == "__main__":
    main(*sys.argv[1:])
