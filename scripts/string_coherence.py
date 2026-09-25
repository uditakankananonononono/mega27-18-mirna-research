"""STRING v12 network-coherence audit of miRNA target sets (tool 21).

PRE-REGISTERED (written before any group comparison was run):
Same 13-miRNA panel, K=200 sets and seeded random draws as
scripts/gprofiler_enrichment.py (CNN top-K, CLIP+ top-K, HEK293 expression
top-K, R=3 random from the per-miRNA scored universe).
Statistic: number of STRING v12 edges (combined score >= 0.700) among the
mapped query proteins, via /api/json/network (query nodes only); density =
edges / C(n_mapped, 2).
Gates:
 G1 positive control: 20 canonical cell-cycle genes give >= 50 edges.
 G2 calibration: median random-set density < 0.01.
 G3 mapping: median mapped fraction >= 0.85 (via /api/json/get_string_ids).
Hypotheses (one-sided sign tests across miRNAs, alpha 0.05, on density):
 H1 CNN top-K density > mean random density.
 H2 CLIP+ top-K density > mean random density.
 H3 CNN top-K density > expression top-K density (confound control).
"""
import json, sys, time, random
import numpy as np, pandas as pd, requests
sys.path.insert(0, "scripts")
from gprofiler_enrichment import PANEL, CELLCYCLE, K, R, signtest

BASE = "https://string-db.org/api/json/"
CALLER = "mega27_laneF"


def post(ep, data, tries=4):
    data = dict(data, species=9606, caller_identity=CALLER)
    for i in range(tries):
        try:
            r = requests.post(BASE + ep, data=data, timeout=180)
            r.raise_for_status(); time.sleep(1.0)
            return r.json()
        except Exception as e:
            err = e; time.sleep(4 * (i + 1))
    raise err


def map_ids(genes):
    js = post("get_string_ids", {"identifiers": "\r".join(genes), "limit": 1, "echo_query": 1})
    m = {}
    for x in js:
        m.setdefault(x["queryItem"], x["stringId"])
    return m


def edges(genes):
    m = map_ids(genes)
    ids = sorted(set(m.values()))
    if len(ids) < 2:
        return {"n_query": len(genes), "n_mapped": len(ids), "edges": 0, "density": 0.0}
    js = post("network", {"identifiers": "\r".join(ids), "required_score": 700})
    pairs = {tuple(sorted((x["stringId_A"], x["stringId_B"]))) for x in js}
    n = len(ids)
    return {"n_query": len(genes), "n_mapped": n, "mapped_frac": n / len(genes),
            "edges": len(pairs), "density": len(pairs) / (n * (n - 1) / 2)}


def main(out="results/string_coherence.json"):
    d = pd.read_csv("results/clip_rows_136.csv", usecols=["mirna", "gene", "min", "label"])
    expr = pd.read_csv("results/hpa_hek293_ntpm.tsv", sep="\t").groupby("Gene name")["nTPM"].max()
    pc = edges(CELLCYCLE)
    per = {}
    for mi in PANEL:
        s = d[d.mirna == mi].drop_duplicates("gene"); U = s.gene.tolist()
        cnn = s.nsmallest(K, "min").gene.tolist()
        clip = s[s.label == 1].nsmallest(K, "min").gene.tolist()
        e = s.assign(x=s.gene.map(expr)).dropna(subset=["x"]); ex = e.nlargest(K, "x").gene.tolist()
        rng = random.Random(f"{mi}-seed"); rnd = [rng.sample(U, K) for _ in range(R)]
        row = {"universe": len(U), "cnn": edges(cnn), "clip": edges(clip), "expr": edges(ex),
               "random": [edges(r) for r in rnd]}
        row["rand_mean_density"] = float(np.mean([r["density"] for r in row["random"]]))
        per[mi] = row
        print(mi, row["cnn"]["edges"], row["clip"]["edges"], row["expr"]["edges"], [r["edges"] for r in row["random"]], flush=True)
    dens_r = [r["density"] for v in per.values() for r in v["random"]]
    mf = [x["mapped_frac"] for v in per.values() for x in [v["cnn"], v["clip"], v["expr"]] + v["random"]]

    def st(f):
        return signtest(sum(f(v) > 0 for v in per.values()), sum(f(v) < 0 for v in per.values()))
    h1 = st(lambda v: v["cnn"]["density"] - v["rand_mean_density"])
    h2 = st(lambda v: v["clip"]["density"] - v["rand_mean_density"])
    h3 = st(lambda v: v["cnn"]["density"] - v["expr"]["density"])
    vd = lambda h: "CONFIRMED" if h["p_one_sided"] < 0.05 else "FALSIFIED"
    ver = requests.get(BASE + "version", timeout=30).json()[0]["string_version"]
    J = {"tool": "STRING API (string-db.org)", "string_version": ver, "required_score": 700,
         "K": K, "R": R, "panel": PANEL,
         "gates": {"G1_cellcycle_edges": pc["edges"], "G1_pass": pc["edges"] >= 50,
                   "G2_median_random_density": float(np.median(dens_r)), "G2_pass": float(np.median(dens_r)) < 0.01,
                   "G3_median_mapped_frac": float(np.median(mf)), "G3_pass": float(np.median(mf)) >= 0.85},
         "H1_cnn_vs_random": dict(h1, verdict=vd(h1)), "H2_clip_vs_random": dict(h2, verdict=vd(h2)),
         "H3_cnn_vs_expression": dict(h3, verdict=vd(h3)), "per_mirna": per}
    json.dump(J, open(out, "w"), indent=1)
    print(json.dumps({k: J[k] for k in ["string_version", "gates", "H1_cnn_vs_random", "H2_clip_vs_random", "H3_cnn_vs_expression"]}))


if __name__ == "__main__":
    main(*sys.argv[1:])
