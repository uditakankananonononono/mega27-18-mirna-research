"""POST-HOC exploratory arm (added after H1-H3 were run and falsified; not pre-registered).
Same CNN top-K and first random set per miRNA, but default whole-genome background
(domain annotated). Tests whether apparent enrichment appears when the matched
universe is dropped - i.e. whether background choice alone manufactures signal."""
import json, random
import numpy as np, pandas as pd
import sys; sys.path.insert(0, "scripts")
from gprofiler_enrichment import gost, summarize, PANEL, K, signtest


def main(path="results/gprofiler_enrichment.json"):
    J = json.load(open(path))
    d = pd.read_csv("results/clip_rows_136.csv", usecols=["mirna", "gene", "min"])
    ph = {}
    for mi in PANEL:
        s = d[d.mirna == mi].drop_duplicates("gene"); U = s.gene.tolist()
        cnn = s.nsmallest(K, "min").gene.tolist()
        rnd = random.Random(f"{mi}-seed").sample(U, K)  # identical to first random set
        a, b = summarize(gost(cnn), K), summarize(gost(rnd), K)
        ph[mi] = {"cnn_nsig": a["n_sig"], "rand_nsig": b["n_sig"], "cnn_top": a["top"][:3], "rand_top": b["top"][:3]}
        print(mi, a["n_sig"], b["n_sig"], flush=True)
    w = sum(v["cnn_nsig"] > v["rand_nsig"] for v in ph.values()); l = sum(v["cnn_nsig"] < v["rand_nsig"] for v in ph.values())
    J["posthoc_genome_background"] = {"note": "exploratory, not pre-registered", "per_mirna": ph,
        "cnn_vs_random": signtest(w, l),
        "median_cnn_nsig": float(np.median([v["cnn_nsig"] for v in ph.values()])),
        "median_rand_nsig": float(np.median([v["rand_nsig"] for v in ph.values()]))}
    json.dump(J, open(path, "w"), indent=1)
    print(json.dumps({k: v for k, v in J["posthoc_genome_background"].items() if k != "per_mirna"}))


if __name__ == "__main__":
    main()
