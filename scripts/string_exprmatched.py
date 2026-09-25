"""POST-HOC confound control for STRING H2 (added after H2 was CONFIRMED; not pre-registered).
For each miRNA, draw R=3 random sets from the scored universe matched to the CLIP+ set's
HEK293-expression decile histogram (deciles over the universe; genes without HPA value form
their own bin). Test: CLIP+ density > mean expression-matched random density (sign test)."""
import json, random, sys
import numpy as np, pandas as pd
sys.path.insert(0, "scripts")
from string_coherence import edges
from gprofiler_enrichment import PANEL, K, R, signtest


def matched(s, clip, rng):
    b = pd.qcut(s.x.rank(method="first"), 10, labels=False).fillna(-1).astype(int)
    s = s.assign(b=b.values); pool = s[~s.gene.isin(clip)]
    need = s[s.gene.isin(clip)].b.value_counts()
    out = []
    for bin_, n in need.items():
        g = pool[pool.b == bin_].gene.tolist()
        out += rng.sample(g, min(n, len(g)))
    return out


def main(path="results/string_coherence.json"):
    J = json.load(open(path))
    d = pd.read_csv("results/clip_rows_136.csv", usecols=["mirna", "gene", "min", "label"])
    expr = pd.read_csv("results/hpa_hek293_ntpm.tsv", sep="\t").groupby("Gene name")["nTPM"].max()
    ph = {}
    for mi in PANEL:
        s = d[d.mirna == mi].drop_duplicates("gene").reset_index(drop=True)
        s = s.assign(x=s.gene.map(expr))
        clip = s[s.label == 1].nsmallest(K, "min").gene.tolist()
        rng = random.Random(f"{mi}-exprmatch")
        sets = [matched(s, clip, rng) for _ in range(R)]
        res = [edges(x) for x in sets]
        md = float(np.mean([r["density"] for r in res]))
        ph[mi] = {"clip_density": J["per_mirna"][mi]["clip"]["density"], "matched_mean_density": md,
                  "matched": res, "set_sizes": [len(x) for x in sets], "clip_size": len(clip)}
        print(mi, round(ph[mi]["clip_density"], 5), round(md, 5), [r["edges"] for r in res], flush=True)
    w = sum(v["clip_density"] > v["matched_mean_density"] for v in ph.values())
    l = sum(v["clip_density"] < v["matched_mean_density"] for v in ph.values())
    t = signtest(w, l)
    J["posthoc_H2_expression_matched_control"] = {"note": "post-hoc confound control, not pre-registered",
        "per_mirna": ph, "clip_vs_matched": dict(t, verdict="SURVIVES" if t["p_one_sided"] < 0.05 else "DOES NOT SURVIVE")}
    json.dump(J, open(path, "w"), indent=1)
    print(json.dumps(J["posthoc_H2_expression_matched_control"]["clip_vs_matched"]))


if __name__ == "__main__":
    main()
