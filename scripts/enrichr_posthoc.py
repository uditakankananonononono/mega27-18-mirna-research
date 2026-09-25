"""POST-HOC confound control for Enrichr H1 (added after H1 was CONFIRMED; not pre-registered).
The per-miRNA feature audit showed the univariate CNN signal is largely site-count confounding.
Arm: site-count top-200 (rank by n8 desc, n7 desc, UTR len desc) from the same universe.
Test: CNN own-term rank < site-count own-term rank (sign test across miRNAs)."""
import json, sys
import pandas as pd
sys.path.insert(0, "scripts")
from enrichr_reverse_lookup import enrich, own_rank, PANEL, K, signtest


def main(path="results/enrichr_reverse_lookup.json"):
    J = json.load(open(path))
    d = pd.read_csv("results/clip_rows_136.csv", usecols=["mirna", "gene", "n8", "n7", "len"])
    ph = {}
    for mi in PANEL:
        s = d[d.mirna == mi].drop_duplicates("gene")
        sc = s.sort_values(["n8", "n7", "len"], ascending=False).head(K).gene.tolist()
        r = own_rank(enrich(sc), mi)
        ph[mi] = {"sitecount": r, "cnn_rank": J["per_mirna"][mi]["cnn"]["rank"]}
        print(mi, ph[mi]["cnn_rank"], r["rank"], flush=True)
    w = sum(v["cnn_rank"] < v["sitecount"]["rank"] for v in ph.values()); l = sum(v["cnn_rank"] > v["sitecount"]["rank"] for v in ph.values())
    t = signtest(w, l)
    J["posthoc_sitecount_control"] = {"note": "post-hoc confound control, not pre-registered", "per_mirna": ph,
        "cnn_vs_sitecount": dict(t, verdict="SURVIVES" if t["p_one_sided"] < 0.05 else "DOES NOT SURVIVE"),
        "sitecount_vs_random_wins": sum(v["sitecount"]["rank"] < J["per_mirna"][m]["rand_mean_rank"] for m, v in ph.items())}
    json.dump(J, open(path, "w"), indent=1)
    print(json.dumps({k: v for k, v in J["posthoc_sitecount_control"].items() if k != "per_mirna"}))


if __name__ == "__main__" and len(sys.argv) == 1:
    main()


def exprmatched(path="results/enrichr_reverse_lookup.json"):
    """POST-HOC (not pre-registered): R=3 random sets matched to the CNN top-200's HEK293
    expression-decile histogram (same matcher as scripts/string_exprmatched.py).
    Test: CNN own-term rank < mean expression-matched rank."""
    import random, numpy as np
    from string_exprmatched import matched
    J = json.load(open(path))
    d = pd.read_csv("results/clip_rows_136.csv", usecols=["mirna", "gene", "min"])
    expr = pd.read_csv("results/hpa_hek293_ntpm.tsv", sep="\t").groupby("Gene name")["nTPM"].max()
    ph = {}
    for mi in PANEL:
        s = d[d.mirna == mi].drop_duplicates("gene").reset_index(drop=True); s = s.assign(x=s.gene.map(expr))
        cnn = s.nsmallest(K, "min").gene.tolist()
        rng = random.Random(f"{mi}-enrichr-exprmatch")
        rk = [own_rank(enrich(matched(s, cnn, rng)), mi)["rank"] for _ in range(3)]
        ph[mi] = {"cnn_rank": J["per_mirna"][mi]["cnn"]["rank"], "matched_ranks": rk, "matched_mean": float(np.mean(rk))}
        print(mi, ph[mi]["cnn_rank"], rk, flush=True)
    w = sum(v["cnn_rank"] < v["matched_mean"] for v in ph.values()); l = sum(v["cnn_rank"] > v["matched_mean"] for v in ph.values())
    t = signtest(w, l)
    J["posthoc_expression_matched_control"] = {"note": "post-hoc confound control, not pre-registered", "per_mirna": ph,
        "cnn_vs_matched": dict(t, verdict="SURVIVES" if t["p_one_sided"] < 0.05 else "DOES NOT SURVIVE")}
    json.dump(J, open(path, "w"), indent=1)
    print(json.dumps(J["posthoc_expression_matched_control"]["cnn_vs_matched"]))


if __name__ == "__main__" and len(sys.argv) > 1:
    exprmatched()
