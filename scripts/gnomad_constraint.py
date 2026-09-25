"""gnomAD loss-of-function constraint of target sets (tool 30).

PRE-REGISTERED (written before any group comparison was run):
Observation to test: conserved miRNA targets are dosage-sensitive, i.e. more
LoF-constrained than random genes. Constraint: gnomAD v2.1.1 LOEUF
(oe_lof_upper; lower = more constrained), canonical-transcript gene table
(gs://gcp-public-data--gnomad/release/2.1.1/constraint/, counted once).
Panel: 13 miRNAs of scripts/gprofiler_enrichment.py, K=200. Genes without LOEUF are
dropped from every universe before any set is drawn.
Arm A (CNN, CLIP universe of non-conserved seed-match genes): CNN top-K;
  site-count top-K (n8, n7, len desc); 3 expression-matched random sets
  (string_exprmatched.matched); 3 uniform random.
Arm B (positive control, TargetScan UTR universe): TargetScan conserved top-K by
  context++ vs 3 uniform random sets from that universe.
Statistic: median LOEUF of the set; a "win" = lower median than the comparator.
Gates:
 G1 positive control: conserved top-K lower median LOEUF than mean random for
    >= 10/13 miRNAs (sign test p < 0.05).
 G2 coverage: >= 85% of CLIP-universe genes (pooled over the panel) have a LOEUF.
Post-hoc block (added after results, labelled not pre-registered): two-sided test of
the reversed direction.
Hypotheses (one-sided sign tests across 13, alpha 0.05, ties excluded):
 H1 CNN more constrained than uniform. H2 than site count. H3 than expression-matched.
"""
import gzip, json, random, sys
import numpy as np, pandas as pd
sys.path.insert(0, "scripts")
from gprofiler_enrichment import PANEL, K, signtest
from string_exprmatched import matched


def main(out="results/gnomad_constraint.json"):
    g = pd.read_csv("data/gnomad/lof_v211.txt.bgz", sep="\t", compression="gzip", usecols=["gene", "oe_lof_upper"]).dropna()
    lo = g.groupby("gene").oe_lof_upper.min()
    med = lambda gs: float(np.median(lo.reindex(gs).values))
    d0 = pd.read_csv("results/clip_rows_136.csv", usecols=["mirna", "gene", "min", "n8", "n7", "len"]).drop_duplicates(["mirna", "gene"])
    d0 = d0[d0.mirna.isin(PANEL)]
    cov = float(d0.gene.isin(lo.index).mean())
    d = d0[d0.gene.isin(lo.index)]
    expr = pd.read_csv("results/hpa_hek293_ntpm.tsv", sep="\t").groupby("Gene name")["nTPM"].max()
    utrg = sorted(set(pd.read_csv("data/human_utrs.tsv", sep="\t", usecols=[1]).iloc[:, 0]) & set(lo.index))
    ts = pd.read_csv("data/human_sites.tsv", sep="\t", usecols=["gene", "mirna", "context_pp"])
    per = {}
    for mi in PANEL:
        s = d[d.mirna == mi].reset_index(drop=True); s = s.assign(x=s.gene.map(expr)); U = s.gene.tolist()
        cnn = s.nsmallest(K, "min").gene.tolist(); sc = s.sort_values(["n8", "n7", "len"], ascending=False).head(K).gene.tolist()
        rm, ru, rb = random.Random(f"{mi}-gn-match"), random.Random(f"{mi}-gn-unif"), random.Random(f"{mi}-gn-utr")
        row = {"universe": len(U), "universe_median": med(U), "cnn": med(cnn), "sitecount": med(sc),
               "matched": [med(matched(s, cnn, rm)) for _ in range(3)], "uniform": [med(ru.sample(U, K)) for _ in range(3)]}
        t = ts[(ts.mirna == mi) & ts.gene.isin(utrg)].groupby("gene").context_pp.min().nsmallest(K)
        row["ts_n"] = len(t); row["ts_median"] = med(list(t.index)); row["ts_rand"] = [med(rb.sample(utrg, len(t))) for _ in range(3)]
        row["matched_mean"] = float(np.mean(row["matched"])); row["uniform_mean"] = float(np.mean(row["uniform"]))
        per[mi] = row
        print(mi, len(U), *(round(row[k], 3) for k in ["cnn", "sitecount", "matched_mean", "uniform_mean", "ts_median"]), round(float(np.mean(row["ts_rand"])), 3), flush=True)
    st = lambda f: signtest(sum(f(v) < 0 for v in per.values()), sum(f(v) > 0 for v in per.values()))
    vd = lambda h: "CONFIRMED" if h["p_one_sided"] < 0.05 else "FALSIFIED"
    g1 = st(lambda v: v["ts_median"] - float(np.mean(v["ts_rand"])))
    h1 = st(lambda v: v["cnn"] - v["uniform_mean"]); h2 = st(lambda v: v["cnn"] - v["sitecount"]); h3 = st(lambda v: v["cnn"] - v["matched_mean"])
    J = {"tool": "gnomAD v2.1.1 constraint (LOEUF)", "n_genes_loeuf": int(len(lo)), "K": K, "panel": PANEL,
         "gates": {"G1_targetscan_vs_random": g1, "G1_pass": bool(g1["wins"] >= 10 and g1["p_one_sided"] < 0.05),
                   "G2_coverage": cov, "G2_pass": bool(cov >= 0.85)},
         "H1_cnn_vs_uniform": dict(h1, verdict=vd(h1)), "H2_cnn_vs_sitecount": dict(h2, verdict=vd(h2)), "H3_cnn_vs_exprmatched": dict(h3, verdict=vd(h3)),
         "post_hoc_not_preregistered": {"note": "direction reversed: CNN top-K LESS constrained; two-sided sign test",
             "cnn_higher_loeuf_than_uniform": int(h1["losses"]), "cnn_higher_loeuf_than_matched": int(h3["losses"]),
             "p_two_sided_uniform": min(1.0, 2 * signtest(h1["losses"], h1["wins"])["p_one_sided"]),
             "p_two_sided_matched": min(1.0, 2 * signtest(h3["losses"], h3["wins"])["p_one_sided"])},
         "per_mirna": per}
    json.dump(J, open(out, "w"), indent=1)
    print(json.dumps({k: J[k] for k in ["gates", "H1_cnn_vs_uniform", "H2_cnn_vs_sitecount", "H3_cnn_vs_exprmatched", "post_hoc_not_preregistered"]}))


if __name__ == "__main__":
    main(*sys.argv[1:])
