"""UniProt transcription-regulator enrichment of target sets (tool 29).

PRE-REGISTERED (written before any group comparison was run):
Classic observation: conserved miRNA targets are enriched for transcription
regulators. Annotation: UniProtKB reviewed human, keyword KW-0805
(Transcription regulation), release header recorded; universe restricted to genes
with a reviewed human UniProt entry.
Panel: 13 miRNAs of scripts/gprofiler_enrichment.py, K=200.
Arm A (CNN, CLIP universe of non-conserved seed-match genes): CNN top-K;
  site-count top-K (n8, n7, len desc); 3 expression-matched random sets
  (string_exprmatched.matched); 3 uniform random.
Arm B (positive control, all-UTR universe = TargetScan UTR genes with UniProt entry):
  TargetScan conserved top-K by context++ vs 3 uniform random from that universe.
Statistic: fraction of set genes carrying KW-0805.
Gates:
 G1 positive control: TargetScan conserved TF fraction > mean random for >= 10/13
    miRNAs with >= 20 conserved-site genes (sign test p < 0.05).
 G2 KW-0805 prevalence in reviewed human proteome between 5% and 20%.
Hypotheses (one-sided sign tests across the 13, alpha 0.05, ties excluded):
 H1 CNN TF fraction > mean uniform.  H2 CNN > site count.  H3 CNN > mean expression-matched.
"""
import json, random, sys
import numpy as np, pandas as pd, requests
sys.path.insert(0, "scripts")
from gprofiler_enrichment import PANEL, K, signtest
from string_exprmatched import matched


def main(out="results/uniprot_tf_enrichment.json"):
    tf = set(pd.read_csv("data/uniprot/kw0805_human.tsv", sep="\t")["Gene Names (primary)"].dropna())
    allp = set(pd.read_csv("data/uniprot/human_reviewed.tsv", sep="\t")["Gene Names (primary)"].dropna())
    rel = requests.head("https://rest.uniprot.org/uniprotkb/search?query=reviewed:true&size=1", timeout=30).headers.get("x-uniprot-release")
    frac = lambda g: float(np.mean([x in tf for x in g])) if len(g) else float("nan")
    d = pd.read_csv("results/clip_rows_136.csv", usecols=["mirna", "gene", "min", "n8", "n7", "len"]).drop_duplicates(["mirna", "gene"])
    d = d[d.gene.isin(allp)]
    expr = pd.read_csv("results/hpa_hek293_ntpm.tsv", sep="\t").groupby("Gene name")["nTPM"].max()
    utrg = sorted(set(pd.read_csv("data/human_utrs.tsv", sep="\t", usecols=[1]).iloc[:, 0]) & allp)
    ts = pd.read_csv("data/human_sites.tsv", sep="\t", usecols=["gene", "mirna", "context_pp"])
    per = {}
    for mi in PANEL:
        s = d[d.mirna == mi].reset_index(drop=True); s = s.assign(x=s.gene.map(expr)); U = s.gene.tolist()
        cnn = s.nsmallest(K, "min").gene.tolist(); sc = s.sort_values(["n8", "n7", "len"], ascending=False).head(K).gene.tolist()
        rm, ru, rb = random.Random(f"{mi}-uni-match"), random.Random(f"{mi}-uni-unif"), random.Random(f"{mi}-uni-utr")
        row = {"universe": len(U), "universe_tf_frac": frac(U), "cnn": frac(cnn), "sitecount": frac(sc),
               "matched": [frac(matched(s, cnn, rm)) for _ in range(3)], "uniform": [frac(ru.sample(U, K)) for _ in range(3)]}
        t = ts[(ts.mirna == mi) & ts.gene.isin(utrg)].groupby("gene").context_pp.min().nsmallest(K)
        row["ts_n"] = len(t); row["ts_frac"] = frac(list(t.index)) if len(t) >= 20 else None
        row["ts_rand"] = [frac(rb.sample(utrg, max(len(t), 1))) for _ in range(3)] if len(t) >= 20 else None
        row["matched_mean"] = float(np.mean(row["matched"])); row["uniform_mean"] = float(np.mean(row["uniform"]))
        per[mi] = row
        print(mi, len(U), round(row["cnn"], 3), round(row["sitecount"], 3), round(row["matched_mean"], 3), round(row["uniform_mean"], 3), "| TS", row["ts_n"], row["ts_frac"] and round(row["ts_frac"], 3), row["ts_rand"] and round(float(np.mean(row["ts_rand"])), 3), flush=True)
    st = lambda f, vs: signtest(sum(f(v) > 0 for v in vs), sum(f(v) < 0 for v in vs))
    vd = lambda h: "CONFIRMED" if h["p_one_sided"] < 0.05 else "FALSIFIED"
    tsv = [v for v in per.values() if v["ts_frac"] is not None]
    g1 = st(lambda v: v["ts_frac"] - float(np.mean(v["ts_rand"])), tsv)
    prev = len(tf & allp) / len(allp)
    h1 = st(lambda v: v["cnn"] - v["uniform_mean"], per.values()); h2 = st(lambda v: v["cnn"] - v["sitecount"], per.values())
    h3 = st(lambda v: v["cnn"] - v["matched_mean"], per.values())
    J = {"tool": "UniProt REST API (rest.uniprot.org, KW-0805)", "uniprot_release": rel, "n_tf": len(tf), "n_reviewed": len(allp), "K": K, "panel": PANEL,
         "gates": {"G1_targetscan_vs_random": g1, "G1_n_eligible": len(tsv), "G1_pass": bool(len(tsv) >= 10 and g1["wins"] >= 10 and g1["p_one_sided"] < 0.05),
                   "G2_kw0805_prevalence": prev, "G2_pass": bool(0.05 <= prev <= 0.20)},
         "H1_cnn_vs_uniform": dict(h1, verdict=vd(h1)), "H2_cnn_vs_sitecount": dict(h2, verdict=vd(h2)), "H3_cnn_vs_exprmatched": dict(h3, verdict=vd(h3)),
         "per_mirna": per}
    json.dump(J, open(out, "w"), indent=1)
    print(json.dumps({k: J[k] for k in ["uniprot_release", "gates", "H1_cnn_vs_uniform", "H2_cnn_vs_sitecount", "H3_cnn_vs_exprmatched"]}))


if __name__ == "__main__":
    main(*sys.argv[1:])
