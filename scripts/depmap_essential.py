"""PRE-REGISTERED DepMap common-essential test (tool 35).

Label: DepMap 24Q4 Public CRISPRInferredCommonEssentials.csv (figshare file 51064916; counted once);
entries "SYMBOL (EntrezID)" parsed to symbols. Question: are CNN top-200 non-conserved
candidates depleted of cell-essential genes, extending the dosage/disease finding (tools 30-33)?
Committed BEFORE running.
Panel: 13 discovery + 30 held-out miRNAs; CLIP universe; K=200. Comparators (3 draws, seeds
f"{mi}-dm-<arm>"): uniform, expression-matched, UTR-length-matched. Essential genes are
highly expressed, so the expression-matched comparison is the key one.
Positive control: TargetScan conserved top-200 vs 3 uniform random UTR genes (>= 20 conserved-site genes).
Gates: G1 positive control direction unrestricted: two-sided sign test p<0.05 (conserved targets
  are known to be enriched in some essential classes and depleted in housekeeping machinery,
  so only a detectable difference is required to show the label is informative);
  G2 essential prevalence in pooled universe between 2% and 20%.
Hypotheses (CNN fewer essential genes), one-sided sign tests over 43 miRNAs, ties excluded,
alpha 0.05: E1 vs uniform; E2 vs expression-matched; E3 vs UTR-length-matched.
If fewer than 15 non-tied miRNAs remain, the hypothesis is INCONCLUSIVE.
"""
import json, random, sys
import numpy as np, pandas as pd
sys.path.insert(0, "scripts")
from gprofiler_enrichment import PANEL, K, signtest
from string_exprmatched import matched


def main(out="results/depmap_essential.json"):
    raw = pd.read_csv("data/depmap/CRISPRInferredCommonEssentials_24Q4.csv").iloc[:, 0].astype(str)
    ess = set(raw.str.split(" ").str[0])
    panel = PANEL + json.load(open("results/gnomad_replication.json"))["panel"]
    d = pd.read_csv("results/clip_rows_136.csv", usecols=["mirna", "gene", "min", "len"]).drop_duplicates(["mirna", "gene"])
    d = d[d.mirna.isin(panel)]
    expr = pd.read_csv("results/hpa_hek293_ntpm.tsv", sep="\t").groupby("Gene name")["nTPM"].max()
    utrg = sorted(set(pd.read_csv("data/human_utrs.tsv", sep="\t", usecols=[1]).iloc[:, 0]))
    ts = pd.read_csv("data/human_sites.tsv", sep="\t", usecols=["gene", "mirna", "context_pp"])
    c = lambda gs: int(sum(g in ess for g in gs))
    per = {}
    for mi in panel:
        s = d[d.mirna == mi].reset_index(drop=True); U = s.gene.tolist(); cnn = s.nsmallest(K, "min").gene.tolist()
        se, sl = s.assign(x=s.gene.map(expr)), s.assign(x=s["len"].astype(float))
        ru, re_, rl, rb = (random.Random(f"{mi}-dm-{a}") for a in ("unif", "expr", "len", "utr"))
        row = {"universe": len(U), "universe_ess": c(U), "cnn": c(cnn), "uniform_mean": float(np.mean([c(ru.sample(U, K)) for _ in range(3)])),
               "expr_mean": float(np.mean([c(matched(se, cnn, re_)) for _ in range(3)])), "len_mean": float(np.mean([c(matched(sl, cnn, rl)) for _ in range(3)]))}
        t = ts[(ts.mirna == mi) & ts.gene.isin(utrg)].groupby("gene").context_pp.min().nsmallest(K); row["ts_n"] = len(t)
        if len(t) >= 20:
            row["ts"] = c(list(t.index)); row["ts_rand_mean"] = float(np.mean([c(rb.sample(utrg, len(t))) for _ in range(3)]))
        per[mi] = row
    tsr = [v for v in per.values() if "ts" in v]
    w, l = sum(v["ts"] > v["ts_rand_mean"] for v in tsr), sum(v["ts"] < v["ts_rand_mean"] for v in tsr)
    g1 = dict(signtest(w, l), p_two_sided=min(1.0, 2 * signtest(max(w, l), min(w, l))["p_one_sided"]))
    prev = len(set(d.gene) & ess) / d.gene.nunique()
    J = {"tool": "DepMap 24Q4 Public CRISPRInferredCommonEssentials (figshare)", "n_essential": len(ess), "panel": panel, "K": K,
         "gates": {"G1_ts_vs_random": g1, "G1_n": len(tsr), "G1_pass": bool(g1["p_two_sided"] < 0.05), "G2_prevalence": prev, "G2_pass": bool(0.02 <= prev <= 0.20)}}
    for name, k in (("E1_vs_uniform", "uniform_mean"), ("E2_vs_expr", "expr_mean"), ("E3_vs_utrlen", "len_mean")):
        r = signtest(sum(v["cnn"] < v[k] for v in per.values()), sum(v["cnn"] > v[k] for v in per.values())); nt = r["wins"] + r["losses"]
        pc, pk = sum(v["cnn"] for v in per.values()), sum(v[k] for v in per.values())
        J[name] = dict(r, non_tied=nt, pooled_cnn=pc, pooled_comparator=pk, pooled_ratio=pc / pk if pk else None,
                       verdict="INCONCLUSIVE" if nt < 15 else ("CONFIRMED" if r["p_one_sided"] < 0.05 else "FALSIFIED"))
    J["per_mirna"] = per
    json.dump(J, open(out, "w"), indent=1)
    print(json.dumps({k: J[k] for k in ["n_essential", "gates", "E1_vs_uniform", "E2_vs_expr", "E3_vs_utrlen"]}))


if __name__ == "__main__":
    main(*sys.argv[1:])
