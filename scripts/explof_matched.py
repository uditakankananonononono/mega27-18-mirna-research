"""PRE-REGISTERED coding-length (expected-LoF) control for the LoF-count labels (follow-up to tools 30, 39).

Explanation under test (genebayes_sec): LOEUF and s_het are weakly informative for genes with few
expected LoF variants (short coding sequences); if CNN candidates are such genes, the LoF-count
labels would show depletion beyond UTR length as an artefact. Proxy for coding length: gnomAD
v2.1.1 exp_lof (expected LoF count, canonical transcript). Committed BEFORE running.
Panel: 13 discovery + 30 held-out miRNAs; CLIP universe restricted to genes with LOEUF, exp_lof and
s_het; K=200. Comparator: 3 exp_lof-matched draws (decile matcher string_exprmatched.matched on
exp_lof; seeds f"{mi}-el-<k>").
Descriptive D1: CNN median exp_lof < universe median (count of miRNAs, sign test).
Hypotheses (the depletion survives coding-length matching), one-sided sign tests over 43, ties
excluded, alpha 0.05:
 L1 CNN median LOEUF > exp_lof-matched mean.
 L2 CNN median log10 s_het < exp_lof-matched mean.
If D1 holds and L1, L2 fail, the LoF-count split is attributed to coding length.
"""
import json, random, sys
import numpy as np, pandas as pd
sys.path.insert(0, "scripts")
from gprofiler_enrichment import PANEL, K, signtest
from string_exprmatched import matched
from genebayes_shet import load_hi as load_shet


def main(out="results/explof_matched.json"):
    g = pd.read_csv("data/gnomad/lof_v211.txt.bgz", sep="\t", compression="gzip", usecols=["gene", "oe_lof_upper", "exp_lof"]).dropna()
    g = g.sort_values("exp_lof", ascending=False).drop_duplicates("gene").set_index("gene")
    sh = load_shet()
    keep = set(g.index) & set(sh.index)
    panel = PANEL + json.load(open("results/gnomad_replication.json"))["panel"]
    d = pd.read_csv("results/clip_rows_136.csv", usecols=["mirna", "gene", "min"]).drop_duplicates(["mirna", "gene"])
    d = d[d.mirna.isin(panel) & d.gene.isin(keep)]
    mlo = lambda gs: float(np.median(g.oe_lof_upper.reindex(gs).values)); msh = lambda gs: float(np.median(sh.reindex(gs).values))
    per = {}
    for mi in panel:
        s = d[d.mirna == mi].reset_index(drop=True); s = s.assign(x=s.gene.map(g.exp_lof)); cnn = s.nsmallest(K, "min").gene.tolist()
        rng = random.Random(f"{mi}-el-0"); M = [matched(s, cnn, rng) for _ in range(3)]
        per[mi] = {"universe": len(s), "cnn_explof": float(g.exp_lof.reindex(cnn).median()), "universe_explof": float(s.x.median()),
                   "cnn_loeuf": mlo(cnn), "m_loeuf": float(np.mean([mlo(x) for x in M])), "m_explof": float(np.mean([g.exp_lof.reindex(x).median() for x in M])),
                   "cnn_shet": msh(cnn), "m_shet": float(np.mean([msh(x) for x in M]))}
    V = list(per.values())
    d1 = signtest(sum(v["cnn_explof"] < v["universe_explof"] for v in V), sum(v["cnn_explof"] > v["universe_explof"] for v in V))
    l1 = signtest(sum(v["cnn_loeuf"] > v["m_loeuf"] for v in V), sum(v["cnn_loeuf"] < v["m_loeuf"] for v in V))
    l2 = signtest(sum(v["cnn_shet"] < v["m_shet"] for v in V), sum(v["cnn_shet"] > v["m_shet"] for v in V))
    vd = lambda h: "CONFIRMED" if h["p_one_sided"] < 0.05 else "FALSIFIED"
    J = {"panel": panel, "K": K, "D1_cnn_shorter_cds": dict(d1, median_ratio=float(np.median([v["cnn_explof"] / v["universe_explof"] for v in V]))),
         "L1_loeuf_vs_explof_matched": dict(l1, median_diff=float(np.median([v["cnn_loeuf"] - v["m_loeuf"] for v in V])), verdict=vd(l1)),
         "L2_shet_vs_explof_matched": dict(l2, median_diff=float(np.median([v["cnn_shet"] - v["m_shet"] for v in V])), verdict=vd(l2)),
         "matching_check_median_explof_ratio": float(np.median([v["m_explof"] / v["cnn_explof"] for v in V])), "per_mirna": per}
    json.dump(J, open(out, "w"), indent=1)
    print(json.dumps({k: J[k] for k in J if k not in ("per_mirna", "panel")}))


if __name__ == "__main__":
    main(*sys.argv[1:])
