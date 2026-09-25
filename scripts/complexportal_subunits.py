"""PRE-REGISTERED Complex Portal subunit test (tool 36).

Label: EBI Complex Portal human complexes (ftp.ebi.ac.uk/pub/databases/intact/complex/current/
complextab/9606.tsv; counted once). UniProt accessions in "Expanded participant list" are
mapped to genes with data/uniprot/human_reviewed.tsv (tool 29); gene in any complex =
"subunit gene". Rationale: the dosage-balance hypothesis says complex subunits are
dosage-sensitive, so the dosage finding (tools 30-32) predicts CNN depletion of subunits.
Committed BEFORE running.
Panel: 13 discovery + 30 held-out miRNAs; CLIP universe; K=200. Comparators (3 draws, seeds
f"{mi}-cp-<arm>"): uniform, expression-matched, UTR-length-matched.
Positive control: TargetScan conserved top-200 vs 3 uniform random UTR genes; direction not
assumed (conserved targets are thought to avoid some housekeeping complexes), so G1 asks
for a two-sided sign test p<0.05.
Gates: G1 as above; G2 accession-to-gene mapping >= 90% of participant accessions;
G3 subunit prevalence in pooled universe between 3% and 30%.
Hypotheses (CNN fewer subunit genes), one-sided sign tests over 43 miRNAs, ties excluded,
alpha 0.05: X1 vs uniform; X2 vs expression-matched; X3 vs UTR-length-matched.
"""
import json, random, re, sys
import numpy as np, pandas as pd
sys.path.insert(0, "scripts")
from gprofiler_enrichment import PANEL, K, signtest
from string_exprmatched import matched

ACC = re.compile(r"\b([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})\b")


def main(out="results/complexportal_subunits.json"):
    cp = pd.read_csv("data/complexportal/9606.tsv", sep="\t", dtype=str)
    accs = set(a for x in cp["Expanded participant list"].fillna("") for a in ACC.findall(x))
    up = pd.read_csv("data/uniprot/human_reviewed.tsv", sep="\t").dropna(subset=["Gene Names (primary)"])
    a2g = dict(zip(up.Entry, up["Gene Names (primary)"]))
    mapped = {a2g[a] for a in accs if a in a2g}; map_rate = sum(a in a2g for a in accs) / len(accs)
    panel = PANEL + json.load(open("results/gnomad_replication.json"))["panel"]
    d = pd.read_csv("results/clip_rows_136.csv", usecols=["mirna", "gene", "min", "len"]).drop_duplicates(["mirna", "gene"])
    d = d[d.mirna.isin(panel)]
    expr = pd.read_csv("results/hpa_hek293_ntpm.tsv", sep="\t").groupby("Gene name")["nTPM"].max()
    utrg = sorted(set(pd.read_csv("data/human_utrs.tsv", sep="\t", usecols=[1]).iloc[:, 0]))
    ts = pd.read_csv("data/human_sites.tsv", sep="\t", usecols=["gene", "mirna", "context_pp"])
    c = lambda gs: int(sum(g in mapped for g in gs))
    per = {}
    for mi in panel:
        s = d[d.mirna == mi].reset_index(drop=True); U = s.gene.tolist(); cnn = s.nsmallest(K, "min").gene.tolist()
        se, sl = s.assign(x=s.gene.map(expr)), s.assign(x=s["len"].astype(float))
        ru, re_, rl, rb = (random.Random(f"{mi}-cp-{a}") for a in ("unif", "expr", "len", "utr"))
        row = {"universe": len(U), "cnn": c(cnn), "uniform_mean": float(np.mean([c(ru.sample(U, K)) for _ in range(3)])),
               "expr_mean": float(np.mean([c(matched(se, cnn, re_)) for _ in range(3)])), "len_mean": float(np.mean([c(matched(sl, cnn, rl)) for _ in range(3)]))}
        t = ts[(ts.mirna == mi) & ts.gene.isin(utrg)].groupby("gene").context_pp.min().nsmallest(K); row["ts_n"] = len(t)
        if len(t) >= 20:
            row["ts"] = c(list(t.index)); row["ts_rand_mean"] = float(np.mean([c(rb.sample(utrg, len(t))) for _ in range(3)]))
        per[mi] = row
    tsr = [v for v in per.values() if "ts" in v]
    w, l = sum(v["ts"] > v["ts_rand_mean"] for v in tsr), sum(v["ts"] < v["ts_rand_mean"] for v in tsr)
    g1 = dict(signtest(w, l), p_two_sided=min(1.0, 2 * signtest(max(w, l), min(w, l))["p_one_sided"]))
    prev = len(set(d.gene) & mapped) / d.gene.nunique()
    J = {"tool": "EBI Complex Portal (complextab 9606.tsv)", "n_complexes": int(len(cp)), "n_participant_accessions": len(accs), "n_subunit_genes": len(mapped),
         "panel": panel, "K": K,
         "gates": {"G1_ts_vs_random": g1, "G1_n": len(tsr), "G1_pass": bool(g1["p_two_sided"] < 0.05), "G2_map_rate": map_rate, "G2_pass": bool(map_rate >= 0.90),
                   "G3_prevalence": prev, "G3_pass": bool(0.03 <= prev <= 0.30)}}
    for name, k in (("X1_vs_uniform", "uniform_mean"), ("X2_vs_expr", "expr_mean"), ("X3_vs_utrlen", "len_mean")):
        r = signtest(sum(v["cnn"] < v[k] for v in per.values()), sum(v["cnn"] > v[k] for v in per.values()))
        pc, pk = sum(v["cnn"] for v in per.values()), sum(v[k] for v in per.values())
        J[name] = dict(r, pooled_cnn=pc, pooled_comparator=pk, pooled_ratio=pc / pk, verdict="CONFIRMED" if r["p_one_sided"] < 0.05 else "FALSIFIED")
    J["per_mirna"] = per
    json.dump(J, open(out, "w"), indent=1)
    print(json.dumps({k: J[k] for k in ["n_complexes", "n_subunit_genes", "gates", "X1_vs_uniform", "X2_vs_expr", "X3_vs_utrlen"]}))


if __name__ == "__main__":
    main(*sys.argv[1:])
