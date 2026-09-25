"""GTEx v8 tissue target-avoidance audit (tool 24).

PRE-REGISTERED (written before any group comparison was run):
Classic prediction: genes targeted by a tissue-specific miRNA are relatively
depleted in that tissue. Data: GTEx v8 gene median TPM per tissue (bulk file
GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz).
Relative expression rel_g(t) = log2(TPM_g,t + 1) - median_t' log2(TPM_g,t' + 1).
Pairs (8): miR-1-3p, miR-133b, miR-206 -> Muscle - Skeletal; miR-122-5p -> Liver;
miR-9-5p, miR-128-3p -> Brain - Cortex; miR-223-3p -> Whole Blood; miR-375 -> Pancreas.
Universe per miRNA = genes scored in results/clip_rows_136.csv that map to GTEx.
Sets: CNN top-200 (lowest 'min'); TargetScan conserved-site genes (data/human_sites.tsv,
top-200 by context++ within the universe; reference arm). Statistic per pair:
AUC_dep = P(rel of a set gene < rel of a non-set gene) (Mann-Whitney), 0.5 = none.
Gates:
 G1 ALB highest in Liver; G2 ACTA1 highest in Muscle - Skeletal;
 G3 median universe mapping to GTEx >= 0.90.
Hypotheses (one-sided sign tests across the 8 pairs, alpha 0.05):
 H1 CNN AUC_dep > 0.5.   H2 TargetScan AUC_dep > 0.5.
 H3 specificity: CNN AUC_dep in the matched tissue > mean CNN AUC_dep over the
    other 4 tissues of the design (controls generic UTR/expression bias).
"""
import json, sys
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
sys.path.insert(0, "scripts")
from gprofiler_enrichment import signtest

PAIRS = {"hsa-miR-1-3p": "Muscle - Skeletal", "hsa-miR-133b": "Muscle - Skeletal", "hsa-miR-206": "Muscle - Skeletal",
         "hsa-miR-122-5p": "Liver", "hsa-miR-9-5p": "Brain - Cortex", "hsa-miR-128-3p": "Brain - Cortex",
         "hsa-miR-223-3p": "Whole Blood", "hsa-miR-375": "Pancreas"}
TISSUES = sorted(set(PAIRS.values()))
K = 200


def load_gtex(path="data/gtex/gtex_median_tpm.gct.gz"):
    g = pd.read_csv(path, sep="\t", skiprows=2)
    g = g.drop(columns=["Name"]).groupby("Description").max()
    L = np.log2(g + 1)
    return g, L.sub(L.median(axis=1), axis=0)


def auc_dep(rel, inset):
    a, b = rel[inset], rel[~inset]
    u = mannwhitneyu(a, b, alternative="less")
    return {"auc_dep": float(1 - u.statistic / (len(a) * len(b))), "p_less": float(u.pvalue), "n_set": int(inset.sum()), "n_rest": int((~inset).sum())}


def main(out="results/gtex_target_avoidance.json"):
    tpm, REL = load_gtex()
    d = pd.read_csv("results/clip_rows_136.csv", usecols=["mirna", "gene", "min"])
    ts = pd.read_csv("data/human_sites.tsv", sep="\t", usecols=["gene", "mirna", "context_pp"])
    g1 = tpm.loc["ALB"].idxmax() == "Liver"; g2 = tpm.loc["ACTA1"].idxmax() == "Muscle - Skeletal"
    per, maps = {}, []
    for mi, tis in PAIRS.items():
        s = d[d.mirna == mi].drop_duplicates("gene")
        maps.append(s.gene.isin(REL.index).mean())
        s = s[s.gene.isin(REL.index)].reset_index(drop=True)
        cnn = set(s.nsmallest(K, "min").gene)
        t = ts[(ts.mirna == mi) & ts.gene.isin(s.gene)].groupby("gene").context_pp.min().nsmallest(K)
        tset = set(t.index)
        row = {"tissue": tis, "universe": len(s), "targetscan_n": len(tset)}
        for tt in TISSUES:
            rel = REL.loc[s.gene, tt].values
            row.setdefault("cnn", {})[tt] = auc_dep(rel, s.gene.isin(cnn).values)
            if len(tset) >= 20:
                row.setdefault("targetscan", {})[tt] = auc_dep(rel, s.gene.isin(tset).values)
        row["cnn_matched"] = row["cnn"][tis]["auc_dep"]
        row["cnn_mismatched_mean"] = float(np.mean([row["cnn"][x]["auc_dep"] for x in TISSUES if x != tis]))
        row["ts_matched"] = row["targetscan"][tis]["auc_dep"] if "targetscan" in row else None
        per[mi] = row
        print(mi, tis, len(s), round(row["cnn_matched"], 3), round(row["cnn_mismatched_mean"], 3), row["ts_matched"] and round(row["ts_matched"], 3), len(tset), flush=True)
    sg = lambda f: signtest(sum(f(v) > 0 for v in per.values() if f(v) is not None), sum(f(v) < 0 for v in per.values() if f(v) is not None))
    h1 = sg(lambda v: v["cnn_matched"] - 0.5)
    h2 = sg(lambda v: None if v["ts_matched"] is None else v["ts_matched"] - 0.5)
    h3 = sg(lambda v: v["cnn_matched"] - v["cnn_mismatched_mean"])
    vd = lambda h: "CONFIRMED" if h["p_one_sided"] < 0.05 else "FALSIFIED"
    J = {"tool": "GTEx Portal v8 bulk gene median TPM (adult-gtex storage)", "file": "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz",
         "K": K, "pairs": PAIRS,
         "gates": {"G1_ALB_liver": bool(g1), "G2_ACTA1_muscle": bool(g2), "G3_median_mapping": float(np.median(maps)), "G3_pass": bool(np.median(maps) >= 0.90)},
         "H1_cnn_depletion": dict(h1, verdict=vd(h1)), "H2_targetscan_depletion": dict(h2, verdict=vd(h2)),
         "H3_cnn_tissue_specificity": dict(h3, verdict=vd(h3)), "per_mirna": per}
    json.dump(J, open(out, "w"), indent=1)
    print(json.dumps({k: J[k] for k in ["gates", "H1_cnn_depletion", "H2_targetscan_depletion", "H3_cnn_tissue_specificity"]}))



def posthoc_targetscan(path="results/gtex_target_avoidance.json"):
    """POST-HOC (not pre-registered): pre-registered H2 is empty by construction because
    the CLIP universe excludes every gene with a conserved TargetScan site for the seed
    (scripts/clip_falsify.py). Reference arm here: universe = genes with a TargetScan
    3'UTR that map to GTEx; set = top-200 TargetScan conserved-site genes by context++."""
    J = json.load(open(path))
    tpm, REL = load_gtex()
    utr_genes = pd.read_csv("data/human_utrs.tsv", sep="\t", usecols=[1]).iloc[:, 0].drop_duplicates()
    U = pd.Series(sorted(set(utr_genes) & set(REL.index)))
    ts = pd.read_csv("data/human_sites.tsv", sep="\t", usecols=["gene", "mirna", "context_pp"])
    per = {}
    for mi, tis in PAIRS.items():
        t = ts[(ts.mirna == mi) & ts.gene.isin(U)].groupby("gene").context_pp.min().nsmallest(K)
        inset = U.isin(set(t.index)).values
        r = {tt: auc_dep(REL.loc[U, tt].values, inset) for tt in TISSUES}
        per[mi] = {"tissue": tis, "n_set": int(inset.sum()), "universe": len(U), "by_tissue": r,
                   "matched": r[tis]["auc_dep"], "mismatched_mean": float(np.mean([r[x]["auc_dep"] for x in TISSUES if x != tis]))}
        print(mi, tis, per[mi]["n_set"], round(per[mi]["matched"], 3), round(per[mi]["mismatched_mean"], 3), "%.2g" % r[tis]["p_less"], flush=True)
    w = lambda f: signtest(sum(f(v) > 0 for v in per.values()), sum(f(v) < 0 for v in per.values()))
    J["posthoc_targetscan_reference"] = {"note": "post-hoc, not pre-registered; replaces structurally empty H2",
        "per_mirna": per, "depletion": w(lambda v: v["matched"] - 0.5), "specificity": w(lambda v: v["matched"] - v["mismatched_mean"])}
    json.dump(J, open(path, "w"), indent=1)
    print(json.dumps({k: J["posthoc_targetscan_reference"][k] for k in ["depletion", "specificity"]}))


if __name__ == "__main__":
    main() if len(sys.argv) < 2 else posthoc_targetscan()
