"""ENCODE microRNA-seq abundance vs per-miRNA predictability (tool 22; Ensembl REST for IDs).

PRE-REGISTERED (written before any correlation was computed):
Data: 21 ENCODE GRCh38 "microRNA quantifications" TSVs (12 cell lines, listed in
data/encode_mirna/files.tsv, each an ENCFF accession). Per file: unstranded counts,
CPM over gene rows, log2(CPM+1). Gene IDs -> symbols via Ensembl REST POST /lookup/id.
Mature name -> precursor symbols (hsa-miR-21-5p -> MIR21; let-7a -> MIRLET7A1..3;
base ending in a digit takes only '-N' suffixes, base ending in a letter takes 'N'
or '-N'); precursor CPMs summed. Precursor counts cannot separate 5p/3p arms
(limitation). Cell-line value = mean over replicate files; abundance = median over lines.
Outcomes: per-miRNA rows of results/clip_per_mirna_delta.json (131 miRNAs:
base_auroc, cnn_auroc, delta, n_pos).
Gates:
 G0 every file has >= 1500 gene rows.
 G1 miR-21-5p abundance is in the panel's top decile.
 G2 miR-122-5p (liver-specific) is highest in HepG2 among the 12 lines.
 G3 >= 90% of the 131 panel miRNAs map to >= 1 precursor.
Hypotheses (one-sided, alpha 0.05):
 H1 Spearman(abundance, cnn_auroc) > 0.
 H2 Spearman(abundance, delta) > 0.
 H3 partial Spearman(abundance, cnn_auroc | log n_pos) > 0 (CLIP+ count confound).
"""
import json, re, sys, os
import numpy as np, pandas as pd, requests
from scipy.stats import spearmanr, rankdata, t as tdist

DIR = "data/encode_mirna"


def load_files():
    f = pd.read_csv(f"{DIR}/files.tsv", sep="\t", header=None, names=["line", "acc", "exp", "href", "size", "reps"])
    mats = {}
    for _, r in f.iterrows():
        d = pd.read_csv(f"{DIR}/{r.acc}.tsv", sep="\t", header=None, usecols=[0, 1], names=["gid", "c"])
        d = d[d.gid.str.startswith("ENSG")]
        d["gid"] = d.gid.str.split(".").str[0]
        mats[r.acc] = d.set_index("gid").c
    return f, mats


def ensembl_symbols(ids, cache="results/encode_mirna_ensg_symbols.json"):
    if os.path.exists(cache):
        return json.load(open(cache))
    out = {}
    ids = sorted(ids)
    import time
    for i in range(0, len(ids), 200):
        for k in range(5):
            r = requests.post("https://rest.ensembl.org/lookup/id", json={"ids": ids[i:i + 200]},
                              headers={"Content-Type": "application/json", "Accept": "application/json"}, timeout=180)
            if r.ok: break
            time.sleep(3 * (k + 1))
        r.raise_for_status()
        for k, v in r.json().items():
            if v and v.get("display_name"):
                out[k] = v["display_name"]
    json.dump(out, open(cache, "w"), indent=0, sort_keys=True)
    return out


def precursors(mature, symbols):
    b = re.sub(r"^hsa-", "", mature)
    b = re.sub(r"-(5p|3p)$", "", b)
    b = ("MIRLET" + b[4:] if b.startswith("let-") else "MIR" + b[4:]).upper()
    pat = re.compile("^" + re.escape(b) + (r"(-\d+)?$" if b[-1].isdigit() else r"(-?\d+)?$"))
    return sorted(s for s in symbols if pat.match(s))


def partial_spearman(x, y, z):
    rx, ry, rz = rankdata(x), rankdata(y), rankdata(z)
    res = lambda a: a - np.polyval(np.polyfit(rz, a, 1), rz)
    r = np.corrcoef(res(rx), res(ry))[0, 1]
    n = len(x); tt = r * np.sqrt((n - 3) / (1 - r * r))
    return float(r), float(tdist.sf(tt, n - 3))


BLOOD = {"GM12878", "K562", "HL-60"}


def posthoc_controls(B):
    """POST-HOC replacement for failed G2 (chosen after G2 failed; not pre-registered):
    MIR142 top-3 lines all hematopoietic; MIR302A, MIR302B, MIR367 each highest in H1 (ESC)."""
    out = {"note": "post-hoc, chosen after pre-registered G2 (miR-122 in HepG2) failed; not pre-registered"}
    top3 = set(B.loc["MIR142"].nlargest(3).index)
    out["MIR142_top3"] = sorted(top3); out["MIR142_pass"] = top3 == BLOOD
    for s in ["MIR302A", "MIR302B", "MIR367"]:
        out[s + "_max_line"] = str(B.loc[s].idxmax()); out[s + "_pass"] = out[s + "_max_line"] == "H1"
    out["MIR122_all_zero"] = bool((B.loc["MIR122"] == 0).all())
    out["all_pass"] = out["MIR142_pass"] and all(out[s + "_pass"] for s in ["MIR302A", "MIR302B", "MIR367"])
    return out


def main(out="results/encode_mirna_abundance.json"):
    f, mats = load_files()
    g0 = all(len(m) >= 1500 for m in mats.values())
    sym = ensembl_symbols(set().union(*[set(m.index) for m in mats.values()]))
    cpm = {a: np.log2(m / m.sum() * 1e6 + 1) for a, m in mats.items()}
    line_vals = {}
    for line, grp in f.groupby("line", sort=False):
        line_vals[line] = pd.concat([cpm[a] for a in grp.acc], axis=1).mean(axis=1)
    L = pd.DataFrame(line_vals)  # gene x line (log2 CPM+1)
    bysym = L.groupby(L.index.map(lambda g: sym.get(g, g))).sum()
    rows = json.load(open("results/clip_per_mirna_delta.json"))["rows"]
    tab = []
    for r in rows:
        pre = precursors(r["mirna"], bysym.index)
        if not pre:
            tab.append(dict(r, precursors=[], abundance=None)); continue
        v = bysym.loc[pre].sum(axis=0)
        tab.append(dict(r, precursors=pre, abundance=float(v.median()), per_line={k: float(x) for k, x in v.items()}))
    T = pd.DataFrame([t for t in tab if t["abundance"] is not None])
    thr = T.abundance.quantile(0.9)
    m21 = T[T.mirna == "hsa-miR-21-5p"]
    g1 = bool(len(m21) and m21.abundance.iloc[0] >= thr)
    m122 = next((t for t in tab if t["mirna"] == "hsa-miR-122-5p" and t["abundance"] is not None), None)
    g2 = bool(m122 and max(m122["per_line"], key=m122["per_line"].get) == "HepG2")
    g3 = len(T) / len(rows) >= 0.90
    s1 = spearmanr(T.abundance, T.cnn_auroc, alternative="greater")
    s2 = spearmanr(T.abundance, T.delta, alternative="greater")
    s0 = spearmanr(T.abundance, T.base_auroc)
    sn = spearmanr(T.abundance, T.n_pos)
    p3 = partial_spearman(T.abundance.values, T.cnn_auroc.values, np.log(T.n_pos.values))
    vd = lambda p: "CONFIRMED" if p < 0.05 else "FALSIFIED"
    J = {"tool": "ENCODE portal microRNA-seq quantifications + Ensembl REST lookup/id",
         "files": f[["line", "acc", "exp"]].to_dict("records"), "n_files": len(f), "n_lines": f.line.nunique(),
         "n_panel": len(rows), "n_mapped": len(T),
         "gates": {"G0_rows_ok": g0, "G1_miR21_top_decile": g1, "G1_miR21_abundance": float(m21.abundance.iloc[0]) if len(m21) else None,
                   "G1_top_decile_threshold": float(thr), "G2_miR122_max_line": max(m122["per_line"], key=m122["per_line"].get) if m122 else None,
                   "G2_miR122_per_line": m122["per_line"] if m122 else None, "G2_pass": g2,
                   "G3_mapped_frac": len(T) / len(rows), "G3_pass": g3},
         "H1_abundance_vs_cnn_auroc": {"rho": float(s1.statistic), "p_one_sided": float(s1.pvalue), "n": len(T), "verdict": vd(s1.pvalue)},
         "H2_abundance_vs_delta": {"rho": float(s2.statistic), "p_one_sided": float(s2.pvalue), "n": len(T), "verdict": vd(s2.pvalue)},
         "H3_partial_given_log_npos": {"rho": p3[0], "p_one_sided": p3[1], "n": len(T), "verdict": vd(p3[1])},
         "descriptive": {"abundance_vs_base_auroc": {"rho": float(s0.statistic), "p_two_sided": float(s0.pvalue)},
                         "abundance_vs_n_pos": {"rho": float(sn.statistic), "p_two_sided": float(sn.pvalue)}},
         "posthoc_G2b_tissue_controls": posthoc_controls(bysym),
         "unmapped": [t["mirna"] for t in tab if t["abundance"] is None],
         "per_mirna": [{k: t.get(k) for k in ["mirna", "precursors", "abundance", "cnn_auroc", "base_auroc", "delta", "n_pos"]} for t in tab]}
    json.dump(J, open(out, "w"), indent=1)
    print(json.dumps({k: J[k] for k in ["n_files", "n_lines", "n_mapped", "gates", "H1_abundance_vs_cnn_auroc", "H2_abundance_vs_delta", "H3_partial_given_log_npos", "descriptive", "unmapped"]}, indent=0))


if __name__ == "__main__":
    main(*sys.argv[1:])
