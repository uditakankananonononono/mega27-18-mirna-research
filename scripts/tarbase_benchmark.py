"""Second independent validation: TarBase v9.0 (DIANA; Homo_sapiens_TarBase-v9.tsv.gz,
4,724,537 human interaction records) as labels on the SAME (miRNA, gene) rows used for
the miRTarBase benchmark (results/mirtarbase_rows.csv; gene universe and CNN scores fixed).

Label sets (row = 1 if the pair has >=1 TarBase record of that evidence class, else 0):
  any        - any experimental method
  low_yield  - targeted methods (luciferase reporter, western blot, qPCR, ELISA, IHC, 3LIFE,
               Biotin-qPCR, 2D-DIGE, SILAC/pSILAC-free targeted assays excluded below)
  clip       - AGO-CLIP / chimeric methods (HITS-CLIP, PAR-CLIP, CLASH, qCLASH,
               Chimeric fragments, AGO-IP, IMPACT-Seq)
Predictors: our DuplexCNN (-min context++), miRDB v6.0 gene score (max over RefSeq
transcripts, 0 if unlisted), seed-site count (8mer + 7mer-A1; same enumeration as
scripts/mirtarbase_sitecount_control.py).
Caveats: TarBase CLIP records overlap ENCORI CLIP data used in our falsification analyses
(not in CNN training); miRDB MirTarget was trained on CLIP/array data (leakage in its favour).
Output: results/tarbase_benchmark.json
"""
import csv, gzip, json, sys
from collections import defaultdict
from pathlib import Path
import numpy as np
from sklearn.metrics import roc_auc_score
from scipy.stats import wilcoxon

R = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(R / "src"))
from mirna.data import load_utrs, load_mirnas
from mirna.cli import rc

LOW = {"Luciferase Reporter Assay", "Western Blot", "qPCR", "ELISA", "Immunohistochemistry",
       "3LIFE", "Biotin-qPCR", "2D-DIGE"}
CLIP = {"HITS-CLIP", "PAR-CLIP", "CLASH", "qCLASH", "Chimeric fragments", "AGO-IP", "IMPACT-Seq"}

rows = list(csv.DictReader(open(R / "results/mirtarbase_rows.csv")))
want = {(r["mirna"], r["gene"]) for r in rows}
wmir = {r["mirna"] for r in rows}
ev = defaultdict(set)
n_rec = 0
with gzip.open(R / "data/tarbase/Homo_sapiens_TarBase-v9.tsv.gz", "rt") as fh:
    next(fh)
    for line in fh:
        n_rec += 1
        f = line.split("\t", 14)
        if f[1] not in wmir:
            continue
        k = (f[1], f[3])
        if k in want:
            m = f[12]
            ev[k].add("any")
            if m in LOW: ev[k].add("low_yield")
            if m in CLIP: ev[k].add("clip")

# miRDB
m2s = json.load(open(R / "data/mirdb/refseq2symbol.json"))
mirdb = defaultdict(dict)
for line in open(R / "data/mirdb/rows_mirs.tsv"):
    mi, nm, sc = line.rstrip("\n").split("\t")
    g = m2s.get(nm)
    if g:
        mirdb[mi][g] = max(mirdb[mi].get(g, 0.0), float(sc))

# site count
utrs = load_utrs(R / "data/human_utrs.tsv")
gene2tids = defaultdict(list)
with open(R / "data/human_utrs.tsv") as fh:
    next(fh)
    for line in fh:
        t, g, _ = line.rstrip("\n").split("\t", 2)
        gene2tids[g].append(t.split(".")[0])
utr_by = {k.split(".")[0]: v for k, v in utrs.items()}
mirs = load_mirnas(R / "data/miR_Family_Info.txt")
def seed_count(seq, g):
    s8, s7 = rc(seq[1:8]), rc(seq[1:7]); n = 0
    for t in gene2tids.get(g, []):
        u = utr_by.get(t, "")
        i = u.find(s8)
        while i >= 0: n += 1; i = u.find(s8, i + 1)
        i = u.find(s7)
        while i >= 0:
            if i + 6 < len(u) and u[i + 6] == "A" and u[i:i + 7] != s8: n += 1
            i = u.find(s7, i + 1)
    return n

rows = [r for r in rows if r["mirna"] in mirdb and r["mirna"] in mirs]
cnn = np.array([-float(r["min"]) for r in rows])
md = np.array([mirdb[r["mirna"]].get(r["gene"], 0.0) for r in rows])
sc = np.array([seed_count(mirs[r["mirna"]], r["gene"]) for r in rows], float)
mtb = np.array([int(r["label"]) for r in rows])
by = defaultdict(list)
for i, r in enumerate(rows): by[r["mirna"]].append(i)
out = {"tarbase_records_scanned": n_rec, "n_rows": len(rows), "n_mirnas": len(by),
       "mirtarbase_pos": int(mtb.sum()),
       "predictors": ["cnn", "mirdb", "site_count"], "labels": {}}
for lab in ["any", "low_yield", "clip"]:
    y = np.array([int(lab in ev.get((r["mirna"], r["gene"]), ())) for r in rows])
    res = {"n_pos": int(y.sum()),
           "overlap_with_mirtarbase_pos": int(((y == 1) & (mtb == 1)).sum()),
           "tarbase_pos_among_mirtarbase_neg": int(((y == 1) & (mtb == 0)).sum()),
           "pooled": {n: float(roc_auc_score(y, v)) for n, v in [("cnn", cnn), ("mirdb", md), ("site_count", sc)]}}
    per = {}
    for mi, idx in by.items():
        yy = y[idx]
        if yy.sum() < 5 or yy.sum() == len(yy): continue
        per[mi] = {n: float(roc_auc_score(yy, v[idx])) for n, v in [("cnn", cnn), ("mirdb", md), ("site_count", sc)]}
    if len(per) >= 5:
        a = {n: np.array([p[n] for p in per.values()]) for n in ["cnn", "mirdb", "site_count"]}
        res["per_mirna"] = {"n_mirnas": len(per),
            **{f"median_{n}": float(np.median(a[n])) for n in a},
            "mirdb_beats_cnn": int((a["mirdb"] > a["cnn"]).sum()),
            "cnn_beats_site_count": int((a["cnn"] > a["site_count"]).sum()),
            "wilcoxon_mirdb_vs_cnn_p": float(wilcoxon(a["mirdb"], a["cnn"]).pvalue),
            "wilcoxon_cnn_vs_site_count_p": float(wilcoxon(a["cnn"], a["site_count"]).pvalue)}
    out["labels"][lab] = res
json.dump(out, open(R / "results/tarbase_benchmark.json", "w"), indent=1)
print(json.dumps(out, indent=1))
