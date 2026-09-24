"""Falsification test for the non-conserved-site discovery (item 18).

Hypothesis H1: the DuplexCNN+context score of NON-CONSERVED seed matches predicts
independent AGO-CLIP support (ENCORI/starBase, Zhou et al. Nat Methods 2026)
beyond naive covariates (8mer/7mer counts, UTR length, miRNA identity).
H0 (falsifier): adding the CNN score does not raise cross-validated AUROC over the
covariate-only model (bootstrap 95% CI of delta-AUROC includes 0, and the
within-miRNA permutation null matches the real CNN columns).

Unit = (miRNA, gene). Genes with ANY TargetScan conserved site for a miRNA sharing
the same seed are excluded, so the CLIP label cannot come from a conserved site.
"""
import csv, json, os, sys, time, urllib.request
from pathlib import Path
import numpy as np
import torch
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from mirna.data import one_hot, load_utrs, site_features, build, gene_split, load_mirnas, WIN, MIR_LEN
from mirna.model import DuplexCNN

D = ROOT / "data"
PANEL = ["hsa-miR-7-5p", "hsa-let-7a-5p", "hsa-miR-21-5p", "hsa-miR-17-5p", "hsa-miR-19a-3p",
         "hsa-miR-92a-3p", "hsa-miR-16-5p", "hsa-miR-24-3p", "hsa-miR-26a-5p", "hsa-miR-27a-3p",
         "hsa-miR-29a-3p", "hsa-miR-30a-5p", "hsa-miR-103a-3p", "hsa-miR-124-3p", "hsa-miR-125b-5p",
         "hsa-miR-130a-3p", "hsa-miR-137", "hsa-miR-142-3p", "hsa-miR-155-5p", "hsa-miR-181a-5p",
         "hsa-miR-199a-5p", "hsa-miR-200c-3p", "hsa-miR-221-3p", "hsa-miR-23a-3p", "hsa-miR-1-3p"]
if os.environ.get("CLIP_PANEL_JSON"):
    _ext = json.load(open(os.environ["CLIP_PANEL_JSON"]))
    PANEL = sorted(set(PANEL) | set(_ext.get("kept_new", {})))
OUT_PATH = os.environ.get("CLIP_OUT", str(ROOT / "results" / "clip_falsification.json"))
COMP = str.maketrans("ACGU", "UGCA")
def rc(s): return s.translate(COMP)[::-1]


def encori(mir):
    f = D / "encori" / f"{mir}.tsv"
    if not f.exists():
        url = ("https://rnasysu.com/encori/api/miRNATarget/?assembly=hg38&geneType=mRNA&miRNA="
               f"{mir}&clipExpNum=1&degraExpNum=0&pancancerNum=0&programNum=0&program=None&target=all&cellType=all")
        for _ in range(3):
            try:
                f.write_bytes(urllib.request.urlopen(url, timeout=90).read()); break
            except Exception as e:
                print("retry", mir, e, flush=True); time.sleep(5)
    genes = set()
    if f.exists():
        for line in f.read_text().splitlines():
            if line.startswith("#") or line.startswith("miRNAid"):
                continue
            p = line.split("\t")
            if len(p) > 11 and p[1] == mir:
                genes.add(p[3])
    return genes


def score(net, utrs, seq, cands, mu, sd):
    sc = np.zeros(len(cands), np.float32)
    with torch.no_grad():
        for lo in range(0, len(cands), 2048):
            ch = cands[lo:lo + 2048]
            Xu = torch.from_numpy(np.stack([one_hot("".join(
                utrs[t][k] if 0 <= k < len(utrs[t]) else "N"
                for k in range(p + 3 - WIN // 2, p + 3 - WIN // 2 + WIN)), WIN) for t, p, _ in ch]))
            Xm = torch.from_numpy(np.stack([one_hot(seq, MIR_LEN)] * len(ch)))
            Xs = torch.from_numpy(np.array([s for *_, s in ch], dtype=np.int64))
            Xf = torch.from_numpy(((np.stack([site_features(utrs[t], p + 1, p + 7)
                                              for t, p, _ in ch]) - mu) / sd).astype(np.float32))
            sc[lo:lo + len(ch)] = net(Xu, Xm, Xs, Xf).numpy()
    return sc


def main():
    fam = {}
    with open(D / "miR_Family_Info.txt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["Species ID"] == "9606":
                fam[row["MiRBase ID"]] = row["Mature sequence"].upper().replace("T", "U")
    mirs = {m: fam[m] for m in PANEL if m in fam}
    print("panel present:", len(mirs), "missing:", [m for m in PANEL if m not in fam], flush=True)
    seed_of = {m: s[1:8] for m, s in fam.items() if len(s) >= 8}
    cons_genes, cons = {}, {}
    with open(D / "human_sites.tsv") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            s = seed_of.get(row["mirna"])
            if s:
                cons_genes.setdefault(s, set()).add(row["gene"])
            cons.setdefault(row["transcript"], []).append((int(row["start"]), int(row["end"])))
    utrs = load_utrs(D / "human_utrs.tsv")
    tid2gene = {}
    with open(D / "human_utrs.tsv") as fh:
        next(fh)
        for line in fh:
            t, g, _ = line.rstrip("\n").split("\t", 2)
            tid2gene[t.split(".")[0]] = g

    def overlaps(t, lo, hi):
        return any(not (hi < a - 1 or lo > b - 1) for a, b in cons.get(t, []))

    all_mirs = load_mirnas(D / "miR_Family_Info.txt")
    U, M, ST, y, genes, meta, F = build(D / "human_sites.tsv", utrs, all_mirs, max_rows=80000, seed=0)
    tr, _ = gene_split(genes)
    mu, sd = F[tr].mean(0), F[tr].std(0) + 1e-6
    del U, M, ST, y, F
    net = DuplexCNN(n_feat=4)
    net.load_state_dict(torch.load(ROOT / "results" / "duplex_cnn.pt"))
    net.eval()

    rows, per_mir = [], {}
    for mi, seq in mirs.items():
        clip = encori(mi)
        excl = cons_genes.get(seq[1:8], set())
        seed8, seed7 = rc(seq[1:8]), rc(seq[1:7])
        cands = []
        for t, u in utrs.items():
            g = tid2gene.get(t)
            if g is None or g in excl:
                continue
            i = u.find(seed8)
            while i >= 0:
                if not overlaps(t, i, i + 6):
                    cands.append((t, i, 3 if (i + 7 < len(u) and u[i + 7] == "A") else 2))
                i = u.find(seed8, i + 1)
            i = u.find(seed7)
            while i >= 0:
                if i + 6 < len(u) and u[i + 6] == "A" and u[i:i + 7] != seed8 and not overlaps(t, i, i + 5):
                    cands.append((t, i, 2))
                i = u.find(seed7, i + 1)
        if not cands or not clip:
            print(mi, "skip", len(cands), len(clip), flush=True)
            continue
        sc = score(net, utrs, seq, cands, mu, sd)
        agg = {}
        for (t, p, s), v in zip(cands, sc):
            a = agg.setdefault(tid2gene[t], {"min": 0.0, "sum": 0.0, "n8": 0, "n7": 0, "len": 0})
            a["min"] = min(a["min"], float(v)); a["sum"] += float(v)
            a["n8" if s == 3 else "n7"] += 1
            a["len"] = max(a["len"], len(utrs[t]))
        yv = np.array([int(g in clip) for g in agg])
        for (g, a), l in zip(agg.items(), yv):
            rows.append((mi, g, a["min"], a["sum"], a["n8"], a["n7"], a["len"], int(l)))
        mins = np.array([a["min"] for a in agg.values()])
        per_mir[mi] = {"genes": len(agg), "clip_pos": int(yv.sum()), "clip_genes_total": len(clip),
                       "auroc_cnn_min": round(float(roc_auc_score(yv, -mins)), 4) if 0 < yv.sum() < len(yv) else None}
        print(mi, per_mir[mi], flush=True)

    mir_idx = {m: k for k, m in enumerate(sorted({r[0] for r in rows}))}
    y = np.array([r[7] for r in rows])
    mids = np.array([mir_idx[r[0]] for r in rows])
    base = np.column_stack([np.log1p([r[4] for r in rows]), np.log1p([r[5] for r in rows]),
                            np.log([max(r[6], 1) for r in rows]), np.eye(len(mir_idx), dtype=np.float32)[mids]]).astype(np.float32)
    cnn = np.column_stack([[r[2] for r in rows], [r[3] for r in rows]]).astype(np.float32)

    def cv_pred(X):
        p = np.zeros(len(y))
        for a, b in StratifiedKFold(5, shuffle=True, random_state=0).split(X, y):
            p[b] = LogisticRegression(max_iter=3000).fit(X[a], y[a]).predict_proba(X[b])[:, 1]
        return p

    pb, pf = cv_pred(base), cv_pred(np.column_stack([base, cnn]))
    rng = np.random.default_rng(0)
    deltas = []
    for _ in range(1000):
        i = rng.integers(0, len(y), len(y))
        deltas.append(roc_auc_score(y[i], pf[i]) - roc_auc_score(y[i], pb[i]))
    perm = []
    for _ in range(20):
        c2 = cnn.copy()
        for m in range(len(mir_idx)):
            ix = np.where(mids == m)[0]
            c2[ix] = cnn[rng.permutation(ix)]
        perm.append(roc_auc_score(y, cv_pred(np.column_stack([base, c2]))))
    a_f = float(roc_auc_score(y, pf))
    out = {"hypothesis": __doc__.strip(), "source": "ENCORI miRNATarget API, clipExpNum>=1, hg38",
           "n_pairs": int(len(y)), "n_clip_pos": int(y.sum()),
           "auroc_base": round(float(roc_auc_score(y, pb)), 4), "auroc_base_plus_cnn": round(a_f, 4),
           "delta_auroc_mean": round(float(np.mean(deltas)), 4),
           "delta_ci95": [round(float(np.percentile(deltas, 2.5)), 4), round(float(np.percentile(deltas, 97.5)), 4)],
           "perm_null_auroc_mean": round(float(np.mean(perm)), 4), "perm_null_auroc_max": round(float(np.max(perm)), 4),
           "perm_p": (1 + sum(p >= a_f for p in perm)) / (1 + len(perm)),
           "per_mirna": per_mir}
    out["verdict"] = ("H0 rejected: CNN adds CLIP signal" if out["delta_ci95"][0] > 0 and a_f > out["perm_null_auroc_max"]
                      else "H0 NOT rejected: CNN adds no CLIP signal beyond covariates (negative kept)")
    json.dump(out, open(OUT_PATH, "w"), indent=1)
    print(json.dumps({k: v for k, v in out.items() if k not in ("per_mirna", "hypothesis")}, indent=1))


if __name__ == "__main__":
    main()
