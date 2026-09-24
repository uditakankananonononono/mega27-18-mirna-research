"""Independent literature-curated validation of the DuplexCNN (item 18).

miRTarBase 10.0 (Huang et al. NAR 2025) strong-evidence set (reporter assay,
miRTarBase_SE_R.csv) provides experimentally validated (miRNA, gene) pairs that
are INDEPENDENT of the ENCORI AGO-CLIP labels used in clip_falsification_136.

Design, per panel miRNA with >= MINPOS validated genes:
  positives = hsa SE validated genes with a 3'UTR in TargetScan UTR_Sequences
  negatives = genes with a UTR, NO miRTarBase 10.0 human MTI of ANY evidence
              class for this miRNA (hsa_MTI.csv full set), sampled <= NEG_MULT x
              positives (max 1500)
  gene score = min DuplexCNN score over 8mer/7mer-A1 seed-match windows across
               all transcripts of the gene (lower = stronger predicted site)
Metrics: per-miRNA AUROC(-min); pooled AUROC; within-miRNA z-scored Mann-Whitney;
label-shuffle permutation p (shuffles within miRNA). Seed-match recall reported:
validated genes with >=1 candidate window / all validated genes with a UTR.
Checkpointed per miRNA (rows CSV + per-miRNA JSONL); resume-safe.
"""
import csv, json, os, sys, gc
from pathlib import Path
import numpy as np
import torch
from sklearn.metrics import roc_auc_score
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from mirna.data import one_hot, load_utrs, site_features, load_mirnas, WIN, MIR_LEN
from mirna.model import DuplexCNN

D = ROOT / "data"
MTB = D / "mirtarbase"
MINPOS = 10
NEG_MULT = 10
MAXNEG = 1500
ROWS_CSV = str(ROOT / "results" / "mirtarbase_rows.csv")
PER_JSONL = str(ROOT / "results" / "mirtarbase_perm.jsonl")
OUT = str(ROOT / "results" / "mirtarbase_validation.json")

PANEL = ["hsa-miR-7-5p", "hsa-let-7a-5p", "hsa-miR-21-5p", "hsa-miR-17-5p", "hsa-miR-19a-3p",
         "hsa-miR-92a-3p", "hsa-miR-16-5p", "hsa-miR-24-3p", "hsa-miR-26a-5p", "hsa-miR-27a-3p",
         "hsa-miR-29a-3p", "hsa-miR-30a-5p", "hsa-miR-103a-3p", "hsa-miR-124-3p", "hsa-miR-125b-5p",
         "hsa-miR-130a-3p", "hsa-miR-137", "hsa-miR-142-3p", "hsa-miR-155-5p", "hsa-miR-181a-5p",
         "hsa-miR-199a-5p", "hsa-miR-200c-3p", "hsa-miR-221-3p", "hsa-miR-23a-3p", "hsa-miR-1-3p"]
_ext = json.load(open(ROOT / "results" / "panel_extension.json"))
PANEL = sorted(set(PANEL) | set(_ext.get("kept_new", {})))
COMP = str.maketrans("ACGU", "UGCA")
def rc(s): return s.translate(COMP)[::-1]

# ---- miRTarBase positives (strong evidence) and any-evidence cleaning set ----
pos = {}          # mirna -> set(gene)
with open(MTB / "miRTarBase_SE_R.csv", encoding="utf-8-sig") as fh:
    for row in csv.DictReader(fh):
        if row["Species (miRNA)"] == "hsa" and row["Species (Target Gene)"] == "hsa":
            pos.setdefault(row["miRNA"], set()).add(row["Target Gene"])
any_ev = {}       # mirna -> set(gene) with ANY human MTI evidence
with open(MTB / "hsa_MTI.csv", encoding="utf-8-sig") as fh:
    for row in csv.DictReader(fh):
        if row["Species (Target Gene)"] == "hsa":
            any_ev.setdefault(row["miRNA"], set()).add(row["Target Gene"])

fam = {}
with open(D / "miR_Family_Info.txt") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        if row["Species ID"] == "9606":
            fam[row["MiRBase ID"]] = row["Mature sequence"].upper().replace("T", "U")
mirs = {m: fam[m] for m in PANEL if m in fam}
print("panel:", len(mirs), "miRTarBase SE miRNAs:", len(pos), "overlap:",
      len(set(mirs) & set(pos)), flush=True)

utrs = load_utrs(D / "human_utrs.tsv")
tid2gene, gene2tids = {}, {}
with open(D / "human_utrs.tsv") as fh:
    next(fh)
    for line in fh:
        t, g, _ = line.rstrip("\n").split("\t", 2)
        tid2gene[t.split(".")[0]] = g
        gene2tids.setdefault(g, []).append(t.split(".")[0])

# feature normalisation: rebuild the training-row stats exactly as clip_falsify
from mirna.data import build, gene_split
all_mirs = load_mirnas(D / "miR_Family_Info.txt")
U, M, ST, y, genes, meta, F = build(D / "human_sites.tsv", utrs, all_mirs, max_rows=80000, seed=0)
tr, _ = gene_split(genes)
mu, sd = F[tr].mean(0), F[tr].std(0) + 1e-6
del U, M, ST, y, meta, F; gc.collect()
net = DuplexCNN(n_feat=4)
net.load_state_dict(torch.load(ROOT / "results" / "duplex_cnn.pt"))
net.eval()

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

rng = np.random.default_rng(0)
all_genes = sorted(gene2tids)
rows, per_mir, done = [], {}, set()
if os.path.exists(PER_JSONL):
    for line in open(PER_JSONL):
        per_mir.update(json.loads(line))
if os.path.exists(ROWS_CSV):
    with open(ROWS_CSV) as fh:
        for line in fh:
            p_ = line.rstrip("\n").split(",")
            if p_[0] == "mirna":
                continue
            rows.append((p_[0], p_[1], float(p_[2]), int(p_[3])))
            done.add(p_[0])
    print("resuming:", len(done), "done", flush=True)
rows_fh = open(ROWS_CSV, "a")
if os.path.getsize(ROWS_CSV) == 0:
    rows_fh.write("mirna,gene,min,label\n"); rows_fh.flush()

for mi, seq in sorted(mirs.items()):
    if mi in done:
        continue
    vp = {g for g in pos.get(mi, set()) if g in gene2tids}
    if len(vp) < MINPOS:
        per_mir[mi] = {"skipped": True, "validated": len(pos.get(mi, set())), "with_utr": len(vp)}
        continue
    banned = any_ev.get(mi, set())
    neg_pool = [g for g in all_genes if g not in banned]
    n_neg = min(len(vp) * NEG_MULT, MAXNEG, len(neg_pool))
    negs = set(rng.choice(neg_pool, size=n_neg, replace=False).tolist())
    eval_genes = vp | negs
    seed8, seed7 = rc(seq[1:8]), rc(seq[1:7])
    cands = []
    for g in eval_genes:
        for t in gene2tids[g]:
            u = utrs[t]
            i = u.find(seed8)
            while i >= 0:
                cands.append((t, i, 3 if (i + 7 < len(u) and u[i + 7] == "A") else 2))
                i = u.find(seed8, i + 1)
            i = u.find(seed7)
            while i >= 0:
                if i + 6 < len(u) and u[i + 6] == "A" and u[i:i + 7] != seed8:
                    cands.append((t, i, 2))
                i = u.find(seed7, i + 1)
    matched_genes = {tid2gene[t] for t, _, _ in cands}
    if cands:
        sc = score(net, utrs, seq, cands, mu, sd)
        agg = {}
        for (t, p_, s), v in zip(cands, sc):
            g = tid2gene[t]
            agg[g] = min(agg.get(g, 0.0), float(v))
        for g, mn in agg.items():
            lab = int(g in vp)
            rows.append((mi, g, mn, lab))
            rows_fh.write(f"{mi},{g},{mn},{lab}\n")
        rows_fh.flush()
        yv = np.array([int(g in vp) for g in agg])
        mins = np.array([agg[g] for g in agg])
        auroc = round(float(roc_auc_score(yv, -mins)), 4) if 0 < yv.sum() < len(yv) else None
    else:
        yv = np.array([]); auroc = None
    per_mir[mi] = {"validated": len(pos.get(mi, set())), "with_utr": len(vp),
                   "neg": len(negs), "seed_match_recall": round(len(matched_genes & vp) / len(vp), 4),
                   "scored_pos": int(sum(1 for g in (matched_genes & vp))),
                   "scored_genes": len(matched_genes & eval_genes), "auroc_cnn_min": auroc}
    with open(PER_JSONL, "a") as pf:
        pf.write(json.dumps({mi: per_mir[mi]}) + "\n")
    del cands; gc.collect()
    print(mi, per_mir[mi], flush=True)

# ---- pooled metrics over scored rows ----
y = np.array([r[3] for r in rows])
mins = np.array([r[2] for r in rows])
mids = np.array([r[0] for r in rows])
keep = np.array([not per_mir.get(m, {}).get("skipped") for m in mids])
y, mins, mids = y[keep], mins[keep], mids[keep]
pooled_auroc = round(float(roc_auc_score(y, -mins)), 4)
# within-miRNA z-score then Mann-Whitney
z = np.zeros_like(mins)
for m in set(mids):
    ix = mids == m
    z[ix] = (mins[ix] - mins[ix].mean()) / (mins[ix].std() + 1e-9)
mw_p = float(mannwhitneyu(-z[y == 1], -z[y == 0], alternative="greater").pvalue)
rng2 = np.random.default_rng(1)
cnt = 0
for _ in range(1000):
    yp = y.copy()
    for m in set(mids):
        ix = np.where(mids == m)[0]
        yp[ix] = rng2.permutation(yp[ix])
    if roc_auc_score(yp, -mins) >= pooled_auroc:
        cnt += 1
perm_p = (cnt + 1) / 1001
out = {"source": "miRTarBase 10.0 strong-evidence (reporter assay), hsa rows; negatives = no human MTI of any evidence in miRTarBase 10.0",
       "n_mirnas_evaluated": int(sum(1 for v in per_mir.values() if not v.get("skipped"))),
       "n_mirnas_skipped_minpos": int(sum(1 for v in per_mir.values() if v.get("skipped"))),
       "n_rows": int(len(y)), "n_pos": int(y.sum()),
       "pooled_auroc_cnn_min": pooled_auroc,
       "within_mirna_z_mannwhitney_p": mw_p, "permutation_p_within_mirna": perm_p,
       "per_mirna": per_mir}
json.dump(out, open(OUT, "w"), indent=1)
print(json.dumps({k: v for k, v in out.items() if k != "per_mirna"}, indent=1))
