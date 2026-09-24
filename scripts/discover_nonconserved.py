"""Discovery scan: strong NON-CONSERVED miRNA sites in human 3'UTRs.

For 221 broadly-conserved human miRNAs (TargetScan Family Conservation == 2),
scan all 25k human 3'UTRs for 8mer / 7mer seed matches that do NOT overlap any
TargetScan conserved site, then score every candidate with the trained
DuplexCNN+context model. Output: ranked named miRNA->transcript candidates
with predicted context++ scores; falsifiable against public CLIP / knockdown
data (e.g. ENCORI/starBase, miRTarBase) - checked in a follow-up step.
"""
import csv, json, sys
from pathlib import Path
import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from mirna.data import one_hot, load_utrs, site_features, build, gene_split, WIN, MIR_LEN
from mirna.model import DuplexCNN

COMP = str.maketrans("ACGU", "UGCA")
def rc(s): return s.translate(COMP)[::-1]

D = ROOT / "data"
# 1. broadly conserved human miRNAs
mirs = {}
with open(D / "miR_Family_Info.txt") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        if row["Species ID"] == "9606" and row["Family Conservation?"] == "2":
            seq = row["Mature sequence"].upper().replace("T", "U")
            if len(seq) >= 9:
                mirs[row["MiRBase ID"]] = seq
print("broadly conserved miRNAs:", len(mirs), flush=True)

# 2. conserved-site intervals per transcript (to exclude)
cons = {}
with open(D / "human_sites.tsv") as fh:
    r = csv.DictReader(fh, delimiter="\t")
    for row in r:
        cons.setdefault(row["transcript"], []).append((int(row["start"]), int(row["end"])))
print("transcripts with conserved sites:", len(cons), flush=True)

utrs = load_utrs(D / "human_utrs.tsv")
tid2gene = {}
with open(D / "human_utrs.tsv") as fh:
    next(fh)
    for line in fh:
        t, g, _ = line.rstrip("\n").split("\t", 2)
        tid2gene[t.split(".")[0]] = g

def overlaps(t, lo, hi):  # 0-based candidate span vs 1-based conserved intervals
    return any(not (hi < a - 1 or lo > b - 1) for a, b in cons.get(t, []))

# 3. scan
cands = []  # (mir, tid, pos, st)
for mi, seq in mirs.items():
    seed8 = rc(seq[1:8])      # nt 2-8 (7nt)
    seed7 = rc(seq[1:7])      # nt 2-7 (6nt)
    for t, u in utrs.items():
        i = u.find(seed8)
        while i >= 0:
            a1 = u[i + 7] == "A" if i + 7 < len(u) else False
            if not overlaps(t, i, i + 6):
                cands.append((mi, t, i, 3 if a1 else 2))
            i = u.find(seed8, i + 1)
        i = u.find(seed7)
        while i >= 0:
            a1 = u[i + 6] == "A" if i + 6 < len(u) else False
            covered = u[i:i + 7] == seed8
            if a1 and not covered and not overlaps(t, i, i + 5):
                cands.append((mi, t, i, 2))  # 7mer-A1
            i = u.find(seed7, i + 1)
print("novel candidates:", len(cands), flush=True)

# 4. rebuild training arrays (same seed -> identical subsample) for F normalisation
print("rebuilding training normalisation...", flush=True)
from mirna.data import load_mirnas
all_mirs = load_mirnas(D / "miR_Family_Info.txt")
U, M, ST, y, genes, meta, F = build(D / "human_sites.tsv", utrs, all_mirs, max_rows=80000, seed=0)
tr, te = gene_split(genes)
mu, sd = F[tr].mean(0), F[tr].std(0) + 1e-6
net = DuplexCNN(n_feat=4)
net.load_state_dict(torch.load(ROOT / "results" / "duplex_cnn.pt"))
net.eval()

# 5. score candidates (subsample if huge)
rng = np.random.default_rng(0)
if len(cands) > 400_000:
    cands = [cands[i] for i in rng.choice(len(cands), 400_000, replace=False)]
    print("subsampled to 400k", flush=True)
Bs = 2048
scores = np.zeros(len(cands), dtype=np.float32)
with torch.no_grad():
    for lo in range(0, len(cands), Bs):
        chunk = cands[lo:lo + Bs]
        Xu = torch.from_numpy(np.stack([
            one_hot("".join(utrs[t][k] if 0 <= k < len(utrs[t]) else "N"
                            for k in range(p + 3 - WIN // 2, p + 3 - WIN // 2 + WIN)), WIN)
            for _, t, p, _ in chunk]))
        Xm = torch.from_numpy(np.stack([one_hot(mirs[m], MIR_LEN) for m, _, _, _ in chunk]))
        Xs = torch.from_numpy(np.array([s for *_, s in chunk], dtype=np.int64))
        Xf = torch.from_numpy(((np.stack([
            site_features(utrs[t], p + 1, p + 7) for _, t, p, _ in chunk
        ]) - mu) / sd).astype(np.float32))
        scores[lo:lo + len(chunk)] = net(Xu, Xm, Xs, Xf).numpy()
        if lo % (Bs * 50) == 0:
            print("scored", lo, flush=True)

order = np.argsort(scores)
def entry(i):
    return {"mirna": cands[i][0], "transcript": cands[i][1],
            "gene": tid2gene.get(cands[i][1]),
            "utr_pos": cands[i][2], "site_class": {2: "7mer", 3: "8mer"}[cands[i][3]],
            "pred_context_pp": round(float(scores[i]), 4)}
out = [entry(i) for i in order[:300]]
seen, diverse = set(), []
for i in order:
    k = (cands[i][0], tid2gene.get(cands[i][1]))
    if k not in seen:
        seen.add(k); diverse.append(entry(i))
    if len(diverse) >= 300:
        break
json.dump({"model": "DuplexCNN+context v2 (held-out-gene pearson 0.811)",
           "n_candidates_scored": len(cands),
           "definition": "8mer/7mer seed match, no overlap with any TargetScan conserved site",
           "top_overall": out, "top_diverse_per_pair": diverse}, open(ROOT / "results" / "nonconserved_discovery.json", "w"), indent=1)
for e in out[:15]:
    print(e, flush=True)
print("saved results/nonconserved_discovery.json")
