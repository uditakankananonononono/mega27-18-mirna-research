"""PRE-REGISTERED falsifier (judge round 1, docs/JUDGE_ROUNDS.md): is the
dosage-sensitivity depletion of CNN top-ranked non-conserved candidates a
SITE-ENUMERATION ARTEFACT?

Locked design (committed BEFORE running, 2026-09-26):
- Candidate enumeration identical to scripts/discover_nonconserved.py
  (8mer/7mer-A1 seed matches in human 3'UTRs not overlapping TargetScan
  conserved sites), subsampled to 200,000 with seed 0 (declared; original
  used 400k - reduced for CPU budget, declared here).
- Universe + statistic identical to scripts/explof_matched.py: CLIP universe
  (results/clip_rows_136.csv), panel = 13 discovery + 30 held-out miRNAs,
  K=200, exp_lof-matched comparators (3 draws, decile matcher), one-sided
  sign tests over the panel: L1 (LOEUF depletion), L2 (s_het depletion).
- NULL A (random weights): 100 DuplexCNN instances, torch.manual_seed(seed),
  seed=0..99, NEVER trained, same normalisation (mu/sd from training rebuild
  seed 0). Statistic per seed: number of panel miRNAs with CNN median LOEUF
  above matched mean (L1 count) and below for s_het (L2 count).
- NULL B (score permutation): within each miRNA, permute the trained model's
  candidate scores across genes, 100 permutations, same statistics.
- VERDICT RULE (locked, alpha 0.01): the depletion survives only if the
  trained model's L1 count exceeds the 99th percentile of BOTH null
  distributions in the depletion direction (and L2 likewise). If either
  null reproduces the counts, the corresponding claim is a site-enumeration
  artefact and is withdrawn (pivot ladder in docs/PREREGISTER.md).
- Compute honesty note: full label-permutation RETRAINING (100x) is
  infeasible on this CPU budget; the permutation null here permutes scores,
  not labels. Up to 5 label-permutation retrains may follow as a spot check
  if Null A/B are beaten; declared now, before outcomes.
"""
import json, random, sys
from pathlib import Path
import numpy as np, pandas as pd, torch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
from mirna.data import one_hot, load_utrs, load_mirnas, site_features, build, gene_split, WIN, MIR_LEN
from mirna.model import DuplexCNN
from gprofiler_enrichment import PANEL, K
from string_exprmatched import matched
from genebayes_shet import load_hi as load_shet

D = ROOT / "data"
COMP = str.maketrans("ACGU", "UGCA")
def rc(s): return s.translate(COMP)[::-1]

def enumerate_candidates():
    mirs = {}
    import csv
    with open(D / "miR_Family_Info.txt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["Species ID"] == "9606" and row["Family Conservation?"] == "2":
                seq = row["Mature sequence"].upper().replace("T", "U")
                if len(seq) >= 9:
                    mirs[row["MiRBase ID"]] = seq
    cons = {}
    with open(D / "human_sites.tsv") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            cons.setdefault(row["transcript"], []).append((int(row["start"]), int(row["end"])))
    utrs = load_utrs(D / "human_utrs.tsv")
    tid2gene = {}
    with open(D / "human_utrs.tsv") as fh:
        next(fh)
        for line in fh:
            t, g, _ = line.rstrip("\n").split("\t", 2)
            tid2gene[t.split(".")[0]] = g
    def overlaps(t, lo, hi):
        return any(lo < e and hi > s - 1 for s, e in cons.get(t, []))
    cands = []
    for mi, seq in mirs.items():
        seed8 = rc(seq[1:8]); seed7 = rc(seq[1:7])
        for t, u in utrs.items():
            i = u.find(seed8)
            while i >= 0:
                if not overlaps(t, i, i + 6):
                    cands.append((mi, t, i, 3))
                i = u.find(seed8, i + 1)
            i = u.find(seed7)
            while i >= 0:
                a1 = i + 6 < len(u) and u[i + 6] == "A"
                if a1 and u[i:i + 7] != seed8 and not overlaps(t, i, i + 5):
                    cands.append((mi, t, i, 2))
                i = u.find(seed7, i + 1)
    return cands, utrs, mirs, tid2gene

def features(cands, utrs, mirs, mu, sd):
    rng = np.random.default_rng(0)
    if len(cands) > 200_000:
        cands = [cands[i] for i in rng.choice(len(cands), 200_000, replace=False)]
    Xu = np.stack([one_hot("".join(utrs[t][k] if 0 <= k < len(utrs[t]) else "N"
                            for k in range(p + 3 - WIN // 2, p + 3 - WIN // 2 + WIN)), WIN)
                   for _, t, p, _ in cands])
    Xm = np.stack([one_hot(mirs[m], MIR_LEN) for m, _, _, _ in cands])
    Xs = np.array([s for *_, s in cands], dtype=np.int64)
    Xf = ((np.stack([site_features(utrs[t], p + 1, p + 7) for _, t, p, _ in cands]) - mu) / sd).astype(np.float32)
    return cands, Xu, Xm, Xs, Xf

def score(net, Xu, Xm, Xs, Xf, bs=512):
    # Pre-outcome memory safety: smaller inference batches on a 1.9GB/no-swap
    # box after three OOM kills. This changes neither candidates nor statistics.
    net.eval(); out = np.zeros(len(Xs), dtype=np.float32)
    with torch.no_grad():
        for lo in range(0, len(Xs), bs):
            out[lo:lo + bs] = net(torch.from_numpy(Xu[lo:lo + bs]), torch.from_numpy(Xm[lo:lo + bs]),
                                  torch.from_numpy(Xs[lo:lo + bs]), torch.from_numpy(Xf[lo:lo + bs])).numpy()
    return out

def depletion_counts(scores, cands, tid2gene, universe, g, sh, tag):
    df = pd.DataFrame({"mirna": [c[0] for c in cands], "gene": [tid2gene.get(c[1]) for c in cands],
                       "score": scores}).dropna()
    df = df.groupby(["mirna", "gene"], as_index=False)["score"].min()
    df = df.merge(universe, on=["mirna", "gene"])
    l1 = l2 = 0
    for mi, s in df.groupby("mirna"):
        if mi not in PANEL_SET:
            continue
        s = s.reset_index(drop=True); s = s.assign(x=s.gene.map(g.exp_lof))
        cnn = s.nsmallest(K, "score").gene.tolist()
        rng = random.Random(f"{mi}-{tag}")
        M = [matched(s.assign(min=s.score), cnn, rng) for _ in range(3)]
        mlo = lambda gs: float(np.median(g.oe_lof_upper.reindex(gs).values))
        msh = lambda gs: float(np.median(sh.reindex(gs).values))
        if mlo(cnn) > float(np.mean([mlo(x) for x in M])): l1 += 1
        if msh(cnn) < float(np.mean([msh(x) for x in M])): l2 += 1
    return l1, l2

def main():
    g = pd.read_csv(D / "gnomad/lof_v211.txt.bgz", sep="\t", compression="gzip",
                    usecols=["gene", "oe_lof_upper", "exp_lof"]).dropna()
    g = g.sort_values("exp_lof", ascending=False).drop_duplicates("gene").set_index("gene")
    sh = load_shet()
    universe = pd.read_csv(ROOT / "results/clip_rows_136.csv", usecols=["mirna", "gene"]).drop_duplicates()
    universe = universe[universe.gene.isin(set(g.index) & set(sh.index))]
    print("rebuilding training normalisation (seed 0)...", flush=True)
    all_mirs = load_mirnas(D / "miR_Family_Info.txt")
    utrs0 = load_utrs(D / "human_utrs.tsv")
    U0, M0, ST0, y0, genes0, meta0, F0 = build(D / "human_sites.tsv", utrs0, all_mirs, max_rows=80000, seed=0)
    tr0, _ = gene_split(genes0)
    mu, sd = F0[tr0].mean(0), F0[tr0].std(0) + 1e-6
    del U0, M0, ST0, y0, F0
    print("enumerating candidates...", flush=True)
    cands, utrs, mirs, tid2gene = enumerate_candidates()
    print("candidates:", len(cands), flush=True)
    cands, Xu, Xm, Xs, Xf = features(cands, utrs, mirs, mu, sd)
    # memory-only fix (no statistic changes): free the full UTR table after
    # feature extraction; the locked design, seeds, subsample and statistics
    # are unchanged. Two OOM kills (pids 4087, 4247) forced this.
    del utrs, mirs
    import gc; gc.collect()
    print("scored feature arrays ready", flush=True)
    net = DuplexCNN(n_feat=4)
    net.load_state_dict(torch.load(ROOT / "results/duplex_cnn.pt", weights_only=True))
    trained_scores = score(net, Xu, Xm, Xs, Xf)
    t_l1, t_l2 = depletion_counts(trained_scores, cands, tid2gene, universe, g, sh, "trained")
    print("trained counts L1 L2:", t_l1, t_l2, flush=True)
    null_a = []
    for seed in range(100):
        torch.manual_seed(seed)
        rn = DuplexCNN(n_feat=4)
        s = score(rn, Xu, Xm, Xs, Xf)
        null_a.append(depletion_counts(s, cands, tid2gene, universe, g, sh, f"rw-{seed}"))
        if seed % 10 == 9: print("null A seed", seed + 1, flush=True)
    null_b = []
    rng = np.random.default_rng(0)
    mi_idx = np.array([c[0] for c in cands])
    for p in range(100):
        sp = trained_scores.copy()
        for mi in np.unique(mi_idx):
            m = mi_idx == mi
            sp[m] = rng.permutation(sp[m])
        null_b.append(depletion_counts(sp, cands, tid2gene, universe, g, sh, f"perm-{p}"))
        if p % 10 == 9: print("null B perm", p + 1, flush=True)
    a1 = [x[0] for x in null_a]; a2 = [x[1] for x in null_a]
    b1 = [x[0] for x in null_b]; b2 = [x[1] for x in null_b]
    J = {"trained": {"L1": t_l1, "L2": t_l2},
         "null_A_random_weight": {"L1": a1, "L2": a2},
         "null_B_permutation": {"L1": b1, "L2": b2},
         "verdict": {
           "L1_survives": bool(t_l1 > np.percentile(a1, 99) and t_l1 > np.percentile(b1, 99)),
           "L2_survives": bool(t_l2 > np.percentile(a2, 99) and t_l2 > np.percentile(b2, 99))},
         "alpha": 0.01, "K": K, "n_candidates": len(cands)}
    json.dump(J, open(ROOT / "results/random_weight_falsify.json", "w"), indent=1)
    print(json.dumps(J["verdict"]), flush=True)

PANEL_SET = set(PANEL + json.load(open(ROOT / "results/gnomad_replication.json"))["panel"])
if __name__ == "__main__":
    main()
