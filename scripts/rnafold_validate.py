"""ViennaRNA duplex-MFE validation of top non-conserved candidates.

Orthogonal biophysical check: do the predicted top sites form stable
miRNA:target duplexes? RNAcofold MFE of (mature miRNA)&(site window +/-25nt)
vs mononucleotide-shuffled miRNA controls (10 shuffles/pair, seeded).
Real sites should sit below the shuffled distribution.
Output: results/rnafold_validation.json + paper/figs/fig_rnafold.pdf
"""
import csv, json, random
import numpy as np
import RNA
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# mature miRNA sequences (human, species 9606) keyed by MiRBase ID
mirseq = {}
with open("data/miR_Family_Info.txt") as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 7 and p[2] == "9606":
            mirseq[p[3]] = p[4].replace("U", "T")

# UTR sequences
utrs = {}
with open("data/human_utrs.tsv") as fh:
    r = csv.reader(fh, delimiter="\t")
    header = next(r)
    for row in r:
        if len(row) >= 3:
            utrs[row[0].split(".")[0]] = row[2].upper().replace("U", "T")
            utrs[row[1]] = row[2].upper().replace("U", "T")  # gene-name alias

disc = json.load(open("results/nonconserved_discovery.json"))
cands = disc["top_diverse_per_pair"][:20]
print(f"candidates: {len(cands)}", flush=True)

def cofold_mfe(a, b):
    fc = RNA.fold_compound(a.replace("T", "U") + "&" + b.replace("T", "U"))
    m = fc.mfe()
    return float(m[1])

rng = random.Random(0)
rows = []
for c in cands:
    mir, tx, pos = c["mirna"], c["transcript"], int(c["utr_pos"])
    ms, utr = mirseq.get(mir), utrs.get(tx)
    if not ms or not utr:
        print("skip", mir, tx, "missing seq", flush=True)
        continue
    win = utr[max(0, pos - 25): pos + 33]
    real = cofold_mfe(ms, win)
    shufs = []
    for _ in range(10):
        s = list(ms); rng.shuffle(s)
        shufs.append(cofold_mfe("".join(s), win))
    rows.append({"mirna": mir, "gene": c["gene"], "site_class": c["site_class"],
                 "pred_context_pp": c["pred_context_pp"],
                 "mfe_real": round(real, 2),
                 "mfe_shuf_mean": round(float(np.mean(shufs)), 2),
                 "mfe_shuf_sd": round(float(np.std(shufs)), 2),
                 "z": round((float(np.mean(shufs)) - real) / (float(np.std(shufs)) + 1e-9), 2)})
    print(f"{mir}->{c['gene']}: real {real:.1f} vs shuf {np.mean(shufs):.1f}+-{np.std(shufs):.1f} (z {rows[-1]['z']})", flush=True)

zs = np.array([r["z"] for r in rows])
from scipy.stats import wilcoxon
w, p = wilcoxon(zs) if len(zs) > 0 and np.any(zs != 0) else (None, None)
res = {"n": len(rows), "z_median": round(float(np.median(zs)), 2),
       "frac_z_gt1": round(float((zs > 1).mean()), 2),
       "wilcoxon_p_z_gt0": float(p) if p is not None else None,
       "note": "z = (shuffled-mean - real)/shuffled-sd; positive = real duplex more stable than shuffled-miRNA controls",
       "pairs": rows}
json.dump(res, open("results/rnafold_validation.json", "w"), indent=1)

fig, ax = plt.subplots(figsize=(4.4, 3.2))
real = [r["mfe_real"] for r in rows]; shuf = [r["mfe_shuf_mean"] for r in rows]
ax.boxplot([real, shuf], labels=["real miRNA", "shuffled x10 mean"])
ax.set_ylabel("RNAcofold MFE (kcal/mol)")
fig.tight_layout(); fig.savefig("paper/figs/fig_rnafold.pdf")
print("saved results/rnafold_validation.json + fig_rnafold.pdf; median z:", res["z_median"], "frac z>1:", res["frac_z_gt1"])
