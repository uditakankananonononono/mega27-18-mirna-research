"""PRE-REGISTERED miRNA-family-held-out transfer test (judge round 2 foldback;
docs/JUDGE_ROUNDS.md). Committed BEFORE running, 2026-09-26, and before the
random-weight falsifier verdict exists.

Question (judge round 2): does DuplexCNN learn transferable pairing rules, or
within-family sequence patterns? Locked design:
- Arrays: identical build to training (human_sites.tsv, vert_80 UTRs,
  miR_Family_Info, max_rows=80000, seed 0).
- Families: TargetScan "miR family" column; the 5 families with the most
  rows in the 80k build are the held-out cases (declared; no outcome peeked).
- Per family f: train DuplexCNN+context (same architecture, epochs=6,
  lr=2e-3, bs=256, seed 0) on all rows NOT in f; ridge+context (alpha=10)
  on the same train rows; one random-weight DuplexCNN (torch seed 0,
  untrained) as architecture-level null.
- Metric: pearson r on held-out family rows for each method.
- VERDICT RULE (locked): "transferable rules" supported only if the CNN
  beats ridge+context on mean pearson across the 5 held-out families AND
  CNN mean pearson exceeds the random-weight null by more than 0.05.
  Otherwise report the boundary honestly (interpolation, not transfer).
- Compute note: 5 retrains x ~7 min CPU each; declared pre-outcome.
"""
import csv, json, sys
from pathlib import Path
import numpy as np, torch
from scipy.stats import pearsonr
from sklearn.linear_model import Ridge

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from mirna.data import load_mirnas, load_utrs, build
from mirna.model import DuplexCNN
from mirna.train import fit_cnn, predict

D = ROOT / "data"

def main():
    fam_of = {}
    with open(D / "miR_Family_Info.txt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["Species ID"] == "9606":
                fam_of[row["MiRBase ID"]] = row["miR family"]
    mirs = load_mirnas(D / "miR_Family_Info.txt")
    utrs = load_utrs(D / "human_utrs.tsv")
    U, M, ST, y, genes, meta, F = build(D / "human_sites.tsv", utrs, mirs, max_rows=80000, seed=0)
    del utrs
    ST = ST.clip(0, 3)
    fams = np.array([fam_of.get(m[2], m[2]) for m in meta])
    top = [f for f, _ in sorted(((f, int((fams == f).sum())) for f in np.unique(fams)),
                                key=lambda kv: -kv[1])[:5]]
    print("held-out families:", top, flush=True)
    out = {"held_out_families": top, "families": {}}
    for f in top:
        te = np.where(fams == f)[0]; tr = np.where(fams != f)[0]
        # Pre-outcome leakage fix: fit normalization on TRAIN families only.
        # Do not allow held-out family context distributions into preprocessing.
        mu, sd = F[tr].mean(0), F[tr].std(0) + 1e-6
        Fn = ((F - mu) / sd).astype(np.float32)
        Xr = np.concatenate([U.reshape(len(y), -1), M.reshape(len(y), -1), np.eye(4)[ST], Fn], 1)
        rg = Ridge(alpha=10.0).fit(Xr[tr], y[tr])
        r_ridge = float(pearsonr(y[te], rg.predict(Xr[te]))[0])
        net = fit_cnn(U, M, ST, y, tr, epochs=6, log=lambda s: None, F=Fn)
        r_cnn = float(pearsonr(y[te], predict(net, U, M, ST, te, F=Fn))[0])
        torch.manual_seed(0)
        rn = DuplexCNN(n_feat=Fn.shape[1])
        r_null = float(pearsonr(y[te], predict(rn, U, M, ST, te, F=Fn))[0])
        out["families"][f] = {"n_test": int(len(te)), "cnn": r_cnn,
                              "ridge_context": r_ridge, "random_weight": r_null}
        print(f, "cnn", round(r_cnn, 4), "ridge", round(r_ridge, 4), "null", round(r_null, 4), flush=True)
        # crash-safe checkpoint after each family (does not alter computation)
        json.dump({"held_out_families": top, "families": out["families"], "partial": True},
                  open(ROOT / "results/family_transfer_checkpoint.json", "w"), indent=1)
    mc = float(np.mean([v["cnn"] for v in out["families"].values()]))
    mr = float(np.mean([v["ridge_context"] for v in out["families"].values()]))
    mn = float(np.mean([v["random_weight"] for v in out["families"].values()]))
    out["mean"] = {"cnn": mc, "ridge_context": mr, "random_weight": mn}
    out["verdict_transfer_supported"] = bool(mc > mr and mc - mn > 0.05)
    json.dump(out, open(ROOT / "results/family_transfer.json", "w"), indent=1)
    print(json.dumps({"mean": out["mean"], "verdict": out["verdict_transfer_supported"]}), flush=True)

if __name__ == "__main__":
    main()
