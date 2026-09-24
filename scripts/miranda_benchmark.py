"""miRanda arm on the committed miRTarBase benchmark rows (fourth head-to-head).

Question: does the classic alignment-based predictor miRanda (3.3a, aug2010
release, canonical source via Wayback capture of cbio.mskcc.org, built with
gcc -fcommon; default cutoffs, -quiet) rank literature-validated targets
better than the DuplexCNN on the same committed rows
(results/mirtarbase_rows.csv)? Pure alignment+thermodynamics, no training
data, no conservation - the same clean biophysical arm class as the RNAhybrid
head-to-head, scored by a different algorithm (Smith-Waterman + Vienna energy).

Per miRNA: one miRanda call over all UTR transcripts of that miRNA's
candidate genes; per-gene score = max alignment score over transcripts;
genes with no hit at default cutoffs get a floor one step below the weakest
hit (coverage reported). Binary path via MIRANDA_BIN (default
/tmp/miranda/bin/miranda; refetch+rebuild per docs/TOOLS.md if /tmp wiped).
Checkpoints per miRNA to results/miranda_progress.jsonl; rerun to resume.
MIRANDA_BUDGET_S caps new work per invocation (default 240s)."""
import csv, json, os, subprocess, sys, tempfile, time
from pathlib import Path
import numpy as np
from sklearn.metrics import roc_auc_score
from scipy.stats import wilcoxon

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from mirna.data import load_mirnas

BIN = os.environ.get("MIRANDA_BIN", "/tmp/miranda/bin/miranda")
ROWS = ROOT / "results" / "mirtarbase_rows.csv"
OUT = ROOT / "results" / "miranda_benchmark.json"
OUTROWS = ROOT / "results" / "miranda_rows.csv"
PROG = ROOT / "results" / "miranda_progress.jsonl"
MIN_PER_MIRNA = 30
BUDGET = float(os.environ.get("MIRANDA_BUDGET_S", "240"))


def run_mirna(mir_seq, gene_utrs):
    """{gene: best miRanda score} for one miRNA over {gene: [(tx, utr)]}."""
    with tempfile.TemporaryDirectory() as td:
        t, q = Path(td) / "t.fa", Path(td) / "q.fa"
        t2g = {}
        with open(t, "w") as fh:
            for g, tl in gene_utrs.items():
                for tx, u in tl:
                    t2g[tx] = g
                    fh.write(f">{tx}\n{u}\n")
        q.write_text(f">q\n{mir_seq}\n")
        r = subprocess.run([BIN, str(q), str(t), "-quiet"],
                           capture_output=True, text=True, timeout=1800)
    best = {}
    for line in r.stdout.splitlines():
        if not line.startswith(">") or line.startswith(">>"):
            continue
        f = line[1:].split("\t")
        if len(f) < 3:
            continue
        try:
            tx, sc = f[1], float(f[2])
        except ValueError:
            continue
        g = t2g.get(tx)
        if g is not None and (g not in best or sc > best[g]):
            best[g] = sc
    return best


def auroc(labels, scores):
    return float(roc_auc_score(labels, scores))


def main():
    t0 = time.time()
    rows = list(csv.DictReader(open(ROWS)))
    gene2t = {}
    with open(ROOT / "data" / "human_utrs.tsv") as fh:
        next(fh)
        for line in fh:
            t, g, u = line.rstrip("\n").split("\t", 2)
            gene2t.setdefault(g, []).append((t, u))
    mirs = load_mirnas(str(ROOT / "data" / "miR_Family_Info.txt"))
    by_mir = {}
    for r in rows:
        if r["gene"] in gene2t:
            by_mir.setdefault(r["mirna"], []).append(r)
    done = {}
    if PROG.exists():
        for line in open(PROG):
            d = json.loads(line)
            done[d["mirna"]] = d
    n_new = 0
    with open(PROG, "a") as pf:
        for mirna, rs in sorted(by_mir.items()):
            if mirna in done or mirna not in mirs or len(rs) < MIN_PER_MIRNA:
                continue
            y = np.array([int(r["label"]) for r in rs])
            if y.sum() == 0 or y.sum() == len(y):
                continue
            if n_new and time.time() - t0 > BUDGET:
                break
            targs = {r["gene"]: gene2t[r["gene"]] for r in rs}
            best = run_mirna(mirs[mirna], targs)
            m_raw = np.array([best[r["gene"]] if r["gene"] in best else np.nan for r in rs])
            floor = np.nanmin(m_raw) - 1.0 if (~np.isnan(m_raw)).any() else 0.0
            m = np.where(np.isnan(m_raw), floor, m_raw)
            cnn = np.array([-float(r["min"]) for r in rs])
            rec = {"mirna": mirna, "n": len(rs), "n_pos": int(y.sum()),
                   "auroc_miranda": auroc(y, m), "auroc_cnn": auroc(y, cnn),
                   "coverage": float((~np.isnan(m_raw)).mean()),
                   "rows": [{"mirna": mirna, "gene": r["gene"], "label": r["label"],
                             "miranda_score": ("" if np.isnan(raw) else round(float(sc), 3)),
                             "cnn_score": round(float(c), 6)}
                            for r, raw, sc, c in zip(rs, m_raw, m, cnn)]}
            pf.write(json.dumps(rec) + "\n"); pf.flush()
            done[mirna] = rec
            n_new += 1
            print(f"{mirna}: auroc {rec['auroc_miranda']:.3f} vs cnn {rec['auroc_cnn']:.3f} "
                  f"cov {rec['coverage']:.2f} ({time.time()-t0:.0f}s)", flush=True)
    per = [v for v in done.values() if "auroc_miranda" in v]
    out_rows = [row for v in per for row in v["rows"]]
    dm = np.array([p["auroc_miranda"] for p in per])
    dc = np.array([p["auroc_cnn"] for p in per])
    y_all = np.array([int(r["label"]) for r in out_rows])
    res = {
        "design": __doc__.strip().splitlines()[0],
        "tool": "miRanda 3.3a (aug2010 release, canonical cbio.mskcc.org source via Wayback, gcc -fcommon build), default cutoffs, -quiet",
        "n_mirnas": len(per), "n_rows": len(out_rows),
        "pooled_auroc_miranda": auroc(y_all, [float(r["miranda_score"]) if r["miranda_score"] != "" else -1e9 for r in out_rows]),
        "pooled_auroc_cnn": auroc(y_all, [r["cnn_score"] for r in out_rows]),
        "median_per_mirna_miranda": float(np.median(dm)),
        "median_per_mirna_cnn": float(np.median(dc)),
        "miranda_wins": int((dm > dc).sum()),
        "wilcoxon_p_miranda_minus_cnn": float(wilcoxon(dm - dc).pvalue),
        "median_coverage": float(np.median([p["coverage"] for p in per])),
        "per_mirna": [{k: v for k, v in p.items() if k != "rows"} for p in per],
    }
    # three-way identical-row join with the RNAhybrid arm
    rh = {(r["mirna"], r["gene"]): r for r in csv.DictReader(open(ROOT / "results" / "rnahybrid_rows.csv"))}
    common = [r for r in out_rows if (r["mirna"], r["gene"]) in rh]
    if common:
        y3 = np.array([int(r["label"]) for r in common])
        m3 = np.array([float(r["miranda_score"]) if r["miranda_score"] != "" else -1e9 for r in common])
        c3 = np.array([float(r["cnn_score"]) for r in common])
        r3 = np.array([float(rh[(r["mirna"], r["gene"])]["rnahybrid_score"]) if rh[(r["mirna"], r["gene"])]["rnahybrid_score"] != "" else -1e9 for r in common])
        res["identical_row_three_way"] = {
            "n_rows": len(common),
            "pooled_auroc_miranda": auroc(y3, m3),
            "pooled_auroc_cnn": auroc(y3, c3),
            "pooled_auroc_rnahybrid": auroc(y3, r3),
        }
    with open(OUT, "w") as f:
        json.dump(res, f, indent=1)
    with open(OUTROWS, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["mirna", "gene", "label", "miranda_score", "cnn_score"])
        w.writeheader(); w.writerows(out_rows)
    print(json.dumps({k: v for k, v in res.items() if k != "per_mirna"}, indent=1))


if __name__ == "__main__":
    main()
