"""RNAhybrid arm on the committed miRTarBase benchmark rows (third head-to-head).

Question: does pure biophysical duplex prediction (RNAhybrid 2.1.2 minimum free
energy, -s 3utr_human) rank literature-validated targets better than the
DuplexCNN on the same committed rows (results/mirtarbase_rows.csv)? No training
data, no conservation, no ML - the cleanest possible biophysical arm.

Binary: apt-get download rnahybrid + libg20, dpkg-deb -x, run with
LD_LIBRARY_PATH=<extract>/usr/lib; path via RNAHYBRID_BIN (default /tmp/rh/...).
Per miRNA: one RNAhybrid call (-c -b 1) over the UTRs of that miRNA's genes.
Rows with no hit get score -inf (weakest possible rank). Coverage reported.
"""
import csv, json, os, subprocess, sys, tempfile
from pathlib import Path
import numpy as np
from sklearn.metrics import roc_auc_score
from scipy.stats import wilcoxon

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from mirna.data import load_mirnas

BIN = os.environ.get("RNAHYBRID_BIN", "/tmp/rh/usr/bin/RNAhybrid")
ROWS = ROOT / "results" / "mirtarbase_rows.csv"
OUT = ROOT / "results" / "rnahybrid_benchmark.json"
MIN_PER_MIRNA = 30


def parse_compact(line):
    """One -c output line -> (target, energy) or None."""
    f = line.rstrip("\n").split(":")
    if len(f) < 6:
        return None
    try:
        return f[0], float(f[4])
    except ValueError:
        return None


def run_mirna(mir_seq, gene_utrs, binpath=BIN, env=None):
    """Return {gene: best_energy} for one miRNA over {gene: utr}."""
    with tempfile.TemporaryDirectory() as td:
        t, q = Path(td) / "t.fa", Path(td) / "q.fa"
        with open(t, "w") as fh:
            for g, u in gene_utrs.items():
                fh.write(f">{g}\n{u}\n")
        q.write_text(f">q\n{mir_seq}\n")
        r = subprocess.run([binpath, "-s", "3utr_human", "-c", "-b", "1",
                            "-t", str(t), "-q", str(q)],
                           capture_output=True, text=True, env=env, timeout=1800)
    best = {}
    for line in r.stdout.splitlines():
        p = parse_compact(line)
        if p and (p[0] not in best or p[1] < best[p[0]]):
            best[p[0]] = p[1]
    return best


def auroc(labels, scores):
    return float(roc_auc_score(labels, scores))


def main():
    rows = list(csv.DictReader(open(ROWS)))
    gene2t = {}
    with open(ROOT / "data" / "human_utrs.tsv") as fh:
        next(fh)
        for line in fh:
            t, g, u = line.rstrip("\n").split("\t", 2)
            gene2t.setdefault(g, []).append((t, u))
    mirs = load_mirnas(str(ROOT / "data" / "miR_Family_Info.txt"))
    env = dict(os.environ, LD_LIBRARY_PATH="/tmp/rh/usr/lib")
    by_mir = {}
    for r in rows:
        if r["gene"] in gene2t:
            by_mir.setdefault(r["mirna"], []).append(r)
    out_rows, per = [], []
    for mirna, rs in sorted(by_mir.items()):
        if mirna not in mirs or len(rs) < MIN_PER_MIRNA:
            continue
        targs = {}
        for r in rs:
            for t, u in gene2t[r["gene"]]:
                targs[t] = u
        best_t = run_mirna(mirs[mirna], targs, env=env)
        t2g = {t: g for g, tl in ((r["gene"], gene2t[r["gene"]]) for r in rs) for t, _ in tl}
        best = {}
        for t, e in best_t.items():
            g = t2g[t]
            if g not in best or e < best[g]:
                best[g] = e
        y = np.array([int(r["label"]) for r in rs])
        if y.sum() == 0 or y.sum() == len(y):
            continue
        rh_raw = np.array([-best[r["gene"]] if r["gene"] in best else np.nan for r in rs])
        floor = np.nanmin(rh_raw) - 1.0  # no hit -> one step below the weakest hit
        rh = np.where(np.isnan(rh_raw), floor, rh_raw)
        cnn = np.array([-float(r["min"]) for r in rs])
        a_rh, a_cnn = auroc(y, rh), auroc(y, cnn)
        per.append({"mirna": mirna, "n": len(rs), "n_pos": int(y.sum()),
                    "auroc_rnahybrid": a_rh, "auroc_cnn": a_cnn,
                    "coverage": float((~np.isnan(rh_raw)).mean())})
        for r, e, hit in zip(rs, rh, ~np.isnan(rh_raw)):
            out_rows.append({"mirna": mirna, "gene": r["gene"], "label": r["label"],
                             "rnahybrid_score": e if hit else "",
                             "cnn_score": -float(r["min"])})
    dr = np.array([p["auroc_rnahybrid"] for p in per])
    dc = np.array([p["auroc_cnn"] for p in per])
    y_all = np.array([int(r["label"]) for r in out_rows])
    res = {
        "design": __doc__.strip().splitlines()[0],
        "tool": "RNAhybrid 2.1.2-6 (Ubuntu jammy .deb), -s 3utr_human -c -b 1",
        "n_mirnas": len(per), "n_rows": len(out_rows),
        "pooled_auroc_rnahybrid": auroc(y_all, [p["auroc_rnahybrid"] for p in per for _ in [0]]) if False else auroc(
            np.array([int(p_["label"]) for p_ in [r for r in out_rows]]),
            [float(r["rnahybrid_score"]) if r["rnahybrid_score"] != "" else -1e9 for r in out_rows]),
        "pooled_auroc_cnn": auroc(y_all, [r["cnn_score"] for r in out_rows]),
        "median_per_mirna_rnahybrid": float(np.median(dr)),
        "median_per_mirna_cnn": float(np.median(dc)),
        "rnahybrid_wins": int((dr > dc).sum()),
        "wilcoxon_p_rnahybrid_minus_cnn": float(wilcoxon(dr - dc).pvalue),
        "median_coverage": float(np.median([p["coverage"] for p in per])),
        "per_mirna": per,
    }
    with open(OUT, "w") as f:
        json.dump(res, f, indent=1)
    with open(ROOT / "results" / "rnahybrid_rows.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["mirna", "gene", "label", "rnahybrid_score", "cnn_score"])
        w.writeheader(); w.writerows(out_rows)
    print(json.dumps({k: v for k, v in res.items() if k != "per_mirna"}, indent=1))


if __name__ == "__main__":
    main()
