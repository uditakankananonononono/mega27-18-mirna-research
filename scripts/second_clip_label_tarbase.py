"""TarBase v9 label derivation (lane-18 item 8), per
docs/PREREG_SECOND_CLIP_20260927.md. PRE-OUTCOME: derives (miRNA, gene)
support labels from DIANA TarBase v9.0 AGO-CLIP interactions; computes no
model metric.

Row universe and feature columns are copied VERBATIM from
results/clip_rows_136.csv (the committed primary-endpoint rows). Only the
label column is new:

  label(miRNA, gene) = 1 iff the pair appears in TarBase v9 human
  interactions with experimental_method in the locked AGO-CLIP whitelist
  (HITS-CLIP, PAR-CLIP, CLASH, qCLASH; see
  data/encode_eclip/ENCODE_ECLIP_SOURCE.md). Gene-level support, mirroring
  the primary endpoint's ENCORI gene-level label.
"""
import csv, gzip, json, time
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / "data"
ROWS = ROOT / "results" / "clip_rows_136.csv"
OUT = ROOT / "results" / "second_clip_rows_tarbase.csv"
AUDIT = ROOT / "results" / "second_clip_label_tarbase_audit.json"
WHITELIST = {"HITS-CLIP", "PAR-CLIP", "CLASH", "qCLASH"}

t0 = time.time()
pairs = set()
with gzip.open(D / "encode_eclip" / "Homo_sapiens_TarBase-v9.tsv.gz", "rt") as fh:
    header = fh.readline().rstrip("\n").split("\t")
    im, ig, ime = header.index("mirna_name"), header.index("gene_name"), header.index("experimental_method")
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if p[ime] in WHITELIST:
            pairs.add((p[im], p[ig]))
print("whitelisted (miRNA,gene) pairs:", len(pairs), f"[{time.time()-t0:.0f}s]", flush=True)

audit = {}
n = 0
with open(ROWS) as fh, open(OUT, "w", newline="") as out:
    r = csv.reader(fh); header = next(r)
    w = csv.writer(out); w.writerow(header)
    for row in r:
        mi, g = row[0], row[1]
        label = int((mi, g) in pairs)
        a = audit.setdefault(mi, {"rows": 0, "pos": 0, "neg": 0})
        a["rows"] += 1; a["pos"] += label; a["neg"] += 1 - label
        w.writerow(row[:7] + [str(label)]); n += 1

elig = sorted(m for m, a in audit.items() if a["pos"] >= 50 and a["neg"] >= 50)
summary = {"rows_written": n, "mirnas_total": len(audit),
           "mirnas_eligible_50_50": len(elig), "eligible": elig,
           "whitelist": sorted(WHITELIST), "seconds": round(time.time() - t0, 1)}
json.dump({"summary": summary, "per_mirna": audit}, open(AUDIT, "w"), indent=1)
print("SUMMARY", summary, flush=True)
