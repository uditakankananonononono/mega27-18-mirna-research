"""Second-CLIP-resource label derivation (lane-18 item 8), per
docs/PREREG_SECOND_CLIP_20260927.md. PRE-OUTCOME: derives (miRNA, gene)
support labels from ENCODE AGO2 eCLIP replicated peaks (ENCFF111TWF, K562,
GRCh38); computes no model metric.

Row universe and feature columns are copied VERBATIM from
results/clip_rows_136.csv (the committed primary-endpoint rows). Only the
label column is new:

  label(miRNA, gene) = 1 iff any candidate non-conserved seed-match site of
  that miRNA in that gene's UTR transcripts overlaps an AGO2 eCLIP peak
  (same chromosome, same strand). The candidate scan replicates
  scripts/clip_falsify.py exactly (same seed8/seed7 search, same A-at-+7/+6
  conditions, same conserved-transcript-region and conserved-seed-gene
  exclusions), so the scanned site set is identical to the one underlying
  the committed features.

Coordinate mapping: TargetScan UTR sequences are sequence-anchored onto
GENCODE v35 spliced transcripts (data/encode_eclip/utr_anchors_v35.json:
exact-prefix anchor; 24,787/28,352 UTRs anchored - 12,783 full-length exact,
11,752 exact to the annotated transcript end with a 3P-seq-style extension,
252 stop at an internal mismatch). Sites inside the anchored prefix map to
the genome through the GTF v35 exon chain; sites outside it are unmappable
and contribute no positive (conservative; counted in the audit).
"""
import csv, gzip, json, re, time
from bisect import bisect_right
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / "data"
ROWS = ROOT / "results" / "clip_rows_136.csv"
OUT = ROOT / "results" / "second_clip_rows.csv"
AUDIT = ROOT / "results" / "second_clip_label_audit.json"
COMP = str.maketrans("ACGU", "UGCA")
def rc(s): return s.translate(COMP)[::-1]

t0 = time.time()
pair_rows = defaultdict(list)
with open(ROWS) as fh:
    r = csv.reader(fh); header = next(r)
    for row in r:
        pair_rows[row[0]].append(row)
print("miRNAs:", len(pair_rows), "rows:", sum(len(v) for v in pair_rows.values()), flush=True)

fam = {}
with open(D / "miR_Family_Info.txt") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        if row["Species ID"] == "9606":
            fam[row["MiRBase ID"]] = row["Mature sequence"].upper().replace("T", "U")
seed_of = {m: s[1:8] for m, s in fam.items() if len(s) >= 8}

cons_genes, cons = {}, {}
with open(D / "human_sites.tsv") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        s = seed_of.get(row["mirna"])
        if s:
            cons_genes.setdefault(s, set()).add(row["gene"])
        cons.setdefault(row["transcript"], []).append((int(row["start"]), int(row["end"])))

def overlaps(t, lo, hi):
    return any(not (hi < a - 1 or lo > b - 1) for a, b in cons.get(t, []))

utrs, tid2gene = {}, {}
with open(D / "human_utrs.tsv") as fh:
    next(fh)
    for line in fh:
        t, g, s = line.rstrip("\n").split("\t", 2)
        tid2gene[t.split(".")[0]] = g
        utrs[t.split(".")[0]] = s.upper().replace("T", "U")
gene2tids = defaultdict(list)
for t, g in tid2gene.items():
    gene2tids[g].append(t)

anchors = {t: tuple(v) for t, v in json.load(open(D / "encode_eclip" / "utr_anchors_v35.json")).items()}

attr_re = re.compile(r'transcript_id "(ENST\d+)')
tmp_ex, meta = defaultdict(list), {}
with gzip.open(D / "encode_eclip" / "gencode.v35.annotation.gtf.gz", "rt") as fh:
    for line in fh:
        if line.startswith("#"): continue
        m = attr_re.search(line)
        if not m: continue
        t = m.group(1)
        if t not in anchors: continue
        p = line.split("\t")
        if p[2] == "transcript":
            meta[t] = (p[0], p[6])
        elif p[2] == "exon":
            tmp_ex[t].append((int(p[3]), int(p[4])))
exon_chain = {}
for t, ex in tmp_ex.items():
    if t in meta:
        chrom, strand = meta[t]
        ex.sort(reverse=(strand == "-"))
        exon_chain[t] = (chrom, strand, ex)
print("exon chains:", len(exon_chain), flush=True)

pk = defaultdict(list)
with gzip.open(D / "encode_eclip" / "ENCFF111TWF.bed.gz", "rt") as fh:
    for line in fh:
        p = line.split("\t")
        pk[(p[0], p[5])].append((int(p[1]), int(p[2])))
peaks = {}
for k, v in pk.items():
    v.sort()
    peaks[k] = (tuple(s for s, _ in v), tuple(e for _, e in v))
print("peaks:", sum(len(v) for v in pk.values()), flush=True)

def peak_hit(chrom, strand, lo, hi):
    e = peaks.get((chrom, strand))
    if not e: return False
    starts, ends = e
    i = bisect_right(ends, lo)
    return i < len(starts) and starts[i] < hi

def site_segments(t, h, pos, length=7):
    chrom, strand, exons = exon_chain[t]
    sp = h + pos
    cum, segs, remaining = 0, [], length
    for (a, b) in exons:
        L = b - a + 1
        if sp < cum + L:
            off = sp - cum
            take = min(L - off, remaining)
            if strand == "+":
                segs.append((chrom, strand, a + off - 1, a + off - 1 + take))
            else:
                segs.append((chrom, strand, b - off - take, b - off))
            remaining -= take
            if remaining == 0: break
        cum += L
    return segs

def gene_label(mi, g, excl, seed8, seed7):
    """(label, n_candidate_sites, n_unmappable) for one (miRNA, gene) pair."""
    if g in excl:
        return 0, 0, 0
    nsite = unmap = 0
    for t in gene2tids.get(g, []):
        u = utrs.get(t)
        if u is None: continue
        cands = []
        i = u.find(seed8)
        while i >= 0:
            if not overlaps(t, i, i + 6): cands.append(i)
            i = u.find(seed8, i + 1)
        i = u.find(seed7)
        while i >= 0:
            if i + 6 < len(u) and u[i + 6] == "A" and u[i:i + 7] != seed8 and not overlaps(t, i, i + 5):
                cands.append(i)
            i = u.find(seed7, i + 1)
        nsite += len(cands)
        anc = anchors.get(t)
        if anc is None or t not in exon_chain:
            unmap += len(cands); continue
        h, plen, lu, tl = anc
        for i in cands:
            if i + 6 >= plen:
                unmap += 1; continue
            for (chrom, strand, lo, hi) in site_segments(t, h, i):
                if peak_hit(chrom, strand, lo, hi):
                    return 1, nsite, unmap
    return 0, nsite, unmap

audit = {}
n_written = 0
with open(OUT, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(header)
    for mi, rows in pair_rows.items():
        seq = fam.get(mi)
        if seq is None:
            audit[mi] = {"skipped": "no mature sequence", "rows": len(rows)}
            continue
        excl = cons_genes.get(seq[1:8], set())
        seed8, seed7 = rc(seq[1:8]), rc(seq[1:7])
        pos = unmap_t = nsite_t = 0
        for row in rows:
            g = row[1]
            label, nsite, unmap = gene_label(mi, g, excl, seed8, seed7)
            nsite_t += nsite; unmap_t += unmap; pos += label
            w.writerow(row[:7] + [str(label)])
            n_written += 1
        audit[mi] = {"rows": len(rows), "pos": pos, "neg": len(rows) - pos,
                     "candidate_sites": nsite_t, "unmappable_sites": unmap_t}
        print(mi, audit[mi], f"[{time.time()-t0:.0f}s]", flush=True)

elig = {m: a for m, a in audit.items() if a.get("pos", 0) >= 50 and a.get("neg", 0) >= 50}
summary = {"rows_written": n_written, "mirnas_total": len(pair_rows),
           "mirnas_eligible_50_50": len(elig),
           "eligible": sorted(elig),
           "seconds": round(time.time() - t0, 1)}
json.dump({"summary": summary, "per_mirna": audit}, open(AUDIT, "w"), indent=1)
print("SUMMARY", summary, flush=True)
