"""mirtarget-score - score miRNA->UTR target sites with the trained DuplexCNN.

Given a miRNA (name resolved via miRBase family info, or a raw mature sequence)
and 3'UTRs (gene symbols from the bundled TargetScan UTR set, or a FASTA),
find 8mer/7mer-A1 seed matches, score each with DuplexCNN + context features
(training normalization from results/feature_norm.json), and emit ranked sites.

Usage:
  mirtarget-score --mirna hsa-miR-7-5p --gene CDR1as --json out.json
  mirtarget-score --seq UGGAAGACUAGUGAUUUUGUUGU --utr-fasta utrs.fa
"""
import argparse, json, pathlib, sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
COMP = str.maketrans("ACGU", "UGCA")


def rc(s):
    return s.translate(COMP)[::-1]


def find_sites(utr, seed8, seed7):
    out = []
    i = utr.find(seed8)
    while i >= 0:
        out.append((i, 3))  # 8mer
        i = utr.find(seed8, i + 1)
    i = utr.find(seed7)
    while i >= 0:
        if i + 6 < len(utr) and utr[i + 6] == "A" and utr[i:i + 7] != seed8:
            out.append((i, 2))  # 7mer-A1
        i = utr.find(seed7, i + 1)
    return out


def score_sites(net, seq, sites_by_name, mu, sd):
    from .data import one_hot, site_features, WIN, MIR_LEN
    import torch
    rows = []
    for name, (utr, sites) in sites_by_name.items():
        for pos, stype in sites:
            win = "".join(utr[k] if 0 <= k < len(utr) else "N"
                          for k in range(pos - WIN // 2 + 3, pos - WIN // 2 + 3 + WIN))
            feat = (site_features(utr, pos + 1, pos + 7) - mu) / sd
            with torch.no_grad():
                v = float(net(torch.from_numpy(one_hot(win, WIN)).unsqueeze(0),
                              torch.from_numpy(one_hot(seq, MIR_LEN)).unsqueeze(0),
                              torch.tensor([stype]), torch.from_numpy(feat.astype(np.float32)).unsqueeze(0)))
            rows.append({"target": name, "pos": pos, "site_type": {3: "8mer", 2: "7mer-A1"}[stype],
                         "score": round(v, 4)})
    rows.sort(key=lambda r: r["score"])
    return rows


def main(argv=None):
    p = argparse.ArgumentParser(prog="mirtarget-score",
        description="Score miRNA->UTR seed-match sites with DuplexCNN (context-normalized).")
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--mirna", help="miRBase ID, e.g. hsa-miR-7-5p (human only)")
    src.add_argument("--seq", help="raw mature miRNA sequence (ACGU)")
    tgt = p.add_mutually_exclusive_group(required=True)
    tgt.add_argument("--gene", nargs="+", help="gene symbols in the bundled TargetScan UTR set")
    tgt.add_argument("--utr-fasta", help="FASTA of 3'UTR sequences")
    p.add_argument("--top", type=int, default=25)
    p.add_argument("--json", help="write ranked sites JSON here")
    p.add_argument("--data-dir", default=str(ROOT / "data"))
    a = p.parse_args(argv)

    from .data import load_mirnas, load_utrs
    data = pathlib.Path(a.data_dir)
    if a.mirna:
        if not (data / "miR_Family_Info.txt").is_file():
            sys.exit("TargetScan miR_Family_Info.txt is absent. Supply the original TargetScan file in --data-dir, or use --seq with --utr-fasta.")
        mirs = load_mirnas(data / "miR_Family_Info.txt")
        if a.mirna not in mirs:
            sys.exit(f"unknown human miRNA {a.mirna!r} (not in miR_Family_Info.txt)")
        seq = mirs[a.mirna]
    else:
        seq = a.seq.upper().replace("T", "U")
    if len(seq) < 8:
        sys.exit("mature sequence must be >= 8 nt")

    sites_by_name = {}
    if a.gene:
        if not (data / "human_utrs.tsv").is_file():
            sys.exit("TargetScan human_utrs.tsv is absent. Extract original UTR_Sequences.txt.zip with scripts/extract_human.py, or use --utr-fasta.")
        utrs = load_utrs(data / "human_utrs.tsv")
        gene2tid = {}
        with open(data / "human_utrs.tsv") as fh:
            next(fh)
            for line in fh:
                t, g, _ = line.rstrip("\n").split("\t", 2)
                gene2tid.setdefault(g, []).append(t.split(".")[0])
        seed8, seed7 = rc(seq[1:8]), rc(seq[1:7])
        for g in a.gene:
            for t in gene2tid.get(g, []):
                u = utrs.get(t)
                if u:
                    s = find_sites(u, seed8, seed7)
                    if s:
                        sites_by_name[f"{g}|{t}"] = (u, s)
    else:
        name, buf = None, []
        for line in open(a.utr_fasta):
            line = line.strip()
            if line.startswith(">"):
                if name:
                    u = "".join(buf).upper().replace("T", "U")
                    s = find_sites(u, rc(seq[1:8]), rc(seq[1:7]))
                    if s:
                        sites_by_name[name] = (u, s)
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line)
        if name:
            u = "".join(buf).upper().replace("T", "U")
            s = find_sites(u, rc(seq[1:8]), rc(seq[1:7]))
            if s:
                sites_by_name[name] = (u, s)

    if not sites_by_name:
        print(json.dumps({"tool": "mirtarget-score", "sites": [],
                          "note": "no canonical seed matches found"}))
        return 0

    import torch
    from .model import DuplexCNN
    net = DuplexCNN(n_feat=4)
    net.load_state_dict(torch.load(ROOT / "results" / "duplex_cnn.pt",
                                   map_location="cpu", weights_only=False))
    net.eval()
    norm = json.loads((ROOT / "results" / "feature_norm.json").read_text())
    mu, sd = np.array(norm["mu"], np.float32), np.array(norm["sd"], np.float32)
    rows = score_sites(net, seq, sites_by_name, mu, sd)[: a.top]
    out = {"tool": "mirtarget-score",
           "mirna": a.mirna or "custom", "mature_seq": seq,
           "model": "DuplexCNN context++ (held-out-gene Pearson 0.8121)",
           "score_note": "lower = stronger predicted repression (context++ scale)",
           "sites": rows}
    if a.json:
        pathlib.Path(a.json).write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
