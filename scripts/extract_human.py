"""Stream-extract human (taxid 9606) data from TargetScan vert_80 archives."""
import csv, io, sys, zipfile
from pathlib import Path

D = Path(__file__).resolve().parents[1] / "data"

def utrs():
    out = open(D / "human_utrs.tsv", "w")
    out.write("transcript\tgene\tutr\n")
    n = 0
    with zipfile.ZipFile(D / "UTR_Sequences.txt.zip") as z, z.open("UTR_Sequences.txt") as fh:
        for line in io.TextIOWrapper(fh):
            p = line.rstrip("\n").split("\t")
            if len(p) < 5 or p[3] != "9606":
                continue
            seq = p[4].replace("-", "")
            if seq:
                out.write(f"{p[0]}\t{p[2]}\t{seq}\n"); n += 1
    out.close(); print("utrs", n, flush=True)

def sites():
    out = open(D / "human_sites.tsv", "w")
    out.write("transcript\tgene\tmirna\tsite_type\tstart\tend\tcontext_pp\n")
    n = 0
    with zipfile.ZipFile(D / "Conserved_Site_Context_Scores.txt.zip") as z:
        name = z.namelist()[0]
        with z.open(name) as fh:
            r = io.TextIOWrapper(fh); next(r)
            for line in r:
                p = line.rstrip("\n").split("\t")
                if p[3] != "9606" or p[8] == "NULL":
                    continue
                out.write(f"{p[2].split('.')[0]}\t{p[1]}\t{p[4]}\t{p[5]}\t{p[6]}\t{p[7]}\t{p[8]}\n"); n += 1
    out.close(); print("sites", n, flush=True)

if __name__ == "__main__":
    sites(); utrs()
