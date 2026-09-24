"""Extend ENCORI CLIP panel from 23 toward ~100+ miRNAs (dataset-count expansion).
Fetches per-miRNA CLIP-supported target TSVs from the ENCORI API (same endpoint
as clip_falsify.py), keeps those with >=20 CLIP+ genes, records honest counts.
"""
import json, pathlib, time, urllib.request

D = pathlib.Path("data/encori"); D.mkdir(exist_ok=True)
CANDIDATES = """hsa-miR-9-5p hsa-miR-10a-5p hsa-miR-10b-5p hsa-miR-15a-5p hsa-miR-15b-5p
hsa-miR-18a-5p hsa-miR-19b-3p hsa-miR-20a-5p hsa-miR-22-3p hsa-miR-25-3p hsa-miR-28-5p
hsa-miR-31-5p hsa-miR-32-5p hsa-miR-34a-5p hsa-miR-93-5p hsa-miR-96-5p hsa-miR-98-5p
hsa-miR-99a-5p hsa-miR-100-5p hsa-miR-101-3p hsa-miR-106a-5p hsa-miR-106b-5p hsa-miR-107
hsa-miR-122-5p hsa-miR-126-3p hsa-miR-127-5p hsa-miR-128-3p hsa-miR-129-5p hsa-miR-130b-3p
hsa-miR-132-3p hsa-miR-133a-3p hsa-miR-133b hsa-miR-134-5p hsa-miR-135b-5p hsa-miR-138-5p
hsa-miR-139-5p hsa-miR-140-5p hsa-miR-141-3p hsa-miR-143-3p hsa-miR-144-3p hsa-miR-145-5p
hsa-miR-146a-5p hsa-miR-146b-5p hsa-miR-148a-3p hsa-miR-149-5p hsa-miR-150-5p hsa-miR-152-3p
hsa-miR-153-3p hsa-miR-182-5p hsa-miR-183-5p hsa-miR-184 hsa-miR-185-5p hsa-miR-186-5p
hsa-miR-187-3p hsa-miR-190a-5p hsa-miR-191-5p hsa-miR-192-5p hsa-miR-193a-5p hsa-miR-194-5p
hsa-miR-195-5p hsa-miR-196a-5p hsa-miR-197-3p hsa-miR-200a-3p hsa-miR-200b-3p hsa-miR-203a-3p
hsa-miR-204-5p hsa-miR-205-5p hsa-miR-206 hsa-miR-210-3p hsa-miR-214-3p hsa-miR-215-5p
hsa-miR-216a-5p hsa-miR-217-5p hsa-miR-218-5p hsa-miR-222-3p hsa-miR-223-3p hsa-miR-224-5p
hsa-miR-296-5p hsa-miR-302a-3p hsa-miR-320a hsa-miR-324-5p hsa-miR-328-3p hsa-miR-335-5p
hsa-miR-338-3p hsa-miR-339-5p hsa-miR-342-3p hsa-miR-345-5p hsa-miR-361-5p hsa-miR-363-3p
hsa-miR-375 hsa-miR-376a-3p hsa-miR-378a-3p hsa-miR-409-3p hsa-miR-424-5p hsa-miR-425-5p
hsa-miR-429 hsa-miR-451a hsa-miR-484 hsa-miR-486-5p hsa-miR-494-3p hsa-miR-495-3p
hsa-miR-500a-3p hsa-miR-503-5p hsa-miR-532-5p hsa-miR-545-3p hsa-miR-574-5p hsa-miR-590-5p
hsa-miR-615-3p hsa-miR-625-5p hsa-miR-629-5p hsa-miR-652-3p hsa-miR-660-5p hsa-miR-671-5p
hsa-miR-708-5p hsa-miR-744-5p hsa-miR-766-3p hsa-miR-885-5p hsa-miR-92b-3p""".split()

def fetch(mir):
    f = D / f"{mir}.tsv"
    if not f.exists():
        url = ("https://rnasysu.com/encori/api/miRNATarget/?assembly=hg38&geneType=mRNA&miRNA="
               f"{mir}&clipExpNum=1&degraExpNum=0&pancancerNum=0&programNum=0&program=None&target=all&cellType=all")
        for _ in range(3):
            try:
                f.write_bytes(urllib.request.urlopen(url, timeout=90).read()); break
            except Exception as e:
                print("retry", mir, repr(e)[:80], flush=True); time.sleep(4)
    genes = set()
    if f.exists():
        for line in f.read_text().splitlines():
            if line.startswith("#") or line.startswith("miRNAid"):
                continue
            p = line.split("\t")
            if len(p) > 11 and p[1] == mir:
                genes.add(p[3])
    return genes

kept, dropped = {}, {}
for mir in CANDIDATES:
    g = fetch(mir)
    if len(g) >= 20:
        kept[mir] = len(g)
    else:
        dropped[mir] = len(g)
    print(mir, len(g), "KEEP" if len(g) >= 20 else "drop", flush=True)
    time.sleep(0.5)

existing = sorted(p.stem for p in D.glob("*.tsv") if p.stem not in kept and p.stem not in dropped)
out = {"min_genes": 20, "kept_new": kept, "dropped": dropped,
       "n_kept_new": len(kept), "n_prior": len(existing),
       "total_panel": len(kept) + len(existing)}
json.dump(out, open("results/panel_extension.json", "w"), indent=1)
print(json.dumps({k: v for k, v in out.items() if k not in ("kept_new", "dropped")}, indent=1))
