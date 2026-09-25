"""RNAcentral sequence-integrity audit of every miRNA used in item 18 (tool 27).

PRE-REGISTERED (written before any comparison was run):
All distinct human miRNAs appearing in results/clip_rows_136.csv or
results/mirtarbase_rows.csv. For each: TargetScan miR_Family_Info mature sequence +
MiRBase Accession (MIMAT) -> RNAcentral API /api/v1/rna/?external_id=<MIMAT>.
Seed = nt 2-8 (the unit every scorer in this paper keys on).
Gates:
 G1 hsa-let-7a-5p (MIMAT0000062) returns UGAGGUAGUAGGUUGUAUAGUU.
 G2 >= 95% of miRNAs resolve to an RNAcentral sequence.
 G3 negative control: seed identity under a seeded derangement of the
    accession->sequence map is < 20% (the test can fail).
Hypotheses:
 H1 seed identity (TargetScan vs RNAcentral) = 100% of resolved miRNAs.
 H2 full-length identity >= 95% of resolved miRNAs.
"""
import csv, json, random, sys, time
from concurrent.futures import ThreadPoolExecutor
import pandas as pd, requests


def fetch(acc, tries=4):
    for i in range(tries):
        try:
            r = requests.get("https://rnacentral.org/api/v1/rna/", params={"external_id": acc, "format": "json"}, headers={"Accept": "application/json"}, timeout=60)
            r.raise_for_status()
            try:
                res = r.json().get("results", [])
            except ValueError:  # transport quirk: some accessions get an HTML page via python-requests; curl gets JSON
                import subprocess
                res = json.loads(subprocess.check_output(["curl", "-s", "-m", "60", r.url])).get("results", [])
            return [(x["rnacentral_id"], x["sequence"].upper().replace("T", "U")) for x in res]
        except Exception as e:
            err = e; time.sleep(2 * (i + 1))
    return None


def main(out="results/rnacentral_seq_audit.json"):
    fam = {}
    for row in csv.DictReader(open("data/miR_Family_Info.txt"), delimiter="\t"):
        if row["Species ID"] == "9606":
            fam[row["MiRBase ID"]] = (row["Mature sequence"].upper().replace("T", "U"), row["MiRBase Accession"])
    used = set(pd.read_csv("results/clip_rows_136.csv", usecols=["mirna"]).mirna) | set(pd.read_csv("results/mirtarbase_rows.csv", usecols=["mirna"]).mirna)
    names = sorted(m for m in used if m in fam)
    with ThreadPoolExecutor(4) as ex:
        got = dict(zip(names, ex.map(lambda m: fetch(fam[m][1]), names)))
    rows = []
    for m in names:
        ts, acc = fam[m]; hits = got[m] or []
        seqs = [s for _, s in hits]
        rows.append({"mirna": m, "mimat": acc, "targetscan_seq": ts, "rnacentral": hits[:3],
                     "resolved": bool(seqs), "seed_match": any(s[1:8] == ts[1:8] for s in seqs), "full_match": any(s == ts for s in seqs)})
    R = [r for r in rows if r["resolved"]]
    let7 = next((r for r in rows if r["mimat"] == "MIMAT0000062"), None)
    g1 = bool(let7 and any(s == "UGAGGUAGUAGGUUGUAUAGUU" for _, s in let7["rnacentral"]))
    perm = R[:]; rng = random.Random("derangement-2026")
    while True:
        idx = list(range(len(R))); rng.shuffle(idx)
        if all(i != j for i, j in enumerate(idx)): break
    ctrl = sum(any(s[1:8] == R[i]["targetscan_seq"][1:8] for _, s in R[j]["rnacentral"]) for i, j in enumerate(idx)) / len(R)
    seed = sum(r["seed_match"] for r in R) / len(R); full = sum(r["full_match"] for r in R) / len(R)
    J = {"tool": "RNAcentral REST API (rnacentral.org/api/v1/rna)", "n_used_mirnas": len(used), "n_with_targetscan_accession": len(names),
         "n_resolved": len(R),
         "gates": {"G1_let7a_exact": g1, "G2_resolved_frac": len(R) / len(names), "G2_pass": len(R) / len(names) >= 0.95,
                   "G3_deranged_seed_identity": ctrl, "G3_pass": ctrl < 0.20},
         "H1_seed_identity": {"frac": seed, "verdict": "CONFIRMED" if seed == 1.0 else "FALSIFIED"},
         "H2_full_identity": {"frac": full, "verdict": "CONFIRMED" if full >= 0.95 else "FALSIFIED"},
         "seed_mismatches": [r["mirna"] for r in R if not r["seed_match"]],
         "full_mismatches": [{"mirna": r["mirna"], "targetscan": r["targetscan_seq"], "rnacentral": [s for _, s in r["rnacentral"]]} for r in R if not r["full_match"]],
         "unresolved": [r["mirna"] for r in rows if not r["resolved"]], "not_in_targetscan_human": sorted(used - set(names)), "rows": rows}
    json.dump(J, open(out, "w"), indent=1)
    print(json.dumps({k: J[k] for k in ["n_used_mirnas", "n_with_targetscan_accession", "n_resolved", "gates", "H1_seed_identity", "H2_full_identity", "seed_mismatches", "full_mismatches", "unresolved", "not_in_targetscan_human"]}, indent=0)[:3000])


if __name__ == "__main__":
    main(*sys.argv[1:])
