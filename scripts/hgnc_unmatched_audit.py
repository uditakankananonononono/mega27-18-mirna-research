"""Post-hoc (not pre-registered) audit of the G1 coverage failure in hgnc_families."""
import json, sys
import pandas as pd
sys.path.insert(0, "scripts")
from gprofiler_enrichment import PANEL
from hgnc_families import load_hgnc

h = load_hgnc(); a = pd.read_csv("data/hgnc/hgnc_complete_set.txt", sep="\t", dtype=str, usecols=["symbol", "prev_symbol", "alias_symbol"])
panel = PANEL + json.load(open("results/gnomad_replication.json"))["panel"]
d = pd.read_csv("results/clip_rows_136.csv", usecols=["mirna", "gene"]); u = pd.Series(d[d.mirna.isin(panel)].gene.unique())
m = u[~u.isin(h.index)]; cur = set(a.symbol)
prev = {x for s in a.prev_symbol.dropna() for x in s.split("|")} - cur
al = {x for s in a.alias_symbol.dropna() for x in s.split("|")} - cur - prev
J = {"universe": int(len(u)), "unmatched": int(len(m)), "approved_non_protein_coding": int(m.isin(cur).sum()),
     "previous_symbol": int(m.isin(prev).sum()), "alias_only": int(m.isin(al).sum())}
J["unresolved"] = J["unmatched"] - J["approved_non_protein_coding"] - J["previous_symbol"] - J["alias_only"]
json.dump(J, open("results/hgnc_unmatched_audit.json", "w"), indent=1); print(J)
