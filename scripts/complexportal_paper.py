#!/usr/bin/env python3
"""LaTeX section for the Complex Portal subunit test (committed JSON only)."""
import json


def fmtp(p):
    if p >= 0.001: return f"{p:.3f}"
    m, e = f"{p:.1e}".split("e"); return f"{m}\\times10^{{{int(e)}}}"


J = json.load(open("results/complexportal_subunits.json")); G = J["gates"]; g1 = G["G1_ts_vs_random"]
H = [("X1 fewer than uniform", J["X1_vs_uniform"]), ("X2 fewer than expression-matched", J["X2_vs_expr"]), ("X3 fewer than UTR-length-matched", J["X3_vs_utrlen"])]
row = lambda n, h: f"{n} & {h['wins']}/{h['losses']} & ${fmtp(h['p_one_sided'])}$ & {h['pooled_cnn']} / {h['pooled_comparator']:.0f} ({h['pooled_ratio']:.2f}) & {h['verdict']} \\\\"
S = [r"\subsection{Protein-complex subunits (Complex Portal): the dosage-balance prediction fails}",
     r"The dosage-balance hypothesis holds that subunits of protein complexes are dosage-sensitive, so the dosage finding predicts that CNN sets"
     f" carry fewer subunit genes. We used EBI Complex Portal ({J['n_complexes']:,} human complexes, {J['n_subunit_genes']:,} subunit genes after mapping"
     r" UniProt accessions; \texttt{complexportal\_subunits}, committed before running; same 43 miRNAs and controls)."
     f" Gates passed: {100*G['G2_map_rate']:.1f}\\% of accessions mapped, prevalence {100*G['G3_prevalence']:.1f}\\%, and conserved top-200 sets differ from"
     f" random UTR genes ({g1['wins']}/{g1['losses']} more subunits, two-sided $p={fmtp(g1['p_two_sided'])}$).",
     r"\begin{table}[h]\centering\small\begin{tabular}{lcccc}\hline Pre-registered test & wins/losses & $p$ & pooled CNN / control (ratio) & verdict\\\hline",
     *[row(n, h) for n, h in H], r"\hline\end{tabular}\caption{Complex Portal subunit genes in CNN top-200 sets, 43 miRNAs.}\label{tab:complexportal}\end{table}",
     r"The depletion against uniform sets disappears once expression or UTR length is matched (X2, X3 falsified; pooled ratios"
     f" {H[1][1]['pooled_ratio']:.2f} and {H[2][1]['pooled_ratio']:.2f}). So the CNN does not avoid complex subunits beyond what expression explains, and"
     r" the dosage-balance route does not account for the LOEUF and rCNV depletion."]
open("paper/complexportal_sec.tex", "w").write("\n".join(S) + "\n"); print("ok")
