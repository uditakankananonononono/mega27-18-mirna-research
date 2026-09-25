#!/usr/bin/env python3
"""LaTeX section for the DECIPHER HI-prediction test (committed JSON only)."""
import json


def fmtp(p):
    if p >= 0.001: return f"{p:.3f}"
    m, e = f"{p:.1e}".split("e"); return f"{m}\\times10^{{{int(e)}}}"


J = json.load(open("results/decipher_hi.json")); G = J["gates"]; g1 = G["G1_ts_vs_random"]
H = [("H1 lower than uniform", "H1_vs_uniform"), ("H2 lower than expression-matched", "H2_vs_expr"), ("H3 lower than UTR-length-matched", "H3_vs_utrlen"),
     ("H4 lower than publication-matched", "H4_vs_pubmatched")]
row = lambda n, k: f"{n} & {J[k]['wins']}/{J[k]['losses']} & ${fmtp(J[k]['p_one_sided'])}$ & {J[k]['median_diff']:+.4f} & {J[k]['verdict']} \\\\"
S = [r"\subsection{DECIPHER haploinsufficiency predictions: depletion beyond expression and study bias, not beyond UTR length}",
     f"DECIPHER's haploinsufficiency predictions (version 3; {J['n_scored_genes']:,} genes) come from a model of genomic, evolutionary and network"
     r" features that does not use gnomAD LOEUF (\texttt{decipher\_hi}, committed before running; 43 miRNAs; publication-matched control added)."
     f" Gates passed: {100*G['G2_coverage']:.1f}\\% coverage, and conserved top-200 sets score higher than random UTR genes"
     f" ({g1['wins']}/{g1['losses']}, $p={fmtp(g1['p_one_sided'])}$).",
     r"\begin{table}[h]\centering\small\begin{tabular}{lcccc}\hline Pre-registered test & wins/losses & $p$ & median diff. & verdict\\\hline",
     *[row(n, k) for n, k in H], r"\hline\end{tabular}\caption{DECIPHER HI scores of CNN top-200 sets, 43 miRNAs.}\label{tab:decipher}\end{table}",
     r"The CNN depletion holds against expression-matched and publication-matched sets but not against UTR-length-matched sets (H3 falsified;"
     f" median difference {J['H3_vs_utrlen']['median_diff']:+.4f}). The tally for ``beyond UTR length'' is now: supported by LOEUF and the rCNV scores,"
     r" not supported by DECIPHER, ClinGen or HPO (ClinVar and DepMap uninformative). We therefore do not claim the depletion goes beyond UTR length."
     r" The claim that survives every informative label is the weaker one: CNN top-ranked non-conserved candidates are less dosage-sensitive than"
     r" expression-matched genes."]
open("paper/decipher_sec.tex", "w").write("\n".join(S) + "\n"); print("ok")
