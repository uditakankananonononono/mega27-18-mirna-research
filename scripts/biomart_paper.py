#!/usr/bin/env python3
"""LaTeX section for the BioMart paralog test (committed JSON only)."""
import json


def fmtp(p):
    if p >= 0.001: return f"{p:.3f}"
    m, e = f"{p:.1e}".split("e"); return f"{m}\\times10^{{{int(e)}}}"


J = json.load(open("results/biomart_paralogs.json")); G = J["gates"]; g1 = G["G1_ts_vs_random"]
H = [("B1 more paralogs than uniform", "B1_more_paralogs_vs_uniform"), ("B2 more paralogs than expression-matched", "B2_vs_expr"),
     ("B3 LOEUF higher than paralog-matched", "B3_loeuf_vs_paralog_matched")]
row = lambda n, k: f"{n} & {J[k]['wins']}/{J[k]['losses']} & ${fmtp(J[k]['p_one_sided'])}$ & {J[k]['median_diff']:+.3f} & {J[k]['verdict']} \\\\"
S = [r"\subsection{Paralog buffering (Ensembl BioMart): no paralog excess; both gates failed}",
     r"Genes with paralogs can be buffered against dosage changes, so one explanation of the depletion is that CNN candidates have more paralogs."
     r" We counted human paralogues per gene from Ensembl BioMart (\texttt{biomart\_paralogs}, committed before running; 43 miRNAs; statistic median"
     r" $\log_2(1+\text{paralogs})$). Both gates failed:"
     f" {100*G['G2_coverage']:.1f}\\% of universe symbols were found (bar 90\\%), and conserved targets did not differ from random genes"
     f" ({g1['wins']}/{g1['losses']}, two-sided $p={fmtp(g1['p_two_sided'])}$).",
     r"\begin{table}[h]\centering\small\begin{tabular}{lcccc}\hline Pre-registered test (G1, G2 failed) & wins/losses & $p$ & median diff. & verdict\\\hline",
     *[row(n, k) for n, k in H], r"\hline\end{tabular}\caption{Paralog counts of CNN top-200 sets, 43 miRNAs.}\label{tab:biomart}\end{table}",
     r"CNN candidates have no paralog excess (B1, B2 falsified), and the LOEUF depletion stays at 43/0 against paralog-matched sets. With both gates"
     r" failed we give these results little weight, but nothing here supports paralog buffering as the explanation."]
open("paper/biomart_sec.tex", "w").write("\n".join(S) + "\n"); print("ok")
