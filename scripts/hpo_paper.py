#!/usr/bin/env python3
"""LaTeX section for the HPO Mendelian disease-gene test (committed JSON only)."""
import json


def fmtp(p):
    if p >= 0.001: return f"{p:.3f}"
    m, e = f"{p:.1e}".split("e"); return f"{m}\\times10^{{{int(e)}}}"


J = json.load(open("results/hpo_disease_genes.json")); G = J["gates"]; g1 = G["G1_ts_vs_random"]
H = [("M1 CNN fewer than uniform", J["M1_vs_uniform"]), ("M2 fewer than expression-matched", J["M2_vs_expr"]), ("M3 fewer than UTR-length-matched", J["M3_vs_utrlen"])]
row = lambda n, h: f"{n} & {h['wins']}/{h['losses']} & ${fmtp(h['p_one_sided'])}$ & {h['pooled_cnn']} / {h['pooled_comparator']:.0f} ({h['pooled_ratio']:.2f}) & {h['verdict']} \\\\"
S = [r"\subsection{Mendelian disease genes (HPO): depletion explained by UTR length}",
     f"HPO gene--disease annotations give a fourth, curated label ({J['n_mendelian_genes']:,} genes with a Mendelian association;"
     r" \texttt{hpo\_disease\_genes}, committed before running; same 43 miRNAs and controls). Gates passed: disease-gene prevalence in the pooled universe is"
     f" {100*G['G2_prevalence']:.1f}\\%, and TargetScan conserved top-200 sets hold more disease genes than random UTR genes ({g1['wins']}/{g1['losses']},"
     f" $p={fmtp(g1['p_one_sided'])}$).",
     r"\begin{table}[h]\centering\small\begin{tabular}{lcccc}\hline Pre-registered test & wins/losses & $p$ & pooled CNN / control (ratio) & verdict\\\hline",
     *[row(n, h) for n, h in H], r"\hline\end{tabular}\caption{HPO Mendelian disease genes in CNN top-200 sets, 43 miRNAs.}\label{tab:hpo}\end{table}",
     f"CNN sets carry fewer disease genes than uniform and expression-matched sets, but against UTR-length-matched sets the per-miRNA split is"
     f" {H[2][1]['wins']}/{H[2][1]['losses']} (M3 falsified), even though the pooled ratio stays at {H[2][1]['pooled_ratio']:.2f}. This agrees with ClinGen:"
     r" both curated labels support depletion relative to expression, and neither supports it beyond UTR length. The two population-derived scores (LOEUF,"
     r" rCNV) do. We state the finding in its supported form: \emph{CNN top-ranked non-conserved candidates are depleted of dosage-sensitive and disease genes"
     r" relative to expression-matched genes}; whether this goes beyond UTR length depends on the label and remains open."]
open("paper/hpo_sec.tex", "w").write("\n".join(S) + "\n"); print("ok")
