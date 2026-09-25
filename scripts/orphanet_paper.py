#!/usr/bin/env python3
"""LaTeX section for the Orphanet disease-gene test (committed JSON only)."""
import json


def fmtp(p):
    if p >= 0.001: return f"{p:.3f}"
    m, e = f"{p:.1e}".split("e"); return f"{m}\\times10^{{{int(e)}}}"


J = json.load(open("results/orphanet_genes.json")); G = J["gates"]; g1 = G["G1_ts_vs_random"]
H = [("O1 disease genes, vs expression-matched", "O1_dis_vs_expr"), ("O2 disease genes, vs publication-matched", "O2_dis_vs_pub"),
     ("O3 LoF genes, vs publication-matched", "O3_lof_vs_pub"), ("O4 LoF genes, vs UTR-length-matched", "O4_lof_vs_utrlen")]
md = lambda x: x if abs(x) >= 5e-4 else 0.0
row = lambda n, k: f"{n} & {J[k]['wins']}/{J[k]['losses']} & ${fmtp(J[k]['p_one_sided'])}$ & {md(J[k]['median_diff']):+.3f} & {J[k]['verdict']} \\\\"
S = [r"\subsection{Orphanet rare-disease genes: depletion vs expression-matched only; gone against publication-matched sets}",
     r"Orphanet curates disorder--gene links by hand and records the mechanism, so it gives a further disease label and an explicit loss-of-function"
     r" label~\cite{orphanet}. From Orphadata \texttt{en\_product6} we took assessed germline disease-causing links"
     f" ({J['n_disease_genes']} genes) and the loss-of-function subset ({J['n_lof_genes']} genes); statistic = fraction of the top-200"
     r" (\texttt{orphanet\_genes}, committed before running; 43 miRNAs). Both gates passed: conserved TargetScan targets carry more disease genes"
     f" than random UTR genes ({g1['wins']}/{g1['losses']}, $p={fmtp(g1['p_one_sided'])}$), and the parse counts met the bars.",
     r"\begin{table}[h]\centering\small\begin{tabular}{lcccc}\hline Pre-registered test (CNN lower) & wins/losses & $p$ & median diff. & verdict\\\hline",
     *[row(n, k) for n, k in H], r"\hline\end{tabular}\caption{Orphanet disease-gene fractions in CNN top-200 sets, 43 miRNAs.}\label{tab:orphanet}\end{table}",
     r"CNN candidates carry fewer disease genes than expression-matched sets (O1), but the gap disappears against publication-matched sets (O2), and"
     r" the explicit loss-of-function label shows no depletion against publication- or UTR-length-matched sets (O3, O4). This repeats the HPO pattern:"
     r" on curated labels the depletion is explained by study bias, while the LoF-count labels (LOEUF, $s_\mathrm{het}$) keep it."]
open("paper/orphanet_sec.tex", "w").write("\n".join(S) + "\n"); print("ok")
