#!/usr/bin/env python3
"""Figure + LaTeX section for the Reactome audit (committed JSON only)."""
import json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def fmtp(p):
    if p >= 0.001: return f"{p:.3f}"
    m, e = f"{p:.1e}".split("e"); return f"{m}\\times10^{{{int(e)}}}"


J = json.load(open("results/reactome_enrichment.json"))
P, ms, G, R = J["per_mirna"], J["panel"], J["gates"], J["posthoc_gate_review"]
lab = [m.replace("hsa-", "") for m in ms]; x = np.arange(len(ms)); w = 0.2
fig, ax = plt.subplots(figsize=(9.0, 3.1))
for i, (v, c, n) in enumerate([([P[m]["cnn"]["n_sig"] for m in ms], "#2980b9", "CNN top-200"), ([P[m]["sitecount"]["n_sig"] for m in ms], "#f39c12", "site-count top-200"),
                               ([P[m]["matched_mean"] for m in ms], "#16a085", "expression-matched random (mean of 3)"), ([P[m]["uniform_mean"] for m in ms], "#7f8c8d", "uniform random (mean of 3)")]):
    ax.bar(x + (i - 1.5) * w, v, w, color=c, label=n)
ax.set_xticks(x); ax.set_xticklabels(lab, rotation=45, ha="right", fontsize=7); ax.set_ylabel("Reactome pathways, FDR $<$ 0.05")
ax.legend(fontsize=6.5, frameon=False); ax.set_title(f"Reactome v{J['reactome_version']} AnalysisService (projection to human)", fontsize=9)
fig.tight_layout(); fig.savefig("paper/figs/fig_reactome.pdf")
h1, h2, h3 = J["H1_cnn_vs_uniform"], J["H2_cnn_vs_sitecount"], J["H3_cnn_vs_exprmatched"]
S = [r"\subsection{Reactome pathway audit with controls fixed in advance: two gates failed, verdicts inconclusive}",
     f"We repeated the pathway question with Reactome v{J['reactome_version']} (AnalysisService projection; \\texttt{{scripts/reactome\\_enrichment.py}};"
     r" \texttt{results/reactome\_enrichment.json}; hermetic tests), this time with the site-count and expression-matched controls pre-registered. Per miRNA we"
     r" count pathways at entities FDR $<0.05$ for CNN top-200, site-count top-200, three expression-matched and three uniform random sets.",
     f"Two of three gates failed as written. G1 required a pathway named ``Cell Cycle'' among the top three for a 20-gene cell-cycle control; the control returned"
     f" {G['G1_cellcycle_nsig']} significant pathways, but the top three are tied at the same FDR and none has that name (a post-hoc look finds {len(R['G1_cellcycle_named_pathways_anywhere'])} such"
     f" pathways lower in the list). G3 required fewer than 10\\% unmapped genes; the median is {100*G['G3_median_not_found']:.1f}\\%, because Reactome curates only part of"
     f" the genome. Both are design errors on our side, and we do not rewrite the gates. G2 passed (median uniform-set count {G['G2_median_uniform_nsig']:.0f}).",
     f"Under failed gates we report the numbers but call the audit inconclusive: CNN versus uniform random {h1['wins']}/{h1['losses']} ($p={fmtp(h1['p_one_sided'])}$), versus site count"
     f" {h2['wins']}/{h2['losses']} ($p={fmtp(h2['p_one_sided'])}$), versus expression-matched {h3['wins']}/{h3['losses']} ($p={fmtp(h3['p_one_sided'])}$). Most sets in every arm have zero"
     r" significant pathways, with a few large random-set outliers (Figure~\ref{fig:reactome}). None of this points toward pathway coherence, which agrees"
     r" with the g:Profiler result.",
     r"\begin{figure}[h]\centering\includegraphics[width=0.95\linewidth]{figs/fig_reactome.pdf}",
     r"\caption{Significant Reactome pathways per miRNA for four arms of 200 genes. Gates G1 and G3 failed, so these counts are descriptive.}\label{fig:reactome}\end{figure}"]
open("paper/reactome_sec.tex", "w").write("\n".join(S) + "\n")
print("ok")
