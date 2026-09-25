#!/usr/bin/env python3
"""Figure + LaTeX section for the STRING audit (reads committed JSON only)."""
import json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np



def fmtp(p):
    if p >= 0.001: return f"{p:.3f}"
    m, e = f"{p:.1e}".split("e"); return f"{m}\\times10^{{{int(e)}}}"


J = json.load(open("results/string_coherence.json"))
P, ms, G = J["per_mirna"], J["panel"], J["gates"]
PH = J["posthoc_H2_expression_matched_control"]
lab = [m.replace("hsa-", "") for m in ms]
ser = [([1e3 * P[m]["cnn"]["density"] for m in ms], "#2980b9", "CNN top-200"),
       ([1e3 * P[m]["clip"]["density"] for m in ms], "#8e44ad", "CLIP+ top-200"),
       ([1e3 * P[m]["rand_mean_density"] for m in ms], "#7f8c8d", "random (mean of 3)"),
       ([1e3 * PH["per_mirna"][m]["matched_mean_density"] for m in ms], "#16a085", "expression-matched random (post-hoc)")]
ex = [1e3 * P[m]["expr"]["density"] for m in ms]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(9.4, 3.4), gridspec_kw={"width_ratios": [3, 1]})
x = np.arange(len(ms)); w = 0.2
for i, (v, c, n) in enumerate(ser):
    a1.bar(x + (i - 1.5) * w, v, w, color=c, label=n)
a1.set_ylim(0, 1.45 * max(max(v) for v, _, _ in ser)); a1.set_xticks(x); a1.set_xticklabels(lab, rotation=45, ha="right", fontsize=7)
a1.set_ylabel("STRING edge density ($\\times10^{-3}$)"); a1.legend(fontsize=6.5, frameon=False, ncol=2)
a1.set_title("Target-set network density (STRING v12, score $\\geq$ 0.7)", fontsize=9)
a2.boxplot([ser[0][0], ser[1][0], ex], tick_labels=["CNN", "CLIP+", "expr"])
a2.set_title("Expression top-200", fontsize=9); a2.tick_params(labelsize=7)
fig.tight_layout(); fig.savefig("paper/figs/fig_string.pdf")

h1, h2, h3, t = J["H1_cnn_vs_random"], J["H2_clip_vs_random"], J["H3_cnn_vs_expression"], PH["clip_vs_matched"]
n = len(ms)
S = [r"\subsection{Network coherence audit: STRING edges in CLIP target sets are explained by expression (negative)}",
     f"The g:Profiler section tests term enrichment; STRING v{J['string_version']} tests whether a target set is more densely connected than chance"
     r" (\texttt{scripts/string\_coherence.py}; \texttt{results/string\_coherence.json}; hermetic tests). Same panel, sets and seeds as above; the statistic is"
     r" the edge density $\rho=E/\binom{n}{2}$ among mapped query proteins with STRING combined score $\geq 0.7$."
     f" Gates passed: the cell-cycle control gives {G['G1_cellcycle_edges']} edges (bar 50); median random density {G['G2_median_random_density']:.5f} (bar 0.01);"
     f" median mapped fraction {100*G['G3_median_mapped_frac']:.1f}\\%.",
     f"CNN top-200 sets are no denser than random (H1: {h1['wins']} wins, {h1['losses']} losses, $p={fmtp(h1['p_one_sided'])}$, falsified) and lose to the"
     f" expression top-200 control on {h3['losses']}/{n} (H3, falsified). CLIP-positive sets are denser than random on {h2['wins']}/{n}"
     f" ($p={fmtp(h2['p_one_sided'])}$), so pre-registered H2 is confirmed. Because CLIP detection depends on expression, we then added a post-hoc control that we"
     f" did not pre-register: random sets drawn from each miRNA's universe to match the CLIP set's HEK293 expression-decile histogram. Against it the CLIP"
     f" excess does not survive ({t['wins']} wins, {t['losses']} losses, one tie excluded, $p={fmtp(t['p_one_sided'])}$; Figure~\\ref{{fig:string}}).",
     r"So the one positive in this audit is explained by expression, the same confound that dominates the CLIP labels and the g:Profiler terms."
     r" We report H2 as confirmed-then-explained, not as evidence that miRNA targets form network modules.",
     r"\begin{figure}[h]\centering\includegraphics[width=0.95\linewidth]{figs/fig_string.pdf}",
     r"\caption{Left: STRING edge density per miRNA for CNN top-200, CLIP+ top-200, seeded random and expression-matched random sets. Right: distribution"
     r" across miRNAs with the expression top-200 control, which is several-fold denser.}\label{fig:string}\end{figure}"]
open("paper/string_sec.tex", "w").write("\n".join(S) + "\n")
print("ok")
