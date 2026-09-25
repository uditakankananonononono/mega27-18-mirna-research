#!/usr/bin/env python3
"""Figure + LaTeX section for the Enrichr reverse-lookup audit (committed JSON only)."""
import json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def fmtp(p):
    if p >= 0.001: return f"{p:.3f}"
    m, e = f"{p:.1e}".split("e"); return f"{m}\\times10^{{{int(e)}}}"


J = json.load(open("results/enrichr_reverse_lookup.json"))
P, ms, G = J["per_mirna"], J["panel"], J["gates"]
SC, EM = J["posthoc_sitecount_control"], J["posthoc_expression_matched_control"]
lab = [m.replace("hsa-", "") for m in ms]; x = np.arange(len(ms)); w = 0.16
arms = [([P[m]["cnn"]["rank"] for m in ms], "#2980b9", "CNN top-200"), ([P[m]["clip"]["rank"] for m in ms], "#8e44ad", "CLIP+ top-200"),
        ([SC["per_mirna"][m]["sitecount"]["rank"] for m in ms], "#f39c12", "site-count top-200 (post-hoc)"),
        ([EM["per_mirna"][m]["matched_mean"] for m in ms], "#16a085", "expression-matched random (post-hoc)"),
        ([P[m]["rand_mean_rank"] for m in ms], "#7f8c8d", "random (mean of 3)")]
fig, ax = plt.subplots(figsize=(9.2, 3.4))
for i, (v, c, n) in enumerate(arms):
    ax.bar(x + (i - 2) * w, np.log10((J["library_terms"] + 1) / np.array(v)), w, color=c, label=n)
ax.set_ylabel("$\\log_{10}$(3241 / own-term rank)"); ax.set_ylim(0, 4.6)
ax.set_xticks(x); ax.set_xticklabels(lab, rotation=45, ha="right", fontsize=7); ax.legend(fontsize=6.5, frameon=False, ncol=3, loc="upper right")
ax.set_title("Enrichr reverse lookup (miRTarBase_2017): does a target set name its own miRNA? (taller = better)", fontsize=9)
fig.tight_layout(); fig.savefig("paper/figs/fig_enrichr.pdf")
h1, h2, h3 = J["H1_cnn_vs_random"], J["H2_clip_vs_random"], J["H3_cnn_vs_expression"]
t1, t2 = SC["cnn_vs_sitecount"], EM["cnn_vs_matched"]
med = lambda k: float(np.median([P[m][k]["rank"] for m in ms]))
S = [r"\subsection{Reverse lookup (Enrichr): CNN non-conserved top-200 sets name their own miRNA, beyond site count and expression (positive)}",
     r"A stricter use of a target list is to ask whether it identifies its miRNA among thousands. We submit each set to the Enrichr API against the"
     f" miRTarBase\\_2017 library ({J['library_terms']:,} miRNA terms; \\texttt{{scripts/enrichr\\_reverse\\_lookup.py}}; \\texttt{{results/enrichr\\_reverse\\_lookup.json}};"
     r" hermetic tests) and record the rank of the query miRNA's own term by Enrichr's Fisher $p$-value. Panel, sets and seeds match the g:Profiler section."
     f" Gates passed: the library's own miR-21-5p set returns miR-21-5p at rank {G['G1_miR21_rank']}; the median random-set own-term rank is {G['G2_median_random_rank']:.0f};"
     r" all 13 miRNAs are library terms.",
     f"H1 is confirmed: CNN top-200 sets rank their own miRNA better than random for {h1['wins']}/13 miRNAs ($p={fmtp(h1['p_one_sided'])}$; median rank {med('cnn'):.0f}),"
     f" and so do CLIP-positive sets (H2, {h2['wins']}/13; partly circular because miRTarBase includes CLIP-derived entries). H3 is falsified: the HEK293 expression"
     f" top-200 control also ranks the own term well, and the CNN beats it on only {h3['wins']}/13 ({h3['losses']} losses, $p={fmtp(h3['p_one_sided'])}$).",
     r"Because the CNN's univariate CLIP signal was earlier shown to be mostly site-count confounding, and expression dominates the other audits, we added two"
     r" controls after H1 came back positive (post-hoc, not pre-registered). The CNN beats a site-count top-200 set drawn from the same universe"
     f" ({t1['wins']} wins, {t1['losses']} losses, 2 ties, $p={fmtp(t1['p_one_sided'])}$) and random sets matched to the CNN set's expression-decile histogram"
     f" ({t2['wins']} wins, {t2['losses']} losses, $p={fmtp(t2['p_one_sided'])}$; Figure~\\ref{{fig:enrichr}}). This is the first audit in the paper where CNN ranking"
     r" carries signal that survives both confounds.",
     r"Limits: miRTarBase 2017 overlaps the miRTarBase 10.0 labels used in the earlier benchmark, so this is a second view of the same kind of evidence rather than"
     r" an independent one; heavily studied miRNAs (miR-21, miR-155) reach rank 1 even from random sets, so the test carries information mainly for the other 11;"
     r" and the controls are post-hoc. We record it as a named, falsifiable positive: on a held-out panel with a pre-registered site-count and expression-matched"
     r" control, CNN top-200 non-conserved sets should keep ranking their own miRNA better than both.",
     r"\begin{figure}[h]\centering\includegraphics[width=0.95\linewidth]{figs/fig_enrichr.pdf}",
     r"\caption{Rank of each miRNA's own miRTarBase\_2017 term when its target set is submitted to Enrichr plotted as $\log_{10}(3241/\mathrm{rank})$, so taller bars mean better ranks."
     r" Site-count and expression-matched arms are post-hoc controls.}\label{fig:enrichr}\end{figure}"]
open("paper/enrichr_sec.tex", "w").write("\n".join(S) + "\n")
print("ok")
