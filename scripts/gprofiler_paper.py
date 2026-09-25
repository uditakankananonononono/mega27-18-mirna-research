#!/usr/bin/env python3
"""Figure + LaTeX section for the g:Profiler audit (reads committed JSON only)."""
import json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

J = json.load(open("results/gprofiler_enrichment.json"))
P = J["per_mirna"]; ms = J["panel"]; G = J["gates"]; PH = J["posthoc_genome_background"]
lab = [m.replace("hsa-", "") for m in ms]
cnn = [P[m]["cnn"]["n_sig"] for m in ms]; clip = [P[m]["clip"]["n_sig"] for m in ms]
ex = [P[m]["expr"]["n_sig"] for m in ms]; rd = [P[m]["rand_mean_nsig"] for m in ms]
fig, ax = plt.subplots(figsize=(9.0, 3.3)); x = np.arange(len(ms)); w = 0.2
for i, (v, c, n) in enumerate([(cnn, "#2980b9", "CNN top-200"), (clip, "#8e44ad", "CLIP+ top-200"),
                               (rd, "#7f8c8d", "random 200 (mean of 3)"), (ex, "#e67e22", "HEK293 expression top-200")]):
    ax.bar(x + (i - 1.5) * w, v, w, color=c, label=n)
ax.set_xticks(x); ax.set_xticklabels(lab, rotation=45, ha="right", fontsize=7)
ax.set_ylabel("significant terms (g:SCS $<$ 0.05)"); ax.legend(fontsize=7, frameon=False)
ax.set_title("g:GOSt enrichment vs matched per-miRNA background (GO:BP, KEGG, Reactome)", fontsize=9)
fig.tight_layout(); fig.savefig("paper/figs/fig_gprofiler.pdf")


def fmtp(p):
    if p >= 0.001: return f"{p:.3f}"
    m, e = f"{p:.1e}".split("e"); return f"{m}\\times10^{{{int(e)}}}"


h1, h2, h3 = J["H1_cnn_vs_random"], J["H2_termJaccard_cnn_vs_random"], J["H3_cnn_vs_expression"]
n = len(ms)
zc = sum(v == 0 for v in cnn); zk = sum(v == 0 for v in clip)
S = []
S.append(r"\subsection{Functional coherence audit: g:Profiler finds no pathway signal in predicted or CLIP target sets (negative)}")
S.append(r"If the CNN ranks real targets, its top-ranked genes for a miRNA might share biology. We tested this with the g:Profiler g:GOSt service"
         f" (version \\texttt{{{J['gprofiler_version'].replace('_', chr(92)+'_')}}}; \\texttt{{scripts/gprofiler\\_enrichment.py}}; \\texttt{{results/gprofiler\\_enrichment.json}}; hermetic tests)"
         f" on {n} well-studied miRNAs. For each miRNA the background is the set of genes scored for that miRNA, and we query four gene sets of size $K={J['K']}$:"
         r" CNN top-$K$, CLIP-positive genes ranked by CNN score, the top-$K$ genes by HEK293 abundance (confound control) and "
         f"{J['R']} seeded random draws. Terms are called significant under g:SCS, which corrects the one-sided hypergeometric tail")
S.append(r"\begin{equation}p=\sum_{i\ge k}\frac{\binom{K_t}{i}\binom{N-K_t}{n-i}}{\binom{N}{n}}\end{equation}")
S.append(r"(query size $n$, background $N$, term size $K_t$, overlap $k$) for the dependence among GO terms. Gates, fixed before any group comparison, all passed:"
         f" a 20-gene canonical cell-cycle positive control recovers GO:0007049 (top hit KEGG cell cycle, $p={fmtp(G['G1_top'][0][2])}$);"
         f" the median random set has {G['G2_median_random_nsig']:.0f} significant terms; the median unmapped fraction is {100*G['G3_median_unmapped']:.2f}\\%.")
S.append(f"All three hypotheses are falsified. CNN top-$K$ sets have zero significant terms for {zc}/{n} miRNAs and CLIP-positive sets for {zk}/{n},"
         f" so neither beats random (H1: {h1['wins']} wins, {h1['losses']} losses, $p={fmtp(h1['p_one_sided'])}$; H2 term overlap with the CLIP set: {h2['wins']} wins, {h2['losses']} losses)."
         f" The expression control returns {min(ex)}--{max(ex)} terms per miRNA and beats the CNN on {h3['losses']}/{n} (H3, Figure~\\ref{{fig:gprofiler}}).")
S.append(f"A post-hoc arm that we did not pre-register swaps to the default whole-genome background: CNN top-$K$ still loses to the first random set"
         f" ({PH['cnn_vs_random']['wins']} wins, {PH['cnn_vs_random']['losses']} losses; median {PH['median_cnn_nsig']:.0f} vs {PH['median_rand_nsig']:.0f} terms),"
         r" consistent with (but not a test of) the CLIP-scored universe being biased toward expressed, well-annotated genes. The reading matches the expression-confound section:"
         r" at $K=200$, target lists from sequence or CLIP carry no detectable pathway coherence, while abundance alone carries a lot. Pathway enrichment of"
         r" predicted miRNA targets should not be read as evidence of targeting unless it beats an expression-matched control.")
S.append(r"\begin{figure}[h]\centering\includegraphics[width=0.95\linewidth]{figs/fig_gprofiler.pdf}")
S.append(r"\caption{Significant g:GOSt terms per miRNA for four gene sets of 200 against each miRNA's scored-gene background. Expression top-200 (orange) is the confound control.}\label{fig:gprofiler}\end{figure}")
open("paper/gprofiler_sec.tex", "w").write("\n".join(S) + "\n")
print("wrote paper/gprofiler_sec.tex and fig")
