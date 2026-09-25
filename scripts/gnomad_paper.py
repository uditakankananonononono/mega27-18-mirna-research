#!/usr/bin/env python3
"""Figure + LaTeX section for the gnomAD constraint audit (committed JSON only)."""
import json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def fmtp(p):
    if p >= 0.001: return f"{p:.3f}"
    m, e = f"{p:.1e}".split("e"); return f"{m}\\times10^{{{int(e)}}}"


J = json.load(open("results/gnomad_constraint.json")); P = J["per_mirna"]; ms = J["panel"]; G = J["gates"]; PH = J["post_hoc_not_preregistered"]
x = np.arange(len(ms)); w = 0.16
fig, ax = plt.subplots(figsize=(9.0, 3.2))
for i, (k, c, n) in enumerate([("cnn", "#2980b9", "CNN top-200"), ("sitecount", "#f39c12", "site-count top-200"), ("matched_mean", "#16a085", "expression-matched"),
                               ("uniform_mean", "#95a5a6", "uniform random (CLIP)"), ("ts_median", "#c0392b", "TargetScan conserved top-200")]):
    ax.bar(x + (i - 2) * w, [P[m][k] for m in ms], w, color=c, label=n)
ax.plot(x, [np.mean(P[m]["ts_rand"]) for m in ms], "k_", ms=12, label="random UTR genes")
ax.set_ylabel("median LOEUF (lower = constrained)"); ax.set_ylim(0.4, 1.2)
ax.set_xticks(x); ax.set_xticklabels([m.replace("hsa-", "") for m in ms], rotation=45, ha="right", fontsize=7); ax.legend(fontsize=6.3, frameon=False, ncol=3, loc="upper left")
ax.set_title("LoF constraint of target sets (gnomAD v2.1.1)", fontsize=9)
fig.tight_layout(); fig.savefig("paper/figs/fig_gnomad.pdf")
h1, h2, h3, g1 = J["H1_cnn_vs_uniform"], J["H2_cnn_vs_sitecount"], J["H3_cnn_vs_exprmatched"], G["G1_targetscan_vs_random"]
S = [r"\subsection{Loss-of-function constraint (gnomAD): conserved targets are constrained, CNN candidates are the least constrained set}",
     r"Conserved miRNA targets tend to be dosage-sensitive. We compared the median LOEUF (gnomAD v2.1.1 upper bound of observed/expected loss-of-function;"
     f" lower means more constrained; {J['n_genes_loeuf']:,} genes; script and results: \\texttt{{gnomad\\_constraint}}) of 200-gene sets, with the same"
     r" controls fixed in advance.",
     f"Gates passed: {100*G['G2_coverage']:.1f}\\% of CLIP-universe genes have a LOEUF, and TargetScan conserved top-200 sets are more constrained than random"
     f" UTR genes for {g1['wins']}/{g1['wins']+g1['losses']} miRNAs ($p={fmtp(g1['p_one_sided'])}$). All three pre-registered CNN hypotheses are falsified, and"
     f" in the opposite direction: the CNN top-200 has a higher median LOEUF than uniform random ({h1['losses']}/13), site count ({h2['losses']}/13) and"
     f" expression-matched sets ({h3['losses']}/13; Figure~\\ref{{fig:gnomad}}). Not pre-registered: the two-sided sign test of this reversal gives"
     f" $p={fmtp(PH['p_two_sided_uniform'])}$ against uniform and $p={fmtp(PH['p_two_sided_matched'])}$ against expression-matched sets.",
     f"Site-count sets are more constrained than uniform random for {sum(P[m]['sitecount'] < P[m]['uniform_mean'] for m in ms)}/13 miRNAs, as expected if they favour genes with long 3$'$ UTRs. The CNN's preference for genes with"
     r" few, strong sites may select short-UTR, less constrained genes. This post-hoc finding needs a held-out, UTR-length-matched replication before it can be"
     r" read as biology, and we make no claim beyond the description.",
     r"\begin{figure}[h]\centering\includegraphics[width=0.95\linewidth]{figs/fig_gnomad.pdf}",
     r"\caption{Median LOEUF of each 200-gene set. Black dashes: random UTR genes, the comparator for the TargetScan positive control.}\label{fig:gnomad}\end{figure}"]
open("paper/gnomad_sec.tex", "w").write("\n".join(S) + "\n"); print("ok")
