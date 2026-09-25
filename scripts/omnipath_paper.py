#!/usr/bin/env python3
"""Figure + LaTeX section for the OmniPath independent-label audit (committed JSON only)."""
import json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def fmtp(p):
    if p >= 0.001: return f"{p:.3f}"
    m, e = f"{p:.1e}".split("e"); return f"{m}\\times10^{{{int(e)}}}"


J = json.load(open("results/omnipath_independent.json")); P = J["per_mirna"]; ms = list(P); G = J["gates"]; M = J["median_auroc"]
x = np.arange(len(ms)); w = 0.27
fig, ax = plt.subplots(figsize=(9.0, 3.2))
for i, (k, c, n) in enumerate([("auroc_cnn", "#2980b9", "CNN"), ("auroc_sitecount", "#f39c12", "site count"), ("auroc_expr", "#e67e22", "HEK293 expression")]):
    ax.bar(x + (i - 1) * w, [P[m][k] for m in ms], w, color=c, label=n)
ax.axhline(0.5, ls=":", c="k", lw=0.8); ax.set_ylim(0.2, 0.9)
ax.set_xticks(x); ax.set_xticklabels([f"{m.replace('hsa-', '')}\n(n+={P[m]['n_pos']})" for m in ms], rotation=45, ha="right", fontsize=6.5)
ax.set_ylabel("AUROC"); ax.legend(fontsize=7, frameon=False, ncol=3); ax.set_title("Non-miRTarBase curated labels (OmniPath) inside the CLIP universe", fontsize=9)
fig.tight_layout(); fig.savefig("paper/figs/fig_omnipath.pdf")
h1, h2, h3 = J["H1_cnn_above_chance"], J["H2_cnn_vs_sitecount"], J["H3_cnn_vs_expr"]
S = [r"\subsection{Independent curated labels (OmniPath, miRTarBase removed): CNN above chance, but not above site count or expression}",
     r"The miRTarBase and Enrichr analyses share one curation source. OmniPath aggregates several (script and results:"
     r" \texttt{omnipath\_independent}; hermetic tests), so we kept only interactions whose sources exclude miRTarBase"
     f" ({J['n_rows_non_mirtarbase']:,} of {J['n_rows_all']:,}; ncRDeathDB, miRecords, miR2Disease, miRDeathDB, SIGNOR) and scored them inside each miRNA's CLIP universe"
     r" of non-conserved seed-match genes. The site-count, expression and label-permutation controls were fixed before scoring.",
     f"Gates passed: {G['G1_n_eligible']} miRNAs have at least five positives (bar 10), and the median permutation-mean CNN AUROC is {G['G2_median_perm_mean']:.3f}."
     f" The CNN is above chance (H1: {h1['wins']}/{h1['losses']}, $p={fmtp(h1['p_one_sided'])}$, median AUROC {M['cnn']:.3f}), but it does not beat site count"
     f" (H2: {h2['wins']}/{h2['losses']}, $p={fmtp(h2['p_one_sided'])}$, median {M['sitecount']:.3f}) and loses to HEK293 expression (H3: {h3['wins']}/{h3['losses']},"
     f" median {M['expr']:.3f}; Figure~\\ref{{fig:omnipath}}). This matches the Enrichr replication on a label set that shares no miRTarBase rows: curated"
     r" low-throughput targets favour expressed genes, and whatever the CNN adds over counting seed sites is not detectable with 5--29 positives per miRNA.",
     r"\begin{figure}[h]\centering\includegraphics[width=0.95\linewidth]{figs/fig_omnipath.pdf}",
     r"\caption{Per-miRNA AUROC on OmniPath interactions with miRTarBase-sourced rows removed, for the CNN, site count and HEK293 expression.}\label{fig:omnipath}\end{figure}"]
open("paper/omnipath_sec.tex", "w").write("\n".join(S) + "\n"); print("ok")
