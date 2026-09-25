#!/usr/bin/env python3
"""Figure + LaTeX section for the UniProt TF-enrichment audit (committed JSON only)."""
import json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def fmtp(p):
    if p >= 0.001: return f"{p:.3f}"
    m, e = f"{p:.1e}".split("e"); return f"{m}\\times10^{{{int(e)}}}"


J = json.load(open("results/uniprot_tf_enrichment.json")); REL = J["uniprot_release"].replace("_", r"\_"); P = J["per_mirna"]; ms = J["panel"]; G = J["gates"]
d_cnn = [100 * (P[m]["cnn"] - P[m]["uniform_mean"]) for m in ms]; d_sc = [100 * (P[m]["sitecount"] - P[m]["uniform_mean"]) for m in ms]
d_em = [100 * (P[m]["matched_mean"] - P[m]["uniform_mean"]) for m in ms]; d_ts = [100 * (P[m]["ts_frac"] - np.mean(P[m]["ts_rand"])) for m in ms]
x = np.arange(len(ms)); w = 0.2
fig, ax = plt.subplots(figsize=(9.0, 3.2))
for i, (v, c, n) in enumerate([(d_cnn, "#2980b9", "CNN top-200"), (d_sc, "#f39c12", "site-count top-200"), (d_em, "#16a085", "expression-matched random"),
                               (d_ts, "#c0392b", "TargetScan conserved top-200 (positive control)")]):
    ax.bar(x + (i - 1.5) * w, v, w, color=c, label=n)
ax.axhline(0, c="k", lw=0.8); ax.set_ylabel("TF fraction minus random (pp)")
ax.set_xticks(x); ax.set_xticklabels([m.replace("hsa-", "") for m in ms], rotation=45, ha="right", fontsize=7); ax.legend(fontsize=6.5, frameon=False, ncol=2)
ax.set_title(f"Transcription regulators (UniProt KW-0805, release {J['uniprot_release']}) among target sets", fontsize=9)
fig.tight_layout(); fig.savefig("paper/figs/fig_uniprot.pdf")
h1, h2, h3, g1 = J["H1_cnn_vs_uniform"], J["H2_cnn_vs_sitecount"], J["H3_cnn_vs_exprmatched"], G["G1_targetscan_vs_random"]
S = [r"\subsection{Transcription-regulator enrichment (UniProt): present in conserved targets, absent in CNN non-conserved targets}",
     r"Conserved miRNA targets are known to be enriched for transcription regulators. We test whether CNN non-conserved targets share this property"
     f" using UniProtKB release {REL} (reviewed human entries with keyword KW-0805, {J['n_tf']:,} of {J['n_reviewed']:,}; script and results:"
     r" \texttt{uniprot\_tf\_enrichment}; hermetic tests). The statistic is the fraction of a 200-gene set carrying KW-0805; site-count and"
     r" expression-matched controls were fixed in advance.",
     f"Gates passed. The positive control recovers the classic result: TargetScan conserved top-200 sets are richer in transcription regulators than random"
     f" UTR genes for {g1['wins']}/{G['G1_n_eligible']} miRNAs ($p={fmtp(g1['p_one_sided'])}$), and KW-0805 prevalence is {100*G['G2_kw0805_prevalence']:.1f}\\%."
     f" CNN top-200 sets show no enrichment: versus uniform random {h1['wins']}/{h1['losses']} ($p={fmtp(h1['p_one_sided'])}$), versus site count {h2['wins']}/{h2['losses']}"
     f" ($p={fmtp(h2['p_one_sided'])}$), versus expression-matched {h3['wins']}/{h3['losses']} ($p={fmtp(h3['p_one_sided'])}$); all three are falsified (Figure~\\ref{{fig:uniprot}}).",
     r"Together with the GTEx avoidance result this is the second functional signature that conserved targets carry and CNN non-conserved candidates lack."
     r" The same caveat applies: conservation and scorer are confounded, and the universes differ between arms.",
     r"\begin{figure}[h]\centering\includegraphics[width=0.95\linewidth]{figs/fig_uniprot.pdf}",
     r"\caption{Transcription-regulator fraction of each 200-gene set minus its random baseline (percentage points). Red: TargetScan conserved targets"
     r" against random UTR genes (positive control).}\label{fig:uniprot}\end{figure}"]
open("paper/uniprot_sec.tex", "w").write("\n".join(S) + "\n"); print("ok")
