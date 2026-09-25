#!/usr/bin/env python3
"""Figure + LaTeX section for the ENCODE miRNA-abundance audit (committed JSON only)."""
import json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def fmtp(p):
    if p >= 0.001: return f"{p:.3f}"
    m, e = f"{p:.1e}".split("e"); return f"{m}\\times10^{{{int(e)}}}"


J = json.load(open("results/encode_mirna_abundance.json"))
T = [r for r in J["per_mirna"] if r["abundance"] is not None]
a = np.array([r["abundance"] for r in T]); c = np.array([r["cnn_auroc"] for r in T]); d = np.array([r["delta"] for r in T])
H1, H2, H3, G, PH = J["H1_abundance_vs_cnn_auroc"], J["H2_abundance_vs_delta"], J["H3_partial_given_log_npos"], J["gates"], J["posthoc_G2b_tissue_controls"]
fig, (x1, x2) = plt.subplots(1, 2, figsize=(9.0, 3.3))
x1.scatter(a, c, s=12, alpha=0.7, c="#2980b9"); x1.set_xlabel("median log$_2$(CPM+1), 12 ENCODE lines"); x1.set_ylabel("per-miRNA CNN AUROC (CLIP)")
x1.set_title(f"H1: Spearman $\\rho$={H1['rho']:.3f}, p={H1['p_one_sided']:.3f}", fontsize=9)
x2.scatter(a, d, s=12, alpha=0.7, c="#8e44ad"); x2.axhline(0, ls=":", c="k", lw=0.8)
x2.set_xlabel("median log$_2$(CPM+1), 12 ENCODE lines"); x2.set_ylabel("CNN gain over covariates ($\\Delta$AUROC)")
x2.set_title(f"H2: Spearman $\\rho$={H2['rho']:.3f}, p={H2['p_one_sided']:.3f}", fontsize=9)
fig.tight_layout(); fig.savefig("paper/figs/fig_encode_mirna.pdf")
S = [r"\subsection{miRNA abundance does not predict which miRNAs the CNN ranks well (ENCODE; negative)}",
     r"CLIP labels need AGO loaded with the miRNA, so abundant miRNAs could give cleaner labels and higher per-miRNA AUROC. We measured abundance from"
     f" {J['n_files']} ENCODE GRCh38 microRNA-seq quantification files covering {J['n_lines']} cell lines (accessions in \\texttt{{results/encode\\_mirna\\_abundance.json}};"
     r" Ensembl REST \texttt{lookup/id} for gene symbols; \texttt{scripts/encode\_mirna\_abundance.py}; hermetic tests). Per file we take $a_g=\log_2(\mathrm{CPM}_g+1)$,"
     r" sum precursor loci per mature miRNA (precursor counts cannot separate the 5p and 3p arms), average replicates and take the median across lines.",
     f"Gates: every file has at least 1{{,}}500 gene rows; miR-21 sits in the panel's top decile ({G['G1_miR21_abundance']:.1f} vs threshold {G['G1_top_decile_threshold']:.1f});"
     f" {J['n_mapped']}/{J['n_panel']} panel miRNAs map to a precursor. The pre-registered tissue gate failed: we expected miR-122 to peak in HepG2, but it is zero in"
     r" all 12 lines. We record that failure and do not retune it. A replacement control chosen after the failure (post-hoc, not pre-registered) passes:"
     f" MIR142's top three lines are the three blood lines ({', '.join(PH['MIR142_top3'])}), and MIR302A, MIR302B and MIR367 each peak in H1 embryonic stem cells.",
     f"All three hypotheses are falsified on $n={H1['n']}$ miRNAs: abundance versus CNN AUROC $\\rho={H1['rho']:.3f}$ ($p={fmtp(H1['p_one_sided'])}$), versus the CNN's gain"
     f" over site-count covariates $\\rho={H2['rho']:.3f}$ ($p={fmtp(H2['p_one_sided'])}$), and the partial correlation given $\\log n_{{pos}}$ $\\rho={H3['rho']:.3f}$"
     f" ($p={fmtp(H3['p_one_sided'])}$; Figure~\\ref{{fig:encode}}). Abundance is also unrelated to the covariate-only AUROC ($\\rho={J['descriptive']['abundance_vs_base_auroc']['rho']:.3f}$)."
     r" A caveat: ENCORI CLIP pools many cell types, and none of the 12 ENCODE lines is HEK293, so a cell-matched abundance could still matter.",
     r"\begin{figure}[h]\centering\includegraphics[width=0.95\linewidth]{figs/fig_encode_mirna.pdf}",
     r"\caption{Per-miRNA ENCODE abundance (median over 12 cell lines) versus CNN AUROC on CLIP labels (left) and versus the CNN's gain over site-count covariates (right).}\label{fig:encode}\end{figure}"]
open("paper/encode_mirna_sec.tex", "w").write("\n".join(S) + "\n")
print("ok")
