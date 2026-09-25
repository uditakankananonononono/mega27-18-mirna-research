#!/usr/bin/env python3
"""LaTeX section for the NCBI Gene study-bias test (committed JSONs only)."""
import json


def fmtp(p):
    if p >= 0.001: return f"{p:.3f}"
    m, e = f"{p:.1e}".split("e"); return f"{m}\\times10^{{{int(e)}}}"


J = json.load(open("results/ncbi_pubmed_bias.json")); Q = json.load(open("results/ncbi_pubmed_posthoc.json")); G = J["gates"]; g1 = G["G1_ts_vs_random"]
names = [("P1 CNN less studied than uniform", "P1_less_studied_vs_uniform"), ("P2 less studied than expression-matched", "P2_vs_expr"),
         ("P3 less studied than UTR-length-matched", "P3_vs_utrlen"), ("P4 LOEUF higher than publication-matched", "P4_loeuf_vs_pubmatched"),
         ("P5 fewer HPO genes than publication-matched", "P5_hpo_vs_pubmatched")]
row = lambda n, k: f"{n} & {J[k]['wins']}/{J[k]['losses']} & ${fmtp(J[k]['p_one_sided'])}$ & {J[k]['verdict']} & {Q[k]['wins']}/{Q[k]['losses']} (${fmtp(Q[k]['p_one_sided'])}$) \\\\"
u = Q["unmapped_cnn_higher"]
S = [r"\subsection{Study bias (NCBI Gene, PubMed counts): it explains the curated-label depletion but not the LOEUF depletion}",
     r"Curated labels (ClinGen, HPO, ClinVar) favour well-studied genes. We counted PubMed records per human gene from NCBI Gene"
     r" (gene\_info and gene2pubmed; \texttt{ncbi\_pubmed\_bias}, committed before running) and added a publication-matched control"
     f" (deciles of $\\log_{{10}}(1+n)$). The positive control passed (conserved targets are better studied, {g1['wins']}/{g1['losses']}, $p={fmtp(g1['p_one_sided'])}$),"
     f" but the mapping gate failed: {100*G['G2_map_rate']:.1f}\\% of universe symbols matched an NCBI symbol against a 95\\% bar. The pre-registered script"
     r" scores unmapped symbols as zero publications, so we add a sensitivity analysis, labelled not pre-registered, that drops them"
     f" (\\texttt{{ncbi\\_pubmed\\_posthoc}}). Unmapped symbols are more common in CNN sets ({u['wins']}/{u['losses']} miRNAs; median {100*u['median_cnn']:.1f}\\%"
     f" vs {100*u['median_universe']:.1f}\\%).",
     r"\begin{table}[h]\centering\footnotesize\begin{tabular}{lcccc}\hline Test & wins/losses & $p$ & verdict (G2 failed) & unmapped dropped$^\dagger$\\\hline",
     *[row(n, k) for n, k in names], r"\hline\end{tabular}\caption{Study bias of CNN top-200 sets, 43 miRNAs. $^\dagger$Not pre-registered.}\label{tab:pubmed}\end{table}",
     r"Both analyses agree. CNN candidates are less studied than every control. Against publication-matched sets the LOEUF depletion stays at 43/0,"
     r" so study bias does not explain it, which is expected because LOEUF comes from population sequencing, not from the literature. The HPO depletion"
     r" disappears (P5 falsified, and reversed in the post-hoc split). This explains why the curated labels (ClinGen, HPO) disagreed with LOEUF: their"
     r" depletion signal is largely study bias. Because gate G2 failed, the verdicts carry that caveat; the post-hoc split is not a substitute for"
     r" a pre-registered test."]
open("paper/pubmed_sec.tex", "w").write("\n".join(S) + "\n"); print("ok")
