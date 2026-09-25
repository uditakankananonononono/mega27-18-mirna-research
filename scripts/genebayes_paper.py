#!/usr/bin/env python3
"""LaTeX section for the GeneBayes s_het test (committed JSON only)."""
import json


def fmtp(p):
    if p >= 0.001: return f"{p:.3f}"
    m, e = f"{p:.1e}".split("e"); return f"{m}\\times10^{{{int(e)}}}"


J = json.load(open("results/genebayes_shet.json")); G = J["gates"]; g1 = G["G1_ts_vs_random"]
H = [("H1 lower than uniform", "H1_vs_uniform"), ("H2 lower than expression-matched", "H2_vs_expr"), ("H3 lower than UTR-length-matched", "H3_vs_utrlen"),
     ("H4 lower than publication-matched", "H4_vs_pubmatched")]
row = lambda n, k: f"{n} & {J[k]['wins']}/{J[k]['losses']} & ${fmtp(J[k]['p_one_sided'])}$ & {J[k]['median_diff']:+.3f} & {J[k]['verdict']} \\\\"
S = [r"\subsection{GeneBayes $s_\mathrm{het}$: population-based labels agree with each other, curated and model-based labels do not}",
     f"GeneBayes $s_\\mathrm{{het}}$ (Zenodo 10403680; {J['n_scored_genes']:,} genes) estimates selection against heterozygous loss of function from gnomAD"
     r" LoF counts plus a gene-feature prior, so it is a second population-based label and not independent of LOEUF (\texttt{genebayes\_shet},"
     f" committed before running; statistic median $\\log_{{10}} s_\\mathrm{{het}}$). Gates passed: {100*G['G2_coverage']:.1f}\\% coverage, and the positive"
     f" control scores higher than random ({g1['wins']}/{g1['losses']}, $p={fmtp(g1['p_one_sided'])}$).",
     r"\begin{table}[h]\centering\small\begin{tabular}{lcccc}\hline Pre-registered test & wins/losses & $p$ & median diff. & verdict\\\hline",
     *[row(n, k) for n, k in H], r"\hline\end{tabular}\caption{GeneBayes $s_\mathrm{het}$ of CNN top-200 sets, 43 miRNAs (differences in $\log_{10}$ units).}\label{tab:genebayes}\end{table}",
     r"All four hypotheses are confirmed, including H3 against UTR-length-matched sets. The split is now clean: the three labels built from"
     r" gnomAD loss-of-function counts (LOEUF, rCNV, $s_\mathrm{het}$) support depletion beyond UTR length, while DECIPHER, ClinGen and HPO do not."
     r" A plausible explanation, not tested here, is that labels built from LoF counts are less informative for genes with few expected LoF variants"
     r" (short coding sequences), and the CNN may favour such genes. Because the split follows label construction, the clause about UTR length stays"
     r" withdrawn; the surviving claim (depletion relative to expression-matched and publication-matched genes) is confirmed here too."]
open("paper/genebayes_sec.tex", "w").write("\n".join(S) + "\n"); print("ok")
