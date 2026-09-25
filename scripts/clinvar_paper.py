#!/usr/bin/env python3
"""LaTeX section for the ClinVar PLP-gene test (committed JSON only)."""
import json


def fmtp(p):
    if p >= 0.001: return f"{p:.3f}"
    m, e = f"{p:.1e}".split("e"); return f"{m}\\times10^{{{int(e)}}}"


J = json.load(open("results/clinvar_pathogenic.json")); G = J["gates"]; g1 = G["G1_ts_vs_random"]
date = J["clinvar_date"].split("dated ")[-1]
H = [("V1 fewer than uniform", J["V1_vs_uniform"]), ("V2 fewer than expression-matched", J["V2_vs_expr"]),
     ("V3 fewer than UTR-length-matched", J["V3_vs_utrlen"]), ("S3 ($\\geq$1 PLP) vs UTR-length-matched", J["S3_any_vs_utrlen"])]
row = lambda n, h: f"{n} & {h['wins']}/{h['losses']} & ${fmtp(h['p_one_sided'])}$ & {h['pooled_ratio']:.2f} & {h['verdict']} \\\\"
S = [r"\subsection{ClinVar pathogenic-variant genes: the positive-control gate failed}",
     f"The last label was ClinVar (gene summary dated {date}): genes with at least five Pathogenic/Likely pathogenic alleles ({J['n_plp5']:,} genes;"
     r" \texttt{clinvar\_pathogenic}, committed before running). We noted in advance that allele counts grow with coding length and testing effort."
     f" Prevalence passed ({100*G['G2_prevalence']:.1f}\\%), but the positive-control gate failed: TargetScan conserved top-200 sets do not hold"
     f" significantly more PLP genes than random UTR genes ({g1['wins']}/{g1['losses']}, $p={fmtp(g1['p_one_sided'])}$). Because the label does not"
     r" recover the known conserved-target signal, we treat the CNN results below as uninterpretable and do not count them for or against the finding.",
     r"\begin{table}[h]\centering\small\begin{tabular}{lcccc}\hline Pre-registered test (gate G1 failed) & wins/losses & $p$ & pooled ratio & verdict\\\hline",
     *[row(n, h) for n, h in H], r"\hline\end{tabular}\caption{ClinVar PLP genes in CNN top-200 sets, 43 miRNAs. Reported for completeness only.}\label{tab:clinvar}\end{table}",
     r"For completeness, the pattern matches ClinGen and HPO (depletion against uniform and expression-matched sets, none against UTR-length-matched sets)."]
open("paper/clinvar_sec.tex", "w").write("\n".join(S) + "\n"); print("ok")
