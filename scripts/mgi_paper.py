#!/usr/bin/env python3
"""LaTeX section for the MGI mouse-mutant phenotype test (committed JSON only)."""
import json


def fmtp(p):
    if p >= 0.001: return f"{p:.3f}"
    m, e = f"{p:.1e}".split("e"); return f"{m}\\times10^{{{int(e)}}}"


J = json.load(open("results/mgi_mouse_ko.json")); G = J["gates"]; g1 = G["G1_ts_vs_random"]
H = [("M1 mortality, vs expression-matched", "M1_mort_vs_expr"), ("M2 mortality, vs publication-matched", "M2_mort_vs_pub"),
     ("M3 mortality, vs UTR-length-matched", "M3_mort_vs_utrlen"), ("M4 embryo phenotype, vs publication-matched", "M4_emb_vs_pub")]
row = lambda n, k: f"{n} & {J[k]['wins']}/{J[k]['losses']} & ${fmtp(J[k]['p_one_sided'])}$ & {J[k]['median_diff']:+.3f} & {J[k]['verdict']} \\\\"
S = [r"\subsection{Mouse mutant phenotypes (MGI): weak, mixed support}",
     r"Mouse knockouts give an experimental label that owes nothing to human population data or clinical curation. From the MGI"
     r" \texttt{HMD\_HumanPhenotype} report~\cite{mgi} we marked human genes whose mouse orthologue mutants show mortality/aging (MP:0010768)"
     f" or embryo (MP:0005380) phenotypes ({J['n_mortality']} and {J['n_embryo']} of {J['n_phenotyped']} phenotyped genes); genes never"
     r" phenotyped were dropped (\texttt{mgi\_mouse\_ko}, committed before running; 43 miRNAs). Both gates passed: conserved TargetScan targets are"
     f" more often lethal than random genes ({g1['wins']}/{g1['losses']}, $p={fmtp(g1['p_one_sided'])}$), and"
     f" {100*G['G2_coverage']:.1f}\\% of universe symbols have an MGI row (bar 80\\%).",
     r"\begin{table}[h]\centering\small\begin{tabular}{lcccc}\hline Pre-registered test (CNN lower) & wins/losses & $p$ & median diff. & verdict\\\hline",
     *[row(n, k) for n, k in H], r"\hline\end{tabular}\caption{Mouse mutant phenotype fractions in CNN top-200 sets, 43 miRNAs.}\label{tab:mgi}\end{table}",
     r"All four comparisons point the same way (CNN sets less often lethal), but only M1 and M4 reach $p<0.05$; mortality is not significantly depleted"
     r" against publication- or UTR-length-matched sets. The plan set no multiple-testing correction; after the fact we note that with Bonferroni over"
     r" the four tests ($\alpha=0.0125$) none would pass. We read this as weak support from an independent experimental label, not as a replication."]
open("paper/mgi_sec.tex", "w").write("\n".join(S) + "\n"); print("ok")
