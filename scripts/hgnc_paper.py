#!/usr/bin/env python3
"""LaTeX section for the HGNC gene-family test (committed JSON only)."""
import json


def fmtp(p):
    if p >= 0.001: return f"{p:.3f}"
    m, e = f"{p:.1e}".split("e"); return f"{m}\\times10^{{{int(e)}}}"


J = json.load(open("results/hgnc_families.json")); G = J["gates"]; A = json.load(open("results/hgnc_unmatched_audit.json"))
H = [("F1 C2H2-ZNF fraction above expression-matched", "F1_znf_vs_expr"), ("F2 family size above expression-matched", "F2_family_size_vs_expr"),
     ("F3 LOEUF above expression-matched, ZNFs removed", "F3_loeuf_noznf_vs_expr")]
row = lambda n, k: f"{n} & {J[k]['wins']}/{J[k]['losses']} & ${fmtp(J[k]['p_one_sided'])}$ & {(J[k]['median_diff'] if abs(J[k]['median_diff'])>=5e-4 else 0.0):+.3f} & {J[k]['verdict']} \\\\"
S = [r"\subsection{Gene families (HGNC): the depletion is not a zinc-finger artefact; coverage gate failed}",
     r"Large, dosage-tolerant families such as the C2H2 zinc fingers (many with long 3$'$UTRs) could make a candidate set look unconstrained."
     f" We used HGNC gene groups~\\cite{{hgnc}} ({J['n_protein_coding']} approved protein-coding symbols, {J['n_znf']} C2H2-ZNF;"
     r" \texttt{hgnc\_families}, committed before running; 43 miRNAs; family size $=\log_2(1+\text{largest group size})$)."
     f" The parse gate passed (C2H2-ZNF genes are {G['G2_chr19_znf_ratio']:.1f}$\\times$ more frequent on chr19, bar 3$\\times$), but the coverage"
     f" gate failed: {100*G['G1_coverage']:.1f}\\% of universe symbols matched current HGNC protein-coding symbols (bar 90\\%). A post-hoc audit (not pre-registered) found that {A['previous_symbol']} of the {A['unmatched']} unmatched symbols are previous HGNC symbols and {A['approved_non_protein_coding']} are approved non-protein-coding loci.",
     r"\begin{table}[h]\centering\small\begin{tabular}{lcccc}\hline Pre-registered test (G1 failed) & wins/losses & $p$ & median diff. & verdict\\\hline",
     *[row(n, k) for n, k in H], r"\hline\end{tabular}\caption{HGNC family composition of CNN top-200 sets, 43 miRNAs.}\label{tab:hgnc}\end{table}",
     r"CNN candidates are not enriched for zinc fingers or large families (F1, F2 falsified), and with every C2H2-ZNF gene removed the LOEUF"
     f" depletion against expression-matched sets stays at {J['F3_loeuf_noznf_vs_expr']['wins']}/{J['F3_loeuf_noznf_vs_expr']['losses']}."
     r" Because the coverage gate failed, this rules out the family explanation only for the matched genes."]
open("paper/hgnc_sec.tex", "w").write("\n".join(S) + "\n"); print("ok")
