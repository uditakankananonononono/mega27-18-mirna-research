#!/usr/bin/env python3
"""LaTeX section for the DepMap common-essential test (committed JSON only)."""
import json


def fmtp(p):
    if p >= 0.001: return f"{p:.3f}"
    m, e = f"{p:.1e}".split("e"); return f"{m}\\times10^{{{int(e)}}}"


J = json.load(open("results/depmap_essential.json")); G = J["gates"]; g1 = G["G1_ts_vs_random"]
H = [("E1 fewer than uniform", J["E1_vs_uniform"]), ("E2 fewer than expression-matched", J["E2_vs_expr"]), ("E3 fewer than UTR-length-matched", J["E3_vs_utrlen"])]
row = lambda n, h: f"{n} & {h['wins']}/{h['losses']} & ${fmtp(h['p_one_sided'])}$ & {h['pooled_cnn']} / {h['pooled_comparator']:.0f} ({h['pooled_ratio']:.2f}) & {h['verdict']} \\\\"
S = [r"\subsection{Cell-essential genes (DepMap): no depletion, and the label is uninformative for conserved targets}",
     f"To ask whether the dosage depletion extends to cell fitness, we used the DepMap 24Q4 inferred common-essential genes ({J['n_essential']:,} genes;"
     r" \texttt{depmap\_essential}, committed before running; same 43 miRNAs and controls). Prevalence passed"
     f" ({100*G['G2_prevalence']:.1f}\\%). The positive-control gate, which only asked for a difference in either direction, failed: conserved top-200 sets"
     f" and random UTR genes split {g1['wins']}/{g1['losses']} (two-sided $p={fmtp(g1['p_two_sided'])}$).",
     r"\begin{table}[h]\centering\small\begin{tabular}{lcccc}\hline Pre-registered test & wins/losses & $p$ & pooled CNN / control (ratio) & verdict\\\hline",
     *[row(n, h) for n, h in H], r"\hline\end{tabular}\caption{DepMap common-essential genes in CNN top-200 sets, 43 miRNAs.}\label{tab:depmap}\end{table}",
     r"All three hypotheses are falsified. CNN sets hold slightly more essential genes than every control, not fewer"
     f" (pooled ratios {H[0][1]['pooled_ratio']:.2f}--{H[1][1]['pooled_ratio']:.2f}). Whatever drives the depletion of dosage-sensitive and disease genes,"
     r" it does not extend to genes needed for cell growth in culture. That fits a dosage reading better than a general avoidance of important genes,"
     r" but with the positive-control gate failed we do not lean on it."]
open("paper/depmap_sec.tex", "w").write("\n".join(S) + "\n"); print("ok")
