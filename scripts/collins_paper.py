#!/usr/bin/env python3
"""Figure + LaTeX for the Collins 2022 dosage test (committed JSON only)."""
import json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def fmtp(p):
    if p >= 0.001: return f"{p:.3f}"
    m, e = f"{p:.1e}".split("e"); return f"{m}\\times10^{{{int(e)}}}"


J = json.load(open("results/collins_dosage.json")); P = J["per_mirna"]; G = J["gates"]
arms = [("uniform_mean", "uniform"), ("expr_mean", "expression-\nmatched"), ("len_mean", "UTR-length-\nmatched")]
fig, ax = plt.subplots(figsize=(6.4, 2.9))
data = [[v["cnn"] - v[k] for v in P.values()] for k, _ in arms] + [[v["cnn_tri"] - v["len_tri_mean"] for v in P.values()]]
ax.boxplot(data, tick_labels=[n for _, n in arms] + ["pTriplo vs\nUTR-length"])
for i, x in enumerate(data): ax.scatter(np.full(len(x), i + 1) + np.random.default_rng(i).uniform(-.12, .12, len(x)), x, s=7, c="#2980b9", zorder=3)
ax.axhline(0, c="k", lw=.8); ax.set_ylabel("CNN minus control median score"); ax.set_title("Collins 2022 dosage scores, 43 miRNAs (pHaplo unless noted)", fontsize=9)
fig.tight_layout(); fig.savefig("paper/figs/fig_collins.pdf")
g1 = G["G1_ts_vs_random"]; H = {k: J[k] for k in ["D1_phaplo_vs_uniform", "D2_phaplo_vs_expr", "D3_phaplo_vs_utrlen", "T3_ptriplo_vs_utrlen"]}
row = lambda n, h: f"{n} & {h['wins']}/{h['losses']} & ${fmtp(h['p_one_sided'])}$ & {h['median_diff']:+.3f} & {h['verdict']} \\\\"
S = [r"\subsection{Third dosage label (Collins et al.\ 2022 rCNV scores): the depletion survives UTR-length matching}",
     f"We then used the pHaplo and pTriplo scores of Collins et al.\\ (Zenodo 6347673; {J['n_scored_genes']:,} genes), which are estimated from rare CNVs"
     r" in about 950,000 people. We noted in advance that the model uses gene features that include gnomAD constraint, so this label is only partly"
     r" independent of LOEUF. Same 43 miRNAs and controls; the script was committed before it was run (\texttt{collins\_dosage}). Gates passed:"
     f" {100*G['G2_coverage']:.1f}\\% coverage, and TargetScan conserved top-200 sets have higher pHaplo than random UTR genes ({g1['wins']}/{g1['losses']},"
     f" $p={fmtp(g1['p_one_sided'])}$).",
     r"\begin{table}[h]\centering\small\begin{tabular}{lcccc}\hline Pre-registered test (CNN lower) & wins/losses & $p$ & median diff. & verdict\\\hline",
     row("D1 pHaplo vs uniform", H["D1_phaplo_vs_uniform"]), row("D2 pHaplo vs expression-matched", H["D2_phaplo_vs_expr"]),
     row("D3 pHaplo vs UTR-length-matched", H["D3_phaplo_vs_utrlen"]), row("T3 pTriplo vs UTR-length-matched", H["T3_ptriplo_vs_utrlen"]),
     r"\hline\end{tabular}\caption{Collins 2022 dosage scores of CNN top-200 sets, 43 miRNAs.}\label{tab:collins}\end{table}",
     r"All four hypotheses are confirmed, including against UTR-length-matched sets for both haploinsufficiency and triplosensitivity (Figure~\ref{fig:collins})."
     r" Across the three labels the picture is: LOEUF and the rCNV scores support depletion beyond UTR length; the sparse ClinGen curation does not."
     r" Because the rCNV model shares features with LOEUF, it is not a clean third vote, and the finding stays at the strength of its weakest independent test."
     r" The remaining falsifier is still the random-weight CNN.",
     r"\begin{figure}[h]\centering\includegraphics[width=0.66\linewidth]{figs/fig_collins.pdf}",
     r"\caption{Per-miRNA difference between CNN top-200 and control median dosage score (negative = CNN less dosage-sensitive).}\label{fig:collins}\end{figure}"]
open("paper/collins_sec.tex", "w").write("\n".join(S) + "\n"); print("ok")
