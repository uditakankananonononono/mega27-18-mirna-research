#!/usr/bin/env python3
"""Figure + LaTeX for the ClinGen dosage-sensitivity test (committed JSON only)."""
import json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def fmtp(p):
    if p >= 0.001: return f"{p:.3f}"
    m, e = f"{p:.1e}".split("e"); return f"{m}\\times10^{{{int(e)}}}"


J = json.load(open("results/clingen_dosage.json")); G = J["gates"]
H = [J["C1_cnn_vs_uniform"], J["C2_cnn_vs_expr"], J["C3_cnn_vs_utrlen"]]
fig, ax = plt.subplots(figsize=(6.0, 2.8))
lab = ["TargetScan conserved\nvs random UTR", "CNN vs\nuniform", "CNN vs\nexpression-matched", "CNN vs\nUTR-length-matched"]
val = [G["G1_pooled_ts_vs_random"][0] / G["G1_pooled_ts_vs_random"][1]] + [h["pooled_ratio"] for h in H]
ax.bar(range(4), val, color=["#c0392b", "#2980b9", "#2980b9", "#2980b9"]); ax.axhline(1, c="k", lw=.8)
for i, v in enumerate(val): ax.text(i, v + .03, f"{v:.2f}", ha="center", fontsize=8)
ax.set_xticks(range(4)); ax.set_xticklabels(lab, fontsize=7); ax.set_ylabel("pooled HI-gene count ratio")
ax.set_title("ClinGen haploinsufficient genes (score 3) in 200-gene sets", fontsize=9)
fig.tight_layout(); fig.savefig("paper/figs/fig_clingen.pdf")
g1 = G["G1_ts_vs_random"]
row = lambda n, h: f"{n} & {h['wins']}/{h['losses']} & ${fmtp(h['p_one_sided'])}$ & {h['pooled_cnn']} / {h['pooled_comparator']:.1f} & {h['verdict']} \\\\"
S = [r"\subsection{Independent dosage label (ClinGen): the depletion beyond UTR length does not replicate}",
     f"gnomAD LOEUF is a population statistic, so we tested the named finding against an independent label: expert-curated ClinGen haploinsufficient genes"
     f" (score 3, {J['n_hi_genes']} genes), on all {len(J['panel'])} miRNAs (13 discovery + 30 held-out). The script was committed before it was run"
     r" (\texttt{clingen\_dosage}). Gates passed: the pooled CLIP universe holds"
     f" {G['G2_hi_in_universe']} HI genes, and TargetScan conserved top-200 sets hold more HI genes than random UTR genes ({g1['wins']}/{g1['losses']} miRNAs,"
     f" $p={fmtp(g1['p_one_sided'])}$; pooled {G['G1_pooled_ts_vs_random'][0]} vs {G['G1_pooled_ts_vs_random'][1]:.0f}).",
     r"\begin{table}[h]\centering\small\begin{tabular}{lcccc}\hline Pre-registered test & wins/losses & $p$ & pooled CNN / control & verdict\\\hline",
     row("C1 CNN $<$ uniform", H[0]), row("C2 CNN $<$ expression-matched", H[1]), row("C3 CNN $<$ UTR-length-matched", H[2]),
     r"\hline\end{tabular}\caption{ClinGen HI genes in CNN top-200 sets. Per-miRNA sign tests (ties excluded) and pooled counts over 43 miRNAs.}\label{tab:clingen}\end{table}",
     f"The CNN sets hold fewer HI genes than expression-matched sets (C2 confirmed; pooled ratio {H[1]['pooled_ratio']:.2f}), but the depletion against"
     f" uniform sets is not significant (C1) and it vanishes against UTR-length-matched sets (pooled ratio {H[2]['pooled_ratio']:.2f}; C3 falsified;"
     r" Figure~\ref{fig:clingen}). On this independent label, UTR length explains the depletion. The gnomAD finding therefore survives only in its weaker form"
     r" (CNN candidates are less constrained than expression-matched genes); its ``beyond UTR length'' clause is supported by LOEUF but not by ClinGen."
     r" HI genes are sparse (about 4 per set), so C3 has little power; we report the disagreement rather than choosing between the labels.",
     r"\begin{figure}[h]\centering\includegraphics[width=0.62\linewidth]{figs/fig_clingen.pdf}",
     r"\caption{Pooled HI-gene count in each set divided by its control. Red: positive control.}\label{fig:clingen}\end{figure}"]
open("paper/clingen_sec.tex", "w").write("\n".join(S) + "\n"); print("ok")
