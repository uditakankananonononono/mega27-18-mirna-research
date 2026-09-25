#!/usr/bin/env python3
"""Figure + LaTeX section for the GTEx target-avoidance audit (committed JSON only)."""
import json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def fmtp(p):
    if p >= 0.001: return f"{p:.3f}"
    m, e = f"{p:.1e}".split("e"); return f"{m}\\times10^{{{int(e)}}}"


J = json.load(open("results/gtex_target_avoidance.json"))
P, PH, G = J["per_mirna"], J["posthoc_targetscan_reference"], J["gates"]
ms = list(J["pairs"]); lab = [m.replace("hsa-", "") + "\n" + J["pairs"][m].split(" - ")[0] for m in ms]
x = np.arange(len(ms)); w = 0.2
fig, ax = plt.subplots(figsize=(9.2, 3.4))
for i, (v, c, n) in enumerate([([P[m]["cnn_matched"] for m in ms], "#2980b9", "CNN top-200 (non-conserved), matched tissue"),
                               ([P[m]["cnn_mismatched_mean"] for m in ms], "#aed6f1", "CNN, mean of other tissues"),
                               ([PH["per_mirna"][m]["matched"] for m in ms], "#c0392b", "TargetScan conserved top-200, matched (post-hoc)"),
                               ([PH["per_mirna"][m]["mismatched_mean"] for m in ms], "#f5b7b1", "TargetScan, mean of other tissues")]):
    ax.bar(x + (i - 1.5) * w, v, w, color=c, label=n)
ax.axhline(0.5, ls=":", c="k", lw=0.8); ax.set_ylim(0.35, 0.8)
ax.set_xticks(x); ax.set_xticklabels(lab, fontsize=7); ax.set_ylabel("AUC$_{dep}$ (0.5 = no depletion)")
ax.legend(fontsize=6.5, frameon=False, ncol=2, loc="upper right"); ax.set_title("GTEx v8 target avoidance in the miRNA's own tissue", fontsize=9)
fig.tight_layout(); fig.savefig("paper/figs/fig_gtex.pdf")
h1, h2, h3 = J["H1_cnn_depletion"], J["H2_targetscan_depletion"], J["H3_cnn_tissue_specificity"]
dp, sp = PH["depletion"], PH["specificity"]
S = [r"\subsection{Tissue target avoidance (GTEx): conserved targets avoid the miRNA's tissue, CNN non-conserved targets do not}",
     r"Targets of a tissue-specific miRNA should be relatively depleted in that tissue. We test this with GTEx v8 gene median TPM"
     r" (\texttt{scripts/gtex\_target\_avoidance.py}; \texttt{results/gtex\_target\_avoidance.json}; hermetic tests) on 8 pre-registered pairs"
     r" (miR-1-3p, miR-133b, miR-206: skeletal muscle; miR-122-5p: liver; miR-9-5p, miR-128-3p: brain cortex; miR-223-3p: whole blood; miR-375: pancreas)."
     r" Relative expression is $r_g(t)=\log_2(\mathrm{TPM}_{g,t}+1)-\mathrm{median}_{t'}\log_2(\mathrm{TPM}_{g,t'}+1)$ and depletion is the Mann--Whitney"
     r" probability $\mathrm{AUC}_{dep}=P(r_{\text{set}}<r_{\text{rest}})$."
     f" Gates passed: ALB peaks in liver, ACTA1 in skeletal muscle, median mapping {100*G['G3_median_mapping']:.1f}\\%.",
     f"CNN top-200 genes are not depleted in the matched tissue (H1: {h1['wins']} wins, {h1['losses']} losses, $p={fmtp(h1['p_one_sided'])}$) and not more depleted there"
     f" than in the other design tissues (H3: {h3['wins']}/{h3['losses']}, $p={fmtp(h3['p_one_sided'])}$); both are falsified. Pre-registered H2 (TargetScan genes inside"
     r" the CLIP universe) turned out to be empty by construction, since that universe excludes every gene with a conserved site for the seed; we record it as not"
     r" evaluable. A post-hoc reference arm that we did not pre-register uses all TargetScan 3'UTR genes as the universe and the top-200 conserved-site genes"
     f" as the set: they are depleted in the matched tissue for {dp['wins']}/8 pairs ($p={fmtp(dp['p_one_sided'])}$) and more so than in other tissues for"
     f" {sp['wins']}/8 ($p={fmtp(sp['p_one_sided'])}$; Figure~\\ref{{fig:gtex}}), which recovers the classic avoidance signal (Farh et al.\\ 2005).",
     r"So the method sees avoidance where it is expected, and the CNN's non-conserved candidates show none. Two caveats: the arms use different universes"
     r" (CNN scores exist only for the CLIP universe), and conservation and scorer are confounded, so this does not separate ``non-conserved sites are not"
     r" functional in vivo'' from ``the CNN ranks them poorly.'' Both readings are negative for the non-conserved discovery claim.",
     r"\begin{figure}[h]\centering\includegraphics[width=0.95\linewidth]{figs/fig_gtex.pdf}",
     r"\caption{Depletion (AUC$_{dep}$) of target sets in each miRNA's own GTEx tissue versus the mean over the other design tissues. Blue: CNN top-200 in the"
     r" CLIP universe (pre-registered). Red: TargetScan conserved top-200 in the all-UTR universe (post-hoc reference).}\label{fig:gtex}\end{figure}"]
open("paper/gtex_sec.tex", "w").write("\n".join(S) + "\n")
print("ok")
