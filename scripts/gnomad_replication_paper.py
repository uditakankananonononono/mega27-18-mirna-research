#!/usr/bin/env python3
"""Figure + LaTeX for the pre-registered gnomAD replication (committed JSON only)."""
import json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def fmtp(p):
    if p >= 0.001: return f"{p:.3f}"
    m, e = f"{p:.1e}".split("e"); return f"{m}\\times10^{{{int(e)}}}"


J = json.load(open("results/gnomad_replication.json")); P = J["per_mirna"]; ms = J["panel"]
d = {k: [P[m]["cnn"] - P[m][k + "_mean"] for m in ms] for k in ("uniform", "expr", "len")}
lenr = [P[m]["cnn_median_len"] / P[m]["universe_median_len"] for m in ms]
fig, ax = plt.subplots(1, 2, figsize=(9.0, 3.0), gridspec_kw={"width_ratios": [2, 1]})
ax[0].boxplot([d["uniform"], d["expr"], d["len"]], tick_labels=["uniform", "expression-\nmatched", "UTR-length-\nmatched"])
for i, k in enumerate(("uniform", "expr", "len")): ax[0].scatter(np.full(len(ms), i + 1) + np.random.default_rng(i).uniform(-.12, .12, len(ms)), d[k], s=8, c="#2980b9", zorder=3)
ax[0].axhline(0, c="k", lw=.8); ax[0].set_ylabel("CNN minus control median LOEUF"); ax[0].set_title("Held-out panel (30 miRNAs)", fontsize=9)
ax[1].hist(lenr, bins=12, color="#7f8c8d"); ax[1].axvline(1, c="k", lw=.8); ax[1].set_xlabel("CNN / universe median UTR length"); ax[1].set_title("UTR length", fontsize=9)
fig.tight_layout(); fig.savefig("paper/figs/fig_gnomad_rep.pdf")
r1, r2, r3, G = J["R1_cnn_vs_uniform"], J["R2_cnn_vs_expr"], J["R3_cnn_vs_utrlen"], J["gates"]
md = lambda k: float(np.median(d[k]))
S = [r"\paragraph{Pre-registered held-out replication.} We committed a replication script before running it (\texttt{gnomad\_replication}):"
     f" {G['G1_n']} CLIP miRNAs outside the discovery panel, a one-sided test in the discovered direction, and a new UTR-length-matched control"
     f" (decile matching on the longest 3$'$ UTR). Gates passed ({100*G['G2_coverage']:.1f}\\% coverage). The reversal replicates: the CNN top-200 is less"
     f" constrained than uniform random in {r1['wins']}/{r1['wins']+r1['losses']} miRNAs ($p={fmtp(r1['p_one_sided'])}$; median LOEUF difference {md('uniform'):+.3f}),"
     f" than expression-matched sets in {r2['wins']}/{r2['wins']+r2['losses']} ($p={fmtp(r2['p_one_sided'])}$; {md('expr'):+.3f}), and than UTR-length-matched sets"
     f" in {r3['wins']}/{r3['wins']+r3['losses']} ($p={fmtp(r3['p_one_sided'])}$; {md('len'):+.3f}; Figure~\\ref{{fig:gnomadrep}}). CNN top-200 genes have shorter"
     f" UTRs than their universe in {J['cnn_shorter_utr_than_universe']}/30 miRNAs, and length matching shrinks the gap, but it does not remove it.",
     r"We name this as a falsifiable finding: \emph{among non-conserved seed-match genes, the CNN's top-ranked candidates are depleted of loss-of-function-constrained"
     r" genes, beyond what HEK293 expression or 3$'$ UTR length explain.} It would be falsified by (i) a finer length control (for example exact length bins or"
     r" site-density matching) that removes the gap, or (ii) the same depletion appearing for a random-weight CNN, which would show it comes from the input encoding"
     r" rather than the learned model. Neither test has been run. One reading is purifying selection against strong non-conserved sites in dosage-sensitive genes,"
     r" but we have not tested it and do not claim it.",
     r"\begin{figure}[h]\centering\includegraphics[width=0.95\linewidth]{figs/fig_gnomad_rep.pdf}",
     r"\caption{Pre-registered replication on 30 held-out miRNAs. Left: CNN top-200 median LOEUF minus each control (positive = CNN less constrained)."
     r" Right: CNN top-200 median UTR length relative to its universe.}\label{fig:gnomadrep}\end{figure}"]
open("paper/gnomad_rep_sec.tex", "w").write("\n".join(S) + "\n"); print("ok")
