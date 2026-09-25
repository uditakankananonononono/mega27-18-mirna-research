#!/usr/bin/env python3
"""LaTeX paragraph for the exp_lof-matched control (committed JSON only)."""
import json


def fmtp(p):
    if p >= 0.001: return f"{p:.3f}"
    m, e = f"{p:.1e}".split("e"); return f"{m}\\times10^{{{int(e)}}}"


J = json.load(open("results/explof_matched.json")); d1, l1, l2 = J["D1_cnn_shorter_cds"], J["L1_loeuf_vs_explof_matched"], J["L2_shet_vs_explof_matched"]
S = [r"\paragraph{Pre-registered coding-length control.} We tested the explanation above directly (\texttt{explof\_matched}, committed before running),"
     r" matching on gnomAD's expected LoF count, a proxy for coding length."
     f" CNN top-200 genes do have fewer expected LoF variants than their universe ({d1['wins']}/{d1['losses']} miRNAs), but only by a median factor of"
     f" {d1['median_ratio']:.2f}. The matching worked (matched/CNN expected-LoF ratio {J['matching_check_median_explof_ratio']:.3f}). The depletion survives:"
     f" CNN LOEUF is higher than in expected-LoF-matched sets for {l1['wins']}/{l1['losses']} miRNAs ($p={fmtp(l1['p_one_sided'])}$) and $s_\\mathrm{{het}}$ lower"
     f" for {l2['wins']}/{l2['losses']} ($p={fmtp(l2['p_one_sided'])}$). So the coding-length explanation is falsified. The disagreement between LoF-count labels and"
     r" the other labels stays unexplained, and the clause about UTR length stays withdrawn."]
open("paper/explof_sec.tex", "w").write("\n".join(S) + "\n"); print("ok")
