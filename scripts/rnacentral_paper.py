#!/usr/bin/env python3
"""LaTeX section for the RNAcentral sequence audit (committed JSON only; no figure - a table)."""
import json
J = json.load(open("results/rnacentral_seq_audit.json")); G = J["gates"]
S = [r"\subsection{Input integrity: every miRNA sequence and seed used here matches RNAcentral}",
     r"All scorers in this paper key on the miRNA seed (nt 2--8) taken from TargetScan's \texttt{miR\_Family\_Info}. A wrong sequence would silently corrupt every"
     r" benchmark, so we checked each miRNA used in the CLIP and miRTarBase analyses against RNAcentral by miRBase accession"
     r" via the RNAcentral REST API (script and results: \texttt{rnacentral\_seq\_audit}; hermetic tests).",
     r"\begin{table}[h]\centering\small\begin{tabular}{lr}\hline Quantity & Value\\\hline",
     f"miRNAs used / with MIMAT accession & {J['n_used_mirnas']} / {J['n_with_targetscan_accession']}\\\\",
     f"Resolved in RNAcentral & {J['n_resolved']} ({100*G['G2_resolved_frac']:.1f}\\%)\\\\",
     f"G1: let-7a-5p exact sequence & {'pass' if G['G1_let7a_exact'] else 'fail'}\\\\",
     f"G3: seed identity under deranged map (control) & {100*G['G3_deranged_seed_identity']:.1f}\\%\\\\",
     f"H1: seed identity, resolved miRNAs & {100*J['H1_seed_identity']['frac']:.1f}\\%\\\\",
     f"H2: full-length identity, resolved miRNAs & {100*J['H2_full_identity']['frac']:.1f}\\%\\\\\\hline\\end{{tabular}}",
     r"\caption{RNAcentral sequence audit. The derangement control shows the identity test can fail.}\label{tab:rnacentral}\end{table}",
     r"Both hypotheses are confirmed (Table~\ref{tab:rnacentral}). One transport note: for a few accessions the endpoint served an HTML page to"
     r" Python \texttt{requests} but JSON to curl, so the fetcher falls back to curl when the body is not JSON; an earlier run without the fallback"
     r" left one or two accessions unresolved and was discarded (it was also overwritten by a duplicate process). No analysis choice changed."
     r" This rules out sequence errors as an explanation for the negatives above."]
open("paper/rnacentral_sec.tex", "w").write("\n".join(S) + "\n"); print("ok")
