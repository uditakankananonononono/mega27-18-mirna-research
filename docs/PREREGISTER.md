# MEGA27-18 Revival Preregistration (locked 2026-09-26, before any new outcome)

Current audited state (README @ delegated HEAD): miRDB v6.0 beats DuplexCNN
on miRTarBase labels; HEK293 expression beats all sequence models on CLIP
labels - no SOTA claim stands. Falsifiable finding: CNN top-ranked
non-conserved candidates are less dosage-sensitive than matched genes on
gnomAD LOEUF / GeneBayes s_het, surviving coding-length, paralog and
C2H2-ZNF controls. Open falsifier declared in README: a random-weight CNN
producing the same depletion would expose a site-enumeration artefact.
28 negatives documented.

## Provenance of this document (honest attribution)
- USER STANDING RULES (verbatim, WhatsApp channel history): 4:11:18 (complete
  all projects except deleted ones; ask CHATGPT for ideas/redirection; minimum
  10 judging rounds on weaknesses/additions; never count a negative as a
  result), 4:11:49 (each project beats benchmarks - improve until it does -
  and produces an actual new discovery), 4:12:25 (ask ChatGPT how to redirect
  when a negative is not moving forward), 4:14:37 (take inspiration from
  previous ISEF winners, e.g. Natasha Kulviwat).
- RESEARCHER-LOCKED METHODOLOGICAL CHOICES (this agent, 2026-09-26, locked
  before inspecting new outcomes): every numeric threshold, alpha level,
  seed count, pivot ladder, gate name and scope framing below. These are the
  lane's own preregistration decisions, NOT user-specified values; they exist
  so results cannot be fished past moving goalposts. Pre-existing gates
  declared by earlier builders in repo history (e.g. the RMSD < 2.0 A redock
  gate already in this repo's README) are inherited, not invented here.

## Locked gates (declared before outcomes)
- G1 falsifier first: run the locked random-weight CNN control on the
  identical site-enumeration + matching pipeline before any further claim.
  Discovery stands only if the trained model's depletion is significantly
  stronger than the random-weight distribution (pre-declared alpha 0.01,
  100 random seeds).
- G2 benchmark beat (fair task): miRTarBase-label ranking is settled against
  us. The beat must come from a task where the comparison is honest: ranked
  options to be ordered by a rule-6 ChatGPT redirection round before
  testing: (a) per-site scoring calibration on held-out CLIP crosslink
  clusters vs the published context++ model on identical sites; (b) the
  dosage-sensitivity screen itself as a method benchmark vs named published
  prioritization baselines on identical gene sets; (c) an independent
  replication of the depletion on an untouched constraint resource.
- G3 discovery: if G1 survives, the dosage-sensitivity depletion of
  predicted non-conserved targets, replicated on a second independent
  constraint source, is the new discovery. If G1 kills it, pivot ladder
  from the redirection round applies; no negative is terminal (rule 4).
- Judge: >= 10 ChatGPT rounds, verbatim docs/JUDGE_ROUNDS.md (file does not
  exist yet at HEAD; created with round 1).
- ISEF archetype: identify + verify biomarker-style discipline applied to
  predicted regulatory targets; independent verification is the headline.


## Judge requirement amendment, 2026-09-27 10:00 IST
The user changed the numeric requirement from ten ChatGPT checks to ONE
round she provides (original WhatsApp wamid.HBgMOTE4MTM0MDk4NTcxFQIAEhgWM0VCMDJCMTZGRTVEMkQwMTFBQzc4MQA=). The older ten-round
text above records the earlier protocol, not the current finish line. All
previous judge transcripts remain intact, but prior agent-initiated rounds
are not assumed to satisfy the new user-provided courier round without a
verified project-specific handoff. Any supplementary Gemini or other LLM
consult is separate and does not count as the user-provided ChatGPT verdict.

## Judge clarification, 2026-09-27 10:01 IST
The user clarified: "EACH PROJECTS NEED ONE FROM ME TO PASS" (WhatsApp
wamid.HBgMOTE4MTM0MDk4NTcxFQIAEhgWM0VCMDMwREI5RDQ0QUNCRDc2MTNDMwA=).
This lane currently has **0 of 1 user-provided ChatGPT verdicts**; the earlier
agent-initiated ChatGPT interactions remain preserved as supplementary history,
not counted toward this one-verdict gate. A courier paste from the user and
its project-specific original wamid must be recorded before declaring pass.
