# MEGA27-18 Revival Preregistration (locked 2026-09-26, before any new outcome)

Current audited state (README @ delegated HEAD): miRDB v6.0 beats DuplexCNN
on miRTarBase labels; HEK293 expression beats all sequence models on CLIP
labels - no SOTA claim stands. Falsifiable finding: CNN top-ranked
non-conserved candidates are less dosage-sensitive than matched genes on
gnomAD LOEUF / GeneBayes s_het, surviving coding-length, paralog and
C2H2-ZNF controls. Open falsifier declared in README: a random-weight CNN
producing the same depletion would expose a site-enumeration artefact.
28 negatives documented.

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
