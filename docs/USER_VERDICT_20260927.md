# User verdict archive, 2026-09-27 (lane 18 = "main (11): miRNA Target Scoring - DuplexCNN")

Provenance: inbound WhatsApp wamid.HBgMOTE4MTM0MDk4NTcxFQIAEhgWM0VCMEFGRDY4MzY4OTkxNzFEQURGRAA=,
received 2026-09-27 12:02:56 IST, author = user (+918134098571). Full-message body SHA256:
0795f5ff4a86b44e7edd301caa8c9c4ce17abbda034f0a0b29ae7b40ff76a6a6.
Header directive (verbatim, first line of the original message): "IGNORE ABOUT ISEF DELIVERABLES,
IMPROVE PAGE COUNT". Effect: page-length weaknesses fleet-wide are superseded; grow pages with
substantive content; 12-slide storyboard deliverables dropped.

## Verbatim section for this lane
3. main (11): miRNA Target Scoring — DuplexCNN
Weaknesses (20) — computational only:

Pair-map ablation changes correlation by +0.0029 only — branch is not load-bearing.

Expression confound shrinks CNN CLIP gain from +0.0048 to +0.0008.

miRDB v6.0 beats CNN on miRTarBase labels (0.787 vs 0.685).

TarBase v9.0 reproduces the miRDB gap.

Univariate CNN min-score wins 0/132 miRNAs vs seed count.

Pair-map costs 20× training time for near-zero gain.

No pathway coherence in CNN top-200 or CLIP sets.

STRING coherence explained by expression.

GTEx target avoidance: CNN nonconserved targets show none.

Enrichr reverse-lookup site-count win fails pre-registered replication.

Reactome gates failed (design errors).

LOEUF depletion survives only in weaker form; UTR-length clause withdrawn.

ClinGen, HPO, DECIPHER do not support "beyond UTR length."

Study bias explains curated-label depletion.

Random-weight CNN falsifier not run.

Tool ledger ~40 pending direct-use audit.

317 entries are record partitions, not studies.

No external benchmark not derived from TargetScan/ENCORI.

Non-conserved candidates unvalidated.

Paper reads as a graveyard of negatives.

Additions (computational):

Run the random-weight CNN falsifier.

Test exact-length-bin control.

Add a truly external miRNA target benchmark.

Make the expression confound the primary finding.

Report per-miRNA effect sizes with cluster-level CIs.

Pre-register a replication on a second CLIP dataset.

Deprecate the pair-map branch unless it earns its keep.

Reduce to one primary endpoint (CLIP AUROC conditional on expression).

Add a decision rule for when CNN scores add value.


## Verbatim cross-cutting themes

Cross-Cutting Computational Themes
Recurring weaknesses:

Dataset/tool count inflation (nested records counted as independent).

Long papers (47-58 pages) — not ISEF-ready.

Negative-heavy narratives that obscure positive contributions.

Single-seed headline numbers.

Ad hoc gates/thresholds rather than theoretically derived.

Homology leakage in random-split benchmarks.

No leave-family-out CV in most projects.

CIs often overlapping — point-estimate wins only.

Universal computational additions:

One primary question per paper.

One locked primary endpoint.

Cluster-level bootstrap CIs everywhere.

Leave-family-out CV as the primary protocol.

12-slide storyboard as the ISEF deliverable.

One-page summary card.
