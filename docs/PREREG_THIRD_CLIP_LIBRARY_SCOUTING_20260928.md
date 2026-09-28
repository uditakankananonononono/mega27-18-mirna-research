# NEW dated preregistration: scouting for a third, library-independent human AGO-CLIP deposit (2026-09-28, 07:29 IST)

Lane 18 (DuplexCNN). Pivot under the 2026-09-28 posture, parent instruction
2026-09-28 07:26 IST: scouting for a usable second-CLIP-library deposit as the
independent replication path; a scouting null is the documented exhaustion
proof. This file is committed BEFORE any candidate is inspected.

## Background
The primary CLIP endpoint is single-source (ENCORI/starBase AGO-CLIP).
The committed second-resource test (TarBase v9.0 labels + ENCODE K562 AGO2
eCLIP, accession ENCSR048YXG / file ENCFF111TWF) passed its locked
within-resource rule but is explicitly NOT independent-library proof, because
resources can share underlying CLIP experiments. The miRTarBase 10.0 SE
benchmark is essentially CLIP-free but is a validated-target benchmark, not a
CLIP library. The queue's remaining NOT DONE is exactly this:
independent-library replication of the expression-adjusted CLIP endpoint.

## Locked scouting question
Does a human AGO-CLIP deposit exist that is (a) independent of every library
already used in this repository and (b) sufficient to run the locked
replication analysis of docs/PREREG_SECOND_CLIP_20260927.md?

## Locked eligibility criteria (ALL must hold)
1. Human AGO CLIP-derived miRNA-gene support calls (HITS-CLIP, PAR-CLIP,
   iCLIP, eCLIP, or CLEAR-CLIP), or raw/processed data sufficient to derive
   them with standard peak calling.
2. Independence: the deposit's underlying experiments must not appear in the
   ENCORI/starBase study list used for the primary endpoint, must not be
   ENCSR048YXG, and must not be TarBase v9/miRTarBase re-exports. Aggregator
   resources (POSTAR3, CLIPdb, starBase itself) are usable only as INDEXES to
   underlying studies; a deposit whose experiments overlap the ENCORI study
   list (checked by accession/PMID) is recorded INELIGIBLE-overlap.
3. Gene symbols mappable to the primary-endpoint universe.
4. After gene filtering, plausibly >= 50 miRNAs with >= 50 positives and
   >= 50 negatives (same sufficiency bar as the second-CLIP prereg; judged
   from study scale/metadata during scouting, verified only under a follow-up
   prereg if a candidate is designated).
5. Publicly accessible raw (SRA) or processed (BED/table) files via
   GEO/ArrayExpress/publisher.

## Locked scouting procedure
- Metadata/header scouting only: GEO/ArrayExpress accession metadata, study
  pages, publication records. NO label construction, NO model scoring, NO
  overlap statistics beyond accession/PMID-level disjointness checks.
- Every candidate inspected is recorded in the scouting log with: accession,
  URL, library type, cell system, study/PMID, and the eligibility verdict
  with the specific criterion that failed when rejected.
- Candidate enumeration sources: GEO DataSets queries for human AGO CLIP,
  ArrayExpress search, POSTAR3/CLIPdb study indexes (as indexes only), and
  targeted literature for post-2020 human AGO CLIP studies. The enumeration
  is recorded so a null is an exhaustion proof, not an assumption.

## Locked decision
- SCOUTING POSITIVE (>=1 deposit passes all five criteria): no analysis is
  run. The deposit is reported to the parent with its scouting-log entry; any
  replication analysis requires its own NEW dated preregistration designating
  the resource before label construction.
- SCOUTING NULL (every inspected candidate rejected with recorded reasons):
  the scouting log is the documented exhaustion proof for the
  independent-library path; the queue records independent-library replication
  as exhausted with evidence; thresholds and prior verdicts do not move.
- Scouting inspects metadata only and therefore has no performance gate;
  nothing in scouting can upgrade or downgrade any existing claim.
