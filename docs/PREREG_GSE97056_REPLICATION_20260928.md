# NEW dated preregistration: independent-library replication of the CLIP endpoint on GSE97056 (2026-09-28, 07:36 IST)

Lane 18 (DuplexCNN). Queue remainder: independent-library replication of the
expression-adjusted CLIP endpoint. Designation follows scouting commit 95d14b9
(docs/CLIP_LIBRARY_SCOUTING_LOG_20260928.md) under
docs/PREREG_THIRD_CLIP_LIBRARY_SCOUTING_20260928.md; parent approval of the
scouting outcome 2026-09-28 07:36 IST. Written and committed BEFORE any fetch
of the designated resource. Thresholds and the decision rule below are locked
and will not move after outcomes.

## Designated resource (locked)
- GSE97056, "Argonaute CLIP defines a deregulated miR-122 bound transcriptome
  that correlates with patient survival in human liver cancer (Human
  CLIP-Seq)", Homo sapiens, 28 samples, PMID 28735896. SuperSeries GSE97061.
  Series supplementary file: GSE97056_HumanLiver_AgoCLIP_Matrix.txt.gz.
- Fallback if the primary file proves unparsable for the locked label
  construction: GSE98670 (human HK2 AGO2 PAR-CLIP, PMID 28925592). The
  fallback may be invoked only for a parse-level failure of the primary, and
  the invocation reason is recorded in the fetch log before any fallback
  outcome is computed. Tertiary GSE41437 is NOT auto-designated; using it
  requires another dated amendment.

## Designation re-verification (integrity gate, before label construction)
1. GSE97056, GSE97061, GSE98670 and PMIDs 28735896/28925592 are re-checked
   against the archived ENCORI AGO study list
   (results/encori_ago_study_list.json, 361 datasets / 37 GEO series / 41
   PMIDs) and against ENCSR048YXG and TarBase v9.0. Any hit aborts the
   replication as INELIGIBLE-overlap; no outcome is computed.
2. Fetch log committed BEFORE parsing: exact download URL, byte size and
   SHA-256 of the fetched file. Only the recorded file is used.

## Label construction (locked, mirrors PREREG_SECOND_CLIP_20260927.md)
1. The candidate row set is the committed primary-endpoint row universe
   verbatim (same (miRNA, gene) pairs and covariates as the locked primary
   endpoint); only the support labels are replaced.
2. A pair is GSE97056-positive iff a deposited AGO-CLIP peak maps to that
   gene's 3'UTR under the committed annotation (same UTR/gene universe as the
   primary endpoint) AND the peak region carries a canonical seed match to
   that miRNA under the pipeline's committed seed definitions. Peaks mapping
   outside 3'UTRs label nothing. If the deposited matrix is not miRNA-resolved
   (expected: it is an AGO-binding matrix), this seed assignment is the
   standard derivation declared in the scouting prereg criterion 1.
3. Same gene filter as the primary endpoint: genes carrying a conserved site
   for any miRNA sharing the seed are excluded.
4. SUFFICIENCY GATE: per-miRNA eligibility >=50 positives and >=50 negatives
   after gene filtering; the resource is decisive only if >=50 miRNAs are
   eligible. Fewer eligible miRNAs records the resource as
   attempted-but-ineligible and no inference is drawn (same convention as the
   ENCODE eCLIP attempt). Eligibility counts are committed with the fetch log
   before any AUROC is computed.

## Locked analysis (mirrors the second-resource test exactly)
Same features (site counts, log UTR length, HEK293 HPA nTPM expression, CNN
min/sum from the frozen duplex_cnn.pt), same models (logistic regression,
gene-grouped 5-fold CV, AE vs AEC), no re-tuning. Primary metric: pooled
delta AUROC (AEC - AE). Secondary: per-miRNA deltas; miRNA-cluster percentile
bootstrap (10,000 replicates, seed 20260928) of the median per-miRNA delta.
Declared limit: features use HEK293 expression while labels derive from liver
tissue; the cross-tissue mismatch is inherent to mirroring the locked primary
endpoint and is stated, not repaired.

## Locked decision rule
- alpha = 0.01, one-sided in the direction of the committed point estimate.
- Replication supports the independent-library condition only if the pooled
  delta is positive AND its within-resource gene-resampled bootstrap CI
  (B=1000) excludes 0, AND the sufficiency gate passed, AND designation
  re-verification passed.
- A null, negative, or ineligible result is reported as-is; thresholds do not
  move; the independent-library condition then remains unmet (documented
  exhaustion if GSE98670 also fails parse or sufficiency).
- No result from this replication can upgrade the pair-map branch, the
  candidate list, or any tool/dataset gate.

## Out of scope
No new architecture, no CNN retraining, no pooling with ENCORI/TarBase labels,
no threshold changes, no additional files beyond the recorded fetch. Any
deviation requires a new dated commit before outcomes are inspected.
