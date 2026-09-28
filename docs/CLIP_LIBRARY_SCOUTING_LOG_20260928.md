# Third CLIP-library scouting log (2026-09-28)

Executed under docs/PREREG_THIRD_CLIP_LIBRARY_SCOUTING_20260928.md (committed
64e17ba before any candidate was inspected). Metadata/header scouting only:
no label construction, no model scoring, no overlap statistics beyond
accession/PMID-level disjointness.

Reference set: ENCORI CLIP reference table (hg38_clipRef via the ENCORI
download endpoint), 1,984 CLIP datasets, 361 AGO datasets spanning 37 GEO
series and 41 PMIDs, archived verbatim at results/encori_ago_study_list.json.
Already-used resources excluded: ENCORI/starBase aggregate (primary endpoint),
TarBase v9.0 (second-resource test), ENCODE eCLIP ENCSR048YXG/ENCFF111TWF.

Enumeration: NCBI GEO DataSets (eutils), four human queries (argonaute CLIP;
AGO CLIP-seq; AGO2 eCLIP; PAR-CLIP argonaute). Every series-level record was
classified; full verdicts and per-accession reasons are in
results/clip_library_scouting.json.

Outcome: SCOUTING POSITIVE.
- Primary candidate GSE97056 (human liver pan-AGO CLIP, 28 samples, PMID
  28735896, processed matrix deposited): passes all five locked criteria at
  metadata level. miRNA resolution requires seed-assignment derivation
  (allowed by criterion 1); the >=50 miRNAs x >=50 positives sufficiency bar
  is plausible at this scale but UNVERIFIED by design - it can only be
  checked under a designated-resource prereg.
- Secondary GSE98670 (human HK2 AGO2 PAR-CLIP, 21 samples, PMID 28925592).
- Tertiary GSE41437 (LCL PAR-CLIP cellular targetome, EBV caveat, raw-only).
- 12 series rejected on accession-level ENCORI overlap; 4 on same-study
  subseries; 2 universe; 11 type; 2 breadth; 3 scale.

Per the locked decision rule, no analysis was run. Next step: NEW dated
prereg designating GSE97056 (fallback GSE98670) for the independent-library
replication question before any label construction.
