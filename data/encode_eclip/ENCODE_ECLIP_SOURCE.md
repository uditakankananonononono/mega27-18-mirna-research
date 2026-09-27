# ENCODE AGO2 eCLIP inputs (lane-18 item 8, second CLIP resource)

Fetched 2026-09-28 per docs/PREREG_SECOND_CLIP_20260927.md eligibility example
("ENCODE eCLIP processed peak calls obtained directly from ENCODE").

- Experiment: ENCSR048YXG (eCLIP, AGO2, K562, human GRCh38)
  https://www.encodeproject.org/experiments/ENCSR048YXG/
- File: ENCFF111TWF.bed.gz - replicated peaks (biological replicates [1,2]),
  output_type "peaks", assembly GRCh38, 12,307 peak intervals.
  Download URL: https://www.encodeproject.org/files/ENCFF111TWF/@@download/ENCFF111TWF.bed.gz
- Per-replicate peak beds exist (ENCFF689OAC rep1, ENCFF578KGU rep2) but the
  replicated set is the analysis input.
- Not ENCORI/starBase and not a re-export of it: ENCODE portal, ENCODE eCLIP
  processing pipeline, independent consortium.

## Second resource attempt: DIANA TarBase v9.0 (2026-09-28)

- File: Homo_sapiens_TarBase-v9.tsv.gz (111,135,987 bytes, 4,724,537 interactions)
  URL: https://dianalab.e-ce.uth.gr/tarbasev9/data/Homo_sapiens_TarBase-v9.tsv.gz
  SHA-256: 172e383add39018eac7c53bb914142adba5aaea4e734a6e5a2c42eca661f244e
  (git-ignored: 106 MB; this recipe + hash regenerates it)
- Independence from ENCORI/starBase: TarBase v9 is DIANA Lab's own uniform
  re-analysis of raw CLIP datasets plus manual curation (their documentation:
  "Raw datasets ... were uniformly analyzed"); it is not an ENCORI re-export.
- Locked AGO-CLIP method whitelist (defined from method semantics, before any
  label or eligibility count): HITS-CLIP, PAR-CLIP, CLASH, qCLASH -
  4,545,851 interactions. All other methods (microarrays, RNA-Seq, luciferase,
  qPCR, etc.) are excluded as not AGO-CLIP-derived.
- Label semantics: (miRNA, gene) support call = pair present in the
  whitelisted interactions. Gene-level, mirroring the primary endpoint's
  ENCORI gene-level support label.
