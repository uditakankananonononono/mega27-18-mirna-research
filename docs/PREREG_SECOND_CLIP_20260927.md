# Pre-registration: replication of the expression-adjusted CLIP endpoint on a second CLIP resource

Date: 2026-09-27 (committed before any second-resource outcome is computed).
Lane: 18 (DuplexCNN). Verdict item: "Pre-register a replication on a second CLIP dataset".

## Background
The locked primary endpoint (paper, commit 92bb788) is pooled AUROC of CLIP
positivity under gene-grouped 5-fold CV for covariates [site counts, UTR length,
HEK293 expression] versus covariates + CNN [min, sum], on ENCORI AGO-CLIP rows.
Committed point estimate: +0.00077 (AE 0.69417 -> AEC 0.69494, 435,394 rows).
ENCORI is one consortium; a single-source result can be pipeline-specific.

## Dataset eligibility (locked)
A replication resource qualifies only if ALL hold:
1. Human AGO CLIP-derived (miRNA, gene) support calls from a resource that is not
   ENCORI/starBase and not a re-export of it (e.g., POSTAR3, or ENCODE eCLIP
   processed peak calls obtained directly from ENCODE).
2. Gene symbols mappable to the same universe used for the primary endpoint.
3. Per-miRNA pair counts sufficient for >= 50 miRNAs with >= 50 positives and
   >= 50 negatives after gene filtering; otherwise the resource is recorded as
   attempted-but-ineligible and no inference is drawn from it.
If no qualifying resource can be obtained, this pre-registration is reported as
not executable and the primary endpoint stays single-source (stated as a limit).

## Locked analysis (mirrors the primary endpoint exactly)
1. Same gene filter: genes carrying a conserved site for any miRNA sharing the
   seed are excluded; same feature definitions (n8, n7, log UTR length, CNN
   min/sum from the frozen duplex_cnn.pt, HEK293 HPA nTPM expression).
2. Same models: logistic regression, gene-grouped 5-fold CV, covariates AE vs
   AEC. No re-tuning, no feature additions.
3. Primary replication metric: pooled delta AUROC (AEC - AE).
4. Secondary: per-miRNA deltas on eligible miRNAs; miRNA-cluster percentile
   bootstrap (10,000 replicates, seed 20260927) of the median per-miRNA delta,
   matching scripts/cluster_interval.py.

## Locked decision rule
- alpha = 0.01, one-sided in the direction of the committed point estimate.
- Replication supports the primary claim only if the pooled delta is positive
  AND its within-resource bootstrap CI (B=1000, rows resampled by gene) excludes 0.
- A null or negative result is reported as-is; thresholds do not move; the
  primary endpoint reverts to "single-source, not replicated".
- No result from this replication can upgrade the pair-map branch, the discovery
  candidate list, or any tool/dataset gate.

## Out of scope
No new architecture, no retraining of the CNN, no threshold changes, no
multiple-resource pooling. Any deviation requires a new dated commit before
outcomes are inspected.
