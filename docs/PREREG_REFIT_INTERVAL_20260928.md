# Preregistration: refit-aware, leakage-free crossed family-and-gene interval (2026-09-28)

Queue item 5 remainder (docs/QUEUE_20260927.md): "A refit-aware, leakage-free
interval ... remain NOT DONE." The committed crossed interval (719db2d,
results/pooled_twoway_interval.json) holds the out-of-fold predictions fixed and
retains global scaling. Written 2026-09-28 before any refit-aware draw was run.
The design below is locked; it will not move after outcomes.

## Interpretation being locked (declared, parent can veto at hash report)
"Model refitting" in the queue remainder refers to the scoring pipeline behind
the endpoint - the per-arm logistic evaluators and the feature standardization -
not to retraining the DuplexCNN feature extractor. CNN retraining inside the
bootstrap is out of scope by cost: one DuplexCNN fit is ~7 min CPU (declared in
scripts/family_transfer.py), so 300 draws x 5 folds would be ~175 h. The CNN
feature columns ("min", "sum") are fixed outputs of a model trained on
TargetScan rows (labels disjoint from the ENCORI CLIP labels evaluated here);
that is stated, not hidden.

## Endpoint (unchanged from the committed analyses)
Delta = AUROC of the A+E+C arm minus AUROC of the A+E arm, where
A = [n8, n7, len], E = [lexpr], C = [min, sum]; rows = the committed
results/clip_rows_136.csv universe attached to HEK293 expression exactly as
scripts/expression_confound.py does (435,394 rows, 16,835 genes, 104 miRNA
families). Integrity gate: the refit pipeline at uniform weights must reproduce
the committed gene-grouped OOF AUROCs (results/expression_confound.json
pooled_auroc_gene_grouped_cv AE and AEC) to 1e-8, matching the check in
scripts/pooled_twoway_interval.py; any mismatch aborts before draw 1 and is
reported as a blocker, not a result.

## Locked procedure per draw (B = 300, seed 20260928, one rng stream, draw order
fixed; checkpointed every 25 draws, resume-safe)
1. Family counts fc = multinomial over the 104 families; gene counts gc =
   multinomial over the 16,835 genes (same construction as the committed
   fixed-prediction interval). Row weight w = fc[family(row)] * gc[gene(row)].
2. Weighted standardization: per-feature weighted mean and weighted sd over all
   rows with weights w (zero sd set to 1). No unweighted refit shortcut.
3. Folds: the SAME deterministic GroupKFold(5) split by gene identity as the
   committed pipeline (no shuffling; folds do not vary across draws - only the
   weights do).
4. Per arm (AE, AEC): fit sklearn LogisticRegression(max_iter=2000) on the
   weighted standardized training rows with sample_weight = w[train]; predict
   out-of-fold probabilities for every row.
5. Draw statistic: delta_b = weighted AUROC(y, p_AEC, w) - weighted AUROC(y,
   p_AE, w), weights applied in roc_auc_score(sample_weight=w).

## Reporting rule (no gate; this closes an uncertainty-reporting remainder)
- Report: point delta at uniform weights (must equal the committed
  +0.000768 by the integrity gate), the 300 draw deltas, percentile 95% CI,
  and n_nonpositive, plus the comparison to the committed fixed-prediction
  crossed CI [-0.001568, +0.002733].
- If the refit-aware CI includes zero, that is reported plainly as the
  strongest dependence-and-refit uncertainty statement; if it excludes zero,
  that is reported plainly too. Neither outcome moves any existing verdict,
  restores the withdrawn trained-model claim, or upgrades the headline; the
  headline remains the small ENCORI primary effect with its uncertainty.
- Degenerate-draw rule: a draw where any arm's effective weighted training
  class count is zero is recorded as skipped and reported (never silently
  replaced); B counts valid draws only if zero skips occur, otherwise the
  skip count is part of the result.

## Deliverables
- scripts/refit_interval.py (committed before any draw),
  results/refit_interval.json (+ checkpoint file).
- Queue item-5 addendum and a paper paragraph reporting the interval as
  measured.
