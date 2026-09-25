# MEGA27-18: miRNA target scoring (DuplexCNN) and a falsification audit of non-conserved candidates

Code, data manifests, results and paper for item 18 of MEGA27.

## What is here
- `src/mirna/` - DuplexCNN (pairing-map 2D CNN plus Grimson context features), trained on TargetScan vert_80 context++ labels.
- `scripts/` - one script per external tool or test. Each pre-registered test script was committed before its results; the commit history shows the order.
- `results/` - committed JSON outputs. Every number in the paper comes from these files.
- `paper/main.tex`, `paper/main.pdf` - full paper (Times, numbered equations, figures, 28 documented negatives).
- `docs/TOOLS.md` - ledger of external tools (43, infrastructure excluded) and accession-level datasets (317).
- `tests/` - hermetic pytest suite (`python3 -m pytest -q tests`).

## Shipped tool: `mirtarget-score`
```
pip install -e .
mirtarget-score --mirna hsa-miR-7-5p --gene SP1 EGFR --top 5
mirtarget-score --seq UGGAAGACUAGUGAUUUUGUUGU --utr-fasta my_utrs.fa --json out.json
```
It scores 8mer/7mer-A1 seed-match sites in 3'UTRs with the frozen DuplexCNN and reports per-site scores on the context++ scale (lower = stronger predicted repression).

## Main findings (honest summary)
- Benchmarks broken against us: miRDB v6.0 beats the CNN on miRTarBase labels, and HEK293 expression alone beats all sequence models on CLIP labels. No state-of-the-art claim is made.
- Named falsifiable finding: CNN top-ranked non-conserved candidates are less dosage-sensitive than expression- and publication-matched genes on LoF-count labels (gnomAD LOEUF, GeneBayes s_het). The finding survives controls for coding length, paralogs and C2H2-ZNF families. On curated labels (HPO, Orphanet, ClinGen, DECIPHER vs UTR length) the depletion is explained by study bias or UTR length. The earlier "beyond UTR length" clause was withdrawn.
- Open falsifier: a random-weight CNN producing the same depletion would show it is a site-enumeration artefact.

## Reproduce
Data sources and accession lists are in `docs/TOOLS.md`. Run a test with e.g. `python3 scripts/mgi_mouse_ko.py`, then rebuild its paper section with the matching `scripts/*_paper.py` and run `pdflatex` three times in `paper/`.
