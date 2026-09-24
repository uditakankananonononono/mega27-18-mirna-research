# Item 18 external tools & datasets (honest ledger)
## Tools/resources used
1. TargetScan vert_80 Conserved_Site_Context_Scores - training labels
2. TargetScan UTR_Sequences - 25,132 human 3'UTRs
3. TargetScan miR_Family_Info - miRNA families/conservation
4. ENCORI/starBase miRNATarget API - AGO-CLIP labels (23 miRNAs)
5. miRBase (via TargetScan IDs)
6. PyTorch, 7. scikit-learn, 8. NumPy, 9. pandas, 10. matplotlib, 11. pytest
Count: 11. Planned: miRTarBase (site currently 404 - retry via browser), TarBase, miRWalk, RNAhybrid/ViennaRNA.
## Datasets (accession-level)
1. TargetScan conserved sites (264,563 human rows)
2. TargetScan 3'UTR sequences
3. TargetScan miR family info
4-26. ENCORI per-miRNA CLIP sets (23 fetched: miR-7-5p, let-7a, miR-21, 17, 19a, 92a, 16, 24, 26a, 27a, 29a, 30a, 103a, 125b, 130a, 137, 155, 181a, 199a, 200c, 221, 23a, 1)
Count: 26. Extending ENCORI panel to ~100 miRNAs gets past 120 honestly.
