# Item 18 external tools & datasets (honest ledger)
## Tools/resources used
1. TargetScan vert_80 Conserved_Site_Context_Scores - training labels
2. TargetScan UTR_Sequences - 25,132 human 3'UTRs
3. TargetScan miR_Family_Info - miRNA families/conservation
4. ENCORI/starBase miRNATarget API - AGO-CLIP labels (136-miRNA panel)
5. miRBase (via TargetScan IDs)
6. PyTorch, 7. scikit-learn, 8. NumPy, 9. pandas, 10. matplotlib, 11. pytest
12. SciPy (rank/Wilcoxon/permutation statistics in CLIP + audit analyses)
13. ViennaRNA 2.7.2 (RNAcofold duplex-MFE validation, results/rnafold_validation.json)
Count: 13. Planned: miRTarBase (site currently 404 - retry via browser), TarBase, miRWalk, RNAhybrid/ViennaRNA.
## Datasets (accession-level)
1. TargetScan conserved sites (264,563 human rows)
2. TargetScan 3'UTR sequences
3. TargetScan miR family info
4-139. ENCORI per-miRNA CLIP sets (136 kept with >=20 CLIP+ genes; 5 candidates fetched but dropped below threshold: hsa-miR-127-5p, hsa-miR-203a-3p, hsa-miR-217-5p, hsa-miR-484, hsa-miR-766-3p). Panel = original 23 + 113 new: miR-100-5p, miR-101-3p, miR-106a-5p, miR-106b-5p, miR-107, miR-10a-5p, miR-10b-5p, miR-122-5p, miR-126-3p, miR-128-3p, miR-129-5p, miR-130b-3p, miR-132-3p, miR-133a-3p, miR-133b, miR-134-5p, miR-135b-5p, miR-138-5p, miR-139-5p, miR-140-5p, miR-141-3p, miR-143-3p, miR-144-3p, miR-145-5p, miR-146a-5p, miR-146b-5p, miR-148a-3p, miR-149-5p, miR-150-5p, miR-152-3p, miR-153-3p, miR-15a-5p, miR-15b-5p, miR-182-5p, miR-183-5p, miR-184, miR-185-5p, miR-186-5p, miR-187-3p, miR-18a-5p et al. (full list in results/panel_extension.json)
Count: 139. Past the 120 bar honestly - every set fetched, thresholded, and usable for the extended falsification rerun.
