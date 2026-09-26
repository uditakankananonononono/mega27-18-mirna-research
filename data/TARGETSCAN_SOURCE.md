# TargetScan vert_80 inputs (not bundled - fetch recipe)

Downloaded 2026-09-26 from the public TargetScan 8.0 release (all HTTP 200,
zip integrity verified with `unzip -t`):

- https://www.targetscan.org/vert_80/vert_80_data_download/miR_Family_Info.txt.zip
  -> data/miR_Family_Info.txt (verbatim)
- https://www.targetscan.org/vert_80/vert_80_data_download/Conserved_Site_Context_Scores.txt.zip
  -> filtered to Gene Tax ID 9606, columns
     transcript,gene,mirna,site_type,start,end,contextpp
  -> data/human_sites.tsv (264,563 rows + header; matches docs/TOOLS.md count)
- https://www.targetscan.org/vert_80/vert_80_data_download/UTR_Sequences.txt.zip
  -> filtered to Species ID 9606, alignment gaps '-' stripped, columns
     transcript,gene,seq
  -> data/human_utrs.tsv (28,351 ENST rows + header; TOOLS.md's 25,132 figure
     refers to the used subset after the inner join in src/mirna/data.py)
- ENST version suffixes stripped from the ID column of both TSVs so the
  versioned sites file joins the unversioned UTR keys (src/mirna/data.py
  strips versions when keying UTRs).

This recipe is the source of truth for regenerating the inputs; the TSVs are
git-ignored like the original pipeline's.
