===========================================================
Multi-omics analysis of thermal stress response in Aiptasia
===========================================================

Link to manuscript: (still under review)

Disclaimer
----------
The title alludes to multiple -omics analysis--but a careful browser of this repository would notice that the stuff in here applies mainly to the **PROTEOMIC** bits of the project, with transcriptomic bits nowhere to be seen.

This is because the transcriptomics pipelines we used were... fairly standard (``kallisto``, ``sleuth``). Our proteomics stuff weren't, therefore the scripts. Robert Frost's "The Road Not Taken" comes to mind here :p

Brief description of folder contents
------------------------------------
Folders are ordered according to when it appears in the manuscript.

1. ``remove_duplicate_AIPGENEs/`` explains how (and why) we excluded ~10% Aiptasia gene models from our analysis.

2. ``proteome_analysis/`` is the folder that contains the giant Excel sheet that we used to carry out our proteome analysis.

3. ``proteome_glm/`` was a GLM that we ran to check whether strain/temperature/a combination of both leads to significant changes in protein expression.

4. ``compare_abs_values/`` contains scripts that led to plotting of Fig. 2.

5. ``compare_fc/`` contains scripts that led to plotting of Fig. 3 and Fig. S1.

6. ``pick_representative_gene_models/``... picks representative gene models to avoid bias in the downstream functional enrichment analysis.
