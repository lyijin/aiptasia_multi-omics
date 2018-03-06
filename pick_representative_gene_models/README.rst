=============================
Removal of ambiguous peptides
=============================
Some detected peptides from the proteomics analysis can originate from multiple (protein) gene models--this is a problem for downstream analyses e.g. functional enrichment.

Take for example a peptide that is deemed differentially expressed, and the peptide maps to AIPGENE123, AIPGENE124 and AIPGENE125. It is very likely that these genes share similar functional annotations. If we didn't reduce ambiguity, the GO terms associated with these three proteins will be far more enriched than what they should have been.

The input file, ``prot.fold_changes.tsv`` comes from the last sheet of the giant proteome analysis Excel sheet.

It gets processed by ``generate_dep.py`` (DEP stands for differentially expressed proteins), to produce the numerous text files seen in this folder.
