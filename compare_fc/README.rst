=========================================================================
Compiling and plotting fold-changes in transcript and protein expressions
=========================================================================

The script ``compute_fc.py`` takes in results from RNA-seq (the ``*.csv`` files are provided in the compressed form here to save space) and the proteome analysis (``prot.fold_changes.tsv``) to produce ``fold_changes.tsv``.

This table was then used by ``plot.per-strain.fc_values.py`` to produce Fig. 3; and by ``plot.per-strain.fc_values.dep.py`` to produce Fig. S1. The pdfs (provided here) produced by the scripts were subsequently prettified with Adobe Illustrator.
