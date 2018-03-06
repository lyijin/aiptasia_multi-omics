==================================================================
Compiling and plotting absolute transcript and protein expressions
==================================================================

Do transcripts and protein levels of the same gene generally correspond to each other?

TL;DR: Yep... ish.

Nitty-gritty
------------
The script ``compile_tx_vs_prot.py`` takes in abundance values from RNA-seq (``kallisto.tpms.tsv.gz`` is compressed to save space, and tpm = transcripts per million) and the proteome analysis (``prot.scpms.tsv.gz``, where scpm = spectral counts per million) to produce the per-strain ``*.tsv`` files in the same folder.

The plotting script ``plot.per-strain.abs_values.py`` uses the ``*.tsv`` files to produce Fig. 2. The pdf produced by the script (``per-strain.abs_values.pdf``) was subsequently prettified with Adobe Illustrator.
