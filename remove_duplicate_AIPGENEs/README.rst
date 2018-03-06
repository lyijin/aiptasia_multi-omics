==================================================
Removing cloned AIPGENEs from Aiptasia gene models
==================================================
Intellectual input/work in this section was shared with Guoxin Cui.

Theoretical considerations
--------------------------
One might wonder, why go through the trouble of modifying the underlying gene models (cDNA / protein sequences) prior to the RNA-seq analysis? What's the problem about the existing gene models in the first place? Seems to be much ado about nothing, no?

The whole hoopla started when we noticed that there were duplicate genes in Aiptasia--there are some cDNA and protein sequences that are **EXACTLY** identical to each other, down to the last nucleotide/amino acid. Biologically, it's possible: very recent gene duplication; technically (and more likely): it's because the assembler finds two fairly distinct haplotypes, and treated both haplotypes as distinct genomic loci. If that loci had genes in it, there'll be two copies of it.

As our RNA-seq quantifier of choice was ``kallisto``, we tested how it behaved when faced with duplicates of the same gene. Let's say AIPGENE123 and AIPGENE124 were identical clones of each other. When we map reads to them, would ``kallisto``:

1. Assign all mappable reads to AIPGENE123, and none to AIPGENE124?
2. Assign half to AIPGENE123, and half to AIPGENE124?
3. Assign none to AIPGENE123, and all to AIPGENE124?

The answer was... none of the above. We couldn't see a pattern in its behaviour. This is a problem when we're trying to assess differential expression across genes--if different samples had had different proportions of reads assigned to the cloned genes, then we might observe differential expression from biases in read assignment.

One might argue, why use ``kallisto`` if this was a problem for us? We unfortunately have a very strong liking of this pipeline (it's bloody fast, and its results are verifiable in the lab), so the next best solution was to remove clones from the sequence file. Bit of a hack job really, but it didn't really take us very long in practice (typing this thing up and explaining our ideas take about as long--there's a reason why people don't like to document stuff!).

Starting material
-----------------
The Aiptasia gff3 file (v1.0) contains a few genes that were... split across multiple scaffolds, which broke some of the scripts I wrote. I would've understood it if half of the gene model was at one end of scaffold X, while the other half were at the start of scaffold Y, but... they weren't like that. I had to manually go into the gff3, and select for the gene model that produced the cDNA or protein sequences in the companion FASTA files, deleting the rest that didn't make sense.

The genes affected were:

- AIPGENE629
- AIPGENE10730
- AIPGENE13859
- AIPGENE13860
- AIPGENE19416
- AIPGENE19417

This compressed file is provided as ``aiptasia_genome.dups_removed.gff3.gz``. There are still 29,269 gene mdoels in this gff file, similar to the v1.0 gff3 file.

Checking for near-clones that occupy the same loci
--------------------------------------------------
While there are examples where distinct genes share similar coding regions, from observation (i.e. genome browser gazing), many of these overlapping gene models tend to be exact copies/minor variants of each other. I wrote ``check_overlapping_genes.py`` to check which genes overlap--and if they do, which "isoform" is the longer one?

A longer explanation of how things work is in the script itself, including how to use its arcane output. When this script is run on ``aiptasia_genome.dups_removed.gff3.gz``, the 4th column of the output file tells me which genes are the clones that I should be removing from the eventual dataset.

1,716 genes fall foul of this cutoff. They were removed using ``pick_reads.py`` with ``--exclude``, producing a file with 27,553 gene models.

Checking for clones at different loci
-------------------------------------
This was done by checking the cDNA sequence of genes: if clones were detected, they were removed with a shell command.

``paste - - < aip.genome_models.no_isoforms.mRNA.fa | sort -u -f -k 2,2 | sed 's/\t/\n/' > aip.genome_models.no_isoforms.no_dups.mRNA.fa``

The 27,504 gene models that passed these filters can be found in ``../ultimate_universe.txt``. The file is in the root folder of the project because it was used by other scripts, e.g. ``../generate_DEP``.

If you'd like to regenerate the gene model files that we use (cDNA or protein), run

``pick_reads.py aip.genome_models.file_of_interest.fa --include ultimate_universe.txt``.

The script ``pick_reads.py`` can be found at https://github.com/lyijin/common.
