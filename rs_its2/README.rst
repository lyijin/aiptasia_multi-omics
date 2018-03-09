============================
Identifying ITS2 types in RS
============================

Trimming reads
--------------
MiSeq reads of the 4 Red Sea replicates were trimmed with the command

``while read a b c d e; do cutadapt_miseq.sh $b $c $d $e *${a}*R1* *${a}*R2*; done < adapter_seqs.tsv``

This produces reads in the form::

  liewy@kw14764:~/kaust/maha/tx_vs_prot/rs_its2/trimmed_reads$ ls -l
  -rw-r--r-- 1 liewy liewy 16M Mar  7 12:08 RS-25-1_R1.trim.fastq.gz
  -rw-r--r-- 1 liewy liewy 24M Mar  7 12:08 RS-25-1_R2.trim.fastq.gz
  -rw-r--r-- 1 liewy liewy 32M Mar  7 12:09 RS-25-2_R1.trim.fastq.gz
  -rw-r--r-- 1 liewy liewy 43M Mar  7 12:09 RS-25-2_R2.trim.fastq.gz
  -rw-r--r-- 1 liewy liewy 20M Mar  7 12:10 RS-25-3_R1.trim.fastq.gz
  -rw-r--r-- 1 liewy liewy 26M Mar  7 12:10 RS-25-3_R2.trim.fastq.gz
  -rw-r--r-- 1 liewy liewy 16M Mar  7 12:10 RS-25-4_R1.trim.fastq.gz
  -rw-r--r-- 1 liewy liewy 22M Mar  7 12:10 RS-25-4_R2.trim.fastq.gz

Merging reads
-------------
R1s and R2s were merged using ``bbmerge`` at default settings.

``for a in 1 2 3 4; do ~/tools/bbmap/bbmerge.sh in1=RS-25-${a}_R1.trim.fastq in2=RS-25-${a}_R2.trim.fastq out=RS-25-${a}.merged.fastq outu1=RS-25-${a}_R1.unmerged.fastq outu2=RS-25-${a}_R2.unmerged.fastq; done``

Unmerged reads were tossed out at this stage.

Clustering reads
----------------
The merged files were concatenated together with ``cat``, then clustered at 97% similarity with ``cdhit-est``.

``cdhit-est -i RS-25.merged.fa -o RS-25.merged.c0.97.fa -T 20 -M 0 -g 1 -c 0.97 -d 100``

This command produces ``RS-25.merged.c0.97.fa.clstr`` (which is provided here as a compressed file).

My script ``parse_clstr.py`` then tallies the number of reads each cluster had, and sorts it in ascending order. Hence, the last 10 lines of ``RS-25.merged.c0.97.fa.clstr.tsv`` are the cluster representatives of the 10 largest clusters.

The read IDs were noted down and extracted from ``RS-25.merged.c0.97.fa``, then subjected to a BLASTN search with

``blastn -db nt -query 0.97_top10.fa -outfmt 5 -max_target_seqs 20 -num_threads 10 -out 0.97_top10_vs_nt.blastn.xml``

The results from the BLASTN search was how I identified which clade/subclade the *Symbiodinium* belonged to--I preferred hits that had the lowest e-value, and in events of ties, I picked the one with the most detailed description (e.g. between "clade A" and "clade A4", the latter is preferred).
