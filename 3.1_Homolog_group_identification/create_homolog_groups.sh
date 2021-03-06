#!/bin/bash

##### blast all-vs-all
makeblastdb -in allProteins.fasta -dbtype prot -out BLAST.DB.allProteins

# this command was used to add additional blast runs later
# all blast runs need to be comparable
blastp -db BLAST.DB.allProteins -query allProteins.fasta -outfmt 6 -out all-vs-all.tsv -num_alignments 3000 -dbsize 499679 -evalue 1e-5 -matrix=BLOSUM62 -num_threads 11

# default command would be:
# blastp -db BLAST.DB.allProteins -query allProteins.fasta -outfmt 6 -out all-vs-all.tsv -evalue 1e-5 -matrix=BLOSUM62 -num_threads 11

#
# reduce blast hits according:
#   alignment length > 60%
#   identity > 50%
#   sequence coverage > 60%
#   e-value < 1e-20
#

filter_blast_file.py -bi all-vs-all.tsv -bo extracted.all-vs-all.tsv

##### clustering

cut -f 1,2,11 extracted.all-vs-all.tsv > extracted.all-vs-all.abc

mcxload -abc extracted.all-vs-all.abc --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o seq.mci -write-tab seq.tab

mcl seq.mci -I 1.2

mcxdump -icl out.seq.mci.I11 -tabr seq.tab -o dump.seq.mci.I11
