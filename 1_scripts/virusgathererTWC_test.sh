#!/bin/sh

# parameters:
#  1) family
#  2) input fastq file
#  3) project name
#  4) data set ID
#  5) seed file
#  6) base directory where to save temporary and result files
#  7) threads
#  8) assembler
#  9) protein refseq host+virus database
# 10) identifiers of viral entries in protein refseq database
# 11) flag for debugging mode [0=no, 1=yes]

./virusgathererTWC.pl \
	RNAviruses \
	/data/tests/SRA_n1/DRR067307.fastq.gz \
	SRA_n1 \
	DRR067307 \
	/data/virushunter/RNAviruses/SRA_n1/results/DRR067307/virushunter/contigs.singlets.fas.gz \
	/data/virushunter \
	24 \
	cap3 \
	/data/db/RefSeq/refseq_protein/refseq_protein \
	/data/db/RefSeq/viral_protein/viral_protein.acc_list \
	1

