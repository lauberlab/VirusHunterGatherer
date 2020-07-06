#!/bin/sh

# parameters:
#  1) family
#  2) input fastq file
#  3) project name
#  4) data set ID
#  5) base directory where to save temporary and result files
#  6) threads
#  7) filter blast database
#  8) protein refseq host+virus database
#  9) nucleotide refseq virus-only database
# 10) identifiers of viral entries in protein refseq database
# 11) path to workflow base directory
# 12) flag whether to filter against human genome [0=no, 1=yes]
# 13) flag for debugging mode [0=no, 1=yes]

./virushunterTWC.pl \
	RNAviruses \
	/data/tests/SRA_n1/DRR067307.fastq.gz \
	SRA_n1 \
	DRR067307 \
	/data/virushunter \
	24 \
	/data/db/RefSeq/filter/filter \
	/data/db/RefSeq/refseq_protein/refseq_protein \
	/data/db/RefSeq/viral_genomic/viral_genomic \
	/data/db/RefSeq/viral_protein/viral_protein.acc_list \
	/home/lauber/workflows/01_VirusHunterGatherer \
	0 \
	1

