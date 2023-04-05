# Summary

This is a two-stage computational workflow for virus discovery from sequencing data from the [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra) or your own data. Stage 1 (**Virushunter**) searches the raw reads using profile Hidden Markov Models. Stage 2 (**Virusgatherer**) perform a seed-based, iterative viral genome assembly that specifically targets the sequences identified in the first stage.

# Software dependencies

 * EMBOSS
 * seqtk
 * fastp
 * NCBI blast
 * NCBI SRA toolkit
 * HMMer
 * Genseed-HMM
 * CAP3
 * SPAdes
 * newbler
 * sfffile, sffinfo
 * Bowtie 2
 * snakemake

# Support

For questions or support, email chris.lauber at twincore.de

# License

[GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html)
