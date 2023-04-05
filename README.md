# Summary

This is a two-stage computational workflow for data-driven virus discovery from sequencing data from the [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra) or your own data. Stage 1 (**Virushunter**) searches the raw reads using profile Hidden Markov Models. Stage 2 (**Virusgatherer**) perform a seed-based, iterative viral genome assembly that specifically targets the sequences identified in the first stage.

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

For questions or support, email chris.lauber *at* twincore.de

# License

[GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html)

# References

Lauber C*, Seitz S*, Mattei S, Suh A, Beck J, Herstein J, Börold J, Salzburger W, Kaderali L, Briggs JAG, Bartenschlager R. Deciphering the Origin and Evolution of Hepatitis B Viruses by Means of a Family of Non-enveloped Fish Viruses. *Cell Host Microbe*. 2017 Sep 13;22(3):387-399.e6. [doi: 10.1016/j.chom.2017.07.019](https://pubmed.ncbi.nlm.nih.gov/28867387/).

Lauber C, Vaas J, Klingler F, Mutz, M, Gorbalenya AE, Bartenschlager F, Seitz S. Deep mining of the Sequence Read Archive reveals bipartite coronavirus genomes and inter-family Spike glycoprotein recombination. *bioRxiv*. [doi: https://doi.org/10.1101/2021.10.20.465146](https://www.biorxiv.org/content/10.1101/2021.10.20.465146v1) 

\* equal contribution
