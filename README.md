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
 * newbler
 * Bowtie 2
 * snakemake
 * vsearch >=2.15.2 <2.20.0

# Blast databases

You need to install the following Blast databases and specify their file paths and names in the config.yaml:
 * refseq_protein (can be downloaded from https://ftp.ncbi.nlm.nih.gov/blast/db/)
 * viral_protein (the sequences can be downloaded from https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/, only accessions are required)
 * viral_genomic (the sequences can be downloaded from https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/; use `makeblastdb` command with `parse_seqids` to create a BLAST database)

NOTE: to download only RdRp-encoding RNA viruses, the following command can be used: `esearch -db nucleotide -query "txid2559587[Organism:exp] AND refseq[filter] NOT txid2732397[Organism:exp]" | efetch -format fasta > riboviria.no_pararnavirae.genomic.fna`

 * filter (see subfolder 4_databases; use `makeblastdb` command to create a BLAST database)

  *Note:* For detailed instructions on downloading Blast databases, please refer to our [GitHub Wiki](https://github.com/lauberlab/VirusHunterGatherer/wiki). 

# Support

For questions or support, email chris.lauber *at* twincore.de

# License

[GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html)

# References

Lauber C*, Seitz S*, Mattei S, Suh A, Beck J, Herstein J, Börold J, Salzburger W, Kaderali L, Briggs JAG, Bartenschlager R. Deciphering the Origin and Evolution of Hepatitis B Viruses by Means of a Family of Non-enveloped Fish Viruses. *Cell Host Microbe*. 2017 Sep 13;22(3):387-399.e6. [doi: 10.1016/j.chom.2017.07.019](https://pubmed.ncbi.nlm.nih.gov/28867387/).

Lauber C, Zhang X, Vaas J, Klingler F, Mutz P, Dubin A, Pietschmann T, Roth O, Neuman BW, Gorbalenya AE, Bartenschlager R, Seitz S. Deep mining of the Sequence Read Archive reveals major genetic innovations in coronaviruses and other nidoviruses of aquatic vertebrates. PLoS Pathog. 2024 Apr 22;20(4):e1012163. [doi: 10.1371/journal.ppat.1012163](https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1012163)

Lauber C, Chong LC. Viroid-like RNA-dependent RNA polymerase-encoding ambiviruses are abundant in complex fungi. *Frontiers Microbiology*. 2023 May 12; Volume 14. [https://doi.org/10.3389/fmicb.2023.1144003](https://doi.org/10.3389/fmicb.2023.1144003)

\* equal contribution
