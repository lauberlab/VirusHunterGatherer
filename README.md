# Summary

This is a two-stage computational workflow for data-driven virus discovery from sequencing data from the [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra) or your own data. Stage 1 (**Virushunter**) searches the raw reads using profile Hidden Markov Models. Stage 2 (**Virusgatherer**) perform a seed-based, iterative viral genome assembly that specifically targets the sequences identified in the first stage.

# Documentation

* [Software dependencies](#software-dependencies)
* [Installation](#installation)
* [Blast databases](#blast-databases)
* [Support](#support)
* [License](#license)
* [References](#references)

## Software dependencies

 * snakemake
 * Perl
 * EMBOSS
 * seqkit
 * fastp
 * NCBI blast
 * NCBI SRA toolkit
 * HMMer
 * vsearch
 * CAP3
 * Bowtie 2 (optional)

## Installation

1. Download the current repository (Click **Code --> Download ZIP**. Extract the ZIP file and navigate into the extracted folder.)
2. Install the requirements. To ease the process, we recommend building a [Conda](https://github.com/conda-forge/miniforge) environment using the provided `vhvg.yaml` file:

```{bash}
conda env create --name environment_name --file vhvg.yaml   # be patient, it takes some time
conda activate environment_name
```

Alternatively, the environment can be created using this command:

```{bash}
conda create --name environment_name -c bioconda blast cap3 emboss fastp hmmer perl seqkit snakemake sra-tools vsearch #bowtie2
conda activate environment_name
```
3. Before running the tool, set up the required [Blast databases](#blast-databases) and [edit the `config.yaml` file](TODO).

## Blast databases

During pipeline execution, Virushunter and Virusgatherer rely on several Blast databases to filter and annotate assembled contigs. For general use and initial test runs, we recommend setting up the following databases:

### Contaminant DB (DBFILTER in `config.yaml`)

A custom set of contaminant sequences is provided and must be indexed prior to use by running:

```{bash}
cd 4_databases/
makeblastdb -in filter.fasta -dbtype nucl -parse_seqids
```

### RefSeq protein DB (DBREFSEQ in `config.yaml`)

The database is available from the NCBI Blast FTP repository (https://ftp.ncbi.nlm.nih.gov/blast/db/). For downloading and updating the database, we recommend using the `update_blastdb.pl` utility provided with NCBI Blast:

```{bash}
cd 4_databases/   # or any preferred location for storing the database files
mkdir refseq_protein; cd refseq_protein/
update_blastdb.pl --decompress refseq_protein
```

**Important**: The full database can exceed 100 GB in size. Even if this filtering step is disabled (by setting "ENABLE_FILTER_3: 0" in `config.yaml`), the database remains required, as Virusgatherer uses it for viral contig annotation. As an alternative, a reduced viral RefSeq protein database can be used:

```
cd 4_databases/   # or any preferred location for storing the database files
mkdir refseq_protein; cd refseq_protein/
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz
gunzip viral.1.protein.faa.gz
makeblastdb -in viral.1.protein.faa -dbtype prot -parse_seqids -out viral_refseq_protein
```

### Viral genome DB (DBVIRAL in `config.yaml`)

The database is available from the NCBI Blast FTP repository (https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/). It can be downloaded and set up by running:

```{bash}
cd 4_databases/   # or any preferred location for storing the database files
mkdir viral_genomic; cd viral_genomic/
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz
gunzip viral.1.1.genomic.fna.gz
makeblastdb -in filter.fasta -dbtype nucl -parse_seqids -out viral_genomic
```

### Viral RefSeq protein DB (ACCSVIRAL in `config.yaml`)

The database is available from the NCBI Blast FTP repository (https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/). However, only a list of accession identifiers is required for subsequent queries against the RefSeq protein DB created above. This list can be generated using:

```
cd 4_databases/   # or any preferred location for storing the database files
mkdir viral_protein; cd viral_protein/
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz
gunzip viral.1.protein.faa.gz
grep ">" viral.1.protein.faa | cut -d " " -f 1,1 | cut -c 2- > viral_protein.acc_list
```

<!--
NOTE: to download only RdRp-encoding RNA viruses, the following command can be used: `esearch -db nucleotide -query "txid2559587[Organism:exp] AND refseq[filter] NOT txid2732397[Organism:exp]" | efetch -format fasta > riboviria.no_pararnavirae.genomic.fna`
-->

## Support

For questions or support, email chris.lauber *at* twincore.de

## License

[GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html)

## References

Neuman BW*, Smart A*, Gilmer O*, Smyth RP*, Vaas J, Böker N, Samborskiy DV, Bartenschlager R, Seitz S, Gorbalenya AE, Caliskan N, Lauber C. Giant RNA genomes: Roles of host, translation elongation, genome architecture, and proteome in nidoviruses. Proc Natl Acad Sci U S A. 2025 Feb 18;122(7):e2413675122. [doi: 10.1073/pnas.2413675122](https://www.pnas.org/doi/10.1073/pnas.2413675122)
 
 * **GenBank TPA accessions**: BK070317-BK070395
 * More data and information about the paper can be found [here](https://github.com/lauberlab/invertebrate_nidovirus_discovery_paper).

Lauber C*, Zhang X*, Vaas J, Klingler F, Mutz P, Dubin A, Pietschmann T, Roth O, Neuman BW, Gorbalenya AE, Bartenschlager R, Seitz S. Deep mining of the Sequence Read Archive reveals major genetic innovations in coronaviruses and other nidoviruses of aquatic vertebrates. PLoS Pathog. 2024 Apr 22;20(4):e1012163. [doi: 10.1371/journal.ppat.1012163](https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1012163)

 * **GenBank TPA accessions**: BK070274-BK070316

Lauber C, Chong LC. Viroid-like RNA-dependent RNA polymerase-encoding ambiviruses are abundant in complex fungi. *Frontiers Microbiology*. 2023 May 12; Volume 14. [https://doi.org/10.3389/fmicb.2023.1144003](https://doi.org/10.3389/fmicb.2023.1144003)

Lauber C*, Seitz S*, Mattei S, Suh A, Beck J, Herstein J, Börold J, Salzburger W, Kaderali L, Briggs JAG, Bartenschlager R. Deciphering the Origin and Evolution of Hepatitis B Viruses by Means of a Family of Non-enveloped Fish Viruses. *Cell Host Microbe*. 2017 Sep 13;22(3):387-399.e6. [doi: 10.1016/j.chom.2017.07.019](https://pubmed.ncbi.nlm.nih.gov/28867387/).

\* equal contribution
