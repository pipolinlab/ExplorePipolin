# PipolinFinder

Pipolins constitute a new group of self-synthesizing or self-replicating 
mobile genetic elements (MGEs). They are widespread among diverse bacterial 
phyla and mitochondria.

> [Redrejo-Rodríguez, Modesto, et al. "Primer-independent DNA synthesis 
>by a family B DNA polymerase from self-replicating Mobile genetic elements." 
>Cell reports 21.6 (2017): 1574-1587.](https://doi.org/10.1016/j.celrep.2017.10.039)

 **PipolinFinder** is a search tool that identifies, extracts and annotates 
 pipolin elements within bacterial genome.

## Installation
**Install the required dependencies:**

python packages:
 * pip
 * click
 * biopython
 * bcbio-gff
 * prefect

other packages:
 * [Entrez Direct](https://www.ncbi.nlm.nih.gov/books/NBK179288/) 
 (not strictly required)
 * [BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK279690/)
 * [ARAGORN](https://github.com/TheSEED/aragorn)
 * [Prokka](https://github.com/tseemann/prokka)

**To install PipolinFinder:**

 1. Click green button on the right "Clone or Download" => "Download zip"
 1. `unzip ExplorePipolin-master.zip && cd ExplorePipolin-master`
 1. `pip install .` (install in user site-package) or
 `sudo pip install .` (requires superuser privileges)
 
#### How to uninstall:

`(sudo) pip uninstall PipolinFinder`

## Quick usage

As input, PipolinFinder takes FASTA file(s) with genome sequence(s). 
A genome sequence can be either a single complete chromosome (preferred) 
or contigs (in a single multiFASTA file).

**Test**

```bash
--> pipolin_finder -h
Usage: pipolin_finder [OPTIONS] GENOMES...

  PipolinFinder is a search tool that identifies, extracts and annotates
  pipolin elements  within bacterial genome(s).

Options:
  --out-dir PATH  [required]
  -h, --help      Show this message and exit.
```

TODO: test run.

## Running using Docker

See https://docs.docker.com/install/ to install Docker.

NOTE: you will need superuser privileges and around 3GB of disk space to build 
the container.

 1. Click green button on the right "Clone or Download" => "Download zip"
 1. `unzip ExplorePipolin-master.zip && cd ExplorePipolin-master/docker`
 1. `sudo docker build -t pipolin_finder .`
 1. `sudo docker run --rm pipolin_finder --help` (test)
 1. `sudo docker run --rm -v $(pwd):/output -w /output pipolin_finder 
 --out-dir output ./input_genomes/*.fa` (example run)

#### Restrictions of the pipeline (to discuss)

 1. The task `find_atts_denovo` works only for complete genomes (SKIPPED for 
 incomplete genomes). It looks for direct repeats of minimal length 6
 around piPolB(s) with a 100 kbp window.
 1. For our dataset of *E. coli* genomes, no "denovo" atts were found. For that
 reason, only "usual" atts will be used as pipolin boundaries. But, if no 
 "usual" atts were found, "denovo" atts will be used instead. For more details, 
 check `are_atts_present` task.
 1. If no "usual" or "denovo" atts were found, the pipeline is not going to
 continue the analysis (FAIL).
 1. The pipeline is expecting a single pipolin per genome, just because 
 we haven't seen genomes with two or more pipolins so far 
 (`NotImplementedError`).
 1. From the above, we are expecting a single piPolB or several piPolBs 
 (the gene might be disrupted or duplicated) within some restricted area
 flanked by att repeats.
 1. TODO: scaffolding restrictions!

## Background

The main modules (scripts): (TODO: some of them not working, REMOVE or FIX)
 * `download_metadata_ncbi.py` -- downloads the metadata for the analysed 
 genomes, such as accessions, organism and strain names
 * `download_genomes_ncbi.py` -- downloads genome (chromosome) sequences 
 given NCBI assembly accession (i.e. for a non-complete genome, it 
 downloads all its contigs)
 * `identify_pipolins_roughly.py` using a reference `pi-polB.fa` sequence
 * `analyse_pipolin_orientation`
 * `extract_pipolin_regions.py` It is possible to extract regions, 
 determined by the leftmost and rightmost atts (`--long` option).
 * `annotate_pipolins.py`
 * `predict_atts_with_hmmer.py`
 * `store_new_att_bounds.py` -- parses HMMER output for atts
 * `include_atts_into_annotation.py` 
 TODO: So far, the script includes atts only into GB and GFF files.
 * `scaffold_gapped_pipolins.py`
 * `easyfig_add_colours.py`
 
Prediction of ATTs:

 1. Prepared ATT sequences with `prepare_atts_for_msa.py`
 ```
The total number of atts is 198
> Maximum att length is 132
> Minimum att length is 121
```
 1. Built a MSA with att sequences using MAFFT, MUSCLE
 and T-Coffee (https://www.ebi.ac.uk/Tools/msa). 
 The output format -- Pearson FASTA, otherwise long
 sequence names might be truncated.
 1. Compared the alignments. Modified them, using 
 Jalview: deleted not conserved regions from both ends.
 1. Created HMM profile with `hmmbuild` and `hmmpress`.


To extract the subsequence from a genome:
 * `get_subsequence.py`
 
 `$ get_subsequence.py genomes/NZ_JNMI01000006.1.fa 80191 82792 pi-polB.fa`
 
 `$ get_subsequence.py genomes/NZ_JNMI01000006.1.fa 64241 64373 attL.fa`

 `$ get_subsequence.py genome/NZ_JNMI010000006.1.fa 90094 90008 tRNA.fa`

 For Saskia's strains: 
 * `edit_contig_names.sh <in-dir>` -- to shorten the long contig names

To get the sequences from roary groups:
 * `extract_roary_groups.py`
