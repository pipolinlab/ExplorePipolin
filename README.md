Pipolins constitute a new group of self-synthesizing or self-replicating 
mobile genetic elements (MGEs). They are widespread among diverse bacterial 
phyla and mitochondria.

> [**Redrejo-Rodríguez, M., *et al.*** Primer-independent DNA synthesis 
>by a family B DNA polymerase from self-replicating Mobile genetic elements. 
>*Cell reports*, 2017](https://doi.org/10.1016/j.celrep.2017.10.039)
>
>[**Flament-Simon, S.C., de Toro, M., Chuprikova, L., *et al.*** High diversity 
>and variability of pipolins among a wide range of pathogenic *Escherichia 
>coli* strains. *bioRxiv*, 2020](https://www.biorxiv.org/content/10.1101/2020.04.24.059261v1)

 **ExplorePipolin** is a search tool that identifies and analyses
 pipolins within bacterial genome.

# Table of contents

* [Requirements](#requirements)
* [Installation](#installation)
    * [Install from source](#install-from-source)
    * [Install using Conda](#install-using-conda)
* [Quick usage](#quick-usage)
    * [Test run](#test-run)
    * [Output files](#output-files)
* [Running with Docker](#running-with-docker)

# Requirements

 * pip
 * [BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK279690/)
 * [ARAGORN](https://github.com/TheSEED/aragorn)
 * [Prokka](https://github.com/tseemann/prokka)

# Installation
### Install from source

 1. Install the requirements (see above).
 1. TODO: `wget <from_release>`
 1. `unzip ExplorePipolin-master.zip && cd ExplorePipolin-master` 
 1. `pip install .` (install in user site-package) or
 `sudo pip install .` (requires superuser privileges)
 
NOTE: before installing, it is possible to run unit tests:
`pytest` or `python setup.py test` (from the source root directory).
 
**How to uninstall:**

`(sudo) pip uninstall ExplorePipolin`

### Install using Conda
 
 * Before installing ExplorePipolin, make sure you'are running the latest 
 version of Conda:
 
 `conda update conda`
 
 `conda install wget`
 
 * Create a new environment that is specific for ExplorePipolin. You can 
 choose whatever name you'd like for the environment.
 
 `wget github/TODO/explore_pipolin-0.0.1.yml`
 
 `conda env create -n ExplorePipolin-0.0.1 --file explore_pipolin-0.0.1.yml`
 
 * Download and install ExplorePipolin into the created environment:
 
 `wget github/TODO/explore-pipolin-0.0.1-py_0.tar.bz2`
 
 `conda install -n ExplorePipolin-0.0.1 explore-pipolin-0.0.1-py_0.tar.bz2`
 
  * Clean up (optional):
 
 `rm explore_pipolin-0.0.1.yml explore-pipolin-0.0.1-py_0.tar.bz2`
 
 * Activate the environment and check the installation:
 
 `conda activate ExplorePipolin-0.0.1`
 
 `explore_pipolin -h`

# Quick usage

### Test run
As input, **ExplorePipolin** takes FASTA file(s) with genome sequence(s). 
A genome sequence can be either a single complete chromosome (preferred) 
or contigs (in a single multiFASTA file).

```bash
--> explore_pipolin -h
Usage: explore_pipolin [OPTIONS] GENOMES...

  ExplorePipolin is a search tool that identifies and analyses
  pipolin elements  within bacterial genome(s).

Options:
  --out-dir PATH  [required]
  -h, --help      Show this message and exit.
```

TODO: test run.

### Output files

The output directory will contain several folders:
 
 | Folder | Content description |
 |--------|---------------------|
 | `pipolbs_search` | BLAST search results for piPolB genes |
 | `atts_search` | BLAST search results for the known *att* sites |
 | `atts_denovo_search` | Results of *de novo* search for *att* sites |
 | `trnas_search` | ARAGORN search results for tRNAs/tmRNAs |
 | `pipolin_sequences` | extracted pipolin sequences in FASTA format |
 | `prokka_results` | Prokka annotation results (check files description [here](https://github.com/tseemann/prokka/blob/master/README.md#output-files))|
 | `results` | GenBank and GFF annotation results with the *att*s included, log files |


# Running with Docker

See https://docs.docker.com/install/ to install Docker.

**NOTE:** superuser privileges are required to run the analysis and around 3GB of disk space for the image.

```
sudo docker pull link/TODO/to/the/image/on/DockerHub
sudo docker run --rm explore_pipolin -h
sudo docker run --rm -v $(pwd):/output -w /output explore_pipolin 
 --out-dir output ./input_genomes/*.fa   #(example run)
```
