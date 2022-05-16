![banner](banner.svg)

Pipolins constitute a new group of self-synthesizing or self-replicating 
mobile genetic elements (MGEs). They are widespread among diverse bacterial 
phyla and mitochondria.

> [**Redrejo-Rodríguez, M., *et al.*** Primer-independent DNA synthesis 
>by a family B DNA polymerase from self-replicating Mobile genetic elements. 
>*Cell reports*, 2017](https://doi.org/10.1016/j.celrep.2017.10.039)
>
>[**Flament-Simon, S.C., de Toro, M., Chuprikova, L., *et al.*** High diversity 
>and variability of pipolins among a wide range of pathogenic *Escherichia 
>coli* strains. *Scientific Reports*, 2020](https://www.nature.com/articles/s41598-020-69356-6#Sec18)

 **ExplorePipolin** is a search tool that identifies and analyses
 pipolins within bacterial genomes.

# Table of contents

* [Requirements](#requirements)
* [Installation](#installation)
    * [Install from source](#install-from-source)
    * [Install using Conda](#install-using-conda)
* [Quick usage](#quick-usage)
    * [Test run](#test-run)
    * [Output files](#output-files)

# Requirements

 * pip
 * [HMMER](http://hmmer.org/)
 * [BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK279690/)
 * [Prokka](https://github.com/tseemann/prokka)
 * [ARAGORN](https://github.com/TheSEED/aragorn)
 * [Prodigal](https://github.com/hyattpd/Prodigal)

# Installation
### Install from source

 1. Install the requirements (see above).
 1. `wget https://github.com/liubovch/ExplorePipolin/archive/0.0.1.zip`
 1. `unzip 0.0.1.zip && cd ExplorePipolin-0.0.1` 
 1. `pip install .` (install in user site-packages) or
 `sudo pip install .` (requires superuser privileges)
 
NOTE: before installing, it is possible to run unit tests:
`pytest` or `python setup.py test` (from the source root directory).
 
**How to uninstall:**

`(sudo) pip uninstall ExplorePipolin`

### Install using Conda

 * Before installing ExplorePipolin, make sure you're running the latest 
 version of Conda:
 
 `conda update conda`

 * Install into a new environment in one-step:

 `conda create -n <new_env_name> -c bioconda -c conda-forge -c defaults -c liubovch explore-pipolin`

NOTE: solving the environment takes time. Be patient.

# Quick usage

### Test run
As input, **ExplorePipolin** takes FASTA file(s) with genome sequence(s). 
A genome sequence can be either a single complete chromosome (preferred) 
or contigs (in a single multiFASTA file).

```bash
--> explore_pipolin -h
Usage: explore_pipolin [OPTIONS] GENOME...

  ExplorePipolin is a search tool for prediction and analysis of pipolins,
  bacterial mobile genetic elements.

Options:
  --out-dir-prefix TEXT       Use this prefix for the output directory,
                              instead of the default "results" prefix.

  --out-dir PATH              Use this output directory instead.
  --pipolb-hmm-profile PATH   piPolB's HMM profile to use as 1st priority.If
                              not provided, the default profile will be used
                              instead.

  --ref-att PATH              Att sequence in FASTA file to use as 1st
                              priority. If not provided, the default file will
                              be used instead.

  --percent-identity INTEGER  Minimum percent identity for direct repeats
                              search [default: 85]

  --max-inflate INTEGER       If no borders of pipolin are found (no ATTs),
                              inflate the analysed region from both sides of
                              piPolB.  [default: 30000]

  --no-annotation             Do not run the annotation step (i.e. Prokka).

  --proteins PATH             Prokka param: FASTA or GBK file to use as 1st
                              priority. If not provided, the default file will
                              be used instead.

  --skip-colours              Do not add an Easyfig-compatible colouring
                              scheme to the final Genbank files.

  --cpus INTEGER              Prokka param: Number of CPUs to use [0=all]
                              [default: 8]

  --keep-tmp                  Preserve intermediate files produced during the
                              run (it might be useful for debugging).

  -h, --help                  Show this message and exit.
```

### Output files

The output directory will contain a separate folder for each analysed genome. In the genome-specific folder,
the following files can be found:
 
 | Folder/File | Content description |
 |--------|---------------------|
 | `<genome>.log` | log file |
 | `pipolins` | extracted pipolin sequences in FASTA format, GenBank and GFF annotations of pipolin genes with *att*s included |

Files in `pipolins` folder have the following scheme: `<genome>_N_vN.<type>.<ext>`,
where:

 * `<genome>` — genome name, retrieved from a provided FASTA file without 
   extension. Should not exceed 16 characters, due to Biopython restrictions!
 * `N` — number of a found pipolin
 * `vN` — number of a reconstructed variant (each pipolin can be reconstructed in
   different ways)
 * `<type>` — can be `complete` (att---pol---att), `truncated` (att---pol---) or 
   `minimal` (---pol---)
 * `<ext>` — file extension (`.fa`, `.gbk` or `.gff`)  
 
Additional files, when `--keep-tmp` option is present:

 | Folder/File | Content description |
 |-------|----------------------|
 | `<genome>.fa` | genome file copied |
 | `pipolbs` | HMMER search results for piPolB genes |
 | `atts` | BLAST search results for the known *att* sites |
 | `atts_denovo` | Results of *de novo* search for *att* sites |
 | `trnas` | ARAGORN search results for tRNAs/tmRNAs |
 | `prokka` | Prokka annotation results (check files description [here](https://github.com/tseemann/prokka/blob/master/README.md#output-files))|
