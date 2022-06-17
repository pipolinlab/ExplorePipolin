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
    * [Command line options](#command-line-options)
    * [Output files](#output-files)
    * [Easyfig-compatible colouring](#easyfig-compatible-colouring)

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
 1. `wget https://github.com/pipolinlab/ExplorePipolin/archive/0.0.1.zip`
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

 `conda create -n <new_env_name> -c bioconda -c conda-forge -c defaults -c pipolinlab explore-pipolin`

NOTE: solving the environment takes time. Be patient.

# Quick usage

As input, **ExplorePipolin** takes FASTA file(s) with genome sequence(s). 
Due to restrictions of some of the libraries used, filename must not exceed 16 characters long.
A genome sequence can be either a single complete chromosome (preferred) 
or contigs (in a single multiFASTA file).



### Command line options

```bash
Usage: explore_pipolin [OPTIONS] GENOMES...

  ExplorePipolin is a search tool for prediction and analysis of pipolins,
  bacterial mobile genetic elements.

Options:
  --out-dir-prefix TEXT       Use this prefix for the output directory,
                              instead of the default "results" prefix.

  --out-dir PATH              Use this output directory instead.
  --pipolb-hmm-profile PATH   piPolB's HMM profile to use as 1st priority.If
                              not provided, the default profile will be used
                              instead.

  --only-find-pipolbs         Only find piPolB genes.

  --ref-att PATH              ATT sequence in FASTA file to use as 1st
                              priority. If not provided, the default ATT will
                              be used instead.

  --percent-identity INTEGER  Minimum percent identity for direct repeats
                              search  [default: 85]

  --max-inflate INTEGER       If no borders of pipolin are found (no ATTs),
                              inflate the analysed region from both sides of
                              piPolB by this length.  [default: 30000]

  --skip-annotation           Do not run the annotation step (i.e. Prokka).
  --proteins PATH             Prokka param: FASTA or GBK file to use as 1st
                              priority. If not provided, the default file will
                              be used instead.

  --colours PATH               A TSV file describing features to colour. Please,
                              refer to
                              https://github.com/pipolinlab/ExplorePipolin for
                              more information about the file structure.

  --skip-colours              Do not add an Easyfig-compatible colouring
                              scheme to the final Genbank files.

  --cpus INTEGER              Prokka param: Number of CPUs to use [0=all]
                              [default: 8]

  --keep-tmp                  Preserve intermediate files produced during the
                              run (it might be useful for debugging).

  --keep-going                Do not stop analysis if an error is raised for
                              one of the genomes.

  -h, --help                  Show this message and exit.
```

### Output files

The output directory will contain a separate folder for each analysed genome. In the genome-specific folder,
the following files can be found:
 
 | Folder/File | Content description |
 |--------|---------------------|
 | `<genome>.log` | log file |
 | `pipolins` | extracted pipolin sequences in FASTA format, GenBank and GFF annotations of pipolin genes with *att*s included |

Files in `pipolins` folder have the following scheme: `<genome>_NvN.<type>.<ext>`,
where:

 * `<genome>` — genome name, retrieved from a provided FASTA file without 
   extension. Should not exceed 16 characters, due to Biopython restrictions!
 * `N` — index of a found pipolin
 * `vN` — index of a reconstructed variant (each pipolin can be reconstructed in
   different ways)

NOTE: When alternative orientations of the piPolB-containing fragment are possible, both reconstructions are generated. The structure denoted as `v0` will include piPolB gene pointing toward the attR (usually overlaps with tRNA). The alternative pipolin with the piPolB gene leftwards referred as `v1`.

 * `<type>` — can be `complete` (att---pol---att), `truncated` (att---pol---) or 
   `minimal` (---pol---)
 * `<ext>` — file extension (`.fa`, `.gff`, `.gbk`)

**NOTE**: A file with `single_record.gbk` ending is produced when separate contigs are 
concatenated into a single record (after the reconstruction step). Contigs are joined 
with an `assembly_gap` feature that can be distinguished by the presence of 
`/inference="ExplorePipolin"` field. Also, by default, these assembly gaps have 
black `/colour`, while all other assembly gaps are pink 
(see [Easyfig-compatible colouring](#easyfig-compatible-colouring)).

Additional files, when `--keep-tmp` option is present:

 | Folder/File | Content description |
 |-------|----------------------|
 | `<genome>.fa` | genome file copied |
 | `pipolbs` | HMMER search results for piPolB genes |
 | `atts` | BLAST search results for the known *att* sites |
 | `atts_denovo` | Results of *de novo* search for *att* sites |
 | `trnas` | ARAGORN search results for tRNAs/tmRNAs |
 | `prokka` | Prokka annotation results (check files description [here](https://github.com/tseemann/prokka/blob/master/README.md#output-files))|

### Easyfig-compatible colouring

By default, a `/colour` field with an RGB colour code will be added to each feature 
in the output Genbank files. The colours will be assigned according to 
[this scheme](explore_pipolin/data/colours.tsv). **To modify the colouring**, download 
and change this TSV file and provide it using `--colours` option.

An example command to run Easyfig:

```bash
# IMPORTANT: run with python2!
python <path>/Easyfig.py -o my_figure.svg -svg -f1 T \
-f CDS arrow -f repeat_region rect -f tRNA rect -f assembly_gap pointer \
-tblastx -min_length 500 -ann_height 100 \
single_record1.gbk single_record2.gbk ...
```

All options can be checked by running `python <path>/Easyfig.py --help`.

NOTE: Use `--skip-colours` if you do not need the colours.
