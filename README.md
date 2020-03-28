# ExplorePipolin

To extract the subsequence from a genome:
 * `get_subsequence.py`
 
 `$ get_subsequence.py genomes/NZ_JNMI01000006.1.fa 80191 82792 pi-polB.fa`
 
 `$ get_subsequence.py genomes/NZ_JNMI01000006.1.fa 64241 64373 attL.fa`

 `$ get_subsequence.py genome/NZ_JNMI010000006.1.fa 90094 90008 tRNA.fa`

The main modules:
 * `download_metadata_ncbi.py` -- downloads the metadata for the analysed 
 genomes, such as accessions, organism and strain names
 * `download_genomes_ncbi.py` -- downloads genome (chromosome) sequences 
 given NCBI assembly accession (i.e. for a not-complete genome, it downloads
  all its contigs)
 * For Saskia's strains: `edit_contig_names.sh <in-dir>` --
 to shorten the long contig names
 * `identify_pipolins_roughly.py` using reference 
 `att.fa` and `pi-polB.fa` sequences
 (saves "pipolins" object into shelve.db)
 * `analyse_pipolin_orientation` (saves "orientations" 
 object into shelve.db) 
 * `extract_pipolin_regions.py` It is possible to 
 extract regions, determined by the leftmost and 
 rightmost atts (`--long` option).
 * `annotate_pipolins.py`
 * Here the steps to predict ATTs (see below)
 * `predict_atts_with_hmmer.py`
 * `store_new_att_bounds.py` 
 (parse HMMER output for atts)
 * `include_atts_into_annotation.py` 
 TODO: So far, the script includes atts only into GB 
 and GFF files.
 * `scaffold_gapped_pipolins.py`
 * `easyfig_add_colours.py`
 
Prediction of ATTs:

TODO: write a script for de-novo search of direct repeats!!!

 1. Prepared ATT sequences with `prepare_atts_for_msa.py`
 ```
The total number of atts is 198
> Maximum att length is 132
> Minimum att length is 121
```
 2. Built a MSA with att sequences using MAFFT, MUSCLE
 and T-Coffee (https://www.ebi.ac.uk/Tools/msa). 
 The output format -- Pearson FASTA, otherwise long
 sequence names might be truncated.
 3. Compared the alignments. Modified them, using 
 Jalview: deleted not conserved regions from both ends.
 4. Created HMM profile with `hmmbuild` and `hmmpress`.

To get the sequences from roary groups:
 * `extract_roary_groups.py`
