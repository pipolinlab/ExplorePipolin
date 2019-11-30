# ExplorePipolin

To extract the subsequence from a genome:
 * `get_subsequence.py`
 
 `$ get_subsequence.py genomes/NZ_JNMI01000006.1.fa 80191 82792 pi-polB.fa`
 
 `$ get_subsequence.py genomes/NZ_JNMI01000006.1.fa 64241 64373 attL.fa`

The whole analysis:
 * `download_genomes_ncbi.py`
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
 * `store_new_att_bounds.py` 
 (parse HMMER output for atts)
 * `include_atts_into_annotation.py` 
 TODO: So far, the script includes atts only into GB 
 and GFF files. 
 
Prediction of ATTs:
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
 5. `predict_atts_with_hmmer.py`

To get the sequences from roary groups:
 * `extract_roary_groups.py`
