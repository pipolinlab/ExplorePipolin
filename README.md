# ExplorePipolin

NCBI strains analysis:
 * `download_genomes_ncbi.py`
 * `extract_ref_att_seq.py`
 * `identify_pipolins.py`
 * `extract_pipolin_regions.py`
 * `annotate_pipolins.py`
 
New strains analysis:
 * `extract_ref_att_seq.py`
 * `identify_new_pipolins.py`
 * `extract_new_pipolin_regions.py`
 * `annotate_pipolins.py`

To get the sequences from roary groups:
 * `extract_roary_groups.py`
 
Atts:
 1. Prepared att sequences with `prepare_atts_for_msa.py`
 2. Built a MSA with att sequences using 
 [MAFFT](https://www.ebi.ac.uk/Tools/msa/mafft/). 
 The output format -- Pearson FASTA, otherwise some
 sequence names will be truncated.
 It looked like there were two different types of atts.
 3. Clustered the atts using the obtained MSA: `cluster_atts.py`
 4. Built a separate MSA for each type of atts (MAFFT).
 5. ...
 