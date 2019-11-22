# ExplorePipolin

To extract the subsequence from a genome 
(for example, ATT region):
 * `get_subsequence.py`

NCBI strains analysis:
 * `download_genomes_ncbi.py`
 * `identify_pipolins.py`
 * `extract_pipolin_regions.py`
 * `annotate_pipolins.py`
 
New strains analysis:
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
 to reverse-complement the sequences in one of the clusters and 
 save as a new FASTA file.
 4. Built a new MSA with unified atts (MAFFT), 
 saved as Pearson FASTA file.
 5. Modified with Jalview: deleted not conserved regions 
 from both ends.
 6. ...
 