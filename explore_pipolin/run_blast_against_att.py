import os
from prefect import task
from explore_pipolin.utilities import blast_genome_against_seq


@task
def run_blast_against_att(genome, root_dir, reference):
    blast_path = os.path.join(root_dir, 'att_blast')
    blast_genome_against_seq(genome=genome, seq=reference, seq_type='nucleotide', output_dir=blast_path)
    return blast_path
