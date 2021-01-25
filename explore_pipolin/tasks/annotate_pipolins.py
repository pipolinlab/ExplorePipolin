import os
from typing import Sequence

from Bio import SeqIO
from prefect import task

from explore_pipolin.common import Genome, Pipolin, PipolinFragment
from explore_pipolin.utilities.external_tools import run_prokka
from explore_pipolin.utilities.io import create_seqio_records_dict
from explore_pipolin.utilities.logging import genome_specific_logging


@task()
@genome_specific_logging
def save_pipolin_sequences(genome: Genome, pipolins: Sequence[Pipolin]):
    pipolins_dir = os.path.join(genome.results_dir, 'pipolins')
    os.makedirs(pipolins_dir, exist_ok=True)

    genome_dict = create_seqio_records_dict(file=genome.file, file_format='fasta')

    for i, pipolin in enumerate(pipolins):
        with open(os.path.join(pipolins_dir, genome.id + f'_{i}.fa'), 'w') as ouf:
            records = [create_fragment_record(f, genome_dict) for f in pipolin.fragments]
            SeqIO.write(records, ouf, 'fasta')

    return pipolins_dir


def create_fragment_record(fragment: PipolinFragment, genome_dict):
    return genome_dict[fragment.contig_id][fragment.start:fragment.end]


@task()
@genome_specific_logging
def annotate_pipolins(genome: Genome, pipolins_dir):
    prokka_results_dir = os.path.join(genome.results_dir, 'prokka')
    os.makedirs(prokka_results_dir, exist_ok=True)

    for fasta_file in os.listdir(pipolins_dir):
        if fasta_file.startswith(genome.id):
            prefix = os.path.splitext(fasta_file)[0]
            input_file = os.path.join(pipolins_dir, fasta_file)
            run_prokka(prefix=prefix, input_file=input_file, prokka_results_dir=prokka_results_dir)

    return prokka_results_dir
