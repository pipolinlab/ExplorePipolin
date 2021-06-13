import os
from typing import Sequence

from Bio import SeqIO
from prefect import task

from explore_pipolin.common import Genome, PipolinFragment, PipolinVariants
from explore_pipolin.utilities.external_tools import run_prokka
from explore_pipolin.utilities.io import create_seqio_records_dict, SeqIORecords
from explore_pipolin.utilities.logging import genome_specific_logging


@task()
@genome_specific_logging
def save_pipolin_sequences(genome: Genome, pipolins: Sequence[PipolinVariants]):
    pipolins_dir = os.path.join(os.path.dirname(genome.file), 'pipolins')
    os.makedirs(pipolins_dir, exist_ok=True)

    genome_dict = create_seqio_records_dict(file=genome.file, file_format='fasta')

    for i, pipolin in enumerate(pipolins):
        for j, variant in enumerate(pipolin.variants):
            with open(os.path.join(pipolins_dir, genome.id + f'_{i}_v{j}.{pipolin.type.to_str()}.fa'), 'w') as ouf:
                records = [create_fragment_record(f, genome_dict) for f in variant.fragments]
                SeqIO.write(records, ouf, 'fasta')

    return pipolins_dir


def create_fragment_record(
        fragment: PipolinFragment,
        genome_dict: SeqIORecords
) -> SeqIO.SeqRecord:
    record = genome_dict[fragment.contig_id][fragment.start:fragment.end]
    record.description = f'{fragment.start}:{fragment.end}'
    return record


@task()
@genome_specific_logging
def annotate_pipolins(genome: Genome, pipolins_dir):
    prokka_dir = os.path.join(os.path.dirname(pipolins_dir), 'prokka')
    os.makedirs(prokka_dir, exist_ok=True)

    for f in os.listdir(pipolins_dir):
        if f.startswith(genome.id) and f.endswith('.fa'):
            input_file = os.path.join(pipolins_dir, f)
            run_prokka(input_file=input_file, prokka_results_dir=prokka_dir)

    return prokka_dir
