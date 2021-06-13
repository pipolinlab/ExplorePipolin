import os
from typing import MutableSequence

from prefect import task

import explore_pipolin.settings as settings
from explore_pipolin.common import Genome, Contig, ContigID
from explore_pipolin.utilities.io import create_seqio_records_dict, SeqIORecords, write_seqio_records


@task()
def prepare_for_the_analysis(original_file: str) -> Genome:
    genome_dict = create_seqio_records_dict(file=original_file, file_format='fasta')
    genome_dict = make_contig_ids_shorter_if_long(genome_dict)

    contigs = create_contigs_for_genome(genome_dict)
    genome_id = define_genome_id(original_file)

    # create genome-specific output directory and put genome fasta file with short headers there
    genome_dir = os.path.join(settings.get_instance().out_dir, genome_id)
    os.makedirs(genome_dir, exist_ok=True)
    internal_use_file = os.path.join(genome_dir, os.path.splitext(os.path.basename(original_file))[0] + '.fa')
    write_seqio_records(genome_dict, internal_use_file, 'fasta')

    return Genome(genome_id, internal_use_file, contigs)


_MAX_LENGTH_ALLOWED = 16


def make_contig_ids_shorter_if_long(genome_dict):
    new_dict = {}
    for key, value in genome_dict.items():
        new_contig_id = key.split(sep='|')[-1]
        if len(new_contig_id) > _MAX_LENGTH_ALLOWED:
            new_contig_id = new_contig_id[-_MAX_LENGTH_ALLOWED:]
        value.id = new_contig_id
        new_dict[value.id] = value
    return new_dict


def create_contigs_for_genome(genome_dict: SeqIORecords) -> MutableSequence[Contig]:
    contigs = []
    for key, value in genome_dict.items():
        contigs.append(Contig(contig_id=ContigID(key), contig_length=len(value.seq)))
    return contigs


def define_genome_id(genome_path: str) -> str:
    genome_id = os.path.splitext(os.path.basename(genome_path))[0]
    _check_genome_id_length(genome_id)
    return genome_id


def _check_genome_id_length(genome_id: str) -> None:
    if len(genome_id) > _MAX_LENGTH_ALLOWED:
        raise AssertionError('Genome file basename is going to be used as an identifier. '
                             f'Due to Biopython restrictions, it cannot be longer than {_MAX_LENGTH_ALLOWED} '
                             f'characters. Please, rename the file, so that its basename does not exceed the limit!')
