import os
from typing import MutableSequence

from prefect import task

from explore_pipolin.common import Genome, Contig, ContigID
from explore_pipolin.utilities.io import create_seqio_records_dict


@task()
def create_genome(genome_file) -> Genome:
    contigs = create_contigs_for_genome(genome_file)
    genome_id = define_genome_id(genome_file)
    return Genome(genome_id, genome_file, contigs)


def create_contigs_for_genome(genome_file: str) -> MutableSequence[Contig]:
    genome_dict = create_seqio_records_dict(file=genome_file, file_format='fasta')
    _check_contig_ids_length(genome_dict)
    contigs = []
    for key, value in genome_dict.items():
        contigs.append(Contig(contig_id=ContigID(key), contig_length=len(value.seq)))
    return contigs


_max_length_allowed = 16


def _check_contig_ids_length(genome_dict):
    for key in genome_dict:
        if len(key) > _max_length_allowed:
            raise AssertionError(f'Due to Biopython restrictions, fasta identifiers cannot be longer than '
                                 f'{_max_length_allowed} characters. They should also be unique. '
                                 f'Please, rename identifiers to obey the restrictions!')


def define_genome_id(genome_path: str) -> str:
    genome_id = os.path.splitext(os.path.basename(genome_path))[0]
    _check_genome_id_length(genome_id)
    return genome_id


def _check_genome_id_length(genome_id: str) -> None:
    if len(genome_id) > _max_length_allowed:
        raise AssertionError('Genome file basename is going to be used as an identifier. '
                             f'Due to Biopython restrictions, it cannot be longer than {_max_length_allowed} '
                             f'characters. Please, rename the file, so that its basename does not exceed the limit!')
