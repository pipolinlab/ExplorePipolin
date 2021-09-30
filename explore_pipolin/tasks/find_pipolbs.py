import os
from typing import Sequence, Tuple, Any, Optional

from Bio import SearchIO
from explore_pipolin.utilities.io import create_seqio_records_dict, write_seqio_records
from prefect import task, context

from explore_pipolin.common import Genome, Feature, Range, Strand, ContigID, FeatureType
from explore_pipolin.utilities.external_tools import run_prodigal, run_hmmsearch
from explore_pipolin.utilities.logging import genome_specific_logging


@task()
@genome_specific_logging
def find_pipolbs(genome: Genome) -> Genome:
    pipolbs_dir = os.path.join(os.path.dirname(genome.file), 'pipolbs')
    os.makedirs(pipolbs_dir, exist_ok=True)

    predicted_proteins = os.path.join(pipolbs_dir, genome.id + '.faa')
    hmm_output_file = os.path.join(pipolbs_dir, genome.id + '.tbl')

    run_prodigal(genome_file=genome.file, output_file=predicted_proteins)
    if os.stat(predicted_proteins).st_size == 0:   # no genes were found
        return genome

    run_hmmsearch(predicted_proteins, hmm_output_file)
    hmm_hits = read_hmm_hits(hmmsearch_table=hmm_output_file)

    entries = create_pipolb_entries(hmm_hits)
    add_pipolb_features(entries, genome)

    save_found_pipolbs(hmm_hits, predicted_proteins, genome)

    return genome


def read_hmm_hits(hmmsearch_table: str) -> Sequence:
    with open(hmmsearch_table) as inf:
        content = list(SearchIO.parse(inf, 'hmmer3-tab'))
        if len(content) == 0:
            return []
        elif len(content) != 1:
            raise AssertionError(f'More than a single query in {hmmsearch_table}! Should be only one.')
        else:
            return content[0]


def create_pipolb_entries(hits: Sequence) -> Sequence[Tuple[str, int, int, int]]:
    entries = []
    for hit in hits:
        name = '_'.join(hit.id.split(sep='_')[:-1])
        description = hit.description.split(sep=' ')
        entries.append((name, int(description[1]), int(description[3]), int(description[5])))

    return entries


def add_pipolb_features(entries, genome):
    for entry in entries:
        pipolb_feature = Feature(location=Range(start=entry[1], end=entry[2]),
                                 strand=Strand.from_pm_one_encoding(entry[3]),
                                 ftype=FeatureType.PIPOLB,
                                 contig_id=ContigID(entry[0]), genome=genome)
        genome.features.add_features(pipolb_feature)


def save_found_pipolbs(hmm_hits, proteins, genome: Genome):
    pipolbs_file = os.path.join(os.path.dirname(genome.file), genome.id + '_piPolBs.faa')

    if len(hmm_hits) > 0:
        proteins_dict = create_seqio_records_dict(proteins, 'fasta')
        pipolbs = dict()
        for hit in hmm_hits:
            pipolbs[hit.id] = proteins_dict[hit.id]
            pipolbs[hit.id].id = genome.id + ' ' + '_'.join(hit.id.split(sep='_')[:-1])
        write_seqio_records(pipolbs, pipolbs_file, 'fasta')


@task()
@genome_specific_logging
def are_pipolbs_present(genome: Genome):
    logger = context.get('logger')

    if len(genome.features.get_features(FeatureType.PIPOLB)) == 0:
        logger.warning('No piPolBs were found!')
        return False

    return True


@task()
def continue_if_true_else_finished(result_to_filter: Any, filter_1: bool, filter_2: bool) -> Optional[Any]:
    if filter_1 and not filter_2:
        return result_to_filter

    return None
