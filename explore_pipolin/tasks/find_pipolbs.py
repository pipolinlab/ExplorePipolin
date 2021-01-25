import os
from typing import Sequence, Tuple, Any, Optional

from Bio import SearchIO
from prefect import task, context

from explore_pipolin.common import Genome, Feature, Range, Strand, ContigID, FeatureType
from explore_pipolin.utilities.external_tools import run_prodigal, run_hmmsearch
from explore_pipolin.utilities.logging import genome_specific_logging


@task()
@genome_specific_logging
def find_pipolbs(genome: Genome) -> Genome:
    results_dir = os.path.join(genome.results_dir, 'pipolbs')
    os.makedirs(results_dir, exist_ok=True)

    prodigal_output_file = os.path.join(results_dir, genome.id + '.faa')
    run_prodigal(genome_file=genome.file, output_file=prodigal_output_file)
    hmm_output_file = os.path.join(results_dir, genome.id + '.tbl')
    run_hmmsearch(prodigal_output_file, hmm_output_file)
    entries_pipolb = create_pipolb_entries(hmmsearch_table=hmm_output_file)

    for entry in entries_pipolb:
        feature = Feature(location=Range(start=entry[1], end=entry[2]),
                          strand=Strand.from_pm_one_encoding(entry[3]),
                          contig_id=ContigID(entry[0]), genome=genome)
        genome.features.add_features(feature, feature_type=FeatureType.PIPOLB)

    return genome


def create_pipolb_entries(hmmsearch_table: str) -> Sequence[Tuple[str, int, int, int]]:
    entries = []
    with open(hmmsearch_table) as inf:
        content = list(SearchIO.parse(inf, 'hmmer3-tab'))
        if len(content) != 1:
            raise AssertionError(f'More than a single query in {hmmsearch_table}! Should be only one.')
        for hit in content[0]:
            name = '_'.join(hit.id.split(sep='_')[:-1])
            description = hit.description.split(sep=' ')
            entries.append((name, int(description[1]), int(description[3]), int(description[5])))

    return entries


@task()
@genome_specific_logging
def are_pipolbs_present(genome: Genome):
    logger = context.get('logger')

    if len(genome.features.get_features(FeatureType.PIPOLB)) == 0:
        logger.warning('No piPolBs were found!')
        return False

    return True


@task()
def return_result_if_true_else_none(result_to_filter: Any, filter_by: bool) -> Optional[Any]:
    if filter_by:
        return result_to_filter

    return None
