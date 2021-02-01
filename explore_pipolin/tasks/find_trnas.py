import os
from collections import defaultdict
from typing import MutableMapping, MutableSequence

from prefect import task

from explore_pipolin.common import Genome, Feature, Range, FeatureType, Strand
from explore_pipolin.utilities.external_tools import run_aragorn
from explore_pipolin.utilities.logging import genome_specific_logging


@task()
@genome_specific_logging
def find_trnas(genome: Genome, out_dir) -> Genome:
    aragorn_results_dir = os.path.join(out_dir, 'trnas')
    os.makedirs(aragorn_results_dir, exist_ok=True)

    output_file = os.path.join(aragorn_results_dir, genome.id + '.batch')
    run_aragorn(genome_file=genome.file, output_file=output_file)
    entries = read_aragorn_batch(aragorn_batch=output_file)

    add_trna_features(entries=entries, genome=genome)

    return genome


def add_trna_features(entries, genome: Genome):
    for contig_id, hits in entries.items():
        for hit in hits:
            # "correct strange coordinates in -l mode" as in Prokka
            start = max(hit[0], 1)
            end = min(hit[1], genome.get_contig_by_id(contig_id=contig_id).length)
            trna_feature = Feature(location=Range(start=start, end=end),
                                   strand=hit[2], contig_id=contig_id, genome=genome)
            genome.features.add_features(trna_feature, feature_type=FeatureType.TRNA)


def read_aragorn_batch(aragorn_batch) -> MutableMapping[str, MutableSequence]:
    entries = defaultdict(list)
    with open(aragorn_batch) as inf:
        for line in inf:
            if line[0] == '>':
                entry = line.strip().split(sep=' ')[0][1:]
            else:
                hit = line.split(sep='\t')
                if len(hit) > 1:
                    coordinates = hit[0].split(sep=' ')[-1]
                    if coordinates[0] == 'c':
                        start, end = (int(i) for i in coordinates[2:-1].split(sep=','))
                        entries[entry].append((start, end, Strand.REVERSE))
                    else:
                        start, end = (int(i) for i in coordinates[1:-1].split(sep=','))
                        entries[entry].append((start, end, Strand.FORWARD))

    return entries
