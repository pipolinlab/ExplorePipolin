import os
from collections import defaultdict
from typing import MutableMapping, MutableSequence, Tuple

from prefect import task

from explore_pipolin.common import Genome, Feature, Range, FeatureType, Strand, ContigID
from explore_pipolin.utilities.external_tools import run_aragorn
from explore_pipolin.utilities.logging import genome_specific_logging


@task()
@genome_specific_logging
def find_trnas(genome: Genome, out_dir, do_not_reuse) -> Genome:
    aragorn_results_dir = os.path.join(out_dir, 'trnas')
    os.makedirs(aragorn_results_dir, exist_ok=True)

    output_file = os.path.join(aragorn_results_dir, genome.id + '.batch')
    if do_not_reuse or not os.path.exists(output_file):
        run_aragorn(genome_file=genome.file, output_file=output_file)
    entries = read_aragorn_batch(aragorn_batch=output_file)

    add_trna_features(entries=entries, genome=genome)

    return genome


def read_aragorn_batch(aragorn_batch) -> MutableMapping[ContigID, MutableSequence[Tuple[int, int, Strand]]]:
    entries = defaultdict(list)
    with open(aragorn_batch) as inf:
        for line in inf:
            if line[0] == '>':
                # >FHEK01000001.1 Staphylococcus aureus strain st2898 ...
                entry = ContigID(line.strip().split(sep=' ')[0][1:])
            else:
                hit = line.split(sep='\t')
                if len(hit) > 1:
                    #                                           \t    \t
                    # 1   tRNA-Arg              c[175742,175814]	34  	(ccg)
                    # 2   tmRNA                  [193736,194094]	94,129	GKSNNNFAVAA*
                    coordinates = hit[0].split(sep=' ')[-1]
                    if coordinates[0] == 'c':
                        start, end = (int(i) for i in coordinates[2:-1].split(sep=','))
                        entries[entry].append((start, end, Strand.REVERSE))
                    else:
                        start, end = (int(i) for i in coordinates[1:-1].split(sep=','))
                        entries[entry].append((start, end, Strand.FORWARD))

    return entries


def add_trna_features(entries: MutableMapping[ContigID, MutableSequence[Tuple[int, int, Strand]]], genome: Genome):
    for contig_id, hits in entries.items():
        for hit in hits:

            # "correct strange coordinates in -l mode" as it's done in Prokka
            start = max(hit[0], 1)
            end = min(hit[1], genome.get_contig_by_id(contig_id=contig_id).length)

            trna_feature = Feature(location=Range(start=start, end=end), strand=hit[2],
                                   ftype=FeatureType.TRNA, contig_id=contig_id, genome=genome)
            genome.features.add_features(trna_feature)
