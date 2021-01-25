import os
from typing import List

from prefect import task

from explore_pipolin.common import Genome, Range, PairedLocation, ContigID, define_genome_id
from explore_pipolin.utilities.io import read_genome_contigs_from_file


@task()
def create_genome(genome_file, output_dir) -> Genome:
    contigs = read_genome_contigs_from_file(genome_file)
    genome_id = define_genome_id(genome_file)
    results_dir = os.path.join(output_dir, genome_id)
    print(results_dir)
    return Genome(genome_id, genome_file, results_dir, contigs)


def get_ranges_around_pipolbs(genome: Genome) -> List[PairedLocation]:
    pipolbs_dict_by_contig = genome.features.pipolbs_dict()

    range_pairs = []
    for contig_id, pipolbs in pipolbs_dict_by_contig.items():
        contig_length = genome.get_contig_by_id(contig_id=contig_id).length

        for i in range(0, len(pipolbs)):
            for j in range(i, len(pipolbs)):
                pipolbs_range = Range(pipolbs[i].start, pipolbs[j].end)
                pipolbs_range = pipolbs_range.inflate(100000, _max=contig_length)
                range_pairs.append(PairedLocation(Range(pipolbs_range.start, pipolbs[i].start),
                                                  Range(pipolbs[j].end, pipolbs_range.end), ContigID(contig_id)))
    return range_pairs
