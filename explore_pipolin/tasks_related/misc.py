from typing import List

from explore_pipolin.common import Genome, FeatureType, Range, PairedLocation, PipolinFragment, ContigID


def is_single_target_trna_per_contig(genome: Genome):
    # TODO: don't like this
    # there was one case with two target trnas per genome, although usually only one
    target_trnas_dict = genome.features.target_trnas_dict()
    target_trnas = genome.features.get_features(feature_type=FeatureType.TARGET_TRNA)
    if len(target_trnas) != len(target_trnas_dict):
        raise AssertionError("We are expecting a single tRNA to overlap with a single att per contig!")


def join_it(iterable, delimiter):
    it = iter(iterable)
    yield next(it)
    for x in it:
        yield delimiter
        yield x


def create_fragment_record(fragment: PipolinFragment, genome_dict):
    return genome_dict[fragment.contig_id][fragment.start:fragment.end]


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
