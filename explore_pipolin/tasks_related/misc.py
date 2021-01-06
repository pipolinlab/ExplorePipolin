from typing import List

from explore_pipolin.common import Strand, Genome, Feature, FeatureType, Range, RangePair, PipolinFragment


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


def feature_from_blasthit(hit, contig_id: str, genome: Genome) -> Feature:
    return Feature(location=Range(start=hit.hit_start, end=hit.hit_end),
                   strand=Strand.from_pm_one_encoding(hit.hit_strand),
                   contig_id=contig_id, genome=genome)


def get_ranges_around_pipolbs(genome: Genome) -> List[RangePair]:
    pipolbs_dict_by_contig = genome.features.pipolbs_dict()

    range_pairs = []
    for contig_id, pipolbs in pipolbs_dict_by_contig.items():
        contig_length = genome.get_contig_by_id(contig_id=contig_id).length

        for i in range(0, len(pipolbs)):
            for j in range(i, len(pipolbs)):
                pipolbs_range = Range(pipolbs[i].start, pipolbs[j].end)
                pipolbs_range = pipolbs_range.inflate(100000, _max=contig_length)
                range_pairs.append(RangePair(Range(pipolbs_range.start, pipolbs[i].start),
                                             Range(pipolbs[j].end, pipolbs_range.end), contig_id))
    return range_pairs


def add_features_from_blast_entries(entries, feature_type: FeatureType, genome: Genome,):
    for entry in entries:
        for hit in entry:
            feature = feature_from_blasthit(hit=hit, contig_id=entry.id, genome=genome)
            genome.features.add_feature(feature=feature, feature_type=feature_type)
