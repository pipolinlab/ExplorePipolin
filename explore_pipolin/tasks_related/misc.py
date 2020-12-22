from typing import List, MutableSequence

from explore_pipolin.common import Orientation, Contig, Genome, Feature, FeatureType, Range


def is_single_target_trna_per_contig(genome: Genome):
    # TODO: don't like this
    # there was one case with two target trnas per genome, although usually only one
    target_trnas_dict = genome.features.get_features_dict_by_contig_normalized(feature_type=FeatureType.TARGET_TRNA)
    target_trnas = genome.features.get_features(feature_type=FeatureType.TARGET_TRNA)
    if len(target_trnas) != len(target_trnas_dict):
        raise AssertionError("We are expecting a single tRNA to overlap with a single att per contig!")


def get_contig_orientation(contig: Contig, genome: Genome) -> Orientation:
    target_trnas = genome.features.get_features_of_contig_normalized(contig_id=contig.contig_id,
                                                                     feature_type=FeatureType.TARGET_TRNA)
    atts = genome.features.get_features_of_contig_normalized(contig_id=contig.contig_id,
                                                             feature_type=FeatureType.ATT)
    atts_strands = [att.strand for att in atts]
    polbs = genome.features.get_features_of_contig_normalized(contig_id=contig.contig_id,
                                                              feature_type=FeatureType.PIPOLB)
    polbs_strands = [polb.strand for polb in polbs]

    if len(target_trnas) != 0:
        if len(set(atts_strands)) != 1:
            raise AssertionError('ATTs are expected to be in the same frame, as they are direct repeats!')
        if set(atts_strands).pop() == target_trnas[0].strand:
            raise AssertionError('ATT and tRNA are expected to be on the different strands!')
        return - target_trnas[0].strand

    elif len(atts) != 0:
        if len(set(atts_strands)) != 1:
            raise AssertionError('ATTs are expected to be in the same frame, as they are direct repeats!')
        return atts[0].strand

    if len(polbs) != 0:
        if len(set(polbs_strands)) != 1:  # an ambiguous case
            return contig.contig_orientation
        return polbs[0].strand


def join_it(iterable, delimiter):
    it = iter(iterable)
    yield next(it)
    for x in it:
        yield delimiter
        yield x


def create_fragment_record(fragment, genome_dict):
    fragment_record = genome_dict[fragment.contig.contig_id][fragment.start:fragment.end]
    if fragment.contig.contig_orientation == Orientation.REVERSE:
        fragment_record = fragment_record.reverse_complement()
    return fragment_record


def feature_from_blasthit(hit, contig_id: str, genome: Genome) -> Feature:
    return Feature(coords=Range(start=hit.hit_start, end=hit.hit_end),
                   strand=Orientation.from_pm_one_encoding(hit.hit_strand),
                   contig_id=contig_id, genome=genome)


class Window:
    def __init__(self, left: Range, right: Range, pipolbs: MutableSequence[Feature]):
        self.left = left
        self.right = right
        self.pipolbs = pipolbs


def get_windows(genome: Genome) -> List[Window]:
    pipolbs = genome.features.get_features(feature_type=FeatureType.PIPOLB)
    pipolbs = sorted(pipolbs, key=lambda x: x.start)

    windows = []
    for pipolb in pipolbs:   # TODO: actually, I need all the combinations here!
        contig_length = genome.get_contig_by_id(contig_id=pipolb.contig_id).contig_length

        left_edge = pipolb.coords.start - 100000
        left_window = Range(start=left_edge if left_edge >= 0 else 0, end=pipolbs[0].coords.start)
        right_edge = pipolb.coords.end + 100000
        right_window = Range(start=pipolbs[-1].coords.end,
                             end=right_edge if right_edge <= contig_length else contig_length)
        windows.append(Window(left_window, right_window, pipolbs=[pipolb]))
    return windows


def add_features_from_blast_entries(entries, feature_type: FeatureType, genome: Genome,):
    for entry in entries:
        for hit in entry:
            feature = feature_from_blasthit(hit=hit, contig_id=entry.id, genome=genome)
            genome.features.add_feature(feature=feature, feature_type=feature_type)
