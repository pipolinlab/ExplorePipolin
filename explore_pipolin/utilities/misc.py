from __future__ import annotations

from collections import defaultdict
from typing import MutableSequence, Optional, Sequence, Mapping, MutableMapping, List

from explore_pipolin.common import Orientation, Contig, Genome, Feature, FeatureType


class GQuery:
    def __init__(self, genome: Genome):
        self.genome = genome
        #self.features: Mapping[FeatureType, MutableSequence[Feature]] = {}
        self.pipolbs: MutableSequence[Feature] = []
        self.atts: MutableSequence[Feature] = []
        self.trnas: MutableSequence[Feature] = []
        self.target_trnas: MutableSequence[Feature] = []
        self.denovo_atts: MutableSequence[Feature] = []
        self.target_trnas_denovo: MutableSequence[Feature] = []

    @staticmethod
    def _dict_by_contig_normalized(features: Sequence[Feature]) -> Mapping[str, Sequence[Feature]]:
        result: MutableMapping[str, List[Feature]] = defaultdict(list)
        for feature in features:
            result[feature.contig_id].append(feature)
        for _, features in result.items():
            features.sort(key=lambda p: p.start)
        return result

    def get_features_dict_by_contig_normalized(self, feature_type: FeatureType) -> Mapping[str, Sequence[Feature]]:
        return self._dict_by_contig_normalized(self.get_features_by_type(feature_type))

    def get_features_of_contig_normalized(self, contig_id: str, feature_type: FeatureType) -> Sequence[Feature]:
        features_dict = self.get_features_dict_by_contig_normalized(feature_type=feature_type)
        return features_dict[contig_id]

    def get_features_by_type(self, feature_type: FeatureType) -> MutableSequence[Feature]:
        if feature_type is FeatureType.PIPOLB:
            return self.pipolbs
        elif feature_type is FeatureType.ATT:
            return self.atts
        elif feature_type is FeatureType.TARGET_TRNA:
            return self.target_trnas
        elif feature_type is FeatureType.TRNA:
            return self.trnas
        else:
            raise AssertionError(f'Feature must be one of: {list(FeatureType)}, not {feature_type}')

    def find_overlapping_feature(self, feature: Feature, feature_type: FeatureType) -> Optional[Feature]:
        feature_dict = self.get_features_dict_by_contig_normalized(feature_type)

        if feature.contig_id in feature_dict:
            features = feature_dict[feature.contig_id]
            for other_feature in features:
                if other_feature.is_overlapping(feature):
                    return other_feature

    def is_on_the_same_contig(self, *feature_types: FeatureType) -> bool:
        target_contigs = []
        for feature_type in feature_types:
            target_contigs.extend(f.contig_id for f in self.get_features_by_type(feature_type))
        return len(set(target_contigs)) == 1

    def get_left_right_windows(self, feature_type) -> [Feature, Feature]:
        features = self.get_features_by_type(feature_type=feature_type)
        features = sorted(features, key=lambda x: x.start)

        if features[-1].start - features[0].start > 10000:   # TODO: This should be changed!
            raise AssertionError(f'You have several piPolBs per genome and they are too far from each other: '
                                 f'within the region ({features[0].start}, {features[-1].end}). It might be, '
                                 f'that you have two or more pipolins per genome, but we are expecting only one.')

        length = self.genome.get_complete_genome_length()
        left_edge = features[0].start - 100000
        left_window = Feature(start=left_edge if left_edge >= 0 else 0, end=features[0].start,
                              strand=Orientation.FORWARD,
                              contig_id=self.genome.get_complete_genome_contig_id(), genome=self.genome)
        right_edge = features[-1].end + 100000
        right_window = Feature(start=features[-1].end, end=right_edge if right_edge <= length else length,
                               strand=Orientation.FORWARD,
                               contig_id=self.genome.get_complete_genome_contig_id(), genome=self.genome)

        return left_window, right_window

    def is_overlapping_with(self, feature: Feature, feature_type: FeatureType):
        for other_feature in self.get_features_by_type(feature_type):
            if feature.is_overlapping(other_feature):
                return True
        return False


def is_single_target_trna_per_contig(gquery: GQuery):
    # TODO: don't like this
    # there was one case with two target trnas per genome, although usually only one
    target_trnas_dict = gquery.get_features_dict_by_contig_normalized(feature_type=FeatureType.TARGET_TRNA)
    target_trnas = gquery.get_features_by_type(feature_type=FeatureType.TARGET_TRNA)
    if len(target_trnas) != len(target_trnas_dict):
        raise AssertionError("We are expecting a single tRNA to overlap with a single att per contig!")


def get_contig_orientation(contig: Contig, gquery: GQuery) -> Orientation:
    target_trnas = gquery.get_features_of_contig_normalized(contig_id=contig.contig_id,
                                                            feature_type=FeatureType.TARGET_TRNA)
    atts = gquery.get_features_of_contig_normalized(contig_id=contig.contig_id, feature_type=FeatureType.ATT)
    atts_strands = [att.strand for att in atts]
    polbs = gquery.get_features_of_contig_normalized(contig_id=contig.contig_id, feature_type=FeatureType.PIPOLB)
    polbs_strands = [polb.strand for polb in polbs]

    if len(target_trnas) != 0:
        if len(set(atts_strands)) != 1:
            raise AssertionError('ATTs are expected to be in the same frame, as they are direct repeats!')
        if set(atts_strands).pop() == target_trnas[0].strand:
            raise AssertionError(f'ATT and tRNA are expected to be on the different strands! '
                                 f'{set(atts_strands).pop()}, {target_trnas[0].strand}, {gquery.genome.genome_id}')
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
    return Feature(start=hit.hit_start, end=hit.hit_end,
                   strand=Orientation.orientation_from_blast(hit.hit_strand),
                   contig_id=contig_id, genome=genome)
