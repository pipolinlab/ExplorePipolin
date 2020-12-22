import os
from collections import defaultdict
from enum import Enum, auto
from typing import MutableSequence, Optional, Sequence, Mapping, MutableMapping, List


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


class Orientation(Enum):
    FORWARD = auto()
    REVERSE = auto()

    @staticmethod
    def from_pm_one_encoding(hit_strand):
        if hit_strand == 1:
            return Orientation.FORWARD
        elif hit_strand == -1:
            return Orientation.REVERSE
        else:
            raise AssertionError(f'Unknown hit strand: {hit_strand}! Should be 1 or -1.')

    def to_pm_one_encoding(self):
        return 1 if self is self.FORWARD else -1

    def __neg__(self):
        return self.REVERSE if self is self.FORWARD else self.FORWARD


class Contig:
    def __init__(self, contig_id: str, contig_length: int, orientation=Orientation.FORWARD):
        self.contig_id = contig_id
        self.contig_length = contig_length
        self.contig_orientation: Orientation = orientation


class Genome:
    def __init__(self, genome_id: str, genome_file: str, contigs: MutableSequence[Contig]):
        self.genome_id = genome_id
        self.genome_file = genome_file
        self.contigs = contigs
        self.features: FeaturesContainer = FeaturesContainer()

    def get_contig_by_id(self, contig_id: str) -> Optional[Contig]:
        for contig in self.contigs:
            if contig.contig_id == contig_id:
                return contig
        raise AssertionError(f'It is unexpected that you ask for a non-existent contig {contig_id}!')

    def is_single_contig(self):
        return len(self.contigs) == 1

    def get_complete_genome_contig_id(self) -> str:
        if not self.is_single_contig():
            raise AssertionError('Unsupported! Not a complete genome!')
        return self.contigs[0].contig_id

    def get_complete_genome_length(self) -> int:
        if not self.is_single_contig():
            raise AssertionError('Unsupported! Not a complete genome!')
        return self.contigs[0].contig_length


class Feature:
    def __init__(self, start: int, end: int, strand: Orientation, contig_id: str, genome: Genome):
        self.start = start
        self.end = end
        self.strand = strand
        self.contig_id = contig_id
        self.genome = genome

        if self.start > self.end:
            raise AssertionError('Feature start cannot be greater than feature end!')

        if self.end > self.contig.contig_length:
            raise AssertionError('Feature end cannot be greater than contig length!')

    def is_overlapping(self, other) -> bool:
        max_start = max(self.start, other.start)
        min_end = min(self.end, other.end)
        return max_start <= min_end

    @property
    def contig(self) -> Contig:
        return self.genome.get_contig_by_id(self.contig_id)

    def shift(self, shift):
        return Feature(self.start + shift, self.end + shift, self.strand, self.contig_id, self.genome)


class FeatureType(Enum):
    PIPOLB = auto()
    ATT = auto()
    TRNA = auto()
    TARGET_TRNA = auto()
    ATT_DENOVO = auto()
    TARGET_TRNA_DENOVO = auto()


class FeaturesContainer:
    def __init__(self):
        self._features: Mapping[FeatureType, MutableSequence[Feature]] = defaultdict(list)

    @staticmethod
    def _dict_by_contig_normalized(features: Sequence[Feature]) -> Mapping[str, Sequence[Feature]]:
        result: MutableMapping[str, List[Feature]] = defaultdict(list)
        for feature in features:
            result[feature.contig_id].append(feature)
        for _, features in result.items():
            features.sort(key=lambda p: p.start)
        return result

    def get_features_dict_by_contig_normalized(self, feature_type: FeatureType) -> Mapping[str, Sequence[Feature]]:
        return self._dict_by_contig_normalized(self.get_features(feature_type))

    def get_features_of_contig_normalized(self, contig_id: str, feature_type: FeatureType) -> Sequence[Feature]:
        features_dict = self.get_features_dict_by_contig_normalized(feature_type=feature_type)
        return features_dict[contig_id]

    def get_features(self, feature_type: FeatureType) -> MutableSequence[Feature]:
        return self._features[feature_type]

    def add_feature(self, feature, feature_type):
        self.get_features(feature_type).append(feature)

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
            target_contigs.extend(f.contig_id for f in self.get_features(feature_type))
        return len(set(target_contigs)) == 1

    def is_overlapping_with(self, feature: Feature, feature_type: FeatureType):
        for other_feature in self.get_features(feature_type):
            if feature.is_overlapping(other_feature):
                return True
        return False


class RepeatPair:
    def __init__(self, left: Feature, right: Feature, left_seq: str, right_seq: str):
        self.left = left
        self.right = right
        self.left_seq = left_seq
        self.right_seq = right_seq

        left_repeat_length = self.left.end - self.left.start
        right_repeat_length = self.right.end - self.right.start
        left_seq_length = len(self.left_seq)
        right_seq_length = len(self.right_seq)
        if left_seq_length > left_repeat_length or right_seq_length > right_repeat_length:
            raise AssertionError(f'Repeat sequence length cannot be greater than repeat ranges:'
                                 f' {left_seq_length} > {left_repeat_length} or '
                                 f'{right_seq_length} > {right_repeat_length}!')

    def shift(self, left_shift: int, right_shift: int):
        if left_shift > right_shift:
            raise AssertionError('Left shift cannot be greater than right shift!')

        return RepeatPair(self.left.shift(left_shift), self.right.shift(right_shift), self.left_seq, self.right_seq)


class PipolinFragment:
    def __init__(self, contig_id: str, genome: Genome, start: int, end: int):
        self.start = start
        self.end = end
        self.atts: MutableSequence[Feature] = []
        self.contig_id = contig_id
        self.genome = genome

        if self.start > self.end:
            raise AssertionError('Fragment start cannot be greater than fragment end!')

        if self.end > self.contig.contig_length:
            raise AssertionError('Fragment end cannot be greater than contig length!')

    @property
    def contig(self):
        return self.genome.get_contig_by_id(self.contig_id)


class Pipolin:
    def __init__(self, *fragments: PipolinFragment):
        self.fragments: Sequence[PipolinFragment] = fragments


def define_genome_id(genome_path: str):
    genome_id = os.path.splitext(os.path.basename(genome_path))[0]
    check_genome_id_length(genome_id)
    return genome_id


def check_genome_id_length(genome_id):
    if len(genome_id) > 16:
        raise AssertionError('Genome file basename is going to be used as a genome identifier. '
                             'Due to Biopython restrictions, it cannot be longer than 16 characters. '
                             'Please, rename the file, so that its basename does not exceed the limit!')
