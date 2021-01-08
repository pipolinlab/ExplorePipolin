import os
from collections import defaultdict
from dataclasses import dataclass
from enum import Enum, auto
from typing import MutableSequence, Optional, Sequence, Mapping, MutableMapping, List, Set, NewType

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


class Strand(Enum):
    FORWARD = auto()
    REVERSE = auto()

    @staticmethod
    def from_pm_one_encoding(hit_strand):
        if hit_strand == 1:
            return Strand.FORWARD
        elif hit_strand == -1:
            return Strand.REVERSE
        else:
            raise AssertionError(f'Unknown hit strand: {hit_strand}! Should be 1 or -1.')

    def to_pm_one_encoding(self):
        return 1 if self is self.FORWARD else -1

    def __neg__(self):
        return self.REVERSE if self is self.FORWARD else self.FORWARD


class Contig:
    def __init__(self, contig_id: str, contig_length: int):
        self.id = contig_id
        self.length = contig_length


class Genome:
    def __init__(self, genome_id: str, genome_file: str, contigs: MutableSequence[Contig]):
        self.id = genome_id
        self.file = genome_file
        self.contigs = contigs
        self.features: FeaturesContainer = FeaturesContainer()

    def get_contig_by_id(self, contig_id: str) -> Contig:
        for contig in self.contigs:
            if contig.id == contig_id:
                return contig
        raise AssertionError(f'Asking for the non-existent contig name {contig_id}!')

    def is_single_contig(self) -> bool:
        return len(self.contigs) == 1


@dataclass(frozen=True)
class Range:
    start: int
    end: int

    def __post_init__(self):
        if self.start > self.end:
            raise AssertionError('Start cannot be greater than end!')
        if self.start < 0:
            raise AssertionError('Start cannot be less than 0!')

    def shift(self, size: int):
        return Range(self.start + size, self.end + size)

    def inflate(self, size: int, _min: int = 0, _max: int = None):
        return Range(max(_min, self.start - size), self.end + size if _max is None else min(_max, self.end + size))

    def clamp(self, min_coord: int, max_coord: int):
        return Range(start=max(min_coord, self.start), end=min(max_coord, self.end))

    def is_overlapping(self, other) -> bool:
        max_start = max(self.start, other.start)
        min_end = min(self.end, other.end)
        return max_start <= min_end

    def is_overlapping_any(self, others) -> bool:
        return any(self.is_overlapping(r) for r in others)


ContigID = NewType('ContigID', str)


@dataclass(frozen=True)
class PairedLocation:
    left: Range
    right: Range
    contig_id: ContigID


@dataclass(frozen=True)
class MultiLocation:
    ranges: MutableSequence[Range]
    contig_id: ContigID


class Feature:
    def __init__(self, location: Range, strand: Strand, contig_id: str, genome: Genome):
        self.location = location
        self.strand = strand
        self.contig_id = contig_id
        self.genome = genome

        if self.location.end > self.contig.length:
            raise AssertionError(f'Feature end cannot be greater than contig length! '
                                 f'{self.location.end} > {self.contig.length}')

    @property
    def start(self) -> int:
        return self.location.start

    @property
    def end(self) -> int:
        return self.location.end

    @property
    # TODO: do I really need it?
    def contig(self) -> Contig:
        return self.genome.get_contig_by_id(self.contig_id)


class FeatureType(Enum):
    PIPOLB = auto()
    ATT = auto()
    TRNA = auto()
    TARGET_TRNA = auto()
    ATT_DENOVO = auto()
    TARGET_TRNA_DENOVO = auto()


class FeatureSet(Set[Feature]):
    def get_dict_by_contig_sorted(self) -> Mapping[str, Sequence[Feature]]:
        result: MutableMapping[str, List[Feature]] = defaultdict(list)
        for feature in self:
            result[feature.contig_id].append(feature)
        for _, features in result.items():
            features.sort(key=lambda p: p.start)
        return result

    def get_list_of_contig_sorted(self, contig_id: str) -> Sequence[Feature]:
        return self.get_dict_by_contig_sorted()[contig_id]

    @property
    def first(self) -> Feature:
        return next(iter(self))

    # TODO: do I really need it?
    def get_overlapping(self, feature: Feature) -> Optional[Feature]:
        try:
            features_list = self.get_list_of_contig_sorted(feature.contig_id)
            return self._get_overlapping_with_feature(features_list, feature)
        except KeyError:
            return None

    @staticmethod
    def _get_overlapping_with_feature(features_list, feature):
        for other_feature in features_list:
            if feature.location.is_overlapping(other_feature.location):
                return other_feature


class FeaturesContainer:
    def __init__(self):
        self._features: Mapping[FeatureType, FeatureSet] = defaultdict(FeatureSet)

    def add_features(self, *features: Feature, feature_type: FeatureType):
        if feature_type is FeatureType.ATT:
            for feature in features:
                if not isinstance(feature, AttFeature):
                    raise AssertionError('ATT should be instance of AttFeature class!')

        for feature in features:
            self._features[feature_type].add(feature)

    def get_features(self, feature_type: FeatureType) -> FeatureSet:
        return self._features[feature_type]

    def pipolbs_dict(self) -> Mapping[str, Sequence[Feature]]:
        return self.get_features(FeatureType.PIPOLB).get_dict_by_contig_sorted()

    def trnas_dict(self) -> Mapping[str, Sequence[Feature]]:
        return self.get_features(FeatureType.TRNA).get_dict_by_contig_sorted()

    def atts_dict(self):
        return self.get_features(FeatureType.ATT).get_dict_by_contig_sorted()

    def target_trnas_dict(self):
        return self.get_features(FeatureType.TARGET_TRNA).get_dict_by_contig_sorted()

    def is_on_the_same_contig(self, *feature_types: FeatureType) -> bool:
        target_contigs = []
        for feature_type in feature_types:
            target_contigs.extend(f.contig_id for f in self.get_features(feature_type))
        return len(set(target_contigs)) == 1


class AttFeature(Feature):
    def __init__(self, location: Range, strand: Strand, contig_id: str, genome: Genome, repeat_id):
        Feature.__init__(self, location, strand, contig_id, genome)
        self.repeat_id = repeat_id


class PipolinFragment(Feature):
    def __init__(self, location: Range, contig_id: str, genome: Genome):
        # TODO: do I need orientation here?
        super(PipolinFragment, self).__init__(location, Strand.FORWARD, contig_id, genome)
        self.atts: MutableSequence[Feature] = []


class Pipolin:
    def __init__(self, *fragments: PipolinFragment):
        self.fragments: Sequence[PipolinFragment] = fragments


def define_genome_id(genome_path: str) -> str:
    genome_id = os.path.splitext(os.path.basename(genome_path))[0]
    _check_genome_id_length(genome_id)
    return genome_id


def _check_genome_id_length(genome_id: str) -> None:
    max_length_allowed = 16
    if len(genome_id) > max_length_allowed:
        raise AssertionError('Genome file basename is going to be used as a genome identifier. '
                             f'Due to Biopython restrictions, it cannot be longer than {max_length_allowed} '
                             f'characters. Please, rename the file, so that its basename does not exceed the limit!')
