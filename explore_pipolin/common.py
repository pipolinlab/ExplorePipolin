from collections import defaultdict
from dataclasses import dataclass
from enum import Enum, auto
from typing import MutableSequence, Optional, Sequence, Mapping, MutableMapping, List, Set, NewType, Tuple

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


ContigID = NewType('ContigID', str)


class Contig:
    def __init__(self, contig_id: ContigID, contig_length: int):
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

    def inflate_within_contig(self, size: int, _min: int = 0, _contig_length: int = None):
        start = max(_min, self.start - size)
        end = (self.end + size) if _contig_length is None else min(_contig_length, self.end + size)
        return Range(start, end)

    def is_overlapping(self, other) -> bool:
        max_start = max(self.start, other.start)
        min_end = min(self.end, other.end)
        return max_start <= min_end

    def is_overlapping_any(self, others) -> bool:
        return any(self.is_overlapping(r) for r in others)


@dataclass(frozen=True)
class PairedLocation:
    left_range: Range
    right_range: Range
    contig_id: ContigID


@dataclass(frozen=True)
class MultiLocation:
    ranges: MutableSequence[Range]
    contig_id: ContigID


class FeatureType(Enum):
    PIPOLB = auto()
    ATT = auto()
    TRNA = auto()
    TARGET_TRNA = auto()


@dataclass(frozen=True)
class Feature:
    location: Range
    strand: Strand
    ftype: FeatureType
    contig_id: ContigID
    genome: Genome

    def __post_init__(self):
        if self.location.end > self.contig.length:
            raise AssertionError(f'Feature end cannot be greater than contig length! '
                                 f'{self.location.end} > {self.contig.length}')

        if self.ftype is FeatureType.ATT:
            if not isinstance(self, AttFeature):
                raise AssertionError('ATT should be instance of AttFeature class!')

    @property
    def start(self) -> int:
        return self.location.start

    @property
    def end(self) -> int:
        return self.location.end

    def is_right_of(self, other) -> bool:
        return other.end < self.start

    def is_left_of(self, other) -> bool:
        return other.is_right_of(self)

    @property
    # TODO: do I really need it?
    def contig(self) -> Contig:
        return self.genome.get_contig_by_id(self.contig_id)

    def get_att_overlapping_ttrna(self):
        assert self.ftype == FeatureType.TARGET_TRNA
        att = self.genome.features.get_features(FeatureType.ATT).get_overlapping(self)
        if not att:
            raise AssertionError
        return att


@dataclass(frozen=True)
class AttFeature(Feature):
    att_id: int


class FeatureSet(Set[Feature]):
    def get_dict_by_contig_sorted(self) -> Mapping[str, MutableSequence[Feature]]:
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

    def get_next_att_id(self) -> int:
        att_ids = [i.att_id for i in self if isinstance(i, AttFeature)]
        return max(att_ids, default=0) + 1

    def get_atts_dict_by_att_id(self: Set[AttFeature]) -> Mapping[int, MutableSequence[AttFeature]]:
        atts_by_att_id = defaultdict(list)
        for att in self:
            atts_by_att_id[att.att_id].append(att)
        return atts_by_att_id

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

    def add_features(self, *features: Feature):
        for feature in features:
            self._features[feature.ftype].add(feature)

    def get_features(self, feature_type: FeatureType) -> FeatureSet:
        return self._features[feature_type]

    def pipolbs_dict(self) -> Mapping[str, MutableSequence[Feature]]:
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


@dataclass(frozen=True)
class PipolinFragment:
    location: Range
    contig_id: ContigID
    genome: Genome

    features: Tuple[Feature, ...] = ()

    orientation: Strand = Strand.FORWARD

    @property
    def start(self):
        return self.location.start

    @property
    def end(self):
        return self.location.end

    def reverse_complement(self):
        return PipolinFragment(self.location, self.contig_id, self.genome, self.features, -self.orientation)

    def is_overlapping(self, other):
        if self.contig_id == other.contig_id:
            return self.location.is_overlapping(other.location)

    def get_fragment_features_sorted(self) -> Sequence[Feature]:
        features = []
        features.extend(self.get_fragment_features_of_type_sorted(FeatureType.PIPOLB))
        features.extend(self.get_fragment_features_of_type_sorted(FeatureType.ATT))
        features.extend(self.get_fragment_features_of_type_sorted(FeatureType.TARGET_TRNA))

        return sorted(features, key=lambda x: x.start)

    def get_fragment_features_of_type_sorted(
            self, feature_type: FeatureType) -> Sequence[Feature]:

        features = self.genome.features.get_features(feature_type)
        contig_features = features.get_dict_by_contig_sorted()[self.contig_id]
        return [f for f in contig_features if f.start >= self.start and f.end <= self.end]

    def get_prime3_ttrnas(self) -> Sequence[Feature]:
        ttrnas = self.get_fragment_features_of_type_sorted(FeatureType.TARGET_TRNA)
        prime3_ttrnas = self._get_prime3_ttrnas(ttrnas)

        if len(prime3_ttrnas) != 0 or self._is_same_direction(prime3_ttrnas):
            return prime3_ttrnas

        return []

    @staticmethod
    def _get_prime3_ttrnas(ttrnas: Sequence[Feature]) -> Sequence[Feature]:
        prime3_ttrnas = []
        for ttrna in ttrnas:
            att = ttrna.get_att_overlapping_ttrna()

            if ttrna.strand == Strand.FORWARD and ttrnas[0].start < att.start:  # 3'-overlap
                prime3_ttrnas.append(ttrna)
            elif ttrnas[0].strand == Strand.REVERSE and ttrnas[0].end > att.end:  # 3'-overlap
                prime3_ttrnas.append(ttrna)

        return prime3_ttrnas

    @staticmethod
    def _is_same_direction(prime3_ttrnas: Sequence[Feature]) -> bool:
        strands = set(ttrna.strand for ttrna in prime3_ttrnas)
        return True if len(strands) == 1 else False


@dataclass(frozen=True)
class Pipolin:
    fragments: Tuple[PipolinFragment, ...]

    @staticmethod
    def from_fragments(*fragments: PipolinFragment):
        contigs = set(i.contig_id for i in fragments)
        if len(contigs) != len(fragments):
            raise AssertionError('Two fragments from the same contig!')
        return Pipolin(fragments)
