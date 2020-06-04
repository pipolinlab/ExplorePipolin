from __future__ import annotations

from collections import defaultdict
from enum import Enum, auto
from typing import MutableSequence, Optional, Tuple, Sequence, Mapping, MutableMapping, List
from random import randrange
import copy

from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

from explore_pipolin.utilities.common import Orientation, Contig, Genome


class Feature:
    def __init__(self, start: int, end: int, frame: Orientation, contig_id: str, genome: Genome):
        self.start = start
        self.end = end
        self.frame = frame
        self.contig_id = contig_id
        self.genome = genome

    @property
    def contig(self):
        return self.genome.get_contig_by_id(self.contig_id)


class PipolinFragment:
    def __init__(self, contig_id: str, genome: Genome, start: int, end: int):
        self.contig_id = contig_id
        self.start = start
        self.end = end
        self.atts: MutableSequence[Feature] = []
        self.genome = genome

    @property
    def contig(self):
        return self.genome.get_contig_by_id(self.contig_id)


class FeatureType(Enum):
    POLB = auto()
    ATT = auto()
    TARGET_TRNA = auto()
    TRNA = auto()


class GQuery:
    def __init__(self, genome: Genome):
        self.genome = genome
        self.polbs: MutableSequence[Feature] = []
        self.atts: MutableSequence[Feature] = []
        self.trnas: MutableSequence[Feature] = []
        self.target_trnas: MutableSequence[Feature] = []
        self.denovo_atts: MutableSequence[Feature] = []
        self.target_trnas_denovo: MutableSequence[Feature] = []
        self.pipolin_fragments: MutableSequence[PipolinFragment] = []

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

    # TODO: replace it with _dict_by_contig_normalized()!
    def get_features_of_contig(self, contig_id, feature_type: FeatureType) -> MutableSequence[Feature]:
        features_to_return = []
        features = self.get_features_by_type(feature_type=feature_type)
        for feature in features:
            if feature.contig_id == contig_id:
                features_to_return.append(feature)
        return features_to_return

    def get_features_by_type(self, feature_type: FeatureType) -> MutableSequence[Feature]:
        if feature_type is FeatureType.POLB:
            return self.polbs
        elif feature_type is FeatureType.ATT:
            return self.atts
        elif feature_type is FeatureType.TARGET_TRNA:
            return self.target_trnas
        elif feature_type is FeatureType.TRNA:
            return self.trnas
        else:
            raise AssertionError(f'Feature must be one of: {list(FeatureType)}, not {feature_type}')

    def feature_from_blasthit(self, hit, contig_id) -> Feature:
        return Feature(start=hit.hit_start, end=hit.hit_end,
                       frame=Orientation.orientation_from_blast(hit.hit_frame),
                       contig_id=contig_id, genome=self.genome)

    # `add_features_from_aragorn` and `add_features_atts_denovo`
    def find_target_trna(self, att: Feature) -> Optional[Feature]:
        trna_dict = self.get_features_dict_by_contig_normalized(FeatureType.TRNA)

        if att.contig_id in trna_dict:
            trnas = trna_dict[att.contig_id]
            for trna in trnas:
                if self._is_overlapping(range1=(att.start, att.end), range2=(trna.start, trna.end)):
                    return trna

    def is_on_the_same_contig(self) -> bool:
        target_contigs = []
        target_contigs.extend(f.contig_id for f in self.polbs)
        target_contigs.extend(f.contig_id for f in self.atts)
        target_contigs.extend(f.contig_id for f in self.target_trnas)
        return len(set(target_contigs)) == 1

    # `find_atts_denovo` and `scaffold_pipolins`
    def get_left_right_windows(self) -> Tuple[Tuple[int, int], Tuple[int, int]]:
        polymerases = sorted((i for i in self.polbs), key=lambda p: p.start)

        if polymerases[-1].start - polymerases[0].start > 10000:   # TODO: is it small/big enough?
            raise AssertionError(f'You have several piPolBs per genome and they are too far from each other: '
                                 f'within the region ({polymerases[0].start}, {polymerases[-1].end}). It might be, '
                                 f'that you have two or more pipolins per genome, but we are expecting only one.')

        length = self.genome.get_complete_genome_length()
        left_edge = polymerases[0].start - 100000
        left_window = (left_edge if left_edge >= 0 else 0, polymerases[0].start)
        right_edge = polymerases[-1].end + 100000
        right_window = (polymerases[-1].end, right_edge if right_edge <= length else length)

        return left_window, right_window

    def is_att_denovo(self, left_repeat: Tuple[int, int], right_repeat: Tuple[int, int]) -> bool:
        if self._is_overlapping_att(left_repeat=left_repeat):
            return False
        return self._is_overlapping_trna(left_repeat=left_repeat, right_repeat=right_repeat)

    def _is_overlapping_att(self, left_repeat):
        for att in self.atts:
            if self._is_overlapping(left_repeat, (att.start, att.end)):
                return True
        return False

    def _is_overlapping_trna(self, left_repeat, right_repeat):
        for trna in self.trnas:
            trna_range = (trna.start, trna.end)
            if self._is_overlapping(left_repeat, trna_range) or self._is_overlapping(right_repeat, trna_range):
                return True
        return False

    # `analyse_pipolin_orientation`
    def is_single_target_trna_per_contig(self):
        # TODO: don't like this
        # there was one case with two target trnas per genome, although usually only one
        targeted_contigs = [trna.contig_id for trna in self.target_trnas]
        if len(self.target_trnas) != len(targeted_contigs):
            raise AssertionError("We are expecting a single tRNA to overlap with a single att per contig!")

    # `analyse_pipolin_orientation`
    def get_contig_orientation(self, contig: Contig) -> Orientation:
        target_trnas = self.get_features_of_contig(contig_id=contig.contig_id, feature_type=FeatureType.TARGET_TRNA)
        atts = self.get_features_of_contig(contig_id=contig.contig_id, feature_type=FeatureType.ATT)
        atts_frames = [att.frame for att in atts]
        polbs = self.get_features_of_contig(contig_id=contig.contig_id, feature_type=FeatureType.POLB)
        polbs_frames = [polb.frame for polb in polbs]

        if len(target_trnas) != 0:
            if len(set(atts_frames)) != 1:
                raise AssertionError('ATTs are expected to be in the same frame, as they are direct repeats!')
            if set(atts_frames).pop() == target_trnas[0].frame:
                raise AssertionError('ATT and tRNA are expected to be on the different strands!')
            return - target_trnas[0].frame

        elif len(atts) != 0:
            if len(set(atts_frames)) != 1:
                raise AssertionError('ATTs are expected to be in the same frame, as they are direct repeats!')
            return atts[0].frame

        if len(polbs) != 0:
            if len(set(polbs_frames)) != 1:  # an ambiguous case
                return contig.contig_orientation
            return polbs[0].frame

    @staticmethod
    def _is_overlapping(range1, range2):
        max_start = max(range1[0], range2[0])
        min_end = min(range1[1], range2[1])
        return max_start <= min_end

    @staticmethod
    def _is_polymerase_inside(atts, polymerases):
        return atts[0].start < polymerases[0].start and polymerases[-1].end < atts[-1].end

    def create_att_feature(self, start: int, end: int, frame: Orientation, records_format: str) -> SeqFeature:
        random_number = randrange(10000, 99999)
        gb_qualifiers = {'inference': ['HMM:custom'], 'locus_tag': [f'{self.genome.id}_{random_number}'],
                         'rpt_family': ['Att'], 'rpt_type': ['direct']}
        gff_qualifiers = {'phase': ['.'], 'source': ['HMM:custom'],
                          'ID': [f'{self.genome.id}_{random_number}'], 'inference': ['HMM:custom'],
                          'locus_tag': [f'{self.genome.id}_{random_number}'],
                          'rpt_family': ['Att'], 'rpt_type': ['direct']}
        att_feature = SeqFeature(type='repeat_region',
                                 location=FeatureLocation(start=start, end=end, strand=frame.to_pm_one_encoding()),
                                 qualifiers=gb_qualifiers if records_format == 'gb' else gff_qualifiers)
        return att_feature


def add_new_gb_feature(new_feature: SeqFeature, record: SeqRecord):
    record.features.append(new_feature)
    record.features.sort(key=lambda x: x.location.start)


def create_new_gb_record(gquery: GQuery, gb_record: SeqRecord) -> SeqRecord:
    new_record = copy.deepcopy(gb_record)
    new_source_features = []

    in_start = 0
    next_start = 0
    for fragment in gquery.pipolin_fragments:
        fragment_shift = fragment.start
        next_start += (fragment.end - fragment.start) + 100

        for i_f, feature in enumerate(gb_record.features):
            if feature.location.start >= in_start and feature.location.end <= next_start:
                old_start, old_end, old_strand = feature.location.start, feature.location.end, feature.location.strand
                new_record.features[i_f].location = FeatureLocation(start=old_start - in_start + fragment_shift,
                                                                    end=old_end - in_start + fragment_shift,
                                                                    strand=old_strand)

        in_start += next_start

        source_location = FeatureLocation(start=fragment.start, end=fragment.end + 100)
        source_feature = SeqFeature(type='source', location=source_location,
                                    qualifiers=copy.deepcopy(gb_record.features[0].qualifiers))
        source_feature.qualifiers.update({'note': [fragment.contig_id,
                                                   fragment.contig.contig_orientation.to_string()]})
        new_source_features.append(source_feature)

    del new_record.features[0]
    for feature in new_source_features:
        add_new_gb_feature(new_feature=feature, record=new_record)

    return new_record


def create_att_seqfeatures(record_format: str, gquery: GQuery) -> MutableSequence[SeqFeature]:
    att_seqfeatures = []
    in_start = 0
    for fragment in gquery.pipolin_fragments:
        fragment_shift = fragment.start if fragment.contig.contig_orientation == Orientation.FORWARD else fragment.end
        for att in fragment.atts:
            att_start, att_end = sorted([abs(att.start - fragment_shift), abs(att.end - fragment_shift)])
            att_feature = gquery.create_att_feature(start=att_start + in_start, end=att_end + in_start,
                                                    frame=att.frame, records_format=record_format)
            att_seqfeatures.append(att_feature)
        in_start += (fragment.end - fragment.start) + 100

    return att_seqfeatures


# def create_assembly_gap_record(record):
#     source_feature = SeqFeature(type='source', location=FeatureLocation(1, 100, strand=+1),
#                                 qualifiers={'mol_type': record.features[0].qualifiers['mol_type'],
#                                             'organism': record.features[0].qualifiers['organism'],
#                                             'strain': record.features[0].qualifiers['strain']})
#     assembly_gap_seq = Seq('N' * 100, alphabet=IUPACAmbiguousDNA())
#     assembly_gap_qualifiers = {'estimated_length': ['unknown'],
#                                'gap_type': ['within_scaffolds'],
#                                'linkage_evidence': ['pipolin_structure']}
#     assembly_gap_feature = SeqFeature(type='assembly_gap',
#                                       location=FeatureLocation(1, 100, strand=+1),
#                                       qualifiers=assembly_gap_qualifiers)
#     assembly_gap_record = SeqRecord(seq=assembly_gap_seq, id=record.id, name=record.name,
#                                     description=record.description, features=[source_feature, assembly_gap_feature],
#                                     annotations=record.annotations)
#
#     return assembly_gap_record
def create_fragment_record(fragment, genome_dict):
    fragment_record = genome_dict[fragment.contig.contig_id][fragment.start:fragment.end]
    if fragment.contig.contig_orientation == Orientation.REVERSE:
        fragment_record = fragment_record.reverse_complement()
    return fragment_record
