import os
from enum import Enum, auto
from typing import MutableSequence, Set, Optional
from itertools import groupby
from random import randrange
import copy

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

from explore_pipolin.utilities.io import read_blastxml


class Orientation(Enum):
    FORWARD = auto()
    REVERSE = auto()

    @staticmethod
    def orientation_from_blast(hit_frame):
        if hit_frame == 1:
            return Orientation.FORWARD
        elif hit_frame == -1:
            return Orientation.REVERSE

    def to_pm_one_encoding(self):
        return +1 if self.FORWARD else -1

    def to_string(self):
        return 'forward' if self.FORWARD else 'reverse'

    def __neg__(self):
        if self is self.FORWARD:
            return self.REVERSE
        else:
            return self.FORWARD


class Contig:
    def __init__(self, contig_id, contig_length, orientation=Orientation.FORWARD):
        self.contig_id: str = contig_id
        self.contig_length: int = contig_length
        self.contig_orientation: Orientation = orientation


class Feature:
    def __init__(self, start: int, end: int, frame: Orientation, contig: Contig):
        self.start = start
        self.end = end
        self.frame = frame
        self.contig = contig


class PipolinFragment:
    def __init__(self, contig: Contig, start: int, end: int):
        self.contig = contig
        self.start = start
        self.end = end
        self.atts: MutableSequence[Feature] = []


class GQuery:
    def __init__(self, gquery_id: str):
        self.gquery_id = gquery_id
        self.contigs: MutableSequence[Contig] = []
        self.polbs: MutableSequence[Feature] = []
        self.atts: MutableSequence[Feature] = []
        self.trnas: MutableSequence[Feature] = []
        self.target_trnas: MutableSequence[Feature] = []
        self.denovo_atts: MutableSequence[Feature] = []
        self.target_trnas_denovo: MutableSequence[Feature] = []
        self.pipolin_fragments: MutableSequence[PipolinFragment] = []

    def get_contig_by_id(self, contig_id: str) -> Contig:
        for contig in self.contigs:
            if contig.contig_id == contig_id:
                return contig

    @staticmethod
    def dict_by_contig_normalized(features):
        return {contig: sorted(list(ps), key=lambda p: p.start) for contig, ps
                in groupby((i for i in features), key=lambda x: x.contig.contig_id)}

    # TODO: replace it with _dict_by_contig_normalized()!
    def get_features_of_contig(self, contig_id, feature_type: str) -> MutableSequence[Feature]:
        features_to_return = []
        features = self.get_features_by_type(feature_type=feature_type)
        for feature in features:
            if feature.contig.contig_id == contig_id:
                features_to_return.append(feature)

        return features_to_return

    def get_features_by_type(self, feature_type: str) -> MutableSequence[Feature]:
        if feature_type == 'polbs':
            return self.polbs
        elif feature_type == 'atts':
            return self.atts
        elif feature_type == 'target_trnas':
            return self.target_trnas
        else:
            raise AssertionError('Feature type can be: "polbs", "atts" or "target_trnas"!')

    def feature_from_blasthit(self, hit, contig_id) -> Feature:
        return Feature(start=hit.hit_start, end=hit.hit_end,
                       frame=Orientation.orientation_from_blast(hit.hit_frame),
                       contig=self.get_contig_by_id(contig_id))

    # `add_features_from_aragorn` and `add_features_atts_denovo`
    def find_target_trna(self, att: Feature) -> Optional[Feature]:
        trna_dict = self.dict_by_contig_normalized(self.trnas)

        if att.contig.contig_id in trna_dict:
            trnas = trna_dict[att.contig.contig_id]
            for trna in trnas:
                if self._is_overlapping(range1=(att.start, att.end), range2=(trna.start, trna.end)):
                    return trna

    def is_single_contig(self):
        return len(self.contigs) == 1

    def is_on_the_same_contig(self):
        target_contigs = []
        target_contigs.extend(i.contig.contig_id for i in self.polbs)
        target_contigs.extend(i.contig.contig_id for i in self.atts)
        target_contigs.extend(i.contig.contig_id for i in self.target_trnas)
        return len(set(target_contigs)) == 1

    # `find_atts_denovo`
    def get_left_right_windows(self):
        polymerases = sorted((i for i in self.polbs), key=lambda p: p.start)

        if polymerases[-1].start - polymerases[0].start > 10000:   # TODO: is it small/big enough?
            raise AssertionError(f'You have several piPolBs per genome and they are too far from each other: '
                                 f'within the region ({polymerases[0].start}, {polymerases[-1].end}). It might be, '
                                 f'that you have two or more pipolins per genome, but we are expecting only one.')

        length = self.contigs[0].contig_length
        left_edge = polymerases[0].start - 100000
        left_window = (left_edge if left_edge >= 0 else 0, polymerases[0].start)
        right_edge = polymerases[-1].end + 100000
        right_window = (polymerases[-1].end, right_edge if right_edge <= length else length)

        return left_window, right_window

    def is_att_denovo(self, left_repeat, right_repeat):
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
        targeted_contigs = [trna.contig.contig_id for trna in self.target_trnas]
        if len(self.target_trnas) != len(targeted_contigs):
            raise AssertionError("We are expecting a single tRNA to overlap with a single att per contig!")

    # `analyse_pipolin_orientation`
    def get_contig_orientation(self, contig: Contig) -> Orientation:
        target_trnas = self.get_features_of_contig(contig_id=contig.contig_id, feature_type='target_trnas')
        atts = self.get_features_of_contig(contig_id=contig.contig_id, feature_type='atts')
        atts_frames = [att.frame for att in atts]
        polbs = self.get_features_of_contig(contig_id=contig.contig_id, feature_type='polbs')
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

    # scaffolding is not required
    def get_pipolin_bounds(self):
        polymerases = sorted((i for i in self.polbs), key=lambda p: p.start)
        atts = sorted((i for i in self.atts), key=lambda p: p.start)

        length = polymerases[0].contig.contig_length
        left_edge = atts[0].start - 50 if atts[0].start < polymerases[0].start else polymerases[0].start - 50
        right_edge = atts[-1].end + 50 if atts[-1].end > polymerases[-1].end else polymerases[-1].end + 50
        return left_edge if left_edge >= 0 else 0, right_edge if right_edge <= length else length

    # scaffolding is required
    def try_creating_single_record(self):
        polbs_dict = self.dict_by_contig_normalized(self.polbs)
        atts_dict = self.dict_by_contig_normalized(self.atts)
        target_trnas_dict = self.dict_by_contig_normalized(self.target_trnas)

        unchangeable_contigs = self._get_unchangeable_contigs()

        if len(unchangeable_contigs) == 1:
            unchangeable_contig = unchangeable_contigs[0]
            print(f'The unchangeable contig is {unchangeable_contig.contig_id}!')

            att_only_contigs = self._order_att_only_contigs()
            if len(att_only_contigs) == 1:
                direction = self._get_direction_of_unchangeable(
                    polbs_sorted=polbs_dict[unchangeable_contig.contig_id],
                    atts_sorted=atts_dict[unchangeable_contig.contig_id])

                orphan_atts = atts_dict[att_only_contigs[0].contig_id]

                if direction == 'none':
                    contig_length = unchangeable_contig.contig_length
                    left_edge = atts_dict[unchangeable_contig.contig_id][0].start - 50
                    right_edge = atts_dict[unchangeable_contig.contig_id][-1].end + 50
                    pipolin = PipolinFragment(contig=self.get_contig_by_id(unchangeable_contig.contig_id),
                                              start=left_edge if left_edge >= 0 else 0,
                                              end=right_edge if right_edge <= contig_length else contig_length)
                    pipolin.atts.extend(atts_dict[unchangeable_contig.contig_id])
                    self.pipolin_fragments.append(pipolin)
                    # TODO: if att_only_contig has a target_trna, it could be added on the right
                elif direction == 'right':
                    left_fragment = self._create_unchangeable_fragment(atts_dict, unchangeable_contig, direction)
                    self.pipolin_fragments.append(left_fragment)

                    right_fragment = self._create_att_contig_fragment(contig_atts=orphan_atts, direction=direction)
                    self.pipolin_fragments.append(right_fragment)
                elif direction == 'left':
                    left_fragment = self._create_att_contig_fragment(contig_atts=orphan_atts, direction=direction)
                    self.pipolin_fragments.append(left_fragment)

                    right_fragment = self._create_unchangeable_fragment(atts_dict, unchangeable_contig, direction)
                    self.pipolin_fragments.append(right_fragment)
            elif len(att_only_contigs) == 2:
                # TODO: the order can be also [middle_fragment, left_fragment, right_fragment]
                middle_fragment = PipolinFragment(contig=self.get_contig_by_id(unchangeable_contig.contig_id),
                                                  start=0, end=unchangeable_contig.contig_length)
                middle_fragment.atts.extend(atts_dict[unchangeable_contig.contig_id])

                left_atts = atts_dict[att_only_contigs[0].contig_id]
                left_fragment = self._create_att_contig_fragment(contig_atts=left_atts, direction='left')

                right_atts = atts_dict[att_only_contigs[1].contig_id]
                right_fragment = self._create_att_contig_fragment(contig_atts=right_atts, direction='right')

                self.pipolin_fragments.extend([left_fragment, middle_fragment, right_fragment])
        elif len(unchangeable_contigs) == 2:
            target_trnas_contig0 = target_trnas_dict[unchangeable_contigs[0].contig_id]
            target_trnas_contig1 = target_trnas_dict[unchangeable_contigs[1].contig_id]

            if len(target_trnas_contig0) != 0:
                right_contig = unchangeable_contigs[0]
                left_contig = unchangeable_contigs[1]
            elif len(target_trnas_contig1) != 0:
                right_contig = unchangeable_contigs[1]
                left_contig = unchangeable_contigs[0]
            else:
                raise NotImplementedError

            left_direction = self._get_direction_of_unchangeable(polbs_sorted=polbs_dict[left_contig.contig_id],
                                                                 atts_sorted=atts_dict[left_contig.contig_id])
            left_fragment = self._create_unchangeable_fragment(atts_dict, left_contig, left_direction)

            right_direction = self._get_direction_of_unchangeable(polbs_sorted=polbs_dict[right_contig.contig_id],
                                                                  atts_sorted=atts_dict[right_contig.contig_id])
            right_fragment = self._create_unchangeable_fragment(atts_dict, right_contig, right_direction)

            self.pipolin_fragments.extend([left_fragment, right_fragment])
        elif len(unchangeable_contigs) > 2:
            raise NotImplementedError
        else:
            self._try_finish_separate(atts_dict)

    def _create_unchangeable_fragment(self, atts_dict, contig, direction):
        if direction == 'right':
            if contig.contig_orientation == Orientation.FORWARD:
                edge = atts_dict[contig.contig_id][0].start - 50
                fragment = PipolinFragment(contig=self.get_contig_by_id(contig.contig_id),
                                           start=edge if edge >= 0 else 0, end=contig.contig_length)
            else:
                edge = atts_dict[contig.contig_id][-1].end + 50
                fragment = PipolinFragment(contig=self.get_contig_by_id(contig.contig_id), start=0,
                                           end=edge if edge <= contig.contig_length else contig.contig_length)
        else:
            if contig.contig_orientation == Orientation.FORWARD:
                edge = atts_dict[contig.contig_id][-1].end + 50
                fragment = PipolinFragment(contig=self.get_contig_by_id(contig.contig_id), start=0,
                                           end=edge if edge <= contig.contig_length else contig.contig_length)
            else:
                edge = atts_dict[contig.contig_id][0].start - 50
                fragment = PipolinFragment(contig=self.get_contig_by_id(contig.contig_id),
                                           start=edge if edge >= 0 else 0, end=contig.contig_length)

        fragment.atts.extend(atts_dict[contig.contig_id])
        return fragment

    @staticmethod
    def _create_att_contig_fragment(contig_atts: MutableSequence[Feature], direction: str):
        contig_length = contig_atts[0].contig.contig_length
        if direction == 'right':
            if contig_atts[0].contig.contig_orientation == Orientation.FORWARD:
                left_edge = 0
                right_edge = contig_atts[-1].end + 50
            else:
                left_edge = contig_atts[0].start - 50
                right_edge = contig_length
        else:
            if contig_atts[0].contig.contig_orientation == Orientation.FORWARD:
                left_edge = contig_atts[0].start - 50
                right_edge = contig_length
            else:
                left_edge = 0
                right_edge = contig_atts[-1].end + 50

        fragment = PipolinFragment(contig=contig_atts[0].contig, start=left_edge if left_edge >= 0 else 0,
                                   end=right_edge if right_edge <= contig_length else contig_length)
        fragment.atts.extend(contig_atts)
        return fragment

    def _try_finish_separate(self, atts_dict):
        polbs_fragments = self._create_polbs_fragments()

        if len(self.target_trnas) > 1 or len(self.target_trnas) == 0:
            raise NotImplementedError

        right_fragment = self._create_att_contig_fragment(contig_atts=atts_dict[self.target_trnas[0].contig.contig_id],
                                                          direction='right')

        atts_contigs = set([i.contig for i in self.atts])
        if len(atts_contigs) == 2:
            print('The single record can be created!!!\n')

            atts_contigs = list(atts_contigs)
            if atts_contigs[0].contig_id == self.target_trnas[0].contig.contig_id:
                left_contig = atts_contigs[1]
            else:
                left_contig = atts_contigs[0]

            left_fragment = self._create_att_contig_fragment(contig_atts=atts_dict[left_contig.contig_id],
                                                             direction='left')
        else:
            raise NotImplementedError

        self.pipolin_fragments.append(left_fragment)
        self.pipolin_fragments.extend(polbs_fragments)
        self.pipolin_fragments.append(right_fragment)

    def _create_polbs_fragments(self) -> MutableSequence[PipolinFragment]:
        polbs_contigs = set([i.contig for i in self.polbs])
        if len(polbs_contigs) == 2:
            print('Two "polb only contigs" were found!')
            polbs_contigs = list(polbs_contigs)
            if polbs_contigs[0].contig_length + polbs_contigs[1].contig_length < 15000:
                polb0_fragment = PipolinFragment(contig=self.get_contig_by_id(polbs_contigs[0].contig_id), start=0,
                                                 end=polbs_contigs[0].contig_length)
                polb1_fragment = PipolinFragment(contig=self.get_contig_by_id(polbs_contigs[1].contig_id), start=0,
                                                 end=polbs_contigs[1].contig_length)

                polb0_length = sum((i.end - i.start) for i in self.get_features_of_contig(
                    contig_id=polbs_contigs[0].contig_id,
                    feature_type='polbs'))
                polb1_length = sum((i.end - i.start) for i in self.get_features_of_contig(
                    contig_id=polbs_contigs[1].contig_id,
                    feature_type='polbs'))
                # TODO: comparing just by length is an unreliable way! REWRITE if possible!
                if polb0_length < polb1_length:
                    polbs_fragments = [polb0_fragment, polb1_fragment]
                else:
                    polbs_fragments = [polb1_fragment, polb0_fragment]
            else:
                raise NotImplementedError

        elif len(polbs_contigs) == 1:
            polbs_fragments = [PipolinFragment(contig=self.get_contig_by_id(self.polbs[0].contig.contig_id),
                                               start=0, end=self.polbs[0].contig.contig_length)]
        else:
            raise NotImplementedError

        return polbs_fragments

    @staticmethod
    def _get_direction_of_unchangeable(polbs_sorted: MutableSequence[Feature], atts_sorted: MutableSequence[Feature]):
        if polbs_sorted[0].start > atts_sorted[-1].end:
            if polbs_sorted[0].contig.contig_orientation == Orientation.FORWARD:
                return 'right'
            else:
                return 'left'
        elif polbs_sorted[-1].end < atts_sorted[0].start:
            if polbs_sorted[-1].contig.contig_orientation == Orientation.FORWARD:
                return 'left'
            else:
                return 'right'
        else:
            return 'none'

    def _get_unchangeable_contigs(self) -> MutableSequence[Contig]:
        polbs_contigs = set([i.contig for i in self.polbs])

        contigs_to_return = []
        for contig in polbs_contigs:
            atts_next_polbs = self.get_features_of_contig(contig_id=contig.contig_id, feature_type='atts')
            if len(atts_next_polbs) != 0:
                contigs_to_return.append(contig)
        return contigs_to_return

    def _order_att_only_contigs(self) -> MutableSequence[Contig]:
        att_only_contigs = self._get_att_only_contigs()

        if len(att_only_contigs) == 1:
            print('The single record can be created!!!\n')
            return [att_only_contigs.pop()]
        elif len(att_only_contigs) == 2:
            contigs = list(att_only_contigs)

            if len(self.dict_by_contig_normalized(features=self.target_trnas)[contigs[0].contig_id]) != 0:
                # TODO: here sometimes fails for LREC241: KeyError: 'NODE_38'
                print('The single record can be created!!!\n')
                return [contigs[1], contigs[0]]
            elif len(self.dict_by_contig_normalized(features=self.target_trnas)[contigs[1].contig_id]) != 0:
                print('The single record can be created!!!\n')
                return [contigs[0], contigs[1]]
            else:
                raise NotImplementedError

        else:
            raise NotImplementedError

    def _get_att_only_contigs(self) -> Set[Contig]:
        att_only_contigs = set()
        for att in self.atts:
            polbs_next_att = self.get_features_of_contig(contig_id=att.contig.contig_id, feature_type='polbs')
            if len(polbs_next_att) == 0:
                att_only_contigs.add(att.contig)

        return att_only_contigs

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
        gb_qualifiers = {'inference': ['HMM:custom'], 'locus_tag': [f'{self.gquery_id}_{random_number}'],
                         'rpt_family': ['Att'], 'rpt_type': ['direct']}
        gff_qualifiers = {'phase': ['.'], 'source': ['HMM:custom'],
                          'ID': [f'{self.gquery_id}_{random_number}'], 'inference': ['HMM:custom'],
                          'locus_tag': [f'{self.gquery_id}_{random_number}'],
                          'rpt_family': ['Att'], 'rpt_type': ['direct']}
        att_feature = SeqFeature(type='repeat_region',
                                 location=FeatureLocation(start=start, end=end, strand=frame.to_pm_one_encoding()),
                                 qualifiers=gb_qualifiers if records_format == 'gb' else gff_qualifiers)
        return att_feature


def define_gquery_id(genome):
    return os.path.splitext(os.path.basename(genome))[0]


def save_left_right_subsequences(genome, left_window, right_window, repeats_dir):
    os.makedirs(repeats_dir, exist_ok=True)
    genome_seq = SeqIO.read(handle=genome, format='fasta')

    left_seq = genome_seq[left_window[0]:left_window[1]]
    right_seq = genome_seq[right_window[0]:right_window[1]]
    SeqIO.write(sequences=left_seq, format='fasta',
                handle=os.path.join(repeats_dir, define_gquery_id(genome=genome) + '.left'))
    SeqIO.write(sequences=right_seq, format='fasta',
                handle=os.path.join(repeats_dir, define_gquery_id(genome=genome) + '.right'))


def extract_repeats(file):
    repeats = read_blastxml(file)
    left_repeats = []
    right_repeats = []
    sequences = []
    for entry in repeats:
        for hit in entry:
            left_repeats.append((hit.query_start, hit.query_end))
            right_repeats.append((hit.hit_start, hit.hit_end))
            sequences.append(str(hit.hit.seq))
    return left_repeats, right_repeats, sequences


def set_proper_location(seq_range, shift):
    return seq_range[0] + shift, seq_range[1] + shift


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
        source_feature.qualifiers.update({'note': [fragment.contig.contig_id,
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
