from enum import Enum, auto
from typing import Sequence, MutableSequence, Optional, Set

from explore_pipolin.common import PipolinFragment, FeatureType, Contig, Feature, Orientation
from explore_pipolin.utilities.misc import GQuery


class Direction(Enum):
    LEFT = auto()
    RIGHT = auto()


class Scaffolder:
    def __init__(self, gquery: GQuery):
        self.gquery = gquery
        self.polbs_dict = gquery.get_features_dict_by_contig_normalized(FeatureType.PIPOLB)
        self.atts_dict = gquery.get_features_dict_by_contig_normalized(FeatureType.ATT)
        self.target_trnas_dict = gquery.get_features_dict_by_contig_normalized(FeatureType.TARGET_TRNA)

    def try_creating_single_record(self) -> Sequence[PipolinFragment]:
        unchangeable_contigs = self._get_unchangeable_contigs()

        if len(unchangeable_contigs) == 1:
            unchangeable_contig = unchangeable_contigs[0]
            print(f'The unchangeable contig is {unchangeable_contig.contig_id}!')

            att_only_contigs = self._order_att_only_contigs()
            if len(att_only_contigs) == 1:
                direction = self._get_direction_of_unchangeable(unchangeable_contig.contig_id)

                orphan_atts = self.atts_dict[att_only_contigs[0].contig_id]

                if direction == 'none':
                    contig_length = unchangeable_contig.contig_length
                    left_edge = self.atts_dict[unchangeable_contig.contig_id][0].start - 50
                    right_edge = self.atts_dict[unchangeable_contig.contig_id][-1].end + 50
                    pipolin = PipolinFragment(contig_id=unchangeable_contig.contig_id,
                                              genome=self.gquery.genome,
                                              start=left_edge if left_edge >= 0 else 0,
                                              end=right_edge if right_edge <= contig_length else contig_length)
                    pipolin.atts.extend(self.atts_dict[unchangeable_contig.contig_id])
                    return [pipolin]
                    # TODO: if att_only_contig has a target_trna, it could be added on the right
                elif direction == 'right':
                    left_fragment = self._create_unchangeable_fragment(unchangeable_contig, direction)
                    right_fragment = self._create_att_contig_fragment(contig_atts=orphan_atts, direction=direction)
                    return [left_fragment, right_fragment]
                elif direction == 'left':
                    left_fragment = self._create_att_contig_fragment(contig_atts=orphan_atts, direction=direction)
                    right_fragment = self._create_unchangeable_fragment(unchangeable_contig, direction)
                    return [left_fragment, right_fragment]
            elif len(att_only_contigs) == 2:
                # TODO: the order can be also [middle_fragment, left_fragment, right_fragment]
                middle_fragment = PipolinFragment(contig_id=unchangeable_contig.contig_id,
                                                  genome=self.gquery.genome,
                                                  start=0, end=unchangeable_contig.contig_length)
                middle_fragment.atts.extend(self.atts_dict[unchangeable_contig.contig_id])

                left_atts = self.atts_dict[att_only_contigs[0].contig_id]
                left_fragment = self._create_att_contig_fragment(contig_atts=left_atts, direction=Direction.LEFT)

                right_atts = self.atts_dict[att_only_contigs[1].contig_id]
                right_fragment = self._create_att_contig_fragment(contig_atts=right_atts, direction=Direction.RIGHT)

                return [left_fragment, middle_fragment, right_fragment]
        elif len(unchangeable_contigs) == 2:
            target_trnas_contig0 = self.target_trnas_dict[unchangeable_contigs[0].contig_id]
            target_trnas_contig1 = self.target_trnas_dict[unchangeable_contigs[1].contig_id]

            if len(target_trnas_contig0) != 0:
                right_contig = unchangeable_contigs[0]
                left_contig = unchangeable_contigs[1]
            elif len(target_trnas_contig1) != 0:
                right_contig = unchangeable_contigs[1]
                left_contig = unchangeable_contigs[0]
            else:
                raise NotImplementedError

            left_direction = self._get_direction_of_unchangeable(left_contig.contig_id)
            left_fragment = self._create_unchangeable_fragment(left_contig, left_direction)

            right_direction = self._get_direction_of_unchangeable(right_contig.contig_id)
            right_fragment = self._create_unchangeable_fragment(right_contig, right_direction)

            return [left_fragment, right_fragment]
        elif len(unchangeable_contigs) > 2:
            raise NotImplementedError
        else:
            return self._try_finish_separate()

    def _create_unchangeable_fragment(self, contig: Contig, direction: Direction) -> PipolinFragment:
        if direction is Direction.RIGHT:
            if contig.contig_orientation == Orientation.FORWARD:
                edge = self.atts_dict[contig.contig_id][0].start - 50
                fragment = PipolinFragment(contig_id=contig.contig_id,
                                           genome=self.gquery.genome,
                                           start=edge if edge >= 0 else 0, end=contig.contig_length)
            else:
                edge = self.atts_dict[contig.contig_id][-1].end + 50
                fragment = PipolinFragment(contig_id=contig.contig_id, start=0,
                                           genome=self.gquery.genome,
                                           end=edge if edge <= contig.contig_length else contig.contig_length)
        else:
            if contig.contig_orientation == Orientation.FORWARD:
                edge = self.atts_dict[contig.contig_id][-1].end + 50
                fragment = PipolinFragment(contig_id=contig.contig_id, start=0,
                                           genome=self.gquery.genome,
                                           end=edge if edge <= contig.contig_length else contig.contig_length)
            else:
                edge = self.atts_dict[contig.contig_id][0].start - 50
                fragment = PipolinFragment(contig_id=contig.contig_id,
                                           genome=self.gquery.genome,
                                           start=edge if edge >= 0 else 0, end=contig.contig_length)

        fragment.atts.extend(self.atts_dict[contig.contig_id])
        return fragment

    def _try_finish_separate(self) -> Sequence[PipolinFragment]:
        polbs_fragments = self._create_polbs_fragments()

        if len(self.gquery.target_trnas) != 1:
            raise NotImplementedError('At the moment I can only handle one target tRNA')

        the_target_trna, = self.gquery.target_trnas

        right_fragment = self._create_att_contig_fragment(contig_atts=self.atts_dict[the_target_trna.contig_id],
                                                          direction=Direction.RIGHT)

        atts_contigs = set([i.contig for i in self.gquery.atts])
        if len(atts_contigs) == 2:
            print('The single record can be created!!!\n')

            atts_contigs = list(atts_contigs)
            if atts_contigs[0].contig_id == the_target_trna.contig_id:
                left_contig = atts_contigs[1]
            else:
                left_contig = atts_contigs[0]

            left_fragment = self._create_att_contig_fragment(contig_atts=self.atts_dict[left_contig.contig_id],
                                                             direction=Direction.LEFT)
        else:
            raise NotImplementedError

        return [left_fragment, polbs_fragments, right_fragment]

    def _create_polbs_fragments(self) -> MutableSequence[PipolinFragment]:
        polbs_contigs = set([i.contig for i in self.gquery.pipolbs])
        if len(polbs_contigs) == 2:
            print('Two "polb only contigs" were found!')
            polbs_contigs = list(polbs_contigs)
            if polbs_contigs[0].contig_length + polbs_contigs[1].contig_length < 15000:
                polb0_fragment = PipolinFragment(contig_id=polbs_contigs[0].contig_id, start=0,
                                                 genome=self.gquery.genome,
                                                 end=polbs_contigs[0].contig_length)
                polb1_fragment = PipolinFragment(contig_id=polbs_contigs[1].contig_id, start=0,
                                                 genome=self.gquery.genome,
                                                 end=polbs_contigs[1].contig_length)

                polb0_length = sum((i.end - i.start) for i in self.gquery.get_features_of_contig_normalized(
                    contig_id=polbs_contigs[0].contig_id,
                    feature_type=FeatureType.PIPOLB))
                polb1_length = sum((i.end - i.start) for i in self.gquery.get_features_of_contig_normalized(
                    contig_id=polbs_contigs[1].contig_id,
                    feature_type=FeatureType.PIPOLB))
                # TODO: comparing just by length is an unreliable way! REWRITE if possible!
                if polb0_length < polb1_length:
                    polbs_fragments = [polb0_fragment, polb1_fragment]
                else:
                    polbs_fragments = [polb1_fragment, polb0_fragment]
            else:
                raise NotImplementedError

        elif len(polbs_contigs) == 1:
            the_polb_contig, = polbs_contigs
            polbs_fragments = [PipolinFragment(contig_id=the_polb_contig.contig_id,
                                               genome=self.gquery.genome,
                                               start=0, end=the_polb_contig.contig_length)]
        else:
            raise NotImplementedError

        return polbs_fragments

    def _get_unchangeable_contigs(self) -> MutableSequence[Contig]:
        polbs_contigs = set([i.contig for i in self.gquery.pipolbs])

        contigs_to_return = []
        for contig in polbs_contigs:
            atts_next_polbs = self.gquery.get_features_of_contig_normalized(contig_id=contig.contig_id,
                                                                            feature_type=FeatureType.ATT)
            if len(atts_next_polbs) != 0:
                contigs_to_return.append(contig)

        return contigs_to_return

    def _create_att_contig_fragment(self, contig_atts: Sequence[Feature],
                                    direction: Direction) -> PipolinFragment:
        contig = contig_atts[0].contig
        contig_length = contig.contig_length
        if direction is Direction.RIGHT:
            if contig.contig_orientation == Orientation.FORWARD:
                left_edge = 0
                right_edge = contig_atts[-1].end + 50
            else:
                left_edge = contig_atts[0].start - 50
                right_edge = contig_length
        else:
            if contig.contig_orientation == Orientation.FORWARD:
                left_edge = contig_atts[0].start - 50
                right_edge = contig_length
            else:
                left_edge = 0
                right_edge = contig_atts[-1].end + 50

        fragment = PipolinFragment(contig_id=contig.contig_id, genome=self.gquery.genome,
                                   start=left_edge if left_edge >= 0 else 0,
                                   end=right_edge if right_edge <= contig_length else contig_length)
        fragment.atts.extend(contig_atts)
        return fragment

    def _order_att_only_contigs(self) -> MutableSequence[Contig]:
        att_only_contigs = list(self._get_att_only_contigs())
        if len(att_only_contigs) == 1:
            print('The single record can be created!!!\n')
            return att_only_contigs
        elif len(att_only_contigs) == 2:
            if self._contig_has_feature_type(att_only_contigs[0].contig_id, FeatureType.TARGET_TRNA):
                # TODO: here sometimes fails for LREC241: KeyError: 'NODE_38'
                print('The single record can be created!!!\n')
                return [att_only_contigs[1], att_only_contigs[0]]
            elif self._contig_has_feature_type(att_only_contigs[1].contig_id, FeatureType.TARGET_TRNA):
                print('The single record can be created!!!\n')
                return [att_only_contigs[0], att_only_contigs[1]]
            else:
                raise NotImplementedError
        else:
            raise NotImplementedError

    def _get_att_only_contigs(self) -> Set[Contig]:
        att_only_contigs = set()
        for att in self.gquery.atts:
            polbs_next_att = self.gquery.get_features_of_contig_normalized(contig_id=att.contig_id,
                                                                           feature_type=FeatureType.PIPOLB)
            if len(polbs_next_att) == 0:
                att_only_contigs.add(att.contig)

        return att_only_contigs

    def _get_direction_of_unchangeable(self, contig_id: str) -> Optional[Direction]:
        polbs_sorted = self.polbs_dict[contig_id]
        atts_sorted = self.atts_dict[contig_id]
        if polbs_sorted[0].start > atts_sorted[-1].end:
            if polbs_sorted[0].contig.contig_orientation == Orientation.FORWARD:
                return Direction.RIGHT
            else:
                return Direction.LEFT
        elif polbs_sorted[-1].end < atts_sorted[0].start:
            if polbs_sorted[-1].contig.contig_orientation == Orientation.FORWARD:
                return Direction.LEFT
            else:
                return Direction.RIGHT
        else:
            return None

    def _contig_has_feature_type(self, contig_id: str, feature_type: FeatureType) -> bool:
        feature_dict = self.gquery.get_features_dict_by_contig_normalized(feature_type)
        return len(feature_dict[contig_id]) != 0


def create_pipolin_fragments_single_contig(gquery: GQuery) -> Sequence[PipolinFragment]:
    if len(gquery.get_features_dict_by_contig_normalized(FeatureType.ATT)) != 0:
        start, end = _get_pipolin_bounds(gquery)
        pipolin = PipolinFragment(contig_id=gquery.pipolbs[0].contig.contig_id, genome=gquery.genome,
                                  start=start, end=end)

        pipolin.atts.extend(gquery.atts)
        return [pipolin]
    else:
        left_window, right_window = gquery.get_left_right_windows()
        pipolin = PipolinFragment(contig_id=gquery.pipolbs[0].contig.contig_id, genome=gquery.genome,
                                  start=left_window[0], end=right_window[1])
        return [pipolin]


def _get_pipolin_bounds(gquery: GQuery):
    polymerases = sorted((i for i in gquery.pipolbs), key=lambda p: p.start)
    atts = sorted((i for i in gquery.atts), key=lambda p: p.start)

    length = polymerases[0].contig.contig_length
    left_edge = atts[0].start - 50 if atts[0].start < polymerases[0].start else polymerases[0].start - 50
    right_edge = atts[-1].end + 50 if atts[-1].end > polymerases[-1].end else polymerases[-1].end + 50
    return left_edge if left_edge >= 0 else 0, right_edge if right_edge <= length else length
