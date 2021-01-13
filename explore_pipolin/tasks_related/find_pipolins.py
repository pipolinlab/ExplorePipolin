from itertools import groupby
from typing import Sequence, Set, List, MutableSequence

from explore_pipolin.common import Genome, FeatureType, Pipolin, AttFeature, PipolinFragment, Range, Feature


class PipolinFinder:
    def __init__(self, genome: Genome):
        self.genome = genome

    @property
    def atts_by_att_id(self):
        return self.genome.features.get_features(FeatureType.ATT).get_atts_dict_by_att_id()

    def find_pipolins(self) -> Sequence[Pipolin]:
        pipolins = []

        pipolbs_dict = self.genome.features.pipolbs_dict()
        for _, pipolbs in pipolbs_dict.items():
            while len(pipolbs) > 0:
                leftmost_pipolb = pipolbs[0]
                atts_around_pipolb = self._find_atts_around_pipolb(leftmost_pipolb)
                orphan_atts = self._get_orphan_atts_same_att_ids(atts_around_pipolb)

                if len(atts_around_pipolb) == 0:
                    pipolins.append(self._pipolin_from_orphan_pipobs(pipolbs))
                else:
                    if len(orphan_atts) == 0:
                        pipolins.append(self._complete_fragment_pipolin(atts_around_pipolb))
                    else:
                        pipolins.append(self._pipolin_with_orphan_atts(atts_around_pipolb, orphan_atts))

                pipolbs_to_remove = []
                for pipolb in pipolbs:
                    if self._is_pipolb_in_pipolin(pipolb, pipolins[-1]):
                        pipolbs_to_remove.append(pipolb)
                pipolbs = [i for i in pipolbs if i not in pipolbs_to_remove]

        return pipolins

    def _pipolin_from_orphan_pipobs(self, pipolbs: MutableSequence[Feature]) -> Pipolin:
        pipolbs_to_include = [pipolbs[0]]
        contig_atts = self.genome.features.atts_dict()[pipolbs[0].contig_id]
        if len(contig_atts) != 0:
            return self._pipolin_oneside_att(pipolbs_to_include[0], contig_atts)
        else:
            if len(pipolbs) > 1:
                for i in range(1, len(pipolbs)):
                    if pipolbs[i].start - pipolbs[i - 1].end < 50000:
                        pipolbs_to_include.append(pipolbs[i])

            all_pipolbs = self.genome.features.get_features(FeatureType.PIPOLB)
            if len(pipolbs_to_include) == len(pipolbs) and \
                    len(pipolbs_to_include) == len(all_pipolbs) and \
                    len(self.atts_by_att_id) == 1:
                return self._single_pipolin_from_fragments(pipolbs)
            else:
                return self._pipolin_from_pipolbs(pipolbs_to_include)

    def _create_pipolin_fragment(self, start_feature: Feature, end_feature: Feature, inflate_size=100000) \
            -> PipolinFragment:
        contig_length = self.genome.get_contig_by_id(start_feature.contig_id).length
        if not isinstance(start_feature, AttFeature):
            loc = Range(0, min(end_feature.end + inflate_size, contig_length))
        elif not isinstance(end_feature, AttFeature):
            loc = Range(max(0, start_feature.start - inflate_size), contig_length)
        else:
            loc = Range(start_feature.start, end_feature.end).inflate(inflate_size, _max=contig_length)
        fragment = PipolinFragment(location=loc, contig_id=start_feature.contig_id, genome=self.genome)
        return fragment

    def _complete_fragment_pipolin(self, atts: Set[AttFeature]) -> Pipolin:
        atts = sorted(atts, key=lambda x: x.start)
        fragment = self._create_pipolin_fragment(atts[0], atts[-1], inflate_size=50)
        return Pipolin(fragment)

    def _pipolin_with_orphan_atts(self, atts: Set[AttFeature], orphan_atts: Set[AttFeature]) -> Pipolin:
        atts = sorted(atts, key=lambda x: x.start)
        fragment = self._create_pipolin_fragment(atts[0], atts[-1])

        orphan_fragments = self._fragments_from_orphan_atts(orphan_atts)

        return Pipolin(fragment, *orphan_fragments)

    def _fragments_from_orphan_atts(self, orphan_atts: Set[AttFeature]) -> Sequence[PipolinFragment]:
        orphan_fragments = []
        grouped_orphan_atts = groupby(orphan_atts, lambda x: x.contig_id)
        for contig_id, o_atts in grouped_orphan_atts:
            o_atts = sorted(o_atts, key=lambda x: x.start)
            orphan_fragments.append(self._create_pipolin_fragment(o_atts[0], o_atts[-1]))
        return orphan_fragments

    def _find_atts_around_pipolb(self, leftmost_pipolb: Feature) -> Set[AttFeature]:
        atts_around_pipolb = set()
        contig_atts = self.genome.features.atts_dict()[leftmost_pipolb.contig_id]
        for att in contig_atts:
            if att.start > leftmost_pipolb.end:
                other_atts = self.atts_by_att_id[att.att_id]
                other_atts.remove(att)
                for other_att in other_atts:
                    if other_att.contig_id == att.contig_id and other_att.end < leftmost_pipolb.start:
                        atts_around_pipolb.add(att)
                        atts_around_pipolb.add(other_att)
        return atts_around_pipolb

    def _get_orphan_atts_same_att_ids(self, atts_around_pipolbs) -> Set[AttFeature]:
        atts_same_att_id = set()
        att_ids = set(i.att_id for i in atts_around_pipolbs)
        for att_id in att_ids:
            for i in self.atts_by_att_id[att_id]:
                atts_same_att_id.add(i)
        return atts_same_att_id - atts_same_att_id.intersection(atts_around_pipolbs)

    def _single_pipolin_from_fragments(self, pipolbs: MutableSequence[Feature]) -> Pipolin:
        """
        contig1 contains all pipolb
        single att repeat on separate contigs
        """
        fragment = self._create_pipolin_fragment(pipolbs[0], pipolbs[-1])
        atts = self.genome.features.get_features(FeatureType.ATT)
        orphan_fragments = self._fragments_from_orphan_atts(atts)
        return Pipolin(fragment, *orphan_fragments)

    def _pipolin_from_pipolbs(self, pipolbs: List[Feature]) -> Pipolin:
        pipolbs = sorted(pipolbs, key=lambda x: x.start)
        fragment = self._create_pipolin_fragment(pipolbs[0], pipolbs[-1])
        return Pipolin(fragment)

    def _pipolin_oneside_att(self, pipolb: Feature, contig_atts: MutableSequence[AttFeature]) -> Pipolin:
        if contig_atts[0].start > pipolb.end:  # atts are on the right
            right_atts, o_atts = self._get_atts_and_orphan_atts(contig_atts, 0)
            fragment = self._create_pipolin_fragment(pipolb, right_atts[-1])
        elif contig_atts[-1].end < pipolb.start:  # atts are on the left
            left_atts, o_atts = self._get_atts_and_orphan_atts(contig_atts, -1)
            fragment = self._create_pipolin_fragment(left_atts[0], pipolb)
        else:
            raise AssertionError
        orphan_fragments = self._fragments_from_orphan_atts(o_atts)
        return Pipolin(fragment, *orphan_fragments)

    def _get_atts_and_orphan_atts(self, contig_atts, anchor_ind):
        atts = [att for att in contig_atts if att.att_id == contig_atts[anchor_ind].att_id]
        atts_same_att_id = self.atts_by_att_id[contig_atts[anchor_ind].att_id]
        o_atts = set(att for att in atts_same_att_id if att.contig_id != contig_atts[anchor_ind].contig_id)
        return atts, o_atts

    @staticmethod
    def _is_pipolb_in_pipolin(pipolb: Feature, pipolin: Pipolin) -> bool:
        for fragment in pipolin.fragments:
            if pipolb.contig_id == fragment.contig_id:
                if pipolb.start >= fragment.start and pipolb.end <= fragment.end:
                    return True
        return False
