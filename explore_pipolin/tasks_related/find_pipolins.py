from itertools import groupby
from typing import Sequence, Set, List, Tuple

from prefect import task

from explore_pipolin.common import Genome, FeatureType, Pipolin, AttFeature, PipolinFragment, Range, Feature


@task()
def find_pipolins(genome: Genome) -> Sequence[Pipolin]:
    finder = PipolinFinder(genome)
    pipolins = []

    pipolbs_dict = genome.features.pipolbs_dict()
    for _, pipolbs in pipolbs_dict.items():
        pipolbs = delete_used_pipolbs(pipolbs, pipolins)
        while len(pipolbs) > 0:
            pipolins.append(finder.pipolin_from_part_of_pipolbs(pipolbs))
            pipolbs = delete_used_pipolbs(pipolbs, pipolins)

    return pipolins


class PipolinFinder:
    def __init__(self, genome: Genome):
        self.genome = genome

    @property
    def atts_by_att_id(self):
        return self.genome.features.get_features(FeatureType.ATT).get_atts_dict_by_att_id()

    def pipolin_from_part_of_pipolbs(self, pipolbs: Sequence[Feature]) -> Pipolin:
        leftmost_pipolb = pipolbs[0]
        atts_around_pipolb = self._find_atts_around_pipolb(leftmost_pipolb)
        orphan_atts = self._get_orphan_atts(atts_around_pipolb)

        if len(atts_around_pipolb) == 0:
            return self._disrupted_pipolin(pipolbs)
        elif len(orphan_atts) == 0:
            return self._complete_fragment_pipolin(atts_around_pipolb)
        else:
            return self._pipolin_with_orphan_atts(atts_around_pipolb, orphan_atts)

    def _find_atts_around_pipolb(self, pipolb: Feature) -> Set[AttFeature]:
        atts_around_pipolb = set()
        contig_atts = self.genome.features.atts_dict()[pipolb.contig_id]
        for att in contig_atts:
            if att.start > pipolb.end:
                for other_att in contig_atts:
                    if other_att.end < pipolb.start and other_att.att_id == att.att_id:
                        atts_around_pipolb.update([att, other_att])
        return atts_around_pipolb

    def _get_orphan_atts(self, atts_around_pipolbs: Set[AttFeature]) -> Set[AttFeature]:
        atts_same_att_id = set()
        att_ids = set(i.att_id for i in atts_around_pipolbs)
        for att_id in att_ids:
            atts_same_att_id.update(i for i in self.atts_by_att_id[att_id])
        return atts_same_att_id - atts_around_pipolbs

    def _disrupted_pipolin(self, pipolbs: Sequence[Feature]) -> Pipolin:
        contig_atts = self.genome.features.atts_dict()[pipolbs[0].contig_id]
        if len(contig_atts) != 0:
            return self._pipolin_pipolb_nextto_att(pipolbs[0], contig_atts)
        else:
            return self._pipolin_from_orphan_pipolbs(pipolbs)

    def _pipolin_pipolb_nextto_att(self, pipolb: Feature, contig_atts: Sequence[AttFeature]) -> Pipolin:

        if contig_atts[0].start > pipolb.end:  # atts are on the right
            right_atts, orphan_atts = self._atts_nextto_pipolb_and_orphan_atts(contig_atts, 0)
            right_atts = sorted(right_atts, key=lambda x: x.start)

            fragment = self._create_pipolin_fragment(pipolb, right_atts[-1])
        elif contig_atts[-1].end < pipolb.start:  # atts are on the left
            left_atts, orphan_atts = self._atts_nextto_pipolb_and_orphan_atts(contig_atts, -1)
            left_atts = sorted(left_atts, key=lambda x: x.start)

            fragment = self._create_pipolin_fragment(left_atts[0], pipolb)
        else:
            raise AssertionError
        other_fragments = self._fragments_from_orphan_atts(orphan_atts)
        return Pipolin.from_fragments(fragment, *other_fragments)

    def _atts_nextto_pipolb_and_orphan_atts(self, contig_atts: Sequence[AttFeature], anchor_ind: int) \
            -> Tuple[Set[AttFeature], Set[AttFeature]]:
        atts_nextto_pipolb = set(att for att in contig_atts if att.att_id == contig_atts[anchor_ind].att_id)
        atts_same_att_id = set(self.atts_by_att_id[contig_atts[anchor_ind].att_id])
        orphan_atts = atts_same_att_id - atts_nextto_pipolb
        return atts_nextto_pipolb, orphan_atts

    def _create_pipolin_fragment(self, start_feature: Feature, end_feature: Feature,
                                 inflate_size: int = 100000) -> PipolinFragment:

        contig_length = self.genome.get_contig_by_id(start_feature.contig_id).length
        is_start_feature_att = isinstance(start_feature, AttFeature)
        is_end_feature_att = isinstance(end_feature, AttFeature)

        if not is_start_feature_att and not is_end_feature_att:
            loc = Range(0, contig_length)

        elif not is_start_feature_att and is_end_feature_att:
            loc = Range(0, min(end_feature.end + inflate_size, contig_length))

        elif is_start_feature_att and not is_end_feature_att:
            loc = Range(max(0, start_feature.start - inflate_size), contig_length)

        else:
            loc = Range(max(0, start_feature.start - inflate_size),
                        min(end_feature.end + inflate_size, contig_length))

        fragment = PipolinFragment(location=loc, contig_id=start_feature.contig_id)
        fragment_features = self._get_fragment_features_sorted(fragment)

        return PipolinFragment(fragment.location, fragment.contig_id, fragment_features)

    def _fragments_from_orphan_atts(self, orphan_atts: Set[AttFeature]) -> Sequence[PipolinFragment]:
        fragments = []

        grouped_orphan_atts = groupby(sorted(orphan_atts, key=lambda x: x.contig_id), lambda x: x.contig_id)

        for contig_id, atts in grouped_orphan_atts:
            atts: List[AttFeature] = sorted(atts, key=lambda x: x.start)
            fragments.append(self._create_pipolin_fragment(atts[0], atts[-1]))

        return fragments

    def _pipolin_from_orphan_pipolbs(self, pipolbs: Sequence[Feature]) -> Pipolin:
        pipolbs_to_include = [pipolbs[0]]
        if len(pipolbs) > 1:
            for i in range(1, len(pipolbs)):
                if pipolbs[i].start - pipolbs[i - 1].end < 50000:
                    pipolbs_to_include.append(pipolbs[i])
        all_pipolbs = self.genome.features.get_features(FeatureType.PIPOLB)

        if len(pipolbs_to_include) == len(pipolbs) and \
                len(pipolbs_to_include) == len(all_pipolbs) and \
                len(self.atts_by_att_id) == 1:
            return self._single_disrupted_pipolin(pipolbs)
        else:
            return self._pipolin_from_just_pipolbs(pipolbs_to_include)

    def _single_disrupted_pipolin(self, pipolbs: Sequence[Feature]) -> Pipolin:
        """
        a contig contains all pipolbs
        att repeats with the same att_id on separate contigs
        """
        fragment = self._create_pipolin_fragment(pipolbs[0], pipolbs[-1])

        atts = self.genome.features.get_features(FeatureType.ATT)
        orphan_fragments = self._fragments_from_orphan_atts(atts)

        return Pipolin.from_fragments(fragment, *orphan_fragments)

    def _pipolin_from_just_pipolbs(self, pipolbs: List[Feature]) -> Pipolin:
        pipolbs = sorted(pipolbs, key=lambda x: x.start)
        fragment = self._create_pipolin_fragment(pipolbs[0], pipolbs[-1])
        return Pipolin.from_fragments(fragment)

    def _complete_fragment_pipolin(self, atts: Set[AttFeature]) -> Pipolin:
        atts = sorted(atts, key=lambda x: x.start)
        fragment = self._create_pipolin_fragment(atts[0], atts[-1], inflate_size=50)
        return Pipolin.from_fragments(fragment)

    def _pipolin_with_orphan_atts(self, atts: Set[AttFeature], orphan_atts: Set[AttFeature]) -> Pipolin:
        atts = sorted(atts, key=lambda x: x.start)
        fragment = self._create_pipolin_fragment(atts[0], atts[-1])

        orphan_fragments = self._fragments_from_orphan_atts(orphan_atts)

        return Pipolin.from_fragments(fragment, *orphan_fragments)

    def _get_fragment_features_sorted(
            self, fragment: PipolinFragment) -> Sequence[Tuple[Feature, FeatureType]]:

        fragment_pipolbs = self._get_fragment_features_of_type(fragment, FeatureType.PIPOLB)
        fragment_atts = self._get_fragment_features_of_type(fragment, FeatureType.ATT)
        fragment_ttrnas = self._get_fragment_features_of_type(fragment, FeatureType.TARGET_TRNA)

        fragment_features = [(i, FeatureType.PIPOLB) for i in fragment_pipolbs]
        fragment_features.extend([(i, FeatureType.ATT) for i in fragment_atts])
        fragment_features.extend([(i, FeatureType.TARGET_TRNA) for i in fragment_ttrnas])

        return tuple(sorted(fragment_features, key=lambda x: x[0].start))

    def _get_fragment_features_of_type(
            self, fragment: PipolinFragment, feature_type: FeatureType) -> Sequence[Feature]:

        features = self.genome.features.get_features(feature_type)
        contig_features = features.get_dict_by_contig_sorted()[fragment.contig_id]
        return [f for f in contig_features if f.start >= fragment.start and f.end <= fragment.end]


def delete_used_pipolbs(pipolbs, pipolins):
    pipolbs_in_pipolins = []
    for pipolb in pipolbs:
        for pipolin in pipolins:
            if _is_pipolb_in_pipolin(pipolb, pipolin):
                pipolbs_in_pipolins.append(pipolb)

    pipolbs = [i for i in pipolbs if i not in pipolbs_in_pipolins]
    return pipolbs


def _is_pipolb_in_pipolin(pipolb: Feature, pipolin: Pipolin) -> bool:
    for fragment in pipolin.fragments:
        if pipolb.contig_id == fragment.contig_id:
            if pipolb.start >= fragment.start and pipolb.end <= fragment.end:
                return True
    return False
