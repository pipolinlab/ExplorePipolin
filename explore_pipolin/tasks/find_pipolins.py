from dataclasses import dataclass
from itertools import groupby
from typing import Sequence, Set, Tuple, MutableSequence, Optional

from prefect import task

from explore_pipolin.common import Genome, FeatureType, Pipolin, AttFeature, PipolinFragment, Range, Feature
from explore_pipolin.utilities.logging import genome_specific_logging


@task()
@genome_specific_logging
def find_pipolins(genome: Genome) -> Sequence[Pipolin]:
    return PipolinFinder(genome).find_pipolins()


PipolinScore = Tuple[int, int]


@dataclass(frozen=True)
class PipolinCandidate:
    pipolin: Pipolin
    score: PipolinScore


class PipolinFinder:
    def __init__(self, genome: Genome):
        self.genome = genome

    @property
    def atts_by_att_id(self):
        return self.genome.features.get_features(FeatureType.ATT).get_atts_dict_by_att_id()

    def find_pipolins(self) -> Sequence[Pipolin]:
        pipolin_fragment_candidates = self._find_pipolin_fragment_candidates()
        pipolin_candidates = self._find_pipolin_candidates(pipolin_fragment_candidates)
        scores = [self._calc_pipolin_score(i) for i in pipolin_candidates]
        candidates = [PipolinCandidate(i, y) for i, y in zip(pipolin_candidates, scores)]

        best_pipolins = self._find_best_non_overlapping(candidates)
        return best_pipolins

    def _find_best_non_overlapping(self, candidates) -> Sequence[Pipolin]:
        candidates = sorted(candidates, key=lambda x: x.score, reverse=True)
        result: MutableSequence[Pipolin] = []
        while candidates:

            best_candidate = candidates[0]
            result.append(best_candidate.pipolin)

            candidates = self._filter_candidates(best_candidate, candidates)

        return result

    def _filter_candidates(self, best_candidate, candidates):
        filtered_candidates = []
        for cdt in candidates:
            if not self._is_candidate_overlapping_best(cdt.pipolin, best_candidate.pipolin):
                filtered_candidates.append(cdt)
        return filtered_candidates

    @staticmethod
    def _is_candidate_overlapping_best(candidate: Pipolin, best: Pipolin):
        for b_f in best.fragments:
            for c_f in candidate.fragments:
                if b_f.is_overlapping(c_f):
                    return True
        return False

    def _calc_pipolin_score(self, pipolin: Pipolin) -> PipolinScore:
        max_score = max(self._calc_pipolin_fragment_score(f) for f in pipolin.fragments)
        sum_score = sum(self._calc_pipolin_fragment_score(f) for f in pipolin.fragments)
        return max_score, sum_score

    @staticmethod
    def _calc_pipolin_fragment_score(fragment: PipolinFragment) -> int:
        atts_and_pipolbs = [f for f in fragment.features if f.ftype != FeatureType.TARGET_TRNA]
        first_feature_type = atts_and_pipolbs[0].ftype
        last_feature_type = atts_and_pipolbs[-1].ftype

        if any(f.ftype == FeatureType.PIPOLB for f in fragment.features):
            if first_feature_type == FeatureType.ATT and last_feature_type == FeatureType.ATT:
                return 100000 * len(atts_and_pipolbs)
            elif first_feature_type == FeatureType.ATT or last_feature_type == FeatureType.ATT:
                return 1000 * len(atts_and_pipolbs)
            else:
                return 10 * len(atts_and_pipolbs)
        else:
            return 1 * len(atts_and_pipolbs)

    def _find_pipolin_fragment_candidates(self):
        pipolin_fragments: MutableSequence[PipolinFragment] = []
        pipolbs_dict = self.genome.features.pipolbs_dict()
        for pipolbs in pipolbs_dict.values():
            for pipolb in pipolbs:
                pipolin_fragments.extend(self._find_pipolin_fragments_with_pipolb(pipolb))
        return pipolin_fragments

    def _find_pipolin_fragments_with_pipolb(self, pipolb: Feature):
        pipolin_fragments: MutableSequence[PipolinFragment] = []

        atts_around_pipolb = sorted(self._find_atts_around_pipolb(pipolb), key=lambda x: x.start)
        if atts_around_pipolb:
            pipolin_fragments.append(self._create_pipolin_fragment(atts_around_pipolb[0], atts_around_pipolb[-1]))

        att_left = self._get_closest_att_left(pipolb)
        if att_left is not None:
            pipolin_fragments.append(self._create_pipolin_fragment(att_left, pipolb))

        att_right = self._get_closest_att_right(pipolb)
        if att_right is not None:
            pipolin_fragments.append(self._create_pipolin_fragment(pipolb, att_right))

        pipolin_fragments.append(self._create_pipolin_fragment(pipolb, pipolb))

        return pipolin_fragments

    def _get_closest_att_left(self, pipolb: Feature) -> Optional[AttFeature]:
        contig_atts = self.genome.features.atts_dict()[pipolb.contig_id]
        atts_to_the_left = [att for att in contig_atts if att.is_left_of(pipolb)]
        for att in atts_to_the_left:
            if att.att_id == atts_to_the_left[-1].att_id:
                return att

    def _get_closest_att_right(self, pipolb: Feature) -> Optional[AttFeature]:
        contig_atts = self.genome.features.atts_dict()[pipolb.contig_id]
        atts_to_the_right = [att for att in contig_atts if att.is_right_of(pipolb)]
        for att in atts_to_the_right[::-1]:
            if att.att_id == atts_to_the_right[0].att_id:
                return att

    def _find_atts_around_pipolb(self, pipolb: Feature):
        atts_around_pipolb: Set[AttFeature] = set()
        contig_atts = self.genome.features.atts_dict()[pipolb.contig_id]
        for att in contig_atts:
            if att.start > pipolb.end:
                for other_att in contig_atts:
                    if other_att.end < pipolb.start and other_att.att_id == att.att_id:
                        atts_around_pipolb.update([att, other_att])
        return atts_around_pipolb

    def _find_pipolin_candidates(self, candidates: Sequence[PipolinFragment]):
        pipolin_candidates: MutableSequence[Pipolin] = []
        for fragment in candidates:
            pipolin_candidates.append(Pipolin.from_fragments(fragment))

            other_fragments = self._find_remaining_fragments(fragment)
            pipolin_candidates.append(Pipolin.from_fragments(fragment, *other_fragments))

        return pipolin_candidates

    def _find_remaining_fragments(self, fragment: PipolinFragment) -> Sequence[PipolinFragment]:
        fragment_atts = set(f for f in fragment.features if f.ftype == FeatureType.ATT)

        if not fragment_atts:
            all_atts = []
            for atts in self.atts_by_att_id.values():
                if all([a.contig_id != fragment.contig_id for a in atts]):
                    all_atts.extend(atts)
            return self._fragments_from_orphan_atts(all_atts)

        atts_on_other_contigs = self._get_orphan_atts(fragment_atts)
        return self._fragments_from_orphan_atts(atts_on_other_contigs)

    def _fragments_from_orphan_atts(self, other_fragment_atts):
        other_fragments: MutableSequence[PipolinFragment] = []

        grouped_by_contig_id = groupby(sorted(other_fragment_atts, key=lambda x: x.contig_id),
                                       lambda x: x.contig_id)
        for _, atts_by_contig in grouped_by_contig_id:
            atts = sorted(atts_by_contig, key=lambda x: x.start)
            new_fragment = self._create_pipolin_fragment(atts[0], atts[-1])
            if all(f.ftype != FeatureType.PIPOLB for f in new_fragment.features):
                other_fragments.append(new_fragment)
        return other_fragments

    def _get_orphan_atts(self, fragment_atts):
        orphan_atts = set()
        for att in fragment_atts:
            for other_att in self.atts_by_att_id[att.att_id]:
                if att.contig_id != other_att.contig_id:
                    orphan_atts.add(other_att)
        return orphan_atts

    def _create_pipolin_fragment(self, start_feature: Feature, end_feature: Feature) -> PipolinFragment:
        contig_length = self.genome.get_contig_by_id(start_feature.contig_id).length

        start_feature, end_feature = self._include_ttrnas_if_there_are(start_feature, end_feature)

        loc = Range(max(0, start_feature.start), min(end_feature.end, contig_length))
        fragment = PipolinFragment(location=loc, contig_id=start_feature.contig_id, genome=self.genome)
        fragment_features = fragment.get_fragment_features_sorted()
        return PipolinFragment(fragment.location, fragment.contig_id, fragment.genome, tuple(fragment_features))

    def _include_ttrnas_if_there_are(self, start_feature, end_feature):
        target_trnas = self.genome.features.get_features(FeatureType.TARGET_TRNA)
        ttrna_o_start = target_trnas.get_overlapping(start_feature)
        ttrna_o_end = target_trnas.get_overlapping(end_feature)

        if ttrna_o_start is not None and ttrna_o_start.start < start_feature.start:
            start_feature = ttrna_o_start

        if ttrna_o_end is not None and ttrna_o_end.end > end_feature.end:
            end_feature = ttrna_o_end

        return start_feature, end_feature
