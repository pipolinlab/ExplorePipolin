from typing import Sequence, Tuple

from prefect import task

from explore_pipolin.common import PipolinFragment, FeatureType, Feature, Pipolin, Genome, Range, Strand
from explore_pipolin.utilities.logging import genome_specific_logging


@task()
@genome_specific_logging
def refine_pipolins(genome: Genome, pipolins: Sequence[Pipolin]) -> Sequence[Pipolin]:
    refined_pipolins = []

    for pipolin in pipolins:
        refiner = PipolinRefiner(genome, pipolin)
        refined_pipolins.append(refiner.refine_pipolin())

    return refined_pipolins


class PipolinRefiner:
    def __init__(self, genome: Genome, pipolin: Pipolin):
        self.genome = genome
        self.pipolin = pipolin

        self.pipolb_only_fragments = []
        self.att_only_fragments = []
        self.att_pipolb_fragments = []
        self.pipolb_att_fragments = []
        self.trna_fragments = []

        self._classify_fragments(pipolin.fragments)

    def refine_pipolin(self) -> Pipolin:
        new_fragments_ordered = self._create_new_fragments_ordered()
        return Pipolin.from_fragments(*new_fragments_ordered)

    def _create_new_fragments_ordered(self) -> Sequence[PipolinFragment]:
        pass

    def _cut_fragment_if_necessary(
            self, fragment: PipolinFragment, cut_atts: bool = False) -> PipolinFragment:
        fragment_features = self._get_fragment_features_sorted(fragment)

        first = fragment_features[0]
        second = fragment_features[1]
        next_to_last = fragment_features[-2]
        last = fragment_features[-1]

        if first[1] == FeatureType.ATT and cut_atts:
            fragment_start = max(fragment.start, first[0].start - 50)
        elif first[1] == FeatureType.TARGET_TRNA and second[1] == FeatureType.ATT:
            fragment_start = max(fragment.start, second[0].start - 50)
        else:
            fragment_start = fragment.start

        if last[1] == FeatureType.ATT and cut_atts:
            fragment_end = min(last[0].end + 50, fragment.end)
        elif last[1] == FeatureType.TARGET_TRNA and next_to_last[1] == FeatureType.ATT:
            fragment_end = min(next_to_last[0].end + 50, fragment.end)
        else:
            fragment_end = fragment.end

        return PipolinFragment(Range(fragment_start, fragment_end), contig_id=fragment.contig_id)

    def _get_fragment_features_sorted(
            self, fragment: PipolinFragment) -> Sequence[Tuple[Feature, FeatureType]]:

        fragment_pipolbs = self._get_fragment_features_of_type(fragment, FeatureType.PIPOLB)
        fragment_atts = self._get_fragment_features_of_type(fragment, FeatureType.ATT)
        fragment_ttrnas = self._get_fragment_features_of_type(fragment, FeatureType.TARGET_TRNA)

        fragment_features = [(i, FeatureType.PIPOLB) for i in fragment_pipolbs]
        fragment_features.extend([(i, FeatureType.ATT) for i in fragment_atts])
        fragment_features.extend([(i, FeatureType.TARGET_TRNA) for i in fragment_ttrnas])

        return sorted(fragment_features, key=lambda x: x[0].start)

    def _get_fragment_features_of_type(self, fragment: PipolinFragment, feature_type: FeatureType):
        features = self.genome.features.get_features(feature_type)
        contig_features = features.get_dict_by_contig_sorted()[fragment.contig_id]
        return [f for f in contig_features if f.start >= fragment.start and f.end <= fragment.end]

    def _classify_fragments(self, fragments: Sequence[PipolinFragment]):
        for fragment in fragments:
            features = self._get_fragment_features_sorted(fragment)

            if all([i[1] == FeatureType.PIPOLB for i in features]):
                self.pipolb_only_fragments.append(fragment)
            elif all([i[1] == FeatureType.ATT for i in features]):
                self.att_only_fragments.append(fragment)
            elif any([i[1] == FeatureType.TARGET_TRNA for i in features]):
                self.trna_fragments.append(fragment)
            elif features[0][1] == FeatureType.ATT:
                if features[0][0].strand == Strand.FORWARD:
                    self.att_pipolb_fragments.append(fragment)
                else:
                    self.pipolb_att_fragments.append(fragment)
            elif features[-1][1] == FeatureType.ATT:
                if features[-1][0].strand == Strand.FORWARD:
                    self.pipolb_only_fragments.append(fragment)
                else:
                    self.att_pipolb_fragments.append(fragment)
            else:
                raise AssertionError
