from typing import Sequence

from prefect import task

from explore_pipolin.common import PipolinFragment, FeatureType, Pipolin, Genome, Range, Strand
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

    @staticmethod
    def _cut_fragment_if_necessary(
            fragment: PipolinFragment, cut_atts: bool = False) -> PipolinFragment:

        first = fragment.features[0]
        second = fragment.features[1]
        next_to_last = fragment.features[-2]
        last = fragment.features[-1]

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

    def _classify_fragments(self, fragments: Sequence[PipolinFragment]):
        for fragment in fragments:
            if all([i[1] == FeatureType.PIPOLB for i in fragment.features]):
                self.pipolb_only_fragments.append(fragment)
            elif all([i[1] == FeatureType.ATT for i in fragment.features]):
                self.att_only_fragments.append(fragment)
            elif any([i[1] == FeatureType.TARGET_TRNA for i in fragment.features]):
                self.trna_fragments.append(fragment)
            elif fragment.features[0][1] == FeatureType.ATT:
                if fragment.features[0][0].strand == Strand.FORWARD:
                    self.att_pipolb_fragments.append(fragment)
                else:
                    self.pipolb_att_fragments.append(fragment)
            elif fragment.features[-1][1] == FeatureType.ATT:
                if fragment.features[-1][0].strand == Strand.FORWARD:
                    self.pipolb_only_fragments.append(fragment)
                else:
                    self.att_pipolb_fragments.append(fragment)
            else:
                raise AssertionError
