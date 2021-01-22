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


_BORDER_INFLATE_SIZE = 0
_NO_BORDER_INFLATE_SIZE = 100000


class PipolinRefiner:
    def __init__(self, genome: Genome, pipolin: Pipolin):
        self.genome = genome
        self.pipolin = pipolin

        self.pipolb_only_fragments = []
        self.att_only_fragments = []
        self.att_pipolb_fragments = []
        self.pipolb_att_fragments = []
        self.trna_fragments = []
        self.att_pipolb_att_fragments = []

        self._classify_fragments(pipolin.fragments)

        if len(self.att_pipolb_fragments) > 1 or len(self.pipolb_att_fragments) > 1:
            raise AssertionError('Too complex pipolin. Cannot refine and scaffold it!')

        if len(self.att_pipolb_att_fragments) > 1:
            raise AssertionError('Too complex pipolin. Cannot refine and scaffold it!')

        atts_num = \
            len(self.att_only_fragments) - len(self.att_pipolb_fragments) - \
            len(self.pipolb_att_fragments) + len(self.att_pipolb_att_fragments) * 2

        if atts_num != 0 and len(self.trna_fragments) > 1:
            raise AssertionError('Too complex pipolin. Cannot refine and scaffold it!')

        if len(self.trna_fragments) == 0 and atts_num > 2:
            raise AssertionError('Too complex pipolin. Cannot refine and scaffold it!')

        if self.att_pipolb_att_fragments and self.pipolb_only_fragments:
            raise AssertionError('Too complex pipolin. Cannot refine and scaffold it!')

        if len(self.trna_fragments) > 2:
            raise AssertionError('Too complex pipolin. Cannot refine and scaffold it!')

    def refine_pipolin(self) -> Pipolin:
        fragments_ordered = self._order_fragments()
        fragments_ordered_inflated = self._inflate_fragments(fragments_ordered)
        return Pipolin.from_fragments(*fragments_ordered_inflated)

    def _order_fragments(self) -> Sequence[PipolinFragment]:

        if self.att_pipolb_att_fragments:
            return self.att_pipolb_att_fragments + self.trna_fragments

        elif len(self.att_pipolb_fragments) == 0 and \
                len(self.pipolb_att_fragments) == 0 and \
                len(self.trna_fragments) == 0:

            if len(self.att_only_fragments) == 1:
                return self.att_only_fragments + self.pipolb_only_fragments
            elif len(self.att_only_fragments) == 2:
                return [self.att_only_fragments[0]] + self.pipolb_only_fragments + [self.att_only_fragments[1]]
            else:
                return self.pipolb_only_fragments

        elif not self.att_pipolb_fragments and not self.pipolb_att_fragments:
            if len(self.att_only_fragments) == 2:
                return [self.att_only_fragments[0]] + self.pipolb_only_fragments + \
                       [self.att_only_fragments[1]] + self.trna_fragments
            else:
                return self.att_only_fragments + self.pipolb_only_fragments + self.trna_fragments

        elif len(self.trna_fragments) == 2:
            return [self.trna_fragments[0]] + self.att_pipolb_fragments + self.pipolb_only_fragments + \
                   self.pipolb_att_fragments + [self.trna_fragments[1]]

        elif len(self.pipolb_att_fragments) == 0:
            return self.att_pipolb_fragments + self.pipolb_only_fragments + \
                   self.att_only_fragments + self.trna_fragments

        elif len(self.att_pipolb_fragments) == 0:
            return self.att_only_fragments + self.pipolb_only_fragments + \
                   self.pipolb_att_fragments + self.trna_fragments

        else:
            raise AssertionError('Something else is wrong here?\n'
                                 f'pipolb_only_fragments: {self.pipolb_only_fragments}\n\n'
                                 f'att_only_fragments: {self.att_only_fragments}\n\n'
                                 f'att_pipolb_fragments: {self.att_pipolb_fragments}\n\n'
                                 f'pipolb_att_fragments: {self.pipolb_att_fragments}\n\n'
                                 f'trna_fragments: {self.trna_fragments}')

    def _inflate_fragments(self, fragments_ordered) -> Sequence[PipolinFragment]:
        left_border_inflate, right_border_inflate = self._define_inflate_for_pipolin_borders(fragments_ordered)

        if len(fragments_ordered) == 1:
            return [self._create_pipolin_fragment(fragments_ordered[0],
                                                  left_border_inflate, right_border_inflate)]

        first_fragment = self._create_pipolin_fragment(fragments_ordered[0],
                                                       left_border_inflate, _NO_BORDER_INFLATE_SIZE)
        last_fragment = self._create_pipolin_fragment(fragments_ordered[-1],
                                                      _NO_BORDER_INFLATE_SIZE, right_border_inflate)
        middle_fragments = self._inflate_middle_fragments(fragments_ordered)

        return [first_fragment] + middle_fragments + [last_fragment]

    def _inflate_middle_fragments(self, fragments_ordered):
        middle_fragments = []
        if len(fragments_ordered) > 2:
            for fragment in fragments_ordered[1:-1]:
                middle_fragments.append(
                    self._create_pipolin_fragment(fragment, _NO_BORDER_INFLATE_SIZE, _NO_BORDER_INFLATE_SIZE)
                )
        return middle_fragments

    def _define_inflate_for_pipolin_borders(self, fragments_ordered):
        left_fragments = \
            self.att_only_fragments + self.att_pipolb_fragments + \
            self.trna_fragments + self.att_pipolb_att_fragments
        right_fragment = \
            self.att_only_fragments + self.pipolb_att_fragments + \
            self.trna_fragments + self.att_pipolb_att_fragments

        if fragments_ordered[0] in left_fragments:
            left_border_inflate = _BORDER_INFLATE_SIZE
        else:
            left_border_inflate = _NO_BORDER_INFLATE_SIZE

        if fragments_ordered[-1] in right_fragment:
            right_border_inflate = _BORDER_INFLATE_SIZE
        else:
            right_border_inflate = _NO_BORDER_INFLATE_SIZE

        return left_border_inflate, right_border_inflate

    def _classify_fragments(self, fragments: Sequence[PipolinFragment]):
        for fragment in fragments:

            if all([i[1] == FeatureType.PIPOLB for i in fragment.features]):
                self.pipolb_only_fragments.append(fragment)
            elif all([i[1] == FeatureType.ATT for i in fragment.features]):
                self.att_only_fragments.append(fragment)
            elif any([i[1] == FeatureType.TARGET_TRNA for i in fragment.features]):
                self.trna_fragments.append(fragment)

            elif fragment.features[0][1] == FeatureType.ATT and fragment.features[-1][1] == FeatureType.ATT:
                self.att_pipolb_att_fragments.append(fragment)

            elif fragment.features[0][1] == FeatureType.ATT:
                if fragment.features[0][0].strand == Strand.FORWARD:
                    self.att_pipolb_fragments.append(fragment)
                else:
                    self.pipolb_att_fragments.append(fragment)

            elif fragment.features[-1][1] == FeatureType.ATT:
                if fragment.features[-1][0].strand == Strand.FORWARD:
                    self.pipolb_att_fragments.append(fragment)
                else:
                    self.att_pipolb_fragments.append(fragment)

            else:
                raise AssertionError

    def _create_pipolin_fragment(self, fragment, left_inflate, right_inflate) -> PipolinFragment:
        contig_length = self.genome.get_contig_by_id(fragment.contig_id).length

        loc = Range(max(0, fragment.start - left_inflate),
                    min(fragment.end + right_inflate, contig_length))

        new_fragment = PipolinFragment(location=loc, contig_id=fragment.contig_id, genome=self.genome)
        new_fragment_features = new_fragment.get_fragment_features_sorted()

        return PipolinFragment(new_fragment.location, new_fragment.contig_id,
                               new_fragment.genome, new_fragment_features)
