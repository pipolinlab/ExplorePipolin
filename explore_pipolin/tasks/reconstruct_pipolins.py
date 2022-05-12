from typing import Sequence, MutableSequence, List, Optional

from prefect import task, context
from prefect.utilities.logging import get_logger

from explore_pipolin.common import Genome, Pipolin, FeatureType, PipolinFragment, Strand, Range, AttFeature, \
    PipolinVariants, PipolinType, Feature
import explore_pipolin.settings as settings
from explore_pipolin.utilities.logging import genome_specific_logging


def draw_pipolin_structure(pipolin: Pipolin) -> Sequence[str]:
    return [_get_fragment_string(fragment) for fragment in pipolin.fragments]


def _get_fragment_string(fragment: PipolinFragment) -> str:
    features = fragment.get_fragment_features_sorted()
    if fragment.orientation == Strand.REVERSE:
        features = features[::-1]

    fragment_string = f'contig {fragment.contig_id}: {_get_feature_string(features[0])}'
    for f1, f2 in zip(features, features[1:]):
        if f1.location.is_overlapping(f2.location):
            fragment_string += _get_feature_string(f2)
        else:
            fragment_string += '---' + _get_feature_string(f2)

    return fragment_string


def _get_feature_string(feature) -> str:
    if feature.ftype == FeatureType.PIPOLB:
        return 'pol'
    elif feature.ftype == FeatureType.ATT:
        return 'att' + str(feature.att_id)
    elif feature.ftype == FeatureType.TARGET_TRNA:
        return '(t)'


@task(name='reconstruct_pipolins')
@genome_specific_logging
def reconstruct_pipolins(
        genome: Genome, pipolins: Sequence[Pipolin]
):

    logger = context.get('logger')

    result: MutableSequence[PipolinVariants] = []

    for pipolin in pipolins:
        if not _is_overlapping_result(pipolin, result):
            # as pipolins are sorted from the longest to the shortest, let's skip the sorted ones
            # which will overlap the longer ones after their inflation (result)
            logger.info('>>> Trying to reconstruct the structure from fragments:')
            for structure in draw_pipolin_structure(pipolin):
                logger.info(f'{structure}')

            reconstructor = Reconstructor(genome=genome, pipolin=pipolin)
            pipolin_variants = reconstructor.reconstruct_pipolin()
            logger.info('>>> Reconstruction is done! The resulting pipolin variants are:')
            for variant in pipolin_variants.variants:
                logger.info('...'.join([i.split(sep=': ')[1] for i in draw_pipolin_structure(variant)]))
            result.append(pipolin_variants)

    return result


_TOO_MANY_ATTS_MESSAGE = 'Unable to reconstruct pipolin: too many orphan ATTs. ' \
                         'Only the fragment with piPolB gene will be included.'


def _is_overlapping_result(pipolin: Pipolin, result: Sequence[PipolinVariants]) -> bool:
    for r in result:
        for v in r.variants:
            if pipolin.is_overlapping(v):
                return True
    return False


class Reconstructor:
    def __init__(self, genome: Genome, pipolin: Pipolin):
        self.genome = genome
        self.pipolin = pipolin

        self.att_pipolb_att_fragments: List[PipolinFragment] = []
        self.att_pipolb_fragments: List[PipolinFragment] = []
        self.pipolb_only_fragments: List[PipolinFragment] = []
        self.att_only_fragments: List[PipolinFragment] = []

        self._classify_fragments(pipolin.fragments)

        if len(self.att_pipolb_att_fragments) > 1 or \
                len(self.att_pipolb_fragments) > 1 or \
                len(self.pipolb_only_fragments) > 1:
            raise AssertionError(f'att---pol---att: {len(self.att_pipolb_att_fragments)} > 1 or '
                                 f'att---pol: {len(self.att_pipolb_fragments)} > 1 or '
                                 f'pol: {len(self.pipolb_only_fragments)} > 1')

    @property
    def _border_inflate(self):
        return settings.get_instance().border_inflate

    @property
    def _no_border_inflate(self):
        return settings.get_instance().no_border_inflate

    @property
    def _max_pipolin_len(self):
        return settings.get_instance().max_pipolin_len

    def reconstruct_pipolin(self) -> PipolinVariants:
        if self.att_pipolb_att_fragments:
            return self._att_pipolb_att_plus_atts()
        elif self.att_pipolb_fragments:
            return self._att_pipolb_plus_atts()
        elif self.pipolb_only_fragments:
            return self._pipolb_plus_atts()
        else:
            raise AssertionError

    def _single_fragment(self, single_fragment: PipolinFragment, pipolin_type: PipolinType) -> PipolinVariants:
        ttrna_fragments = self._get_ttrna_fragments(single_fragment)
        if ttrna_fragments:

            if self.is_too_long(ttrna_fragments[0]):
                fragment = self.shorten_the_fragment(ttrna_fragments[0])
                pipolin = Pipolin.from_fragments(fragment)
                return Reconstructor(genome=fragment.genome, pipolin=pipolin).reconstruct_pipolin()

            (fragment,) = self._orient_fragments_according_ttrna(ttrna_fragments[0])
            return PipolinVariants.from_variants(self._create_pipolin(single=fragment, is_ttrna=True),
                                                 pipolin_type=pipolin_type)
        else:

            if self.is_too_long(single_fragment):
                fragment = self.shorten_the_fragment(single_fragment)
                pipolin = Pipolin.from_fragments(fragment)
                return Reconstructor(genome=fragment.genome, pipolin=pipolin).reconstruct_pipolin()

            fragment1 = self._orient_according_pipolb(single_fragment)
            fragment2 = fragment1.reverse_complement()
            variant1 = self._create_pipolin(single=fragment1)
            variant2 = self._create_pipolin(single=fragment2)
            return PipolinVariants.from_variants(variant1, variant2, pipolin_type=pipolin_type)

    def _att_pipolb_att_plus_atts(self) -> PipolinVariants:
        att_pipolb_att_fragment = self.att_pipolb_att_fragments[0]
        ttrna_fragments = self._get_ttrna_fragments(att_pipolb_att_fragment, *self.att_only_fragments)

        if len(ttrna_fragments) != 0 and ttrna_fragments[0] == att_pipolb_att_fragment:
            # variant: ---att---pol---att(t)          # skipping other att(t) fragments if present

            if self.is_too_long(att_pipolb_att_fragment):
                fragment = self.shorten_the_fragment(att_pipolb_att_fragment)
                pipolin = Pipolin.from_fragments(fragment)
                return Reconstructor(genome=self.genome, pipolin=pipolin).reconstruct_pipolin()

            (att_pipolb_att_fragment,) = self._orient_fragments_according_ttrna(ttrna_fragments[0])
            pipolin = self._create_pipolin(complete=att_pipolb_att_fragment, is_ttrna=True)
            return PipolinVariants.from_variants(pipolin, pipolin_type=PipolinType.COMPLETE)
        elif len(ttrna_fragments) == 1:
            # variant 1: ---att---pol---att---...---att(t)---   tRNA is required
            (ttrna_fragment, att_pipolb_att_fragment) = self._orient_fragments_according_ttrna(
                ttrna_fragments[0], att_pipolb_att_fragment
            )
            variant1 = self._create_pipolin(
                left=att_pipolb_att_fragment, right=ttrna_fragment, is_ttrna=True
            )
            return PipolinVariants.from_variants(variant1, pipolin_type=PipolinType.COMPLETE)
        else:

            if self.is_too_long(att_pipolb_att_fragment):
                fragment = self.shorten_the_fragment(att_pipolb_att_fragment)
                pipolin = Pipolin.from_fragments(fragment)
                return Reconstructor(genome=fragment.genome, pipolin=pipolin).reconstruct_pipolin()

            # variant 1: ---att---pol---att---          # skip additional att fragments if present
            # variant 2: variant 1 reverse-complement
            fr1 = self._orient_according_pipolb(att_pipolb_att_fragment)
            variant1 = self._create_pipolin(complete=fr1)
            fr2 = fr1.reverse_complement()
            variant2 = self._create_pipolin(complete=fr2)

            return PipolinVariants.from_variants(variant1, variant2, pipolin_type=PipolinType.COMPLETE)

    def _att_pipolb_plus_atts(self) -> PipolinVariants:
        # we can reconstruct the cases:
        # 1) ---att---pol---...---att---                  tRNA is not required
        # 2) ---att---...---pol---att---...---att(t)---   tRNA is required
        # 3) ---pol---att/att(t)--- return correctly truncated
        if len(self.att_only_fragments) == 1:
            return self._att_pipolb_plus_one_att()
        elif len(self.att_only_fragments) == 2:
            return self._att_pipolb_plus_two_atts()
        elif len(self.att_only_fragments) == 0:

            if self.is_too_long(self.att_pipolb_fragments[0]):
                fragment = self.shorten_the_fragment(self.att_pipolb_fragments[0])
                pipolin = Pipolin.from_fragments(fragment)
                return Reconstructor(genome=fragment.genome, pipolin=pipolin).reconstruct_pipolin()

            ttrna_fragments = self._get_ttrna_fragments(self.att_pipolb_fragments[0])
            if ttrna_fragments:
                pipolb_att_ttrna_fragment = self._orient_fragments_according_ttrna(self.att_pipolb_fragments[0])
                return PipolinVariants.from_variants(self._create_pipolin(right=pipolb_att_ttrna_fragment[0],
                                                                          is_ttrna=True),
                                                     pipolin_type=PipolinType.TRUNCATED)
            else:
                fragment1 = self._orient_according_pipolb(self.att_pipolb_fragments[0])
                fragment2 = fragment1.reverse_complement()
                variant1 = self._create_pipolin(single=fragment1)
                variant2 = self._create_pipolin(single=fragment2)
                return PipolinVariants.from_variants(variant1, variant2, pipolin_type=PipolinType.TRUNCATED)
        else:
            logger = get_logger(name='reconstruct_pipolins')
            logger.warning(_TOO_MANY_ATTS_MESSAGE)
            return self._single_fragment(self.att_pipolb_fragments[0], PipolinType.TRUNCATED)

    def _att_pipolb_plus_one_att(self) -> PipolinVariants:
        att_pipolb_fragment = self.att_pipolb_fragments[0]
        att_fragment = self.att_only_fragments[0]
        ttrna_fragments = self._get_ttrna_fragments(att_pipolb_fragment, att_fragment)

        if len(ttrna_fragments) == 1:
            if ttrna_fragments[0] == att_pipolb_fragment:
                # variant: ---att---...---pol---att(t)---
                att_pipolb_fragment, att_fragment = self._orient_fragments_according_ttrna(
                    att_pipolb_fragment, att_fragment
                )
                return PipolinVariants.from_variants(self._create_pipolin(left=att_fragment, right=att_pipolb_fragment,
                                                                          is_ttrna=True),
                                                     pipolin_type=PipolinType.COMPLETE)
            else:
                # variant: ---att---pol---...---att(t)---
                att_fragment, att_pipolb_fragment = self._orient_fragments_according_ttrna(
                    att_fragment, att_pipolb_fragment
                )
                return PipolinVariants.from_variants(self._create_pipolin(left=att_pipolb_fragment, right=att_fragment,
                                                                          is_ttrna=True),
                                                     pipolin_type=PipolinType.COMPLETE)

        elif len(ttrna_fragments) == 2:
            # ambiguous case, drop att fragment
            return self._single_fragment(att_pipolb_fragment, PipolinType.TRUNCATED)

        # variant 1: ---att---pol---...---att---
        # variant 2: ---att---...---pol---att---
        att_pipolb_fragment = self._orient_according_pipolb(att_pipolb_fragment)
        feature_0 = att_pipolb_fragment.features[0]
        feature_m1 = att_pipolb_fragment.features[-1]

        is_0_forward = (feature_0.ftype == FeatureType.PIPOLB) and (att_pipolb_fragment.orientation == Strand.FORWARD)
        is_0_reverse = (feature_0.ftype == FeatureType.PIPOLB) and (att_pipolb_fragment.orientation == Strand.REVERSE)
        is_m1_forward = (feature_m1.ftype == FeatureType.PIPOLB) and (att_pipolb_fragment.orientation == Strand.FORWARD)
        is_m1_reverse = (feature_m1.ftype == FeatureType.PIPOLB) and (att_pipolb_fragment.orientation == Strand.REVERSE)

        if is_0_forward or is_m1_reverse:
            left_fragment = att_fragment
            right_fragment = att_pipolb_fragment
            (left_fragment,) = self._orient_fragment_according_main(right_fragment, left_fragment)
        elif is_0_reverse or is_m1_forward:
            left_fragment = att_pipolb_fragment
            right_fragment = att_fragment
            (right_fragment,) = self._orient_fragment_according_main(left_fragment, right_fragment)
        else:
            raise AssertionError
        variant1 = self._create_pipolin(left=left_fragment, right=right_fragment)
        variant2 = self._create_pipolin(left=right_fragment.reverse_complement(),
                                        right=left_fragment.reverse_complement())
        return PipolinVariants.from_variants(variant1, variant2, pipolin_type=PipolinType.COMPLETE)

    def _att_pipolb_plus_two_atts(self) -> PipolinVariants:
        att_pipolb_fragment = self.att_pipolb_fragments[0]
        ttrna_fragments = self._get_ttrna_fragments(att_pipolb_fragment, *self.att_only_fragments)
        if len(ttrna_fragments) == 1:
            if ttrna_fragments[0] == att_pipolb_fragment:
                # ---att---...---att---...---pol---att(t)--- +different order of atts
                (att_pipolb_fragment, att1_fragment, att2_fragment) = self._orient_fragments_according_ttrna(
                    att_pipolb_fragment, self.att_only_fragments[0], self.att_only_fragments[1]
                )
                v1 = self._create_pipolin(left=att1_fragment, middle=att2_fragment, right=att_pipolb_fragment,
                                          is_ttrna=True)
                v2 = self._create_pipolin(left=att2_fragment, middle=att1_fragment, right=att_pipolb_fragment,
                                          is_ttrna=True)
                return PipolinVariants.from_variants(v1, v2, pipolin_type=PipolinType.COMPLETE)

            elif ttrna_fragments[0] == self.att_only_fragments[0]:
                (right_fragment, middle_fragment, left_fragment) = self._orient_fragments_according_ttrna(
                    self.att_only_fragments[0], self.att_pipolb_fragments[0], self.att_only_fragments[1]
                )
            else:
                (right_fragment, middle_fragment, left_fragment) = self._orient_fragments_according_ttrna(
                    self.att_only_fragments[1], self.att_pipolb_fragments[0], self.att_only_fragments[0]
                )

            if middle_fragment.orientation == Strand.FORWARD:
                pipolb_is_left = middle_fragment.features[0].ftype == FeatureType.PIPOLB
            else:
                pipolb_is_left = middle_fragment.features[-1].ftype == FeatureType.PIPOLB

            if pipolb_is_left:
                # ---att---...---pol---att---...---att(t)---
                return PipolinVariants.from_variants(
                    self._create_pipolin(left=left_fragment, middle=middle_fragment, right=right_fragment,
                                         is_ttrna=True),
                    pipolin_type=PipolinType.COMPLETE
                )
            else:
                # ---att---pol---...---att---...---att(t)---
                # ---att---...---att---pol---...---att(t)---
                return PipolinVariants.from_variants(
                    self._create_pipolin(left=middle_fragment, middle=left_fragment, right=right_fragment,
                                         is_ttrna=True),
                    self._create_pipolin(left=left_fragment, middle=middle_fragment, right=right_fragment,
                                         is_ttrna=True),
                    pipolin_type=PipolinType.COMPLETE
                )
        elif len(ttrna_fragments) >= 2:
            # an ambiguous case, leave the main fragment only
            return self._single_fragment(self.att_pipolb_fragments[0], PipolinType.TRUNCATED)
        else:
            (att1_fragment, att2_fragment) = self._orient_fragment_according_main(
                att_pipolb_fragment, *self.att_only_fragments
            )
            if att_pipolb_fragment.orientation == Strand.FORWARD:
                pipolb_is_left = att_pipolb_fragment.features[0].ftype == FeatureType.PIPOLB
            else:
                pipolb_is_left = att_pipolb_fragment.features[-1].ftype == FeatureType.PIPOLB

            if pipolb_is_left:
                # ---att---...---pol---att---...---att---   +different order of atts
                # ---att---...---att---...---pol---att---   +different order of atts
                variant1 = self._create_pipolin(left=att1_fragment, middle=att_pipolb_fragment, right=att2_fragment)
                variant2 = self._create_pipolin(left=att2_fragment, middle=att_pipolb_fragment, right=att1_fragment)
                variant3 = self._create_pipolin(left=att1_fragment, middle=att2_fragment, right=att_pipolb_fragment)
                variant4 = self._create_pipolin(left=att2_fragment, middle=att1_fragment, right=att_pipolb_fragment)
                return PipolinVariants.from_variants(variant1, variant2, variant3, variant4,
                                                     pipolin_type=PipolinType.COMPLETE)
            else:
                # ---att---pol---...---att---...---att---   +different order of atts
                # ---att---...---att---pol---...---att---   +different order of atts
                variant1 = self._create_pipolin(left=att_pipolb_fragment, middle=att1_fragment, right=att2_fragment)
                variant2 = self._create_pipolin(left=att_pipolb_fragment, middle=att2_fragment, right=att1_fragment)
                variant3 = self._create_pipolin(left=att1_fragment, middle=att_pipolb_fragment, right=att2_fragment)
                variant4 = self._create_pipolin(left=att2_fragment, middle=att_pipolb_fragment, right=att1_fragment)
                return PipolinVariants.from_variants(variant1, variant2, variant3, variant4,
                                                     pipolin_type=PipolinType.COMPLETE)

    def _pipolb_plus_atts(self) -> PipolinVariants:
        # we can try to reconstruct the cases, assuming that pipolb is on the plus strand:
        # 1) ---pol---...---att(t)---                    tRNA is required
        # 2) ---att---...---pol---...---att(t)---        tRNA is required
        if len(self.att_only_fragments) == 1:
            return self._pipolb_plus_one_att()

        elif len(self.att_only_fragments) == 2:
            return self._pipolb_plus_two_atts()

        else:
            if len(self.att_only_fragments) > 2:
                logger = get_logger(name='reconstruct_pipolins')
                logger.warning(_TOO_MANY_ATTS_MESSAGE)
            return self._single_fragment(self.pipolb_only_fragments[0], PipolinType.MINIMAL)

    def _pipolb_plus_one_att(self) -> PipolinVariants:
        # ---pol---...---att(t)---                    tRNA is required
        ttrna_fragments = self._get_ttrna_fragments(self.att_only_fragments[0])

        pipolb_fragment1 = self._orient_according_pipolb(self.pipolb_only_fragments[0])
        pipolb_fragment2 = pipolb_fragment1.reverse_complement()

        if len(ttrna_fragments) == 1:
            variant1 = self._create_pipolin(middle=pipolb_fragment1, right=ttrna_fragments[0], is_ttrna=True)
            variant2 = self._create_pipolin(middle=pipolb_fragment2, right=ttrna_fragments[0], is_ttrna=True)
            return PipolinVariants.from_variants(variant1, variant2, pipolin_type=PipolinType.TRUNCATED)
        else:
            return PipolinVariants.from_variants(
                self._create_pipolin(left=self.att_only_fragments[0], middle=pipolb_fragment1),
                self._create_pipolin(middle=pipolb_fragment1, right=self.att_only_fragments[0]),
                self._create_pipolin(left=self.att_only_fragments[0], middle=pipolb_fragment2),
                self._create_pipolin(middle=pipolb_fragment2, right=self.att_only_fragments[0]),
                pipolin_type=PipolinType.TRUNCATED
            )

    def _pipolb_plus_two_atts(self) -> PipolinVariants:
        # ---att---...---pol---...---att(t)---        tRNA is required
        ttrna_fragments = self._get_ttrna_fragments(*self.att_only_fragments)
        pipolb_fragment1 = self._orient_according_pipolb(self.pipolb_only_fragments[0])
        pipolb_fragment2 = pipolb_fragment1.reverse_complement()

        if len(ttrna_fragments) == 1:
            if ttrna_fragments[0] == self.att_only_fragments[0]:
                right_fragment, left_fragment = self._orient_fragments_according_ttrna(
                    ttrna_fragments[0], self.att_only_fragments[1]
                )
            else:
                right_fragment, left_fragment = self._orient_fragments_according_ttrna(
                    ttrna_fragments[0], self.att_only_fragments[0]
                )
            variant1 = self._create_pipolin(left=left_fragment, middle=pipolb_fragment1, right=right_fragment,
                                            is_ttrna=True)
            variant2 = self._create_pipolin(left=left_fragment, middle=pipolb_fragment2, right=right_fragment,
                                            is_ttrna=True)
            return PipolinVariants.from_variants(variant1, variant2, pipolin_type=PipolinType.COMPLETE)
        elif len(ttrna_fragments) == 2:
            # an ambiguous case, drop att fragments
            return self._single_fragment(self.pipolb_only_fragments[0], PipolinType.MINIMAL)
        else:
            main_att_fragment = self.att_only_fragments[0]
            (dep_att_fragment,) = self._orient_fragment_according_main(main_att_fragment, self.att_only_fragments[1])
            return PipolinVariants.from_variants(
                self._create_pipolin(left=main_att_fragment, middle=pipolb_fragment1, right=dep_att_fragment),
                self._create_pipolin(left=dep_att_fragment, middle=pipolb_fragment1, right=main_att_fragment),
                self._create_pipolin(left=main_att_fragment, middle=pipolb_fragment2, right=dep_att_fragment),
                self._create_pipolin(left=dep_att_fragment, middle=pipolb_fragment2, right=main_att_fragment),
                pipolin_type=PipolinType.COMPLETE
            )

    # -----------------------------------
    # auxiliary functions
    # -----------------------------------

    @staticmethod
    def _get_ttrna_fragments(*fragments: PipolinFragment) -> Sequence[PipolinFragment]:
        return [f for f in fragments if len(f.get_ttrnas_outside_fragment()) != 0]

    def is_too_long(self, *fragment: PipolinFragment) -> bool:
        total_len = sum([abs(f.location.end - f.location.start) for f in fragment])
        return total_len > self._max_pipolin_len

    def _orient_fragments_according_ttrna(
            self, ttrna_fragment: PipolinFragment, *other_fragment: PipolinFragment
    ) -> List[PipolinFragment]:
        prime3_ttrna = ttrna_fragment.get_ttrnas_outside_fragment()[0]

        # check one is enough as of the same direction
        if prime3_ttrna.strand == Strand.FORWARD:
            ttrna_fragment = ttrna_fragment.reverse_complement()

        new_other_fragments = []
        att = prime3_ttrna.get_att_overlapping_ttrna()
        other_fragment_strands = [self._get_fragment_atts_strand(f, att.att_id) for f in other_fragment]

        for f, s in zip(other_fragment, other_fragment_strands):
            if s:
                if ttrna_fragment.orientation == f.orientation and att.strand != s:
                    new_other_fragments.append(f.reverse_complement())
                elif ttrna_fragment.orientation != f.orientation and att.strand == s:
                    new_other_fragments.append(f.reverse_complement())
                else:
                    new_other_fragments.append(f)
            else:
                new_other_fragments.append(f)

        return [ttrna_fragment] + new_other_fragments

    def _inflate_fragment(self, fragment: PipolinFragment, left: int, right: int) -> PipolinFragment:
        contig_length = self.genome.get_contig_by_id(fragment.contig_id).length

        if fragment.orientation == Strand.FORWARD:
            loc = Range(max(0, fragment.start - left), min(fragment.end + right, contig_length))
        else:
            loc = Range(max(0, fragment.start - right), min(fragment.end + left, contig_length))

        new = PipolinFragment(loc, fragment.contig_id, fragment.genome, orientation=fragment.orientation)
        new_features = new.get_fragment_features_sorted()
        return PipolinFragment(new.location, new.contig_id, new.genome, tuple(new_features), new.orientation)

    def _inflate_single(self, fragment: PipolinFragment, is_ttrna: bool) -> PipolinFragment:
        if is_ttrna:
            fragment = self.ensure_ttrna_edge(fragment)

        atts_and_pipolbs = _get_fragment_atts_and_pipolbs(fragment)

        is_0_forward = atts_and_pipolbs[0].ftype == FeatureType.ATT and fragment.orientation == Strand.FORWARD
        is_0_reverse = atts_and_pipolbs[0].ftype == FeatureType.ATT and fragment.orientation == Strand.REVERSE
        is_m1_forward = atts_and_pipolbs[-1].ftype == FeatureType.ATT and fragment.orientation == Strand.FORWARD
        is_m1_reverse = atts_and_pipolbs[-1].ftype == FeatureType.ATT and fragment.orientation == Strand.REVERSE
        left = self._border_inflate if (is_0_forward or is_m1_reverse) else self._no_border_inflate
        right = self._border_inflate if (is_0_reverse or is_m1_forward) else self._no_border_inflate
        return self._inflate_fragment(fragment, left, right)

    def _create_pipolin(
            self, left=None, middle=None, right=None, complete=None, single=None, is_ttrna=False
    ) -> Pipolin:
        fragments = []

        if left:
            fragments.append(self._inflate_fragment(left, self._border_inflate, self._no_border_inflate))
        if middle:
            fragments.append(self._inflate_fragment(middle, self._no_border_inflate, self._no_border_inflate))
        if right:
            if is_ttrna:
                right = self.ensure_ttrna_edge(right)
            fragments.append(self._inflate_fragment(right, self._no_border_inflate, self._border_inflate))
        if complete:
            if is_ttrna:
                complete = self.ensure_ttrna_edge(complete)
            fragments.append(self._inflate_fragment(complete, self._border_inflate, self._border_inflate))
        if single:
            fragments.append(self._inflate_single(single, is_ttrna))
        return Pipolin.from_fragments(*fragments)

    @staticmethod
    def create_new_features_fragment(
            old_fragment: PipolinFragment, new_features: Sequence[Feature],
    ) -> PipolinFragment:
        new_genome = Genome(genome_id=old_fragment.genome.id,
                            genome_file=old_fragment.genome.file,
                            contigs=old_fragment.genome.contigs)
        new_genome.features.add_features(*new_features)

        contig_length = old_fragment.genome.get_contig_by_id(old_fragment.contig_id).length
        new_loc = Range(
            max(0, new_features[0].location.start),
            min(new_features[-1].location.end, contig_length)
        )
        new_fragment = PipolinFragment(location=new_loc,
                                       contig_id=old_fragment.contig_id,
                                       genome=new_genome,
                                       orientation=old_fragment.orientation)
        new_fragment_features = new_fragment.get_fragment_features_sorted()
        new_fragment = PipolinFragment(new_fragment.location,
                                       new_fragment.contig_id,
                                       new_fragment.genome,
                                       tuple(new_fragment_features),
                                       new_fragment.orientation)
        return new_fragment

    def ensure_ttrna_edge(self, fragment: PipolinFragment) -> PipolinFragment:
        ttrnas = [f for f in fragment.features if f.ftype == FeatureType.TARGET_TRNA]
        if fragment.orientation == Strand.FORWARD:
            new_features = [f for f in fragment.features if f.location.end <= ttrnas[-1].location.end]
        else:
            new_features = [f for f in fragment.features if f.location.start >= ttrnas[0].location.start]

        return self.create_new_features_fragment(fragment, new_features)

    def shorten_the_fragment(self, fragment: PipolinFragment) -> PipolinFragment:
        if fragment.features[0].location.is_overlapping(fragment.features[1].location):
            start = 1
        else:
            start = 0

        if fragment.features[-1].location.is_overlapping(fragment.features[-2].location):
            end = -2
        else:
            end = -1

        diff_left = fragment.features[start+1].location.start - fragment.features[start].location.start
        diff_right = fragment.features[end].location.end - fragment.features[end-1].location.end

        # NEVER cut piPolB!
        if fragment.features[start].ftype == FeatureType.PIPOLB and fragment.features[end].ftype == FeatureType.PIPOLB:
            raise AssertionError

        if diff_left > diff_right and fragment.features[start].ftype != FeatureType.PIPOLB:
            new_features = fragment.features[start+1:]
        else:
            if fragment.features[end].ftype != FeatureType.PIPOLB:
                new_features = fragment.features[:end]
            else:
                new_features = fragment.features[start:]

        return self.create_new_features_fragment(fragment, new_features)

    @staticmethod
    def _orient_according_pipolb(pipolb_containing_fragment: PipolinFragment) -> PipolinFragment:
        pipolbs = pipolb_containing_fragment.get_fragment_features_of_type_sorted(FeatureType.PIPOLB)
        pipolbs_strands = set(i.strand for i in pipolbs)
        if len(pipolbs_strands) == 1:
            if pipolbs_strands.pop() == Strand.REVERSE:
                return pipolb_containing_fragment.reverse_complement()
        return pipolb_containing_fragment

    def _orient_fragment_according_main(
            self, main_fragment: PipolinFragment, *dependent_fragment: PipolinFragment) -> List[PipolinFragment]:
        att_ids_sets = [set(m.att_id for m in main_fragment.features if isinstance(m, AttFeature))]
        for fragment in dependent_fragment:
            att_ids_sets.append(set(d.att_id for d in fragment.features if isinstance(d, AttFeature)))
        att_ids_on_all_fragments = set()
        for i, s in enumerate(att_ids_sets[:-1]):
            att_ids_on_all_fragments = s.intersection(att_ids_sets[i + 1])
        if len(att_ids_on_all_fragments) == 0:
            raise AssertionError
        att_id = att_ids_on_all_fragments.pop()   # one is enough

        main_fragment_strand = self._get_fragment_atts_strand(main_fragment, att_id)
        dependent_fragment_strands = [self._get_fragment_atts_strand(d, att_id) for d in dependent_fragment]
        new_dependent_fragments = []
        for i, strand in enumerate(dependent_fragment_strands):
            if main_fragment_strand != strand:
                new_dependent_fragments.append(dependent_fragment[i].reverse_complement())
            else:
                new_dependent_fragments.append(dependent_fragment[i])
        return new_dependent_fragments

    @staticmethod
    def _get_fragment_atts_strand(fragment: PipolinFragment, att_id) -> Optional[Strand]:
        strands = set(f.strand for f in fragment.features if isinstance(f, AttFeature) and f.att_id == att_id)
        if len(strands) != 1:
            raise AssertionError
        return strands.pop()

    def _classify_fragments(self, fragments: Sequence[PipolinFragment]):
        for fragment in fragments:
            fragment_atts_and_pipolbs = _get_fragment_atts_and_pipolbs(fragment)

            if all([i.ftype == FeatureType.PIPOLB for i in fragment_atts_and_pipolbs]):
                self.pipolb_only_fragments.append(fragment)
            elif all([i.ftype == FeatureType.ATT for i in fragment_atts_and_pipolbs]):
                self.att_only_fragments.append(fragment)

            elif fragment_atts_and_pipolbs[0].ftype == FeatureType.ATT and \
                    fragment_atts_and_pipolbs[-1].ftype == FeatureType.ATT:
                self.att_pipolb_att_fragments.append(fragment)
            else:
                self.att_pipolb_fragments.append(fragment)


def _get_fragment_atts_and_pipolbs(fragment: PipolinFragment):
    features = fragment.get_fragment_features_sorted()
    return [f for f in features if f.ftype != FeatureType.TARGET_TRNA]
