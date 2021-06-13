from typing import Sequence, MutableSequence, List, Optional

from prefect import task, context
from prefect.utilities.logging import get_logger

from explore_pipolin.common import Genome, Pipolin, FeatureType, PipolinFragment, Strand, Range, AttFeature, \
    PipolinVariants
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
            # as pipolins are sorted from the longest to the shortest, let's skip the sorter ones
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

    def reconstruct_pipolin(self) -> PipolinVariants:
        if len(self.pipolin.fragments) == 1:
            return self._single_fragment(self.pipolin.fragments[0])
        else:
            if self.att_pipolb_att_fragments:
                return self._att_pipolb_att_plus_atts()
            elif self.att_pipolb_fragments:
                return self._att_pipolb_plus_atts()
            elif self.pipolb_only_fragments:
                return self._pipolb_plus_atts()
            else:
                raise AssertionError

    def _single_fragment(self, single_fragment: PipolinFragment) -> PipolinVariants:
        ttrna_fragments = self._get_ttrna_fragments(single_fragment)
        if ttrna_fragments:
            (fragment,) = self._orient_fragments_according_ttrna(ttrna_fragments[0])
            return PipolinVariants.from_variants(self._create_pipolin(single=fragment))
        else:
            fragment1 = self._orient_according_pipolb(single_fragment)
            fragment2 = fragment1.reverse_complement()
            variant1 = self._create_pipolin(single=fragment1)
            variant2 = self._create_pipolin(single=fragment2)
            return PipolinVariants.from_variants(variant1, variant2)

    def _att_pipolb_att_plus_atts(self) -> PipolinVariants:
        att_pipolb_att_fragment = self.att_pipolb_att_fragments[0]
        ttrna_fragments = self._get_ttrna_fragments(att_pipolb_att_fragment, *self.att_only_fragments)

        if len(ttrna_fragments) == 1:
            if ttrna_fragments[0] == att_pipolb_att_fragment:
                # variant: ---att---pol---att(t)
                (att_pipolb_att_fragment,) = self._orient_fragments_according_ttrna(ttrna_fragments[0])
                return PipolinVariants.from_variants(self._create_pipolin(complete=att_pipolb_att_fragment))
            else:
                # variant 1: ---att---pol---att---...---att(t)---   tRNA is required
                (ttrna_fragment, att_pipolb_att_fragment) = self._orient_fragments_according_ttrna(
                    ttrna_fragments[0], att_pipolb_att_fragment
                )
                variant1 = self._create_pipolin(
                    left=att_pipolb_att_fragment, right=ttrna_fragment
                )
                return PipolinVariants.from_variants(variant1)

        # variant 2: ---att---pol---att---
        # variant 3: variant 2 reverse-complement
        fr2 = self._orient_according_pipolb(att_pipolb_att_fragment)
        variant2 = self._create_pipolin(complete=fr2)
        fr3 = fr2.reverse_complement()
        variant3 = self._create_pipolin(complete=fr3)

        return PipolinVariants.from_variants(variant2, variant3)

    def _att_pipolb_plus_atts(self) -> PipolinVariants:
        # we can reconstruct the cases:
        # 1) ---att---pol---...---att---                  tRNA is not required
        # 2) ---att---...---pol---att---...---att(t)---   tRNA is required
        if len(self.att_only_fragments) == 1:
            return self._att_pipolb_plus_one_att()
        elif len(self.att_only_fragments) == 2:
            return self._att_pipolb_plus_two_atts()
        else:
            logger = get_logger(name='reconstruct_pipolins')
            logger.warning(_TOO_MANY_ATTS_MESSAGE)
            return self._single_fragment(self.att_pipolb_fragments[0])

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
                return PipolinVariants.from_variants(self._create_pipolin(left=att_fragment, right=att_pipolb_fragment))
            else:
                # variant: ---att---pol---...---att(t)---
                att_fragment, att_pipolb_fragment = self._orient_fragments_according_ttrna(
                    att_fragment, att_pipolb_fragment
                )
                return PipolinVariants.from_variants(self._create_pipolin(left=att_pipolb_fragment, right=att_fragment))

        elif len(ttrna_fragments) == 2:
            return self._single_fragment(att_pipolb_fragment)

        # variant 1: ---att---pol---...---att---
        # variant 2: ---att---...---pol---att---
        att_pipolb_fragment = self._orient_according_pipolb(att_pipolb_fragment)
        if att_pipolb_fragment.features[0].ftype == FeatureType.PIPOLB:
            left_fragment = att_fragment
            right_fragment = att_pipolb_fragment
            (left_fragment,) = self._orient_fragment_according_main(right_fragment, left_fragment)
        elif att_pipolb_fragment.features[-1].ftype == FeatureType.PIPOLB:
            left_fragment = att_pipolb_fragment
            right_fragment = att_fragment
            (right_fragment,) = self._orient_fragment_according_main(left_fragment, right_fragment)
        else:
            raise AssertionError
        variant1 = self._create_pipolin(left=left_fragment, right=right_fragment)
        variant2 = self._create_pipolin(left=right_fragment.reverse_complement(),
                                        right=left_fragment.reverse_complement())
        return PipolinVariants.from_variants(variant1, variant2)

    def _att_pipolb_plus_two_atts(self) -> PipolinVariants:
        att_pipolb_fragment = self.att_pipolb_fragments[0]
        ttrna_fragments = self._get_ttrna_fragments(att_pipolb_fragment, *self.att_only_fragments)
        if len(ttrna_fragments) == 1:
            if ttrna_fragments[0] == att_pipolb_fragment:
                # ---att---...---att---...---pol---att(t)--- +different order of atts
                (att_pipolb_fragment, att1_fragment, att2_fragment) = self._orient_fragments_according_ttrna(
                    att_pipolb_fragment, self.att_only_fragments[0], self.att_only_fragments[1]
                )
                v1 = self._create_pipolin(left=att1_fragment, middle=att2_fragment, right=att_pipolb_fragment)
                v2 = self._create_pipolin(left=att2_fragment, middle=att1_fragment, right=att_pipolb_fragment)
                return PipolinVariants.from_variants(v1, v2)

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
                    self._create_pipolin(left=left_fragment, middle=middle_fragment, right=right_fragment)
                )
            else:
                # ---att---pol---...---att---...---att(t)---
                # ---att---...---att---pol---...---att(t)---
                return PipolinVariants.from_variants(
                    self._create_pipolin(left=middle_fragment, middle=left_fragment, right=right_fragment),
                    self._create_pipolin(left=left_fragment, middle=middle_fragment, right=right_fragment)
                )
        elif len(ttrna_fragments) >= 2:
            return self._single_fragment(self.att_pipolb_fragments[0])
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
                return PipolinVariants.from_variants(variant1, variant2, variant3, variant4)
            else:
                # ---att---pol---...---att---...---att---   +different order of atts
                # ---att---...---att---pol---...---att---   +different order of atts
                variant1 = self._create_pipolin(left=att_pipolb_fragment, middle=att1_fragment, right=att2_fragment)
                variant2 = self._create_pipolin(left=att_pipolb_fragment, middle=att2_fragment, right=att1_fragment)
                variant3 = self._create_pipolin(left=att1_fragment, middle=att_pipolb_fragment, right=att2_fragment)
                variant4 = self._create_pipolin(left=att2_fragment, middle=att_pipolb_fragment, right=att1_fragment)
                return PipolinVariants.from_variants(variant1, variant2, variant3, variant4)

    def _pipolb_plus_atts(self) -> PipolinVariants:
        # we can try to reconstruct the cases, assuming that pipolb is on the plus strand:
        # 1) ---pol---...---att(t)---                    tRNA is required
        # 2) ---att---...---pol---...---att(t)---        tRNA is required
        if len(self.att_only_fragments) == 1:
            return self._pipolb_plus_one_att()

        elif len(self.att_only_fragments) == 2:
            return self._pipolb_plus_two_atts()

        else:
            logger = get_logger(name='reconstruct_pipolins')
            logger.warning(_TOO_MANY_ATTS_MESSAGE)
            return self._single_fragment(self.pipolb_only_fragments[0])

    def _pipolb_plus_one_att(self) -> PipolinVariants:
        # ---pol---...---att(t)---                    tRNA is required
        ttrna_fragments = self._get_ttrna_fragments(self.att_only_fragments[0])

        pipolb_fragment1 = self._orient_according_pipolb(self.pipolb_only_fragments[0])
        pipolb_fragment2 = pipolb_fragment1.reverse_complement()

        if len(ttrna_fragments) == 1:
            variant1 = self._create_pipolin(middle=pipolb_fragment1, right=ttrna_fragments[0])
            variant2 = self._create_pipolin(middle=pipolb_fragment2, right=ttrna_fragments[0])
            return PipolinVariants.from_variants(variant1, variant2)
        else:
            return PipolinVariants.from_variants(
                self._create_pipolin(left=self.att_only_fragments[0], middle=pipolb_fragment1),
                self._create_pipolin(middle=pipolb_fragment1, right=self.att_only_fragments[0]),
                self._create_pipolin(left=self.att_only_fragments[0], middle=pipolb_fragment2),
                self._create_pipolin(middle=pipolb_fragment2, right=self.att_only_fragments[0])
            )

    def _pipolb_plus_two_atts(self) -> PipolinVariants:
        # ---att---...---pol---...---att(t)---        tRNA is required
        ttrna_fragments = self._get_ttrna_fragments(*self.att_only_fragments)
        pipolb_fragment1 = self._orient_according_pipolb(self.pipolb_only_fragments[0])
        print(pipolb_fragment1.orientation)
        pipolb_fragment2 = pipolb_fragment1.reverse_complement()
        print(pipolb_fragment2.orientation)

        if len(ttrna_fragments) == 1:
            if ttrna_fragments[0] == self.att_only_fragments[0]:
                right_fragment, left_fragment = self._orient_fragments_according_ttrna(
                    ttrna_fragments[0], self.att_only_fragments[1]
                )
            else:
                right_fragment, left_fragment = self._orient_fragments_according_ttrna(
                    ttrna_fragments[0], self.att_only_fragments[0]
                )
            variant1 = self._create_pipolin(left=left_fragment, middle=pipolb_fragment1, right=right_fragment)
            variant2 = self._create_pipolin(left=left_fragment, middle=pipolb_fragment2, right=right_fragment)
            return PipolinVariants.from_variants(variant1, variant2)
        elif len(ttrna_fragments) == 2:
            return self._single_fragment(self.pipolb_only_fragments[0])
        else:
            main_att_fragment = self.att_only_fragments[0]
            (dep_att_fragment,) = self._orient_fragment_according_main(main_att_fragment, self.att_only_fragments[1])
            return PipolinVariants.from_variants(
                self._create_pipolin(left=main_att_fragment, middle=pipolb_fragment1, right=dep_att_fragment),
                self._create_pipolin(left=dep_att_fragment, middle=pipolb_fragment1, right=main_att_fragment),
                self._create_pipolin(left=main_att_fragment, middle=pipolb_fragment2, right=dep_att_fragment),
                self._create_pipolin(left=dep_att_fragment, middle=pipolb_fragment2, right=main_att_fragment)
            )

    # -----------------------------------
    # auxiliary functions
    # -----------------------------------

    @staticmethod
    def _get_ttrna_fragments(*fragments: PipolinFragment) -> Sequence[PipolinFragment]:
        return [f for f in fragments if len(f.get_prime3_ttrnas()) != 0]

    def _orient_fragments_according_ttrna(
            self, ttrna_fragment: PipolinFragment, *other_fragment: PipolinFragment
    ) -> List[PipolinFragment]:
        prime3_ttrnas = ttrna_fragment.get_prime3_ttrnas()
        if len(prime3_ttrnas) == 0:
            raise AssertionError

        # check one is enough as of the same direction
        if prime3_ttrnas[0].strand == Strand.FORWARD:
            ttrna_fragment = ttrna_fragment.reverse_complement()

        new_other_fragments = []
        for ttrna in prime3_ttrnas:
            att = ttrna.get_att_overlapping_ttrna()
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

    def _inflate_single(self, fragment: PipolinFragment) -> PipolinFragment:
        atts_and_pipolbs = _get_fragment_atts_and_pipolbs(fragment)

        left = self._border_inflate if atts_and_pipolbs[0].ftype == FeatureType.ATT else self._no_border_inflate
        right = self._border_inflate if atts_and_pipolbs[-1].ftype == FeatureType.ATT else self._no_border_inflate
        return self._inflate_fragment(fragment, left, right)

    def _create_pipolin(self, left=None, middle=None, right=None, complete=None, single=None) -> Pipolin:
        fragments = []

        if left:
            fragments.append(self._inflate_fragment(left, self._border_inflate, self._no_border_inflate))
        if middle:
            fragments.append(self._inflate_fragment(middle, self._no_border_inflate, self._no_border_inflate))
        if right:
            fragments.append(self._inflate_fragment(right, self._no_border_inflate, self._border_inflate))
        if complete:
            fragments.append(self._inflate_fragment(complete, self._border_inflate, self._border_inflate))
        if single:
            fragments.append(self._inflate_single(single))
        return Pipolin.from_fragments(*fragments)

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
