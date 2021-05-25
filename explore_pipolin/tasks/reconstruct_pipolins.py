from itertools import chain, groupby
from typing import Sequence, MutableSequence, List, Optional

from prefect import task, context

from explore_pipolin.common import Genome, Pipolin, FeatureType, PipolinFragment, Strand, Feature, Range, AttFeature
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


def _rev_comp_fragment(f: PipolinFragment) -> PipolinFragment:
    return PipolinFragment(f.location, f.contig_id, f.genome, f.features, -f.orientation)


@task()
@genome_specific_logging
def reconstruct_pipolins(genome: Genome, pipolins: Sequence[Pipolin]):

    logger = context.get('logger')

    result: MutableSequence[Sequence[Pipolin]] = []

    for pipolin in pipolins:
        logger.info('>>> Trying to reconstruct the structure from fragments:')
        logger.info(''.join(f'\n\t{i}' for i in draw_pipolin_structure(pipolin)))

        try:
            reconstructor = Reconstructor(genome=genome, pipolin=pipolin)
            pipolin_variants = reconstructor.reconstruct_pipolin()
            logger.info('>>> Reconstruction is done! The resulting pipolin variants are:')
            for variant in pipolin_variants:
                logger.info('...'.join([i.split(sep=': ')[1] for i in draw_pipolin_structure(variant)]))
            result.append(pipolin_variants)

        except CannotReconstructError as e:
            logger.warning(f'>>> Cannot reconstruct the structure! {e}')

    filtered = _filter_redundant(result) if len(result) > 1 else result
    logger.info('Pipolins after filtering:')
    for pipolin_variants in filtered:
        for variant in pipolin_variants:
            logger.info('...'.join([i.split(sep=': ')[1] for i in draw_pipolin_structure(variant)]))
    return filtered


class CannotReconstructError(Exception):
    pass


_BORDER_INFLATE = 0
_NO_BORDER_INFLATE = 100_000


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
            raise CannotReconstructError(f'att---pol---att: {len(self.att_pipolb_att_fragments)} > 1 or '
                                         f'att---pol: {len(self.att_pipolb_fragments)} > 1 or '
                                         f'pol: {len(self.pipolb_only_fragments)} > 1')
        # TODO: check this condition if it still holds
        # If >=2 contigs with tRNAs or a contig with >=2 tRNAs on different strands -> raise CannotReconstructError
        self._check_ttrnas(pipolin.fragments)

    def reconstruct_pipolin(self) -> Sequence[Pipolin]:
        if len(self.pipolin.fragments) == 1:
            return self._single_fragment()
        else:
            if self.att_pipolb_att_fragments:
                return self._att_pipolb_att_plus_atts()
            elif self.att_pipolb_fragments:
                return self._att_pipolb_plus_atts()
            elif self.pipolb_only_fragments:
                return self._pipolb_plus_atts()
            else:
                raise AssertionError

    def _single_fragment(self) -> Sequence[Pipolin]:
        ttrna_fragment = self._get_single_ttrna_containing_fragment([self.pipolin.fragments[0]])
        if ttrna_fragment:
            (fragment,) = self._orient_fragments_according_ttrnas(ttrna_fragment)
            return [self._inflate_single_fragment_pipolin(fragment)]
        else:
            fragment1 = self._orient_according_pipolb(self.pipolin.fragments[0])
            fragment2 = _rev_comp_fragment(fragment1)
            variant1 = self._inflate_single_fragment_pipolin(fragment1)
            variant2 = self._inflate_single_fragment_pipolin(fragment2)
            return [variant1, variant2]

    def _att_pipolb_att_plus_atts(self) -> Sequence[Pipolin]:
        att_pipolb_att_fragment = self.att_pipolb_att_fragments[0]
        ttrna_fragment = self._get_single_ttrna_containing_fragment(
            [att_pipolb_att_fragment] + self.att_only_fragments
        )
        if ttrna_fragment and ttrna_fragment == att_pipolb_att_fragment:
            # variant: ---att---pol---att(t)
            ttrna_fragment = self._orient_fragments_according_ttrnas(ttrna_fragment)[0]
            return [self._create_pipolin(complete=ttrna_fragment)]
        else:
            pipolin_variants = []
            # variant 1: ---att---pol---att---
            # variant 2: variant 1 reverse-complement
            # variant 3: ---att---pol---att---...---att(t)---   tRNA is required
            fr1 = self._orient_according_pipolb(att_pipolb_att_fragment)
            variant1 = self._create_pipolin(complete=fr1)
            fr2 = _rev_comp_fragment(fr1)
            variant2 = self._create_pipolin(complete=fr2)

            pipolin_variants.extend([variant1, variant2])

        if ttrna_fragment:
            try:
                (new_ttrna_fragment, new_att_pipolb_att_fragment) = self._orient_fragments_according_ttrnas(
                    ttrna_fragment, att_pipolb_att_fragment
                )
                variant3 = self._create_pipolin(
                    left=new_att_pipolb_att_fragment, right=new_ttrna_fragment
                )
                pipolin_variants.append(variant3)
            except CannotReconstructError:
                pass

        return pipolin_variants

    def _att_pipolb_plus_atts(self) -> Sequence[Pipolin]:
        # we can reconstruct the cases:
        # 1) ---att---pol---...---att---                  tRNA is not required
        # 2) ---att---...---pol---att---...---att(t)---   tRNA is required
        if len(self.att_only_fragments) == 1:
            return self._att_pipolb_plus_one_att()
        elif len(self.att_only_fragments) == 2:
            return self._att_pipolb_plus_two_atts()
        else:
            raise CannotReconstructError

    def _att_pipolb_plus_one_att(self) -> Sequence[Pipolin]:
        att_pipolb_fragment = self.att_pipolb_fragments[0]
        att_fragment = self.att_only_fragments[0]
        ttrna_fragment = self._get_single_ttrna_containing_fragment([att_pipolb_fragment, att_fragment])

        if ttrna_fragment:
            if ttrna_fragment == att_pipolb_fragment:
                # variant: ---att---...---pol---att(t)---
                att_pipolb_fragment, att_fragment = self._orient_fragments_according_ttrnas(
                    att_pipolb_fragment, att_fragment
                )
                return [self._create_pipolin(left=att_fragment, right=att_pipolb_fragment)]
            else:
                # variant: ---att---pol---...---att(t)---
                att_fragment, att_pipolb_fragment = self._orient_fragments_according_ttrnas(
                    att_fragment, att_pipolb_fragment
                )
                return [self._create_pipolin(left=att_pipolb_fragment, right=att_fragment)]

        # variant 1: ---att---pol---...---att---
        # variant 2: ---att---...---pol---att---
        att_pipolb_fragment = self._orient_according_pipolb(att_pipolb_fragment)
        if att_pipolb_fragment.features[0].ftype == FeatureType.PIPOLB:
            left_fragment = att_fragment
            right_fragment = att_pipolb_fragment
            left_fragment = self._orient_fragment_according_main(right_fragment, left_fragment)
        elif att_pipolb_fragment.features[-1].ftype == FeatureType.PIPOLB:
            left_fragment = att_pipolb_fragment
            right_fragment = att_fragment
            right_fragment = self._orient_fragment_according_main(left_fragment, right_fragment)
        else:
            raise AssertionError
        variant1 = self._create_pipolin(left=left_fragment, right=right_fragment)
        variant2 = self._create_pipolin(left=_rev_comp_fragment(right_fragment),
                                        right=_rev_comp_fragment(left_fragment))
        return [variant1, variant2]

    def _att_pipolb_plus_two_atts(self) -> Sequence[Pipolin]:
        # 2) ---att---...---pol---att---...---att(t)---   tRNA is required
        ttrna_fragment = self._get_single_ttrna_containing_fragment(
            [self.att_only_fragments[0], self.att_only_fragments[1]]
        )
        if ttrna_fragment:
            if ttrna_fragment == self.att_only_fragments[0]:
                (right_fragment, middle_fragment, left_fragment) = self._orient_fragments_according_ttrnas(
                    ttrna_fragment, self.att_pipolb_fragments[0], self.att_only_fragments[1]
                )
            else:
                (right_fragment, middle_fragment, left_fragment) = self._orient_fragments_according_ttrnas(
                    ttrna_fragment, self.att_pipolb_fragments[0], self.att_only_fragments[0]
                )

            if middle_fragment.orientation == Strand.FORWARD:
                pipolb_is_left = middle_fragment.features[0].ftype == FeatureType.PIPOLB
            else:
                pipolb_is_left = middle_fragment.features[-1].ftype == FeatureType.PIPOLB

            if pipolb_is_left:
                return [self._create_pipolin(left=left_fragment, middle=middle_fragment, right=right_fragment)]
            else:
                raise CannotReconstructError('---att---...---att---pol---...---att(t)---')
        else:
            raise CannotReconstructError

    def _pipolb_plus_atts(self) -> Sequence[Pipolin]:
        # we can try to reconstruct the cases, assuming that pipolb is on the plus strand:
        # 1) ---pol---...---att(t)---                    tRNA is required
        # 2) ---att---...---pol---...---att(t)---        tRNA is required
        if len(self.att_only_fragments) == 1:
            return self._pipolb_plus_one_att()

        elif len(self.att_only_fragments) == 2:
            return self._pipolb_plus_two_atts()

        else:
            raise CannotReconstructError

    def _pipolb_plus_one_att(self) -> Sequence[Pipolin]:
        # ---pol---...---att(t)---                    tRNA is required
        ttrna_fragment = self._get_single_ttrna_containing_fragment([self.att_only_fragments[0]])

        pipolb_fragment1 = self._orient_according_pipolb(self.pipolb_only_fragments[0])
        pipolb_fragment2 = _rev_comp_fragment(pipolb_fragment1)

        if ttrna_fragment:
            variant1 = self._create_pipolin(middle=pipolb_fragment1, right=ttrna_fragment)
            variant2 = self._create_pipolin(middle=pipolb_fragment2, right=ttrna_fragment)
            return [variant1, variant2]
        else:
            return [self._create_pipolin(middle=pipolb_fragment1),
                    self._create_pipolin(middle=pipolb_fragment2)]

    def _pipolb_plus_two_atts(self) -> Sequence[Pipolin]:
        # ---att---...---pol---...---att(t)---        tRNA is required
        ttrna_fragment = self._get_single_ttrna_containing_fragment(
            [self.att_only_fragments[0], self.att_only_fragments[1]]
        )
        pipolb_fragment1 = self._orient_according_pipolb(self.pipolb_only_fragments[0])
        pipolb_fragment2 = _rev_comp_fragment(pipolb_fragment1)

        if ttrna_fragment:
            if ttrna_fragment == self.att_only_fragments[0]:
                right_fragment, left_fragment = self._orient_fragments_according_ttrnas(
                    ttrna_fragment, self.att_only_fragments[1]
                )
            else:
                right_fragment, left_fragment = self._orient_fragments_according_ttrnas(
                    ttrna_fragment, self.att_only_fragments[0]
                )
            variant1 = self._create_pipolin(left=left_fragment, middle=pipolb_fragment1, right=right_fragment)
            variant2 = self._create_pipolin(left=left_fragment, middle=pipolb_fragment2, right=right_fragment)
            return [variant1, variant2]
        else:
            # TODO: not just pipolb!
            return [self._create_pipolin(middle=pipolb_fragment1),
                    self._create_pipolin(middle=pipolb_fragment2)]

    # -----------------------------------
    # auxiliary functions
    # -----------------------------------

    def _get_single_ttrna_containing_fragment(
            self, fragments: Sequence[PipolinFragment]
    ) -> Optional[PipolinFragment]:
        ttrna_containing_fragment = None
        for fragment in fragments:
            try:
                self._get_ttrna_of_fragment(fragment)
                if ttrna_containing_fragment:   # more than one contain ttrna
                    return None
                else:
                    ttrna_containing_fragment = fragment
            except CannotReconstructError:
                continue
        return ttrna_containing_fragment

    def _inflate_fragment(self, fragment: PipolinFragment, left: int, right: int) -> PipolinFragment:
        contig_length = self.genome.get_contig_by_id(fragment.contig_id).length

        if fragment.orientation == Strand.FORWARD:
            loc = Range(max(0, fragment.start - left), min(fragment.end + right, contig_length))
        else:
            loc = Range(max(0, fragment.start - right), min(fragment.end + left, contig_length))

        new = PipolinFragment(loc, fragment.contig_id, fragment.genome, orientation=fragment.orientation)
        new_features = new.get_fragment_features_sorted()
        return PipolinFragment(new.location, new.contig_id, new.genome, tuple(new_features), new.orientation)

    def _inflate_single_fragment_pipolin(self, fragment: PipolinFragment) -> Pipolin:
        atts_and_pipolbs = self._get_fragment_atts_and_pipolbs(fragment)

        left = _BORDER_INFLATE if atts_and_pipolbs[0].ftype == FeatureType.ATT else _NO_BORDER_INFLATE
        right = _BORDER_INFLATE if atts_and_pipolbs[-1].ftype == FeatureType.ATT else _NO_BORDER_INFLATE
        fragment = self._inflate_fragment(fragment, left, right)
        return Pipolin.from_fragments(fragment)

    def _create_pipolin(self, left=None, middle=None, right=None, complete=None) -> Pipolin:
        fragments = []
        if left is not None:
            fragments.append(self._inflate_fragment(left, _BORDER_INFLATE, _NO_BORDER_INFLATE))
        if middle is not None:
            fragments.append(self._inflate_fragment(middle, _NO_BORDER_INFLATE, _NO_BORDER_INFLATE))
        if right is not None:
            fragments.append(self._inflate_fragment(right, _NO_BORDER_INFLATE, _BORDER_INFLATE))
        if complete:
            fragments.append(self._inflate_fragment(complete, _BORDER_INFLATE, _BORDER_INFLATE))
        return Pipolin.from_fragments(*fragments)

    @staticmethod
    def _orient_according_pipolb(pipolb_containing_fragment: PipolinFragment) -> PipolinFragment:
        pipolbs = pipolb_containing_fragment.get_fragment_features_of_type_sorted(FeatureType.PIPOLB)
        pipolbs_strands = set(i.strand for i in pipolbs)
        if len(pipolbs_strands) == 1:
            if pipolbs_strands.pop() == Strand.REVERSE:
                return _rev_comp_fragment(pipolb_containing_fragment)
        return pipolb_containing_fragment

    def _orient_fragment_according_main(
            self, main_fragment: PipolinFragment, dependent_fragment: PipolinFragment) -> PipolinFragment:
        att_ids = set(f.att_id for f in main_fragment.features if isinstance(f, AttFeature))
        if len(att_ids) != 1:
            raise CannotReconstructError
        att_id = att_ids.pop()

        main_fragment_strand = self._get_fragment_atts_strand(main_fragment, att_id)
        dependent_fragment_strand = self._get_fragment_atts_strand(dependent_fragment, att_id)

        if main_fragment_strand != dependent_fragment_strand:
            return _rev_comp_fragment(dependent_fragment)
        else:
            return dependent_fragment

    def _orient_fragments_according_ttrnas(
            self, ttrna_containing_fragment: PipolinFragment, *other_fragment: PipolinFragment
    ) -> List[PipolinFragment]:
        ttrnas = [f for f in ttrna_containing_fragment.features if f.ftype == FeatureType.TARGET_TRNA]
        if len(ttrnas) == 0:
            raise CannotReconstructError

        att = self._get_att_overlapping_ttrna(ttrnas[0])

        if ttrnas[0].strand == Strand.FORWARD and ttrnas[0].start < att.start:  # 3'-overlap
            ttrna_containing_fragment = _rev_comp_fragment(ttrna_containing_fragment)

        elif ttrnas[0].strand == Strand.REVERSE and ttrnas[0].end > att.end:  # 3'-overlap
            pass
        else:
            raise CannotReconstructError

        new_other_fragments = []
        for ttrna in ttrnas:
            att = self._get_att_overlapping_ttrna(ttrna)
            other_fragment_strands = [self._get_fragment_atts_strand(f, att.att_id) for f in other_fragment]

            for f, s in zip(other_fragment, other_fragment_strands):
                if s:
                    if ttrna_containing_fragment.orientation == f.orientation and att.strand != s:
                        new_other_fragments.append(_rev_comp_fragment(f))
                    elif ttrna_containing_fragment.orientation != f.orientation and att.strand == s:
                        new_other_fragments.append(_rev_comp_fragment(f))
                    else:
                        new_other_fragments.append(f)
                else:
                    new_other_fragments.append(f)
        return [ttrna_containing_fragment] + new_other_fragments

    @staticmethod
    def _get_ttrna_of_fragment(fragment: PipolinFragment) -> Feature:
        ttrnas = [f for f in fragment.features if f.ftype == FeatureType.TARGET_TRNA]
        if len(ttrnas) == 0:
            raise CannotReconstructError
        # one tRNA is enough as all are of the same direction
        return ttrnas[0]

    def _get_att_overlapping_ttrna(self, ttrna) -> Feature:
        att = self.genome.features.get_features(FeatureType.ATT).get_overlapping(ttrna)
        if not att:
            raise AssertionError
        return att

    @staticmethod
    def _get_fragment_atts_strand(fragment: PipolinFragment, att_id) -> Optional[Strand]:
        atts = set(f.strand for f in fragment.features if isinstance(f, AttFeature) and f.att_id == att_id)
        if len(atts) == 1:
            return atts.pop()

    def _check_ttrnas(self, fragments: Sequence[PipolinFragment]) -> None:
        num_fragments_with_ttrnas = 0
        the_fragment = None
        for fragment in fragments:
            ttrnas = [f for f in fragment.features if f.ftype == FeatureType.TARGET_TRNA]
            if len(ttrnas) != 0:
                num_fragments_with_ttrnas += 1
                the_fragment = fragment

        if num_fragments_with_ttrnas > 1:
            raise CannotReconstructError(f'Number of fragments with tRNAs: {num_fragments_with_ttrnas} > 1')

        if the_fragment is not None:
            self._check_ttrnas_directions(the_fragment)

    @staticmethod
    def _check_ttrnas_directions(fragment: PipolinFragment) -> None:
        ttrnas = [f for f in fragment.features if f.ftype == FeatureType.TARGET_TRNA]
        strands = set(ttrna.strand for ttrna in ttrnas)
        if len(strands) != 1:
            raise CannotReconstructError(f'tRNAs on different strands on the fragment {fragment}')

    def _classify_fragments(self, fragments: Sequence[PipolinFragment]):
        for fragment in fragments:
            fragment_atts_and_pipolbs = self._get_fragment_atts_and_pipolbs(fragment)

            if all([i.ftype == FeatureType.PIPOLB for i in fragment_atts_and_pipolbs]):
                self.pipolb_only_fragments.append(fragment)
            elif all([i.ftype == FeatureType.ATT for i in fragment_atts_and_pipolbs]):
                self.att_only_fragments.append(fragment)

            elif fragment_atts_and_pipolbs[0].ftype == FeatureType.ATT and \
                    fragment_atts_and_pipolbs[-1].ftype == FeatureType.ATT:
                self.att_pipolb_att_fragments.append(fragment)
            else:
                self.att_pipolb_fragments.append(fragment)

    @staticmethod
    def _get_fragment_atts_and_pipolbs(fragment: PipolinFragment):
        features = fragment.get_fragment_features_sorted()
        return [f for f in features if f.ftype != FeatureType.TARGET_TRNA]


def _filter_redundant(
        pipolins_to_filter: MutableSequence[Sequence[Pipolin]]
) -> Sequence[Sequence[Pipolin]]:

    filtered = []

    for i, pipolin_variants in enumerate(pipolins_to_filter[:-1]):
        variants_to_leave = []
        for variant in pipolin_variants:
            if not _is_subset_of_any(variant, pipolins_to_filter[i + 1:]):
                variants_to_leave.append(variant)
        filtered.append(variants_to_leave)
    filtered.append(pipolins_to_filter[-1])

    return filtered


def _is_subset_of_any(pipolin: Pipolin, pipolins_to_filter: Sequence[Sequence[Pipolin]]) -> bool:
    pipolin_features = list(chain.from_iterable(f.features for f in pipolin.fragments))

    for pipolin_variants in pipolins_to_filter:
        for variant in pipolin_variants:
            if _are_comparable(pipolin, variant):
                variant_features = list(chain.from_iterable(f.features for f in variant.fragments))
                if set(pipolin_features).issubset(set(variant_features)):
                    return True
    return False


def _are_comparable(pipolin1: Pipolin, pipolin2: Pipolin) -> bool:
    fragments_by_contig_id = groupby(
        sorted(pipolin1.fragments + pipolin2.fragments, key=lambda x: x.contig_id),
        key=lambda x: x.contig_id
    )
    for _, group in fragments_by_contig_id:
        group = list(group)
        if len(group) == 2:
            if group[0].orientation != group[1].orientation:
                return False
    return True
