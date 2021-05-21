from itertools import combinations, chain
from typing import Sequence, MutableSequence, List, Tuple, Optional

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


@task()
@genome_specific_logging
def reconstruct_pipolins(genome: Genome, pipolins: Sequence[Pipolin]):
    logger = context.get('logger')

    single_fragment_pipolins: MutableSequence[Pipolin] = []
    other_pipolins: MutableSequence[Pipolin] = []

    result: MutableSequence[Pipolin] = []

    for i_p, pipolin in enumerate(pipolins):
        logger.info(f'Reconstructing pipolin {i_p + 1}:')
        if len(pipolin.fragments) == 1:
            logger.info('>>> Reconstruction is not required for pipolin:')

            features = pipolin.fragments[0].get_fragment_features_sorted()
            if features[0].ftype == FeatureType.TARGET_TRNA:
                pipolin.fragments[0].change_orientation()

            logger.warning(f'{draw_pipolin_structure(pipolin)[0]}')
            single_fragment_pipolins.append(pipolin)
        else:
            logger.info('>>> Trying to reconstruct the structure from fragments:')
            logger.info(''.join(f'\n\t{i}' for i in draw_pipolin_structure(pipolin)))
            try:
                reconstructor = Reconstructor(genome=genome, pipolin=pipolin)
                result.append(reconstructor.reconstruct_pipolin())
                logger.info('>>> Reconstruction is done! The resulting pipolin is:')
                logger.info('...'.join([i.split(sep=': ')[1] for i in draw_pipolin_structure(result[-1])]))
                logger.warning('PLEASE, DOUBLE CHECK THE FINAL STRUCTURE MANUALLY!')
            except CannotReconstructError as e:
                logger.warning(f'>>> Cannot reconstruct the structure! {e}')
                other_pipolins.append(pipolin)

    for pipolin in single_fragment_pipolins:
        reconstructor = Reconstructor(genome=genome, pipolin=pipolin)
        result.append(reconstructor.inflate_single_fragment_pipolin())

    for pipolin in other_pipolins:
        reconstructor = Reconstructor(genome=genome, pipolin=pipolin)
        result.append(reconstructor.inflate_unreconstructed_pipolin())

    filtered = _filter_redundant(result) if len(result) > 1 else result
    logger.info('Pipolins after filtering:')
    for f in filtered:
        logger.info('...'.join([i.split(sep=': ')[1] for i in draw_pipolin_structure(f)]))
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
            raise CannotReconstructError

        # If >=2 contigs with tRNAs or a contig with >=2 tRNAs on different strands -> raise CannotReconstructError
        self._check_ttrnas(pipolin.fragments)

    def reconstruct_pipolin(self) -> Pipolin:
        if self.att_pipolb_att_fragments:
            return self._att_pipolb_att_plus_atts()
        elif self.att_pipolb_fragments:
            return self._att_pipolb_plus_atts()
        elif self.pipolb_only_fragments:
            return self._pipolb_plus_atts()
        else:
            raise AssertionError

    def _att_pipolb_att_plus_atts(self) -> Pipolin:
        if len(self.att_only_fragments) != 1:
            raise CannotReconstructError

        # we can reconstruct the case: ---att---pol---att---...---att(t)---   tRNA is required
        try:
            self._orient_fragments_according_ttrnas(
                self.att_only_fragments[0], self.att_pipolb_att_fragments[0]
            )
            return self._create_pipolin(
                left=self.att_pipolb_att_fragments[0], right=self.att_only_fragments[0]
            )

        except CannotReconstructError:
            self._orient_fragments_according_ttrnas(
                self.att_pipolb_att_fragments[0], self.att_only_fragments[0]
            )
            return self._create_pipolin(
                left=self.att_only_fragments[0], right=self.att_pipolb_att_fragments[0]
            )

    def _att_pipolb_plus_atts(self) -> Pipolin:
        # we can reconstruct the cases:
        # 1) ---att---pol---...---att---                  tRNA is not required
        # 2) ---att---...---pol---att---...---att(t)---   tRNA is required
        if len(self.att_only_fragments) == 1:
            return self._att_pipolb_plus_one_att()
        elif len(self.att_only_fragments) == 2:
            return self._att_pipolb_plus_two_atts()
        else:
            raise CannotReconstructError

    def _att_pipolb_plus_one_att(self) -> Pipolin:
        # 1) ---att---pol---...---att---                  tRNA is not required (but helps with orientation)
        if self.att_pipolb_fragments[0].features[0].ftype == FeatureType.PIPOLB:
            left_fragment = self.att_only_fragments[0]
            right_fragment = self.att_pipolb_fragments[0]
            left_fragment = self._orient_fragment_according_main(right_fragment, left_fragment)
        elif self.att_pipolb_fragments[0].features[-1].ftype == FeatureType.PIPOLB:
            left_fragment = self.att_pipolb_fragments[0]
            right_fragment = self.att_only_fragments[0]
            right_fragment = self._orient_fragment_according_main(left_fragment, right_fragment)
        else:
            raise AssertionError

        try:
            left_fragment, right_fragment, ttrna = self._define_left_right_fragments_and_ttrna(
                left_fragment, right_fragment
            )
            self._orient_fragments_according_ttrnas(right_fragment, left_fragment)

        except CannotReconstructError:
            pass

        return self._create_pipolin(left=left_fragment, right=right_fragment)

    def _att_pipolb_plus_two_atts(self) -> Pipolin:
        # 2) ---att---...---pol---att---...---att(t)---   tRNA is required
        middle_fragment = self.att_pipolb_fragments[0]
        left_fragment, right_fragment, ttrna = self._define_left_right_fragments_and_ttrna(
            self.att_only_fragments[0], self.att_only_fragments[1]
        )
        self._orient_fragments_according_ttrnas(right_fragment, left_fragment, middle_fragment)

        if middle_fragment.orientation == Strand.FORWARD:
            pipolb_is_left = middle_fragment.features[0].ftype == FeatureType.PIPOLB
        else:
            pipolb_is_left = middle_fragment.features[-1].ftype == FeatureType.PIPOLB

        if pipolb_is_left:
            return self._create_pipolin(left=left_fragment, middle=middle_fragment, right=right_fragment)
        else:
            raise CannotReconstructError('---att---...---att---pol---...---att(t)---')

    def _pipolb_plus_atts(self) -> Pipolin:
        # we can try to reconstruct the cases, assuming that pipolb is on the plus strand:
        # 1) ---pol---...---att(t)---                    tRNA is required
        # 2) ---att---...---pol---...---att(t)---        tRNA is required
        if len(self.att_only_fragments) == 1:
            return self._pipolb_plus_one_att()

        elif len(self.att_only_fragments) == 2:
            return self._pipolb_plus_two_atts()

        else:
            raise CannotReconstructError

    def _pipolb_plus_one_att(self) -> Pipolin:
        # 1) ---pol---...---att(t)---                    tRNA is required
        self._orient_pipolb_only_fragment()
        self._orient_fragments_according_ttrnas(self.att_only_fragments[0])
        return self._create_pipolin(middle=self.pipolb_only_fragments[0], right=self.att_only_fragments[0])

    def _pipolb_plus_two_atts(self) -> Pipolin:
        # 2) ---att---...---pol---...---att(t)---        tRNA is required
        self._orient_pipolb_only_fragment()
        left_fragment, right_fragment, ttrna = self._define_left_right_fragments_and_ttrna(
            self.att_only_fragments[0], self.att_only_fragments[1]
        )
        self._orient_fragments_according_ttrnas(right_fragment, left_fragment)
        return self._create_pipolin(
            left=left_fragment, middle=self.pipolb_only_fragments[0], right=right_fragment
        )

    # -----------------------------------
    # auxiliary functions
    # -----------------------------------

    def inflate_single_fragment_pipolin(self):
        atts_and_pipolbs = self._get_fragment_atts_and_pipolbs(self.pipolin.fragments[0])

        left = _BORDER_INFLATE if atts_and_pipolbs[0].ftype == FeatureType.ATT else _NO_BORDER_INFLATE
        right = _BORDER_INFLATE if atts_and_pipolbs[-1].ftype == FeatureType.ATT else _NO_BORDER_INFLATE
        fragment = self._inflate_fragment(self.pipolin.fragments[0], left, right)
        return Pipolin.from_fragments(fragment)

    def inflate_unreconstructed_pipolin(self):
        inflated_fragments = []
        for fragment in self.pipolin.fragments:
            inflated_fragments.append(self._inflate_fragment(fragment, _NO_BORDER_INFLATE, _NO_BORDER_INFLATE))
        return Pipolin.from_fragments(*inflated_fragments)

    def _inflate_fragment(self, fragment: PipolinFragment, left: int, right: int) -> PipolinFragment:
        contig_length = self.genome.get_contig_by_id(fragment.contig_id).length

        if fragment.orientation == Strand.FORWARD:
            loc = Range(max(0, fragment.start - left), min(fragment.end + right, contig_length))
        else:
            loc = Range(max(0, fragment.start - right), min(fragment.end + left, contig_length))

        new = PipolinFragment(loc, fragment.contig_id, fragment.genome, orientation=fragment.orientation)
        new_features = new.get_fragment_features_sorted()
        return PipolinFragment(new.location, new.contig_id, new.genome, tuple(new_features), new.orientation)

    def _create_pipolin(self, left=None, middle=None, right=None) -> Pipolin:
        fragments = []
        if left is not None:
            fragments.append(self._inflate_fragment(left, _BORDER_INFLATE, _NO_BORDER_INFLATE))
        if middle is not None:
            fragments.append(self._inflate_fragment(middle, _NO_BORDER_INFLATE, _NO_BORDER_INFLATE))
        if right is not None:
            fragments.append(self._inflate_fragment(right, _NO_BORDER_INFLATE, _BORDER_INFLATE))
        return Pipolin.from_fragments(*fragments)

    def _orient_pipolb_only_fragment(self) -> None:
        pipolb_only_fragment = self.pipolb_only_fragments[0]
        if pipolb_only_fragment.features[0].ftype == FeatureType.PIPOLB:
            if pipolb_only_fragment.features[0].strand == Strand.REVERSE:
                pipolb_only_fragment.change_orientation()
        else:
            raise AssertionError

    def _define_left_right_fragments_and_ttrna(
            self, fragment1, fragment2) -> Tuple[PipolinFragment, PipolinFragment, Feature]:
        try:
            ttrna = self._get_ttrna_of_fragment(fragment1)
            right_fragment, left_fragment = fragment1, fragment2
        except CannotReconstructError:
            ttrna = self._get_ttrna_of_fragment(fragment2)
            right_fragment, left_fragment = fragment2, fragment1

        return left_fragment, right_fragment, ttrna

    def _orient_fragment_according_main(
            self, main_fragment: PipolinFragment, dependent_fragment: PipolinFragment) -> PipolinFragment:
        att_ids = set(f.att_id for f in main_fragment.features if isinstance(f, AttFeature))
        if len(att_ids) != 1:
            raise CannotReconstructError
        att_id = att_ids.pop()

        main_fragment_strand = self._get_fragment_atts_strand(main_fragment, att_id)
        dependent_fragment_strand = self._get_fragment_atts_strand(dependent_fragment, att_id)

        if main_fragment_strand != dependent_fragment_strand:
            dependent_fragment.change_orientation()
        return dependent_fragment

    def _orient_fragments_according_ttrnas(
            self, ttrna_containing_fragment: PipolinFragment, *other_fragment: PipolinFragment
    ) -> None:
        ttrnas = [f for f in ttrna_containing_fragment.features if f.ftype == FeatureType.TARGET_TRNA]
        if len(ttrnas) == 0:
            raise CannotReconstructError

        att = self._get_att_overlapping_ttrna(ttrnas[0])

        if ttrnas[0].strand == Strand.FORWARD and ttrnas[0].start < att.start:  # 3'-overlap
            ttrna_containing_fragment.change_orientation()

        elif ttrnas[0].strand == Strand.REVERSE and ttrnas[0].end > att.end:  # 3'-overlap
            pass
        else:
            raise CannotReconstructError

        for ttrna in ttrnas:
            att = self._get_att_overlapping_ttrna(ttrna)
            other_fragment_strands = [self._get_fragment_atts_strand(f, att.att_id) for f in other_fragment]

            for f, s in zip(other_fragment, other_fragment_strands):
                if s:
                    if ttrna_containing_fragment.orientation == f.orientation and att.strand != s:
                        f.change_orientation()
                    elif ttrna_containing_fragment.orientation != f.orientation and att.strand == s:
                        f.change_orientation()

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
            raise CannotReconstructError

        if the_fragment is not None:
            self._check_ttrnas_directions(the_fragment)

    @staticmethod
    def _check_ttrnas_directions(fragment: PipolinFragment) -> None:
        ttrnas = [f for f in fragment.features if f.ftype == FeatureType.TARGET_TRNA]
        strands = set(ttrna.strand for ttrna in ttrnas)
        if len(strands) != 1:
            raise CannotReconstructError

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


def _filter_redundant(pipolins_to_filter: Sequence[Pipolin]) -> Sequence[Pipolin]:
    filtered = set()
    for comb in combinations(pipolins_to_filter, 2):
        p1_features = list(chain.from_iterable(f.features for f in comb[0].fragments))
        p2_features = list(chain.from_iterable(f.features for f in comb[1].fragments))
        if set(p1_features).issubset(set(p2_features)):
            filtered.add(comb[1])
        elif set(p2_features).issubset(set(p1_features)):
            filtered.add(comb[0])
        else:
            filtered.update(comb)
    return list(filtered)
