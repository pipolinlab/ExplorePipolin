from itertools import combinations, chain
from typing import Sequence, MutableSequence, List, Tuple, Iterable

from prefect import task
from prefect import context

from explore_pipolin.common import Genome, Pipolin, FeatureType, PipolinFragment, Strand, Feature, Range
from explore_pipolin.utilities.logging import genome_specific_logging


@task()
@genome_specific_logging
def scaffold_pipolins(genome: Genome, pipolins: Sequence[Pipolin]):
    logger = context.get('logger')

    single_fragment_pipolins: MutableSequence[Pipolin] = []
    other_pipolins: MutableSequence[Pipolin] = []

    result: MutableSequence[Pipolin] = []

    for pipolin in pipolins:
        if len(pipolin.fragments) == 1:
            logger.warning('>>> Scaffolding is not required!')
            single_fragment_pipolins.append(pipolin)
        else:
            logger.warning('>>> Trying to scaffold...')
            scaffolder = Scaffolder(genome=genome, pipolin=pipolin)
            try:
                result.append(scaffolder.scaffold_pipolin())
                logger.warning('>>> Scaffolding is done!')
            except CannotScaffoldError as e:
                logger.warning(f'>>> Cannot scaffold! {e}')
                other_pipolins.append(pipolin)

    for pipolin in single_fragment_pipolins:
        scaffolder = Scaffolder(genome=genome, pipolin=pipolin)
        result.append(scaffolder.inflate_single_fragment_pipolin())

    for pipolin in other_pipolins:
        scaffolder = Scaffolder(genome=genome, pipolin=pipolin)
        result.append(scaffolder.inflate_unscaffolded_pipolin())
    return _filter_overlapping(result) if len(result) > 1 else result


class CannotScaffoldError(Exception):
    pass


_BORDER_INFLATE = 0
_NO_BORDER_INFLATE = 100000


class Scaffolder:
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
            raise CannotScaffoldError

    def scaffold_pipolin(self) -> Pipolin:
        if self.att_pipolb_att_fragments:
            return self._att_pipolb_att_plus_atts()
        elif self.att_pipolb_fragments:
            return self._att_pipolb_plus_atts()
        elif self.pipolb_only_fragments:
            return self._pipolb_plus_atts()
        else:
            raise AssertionError

    def _att_pipolb_att_plus_atts(self) -> Pipolin:
        if len(self.att_only_fragments) == 1:
            # we can scaffold the case: ---att---pol---att---...---att(t)---   tRNA is required
            right_fragment = self.att_only_fragments[0]
            left_fragment = self.att_pipolb_att_fragments[0]

            try:
                ttrna = self._get_ttrna_of_fragment(right_fragment)
            except CannotScaffoldError:   # let's leave it as ---att---pol---att---
                return Pipolin.from_fragments(self._inflate_fragment(left_fragment, _BORDER_INFLATE, _BORDER_INFLATE))
            right_fragment, (left_fragment,) = self._orient_fragments_according_ttrna(ttrna,
                                                                                      right_fragment, left_fragment)
            return self._create_pipolin(left=left_fragment, right=right_fragment)
        else:
            raise CannotScaffoldError('---att---pol---att---...---att---...---att---...')

    def _att_pipolb_plus_atts(self) -> Pipolin:
        # we can scaffold the cases:
        # 1) ---att---pol---...---att---                  tRNA is not required
        # 2) ---att---...---pol---att---...---att(t)---   tRNA is required
        if len(self.att_only_fragments) == 1:
            return self._att_pipolb_plus_one_att()
        elif len(self.att_only_fragments) == 2:
            return self._att_pipolb_plus_two_atts()
        else:
            raise CannotScaffoldError

    def _att_pipolb_plus_one_att(self) -> Pipolin:
        # 1) ---att---pol---...---att---                  tRNA is not required
        if self.att_pipolb_fragments[0].features[0][1] == FeatureType.PIPOLB:
            left_fragment = self.att_only_fragments[0]
            right_fragment = self.att_pipolb_fragments[0]
            left_fragment = self._orient_fragment_according_main(right_fragment, left_fragment)
        elif self.att_pipolb_fragments[0].features[-1][1] == FeatureType.PIPOLB:
            left_fragment = self.att_pipolb_fragments[0]
            right_fragment = self.att_only_fragments[0]
            right_fragment = self._orient_fragment_according_main(left_fragment, right_fragment)
        else:
            raise AssertionError

        return self._create_pipolin(left=left_fragment, right=right_fragment)

    def _att_pipolb_plus_two_atts(self) -> Pipolin:
        # 2) ---att---...---pol---att---...---att(t)---   tRNA is required
        middle_fragment = self.att_pipolb_fragments[0]
        left_fragment, right_fragment, ttrna = self._define_left_right_fragments_and_ttrna()
        right_fragment, (left_fragment, middle_fragment) = \
            self._orient_fragments_according_ttrna(ttrna, right_fragment, left_fragment, middle_fragment)

        if middle_fragment.orientation == Strand.FORWARD:
            _is_pipolb_left = middle_fragment.features[0][1] == FeatureType.PIPOLB
        else:
            _is_pipolb_left = middle_fragment.features[-1][1] == FeatureType.PIPOLB

        if _is_pipolb_left:
            return self._create_pipolin(left=left_fragment, middle=middle_fragment, right=right_fragment)
        else:
            raise CannotScaffoldError('---att---...---att---pol---...---att(t)---')

    def _pipolb_plus_atts(self) -> Pipolin:
        # we can scaffold the cases, assuming that pipolb is on the plus strand:
        # 1) ---pol---...---att(t)---                    tRNA is required
        # 2) ---att---...---pol---...---att(t)---        tRNA is required
        if len(self.att_only_fragments) == 1:
            return self._pipolb_plus_one_att()

        elif len(self.att_only_fragments) == 2:
            return self._pipolb_plus_two_atts()

        else:
            raise CannotScaffoldError

    def _pipolb_plus_one_att(self):
        # 1) ---pol---...---att(t)---                    tRNA is required
        middle_fragment = self._orient_pipolb_only_fragment()
        right_fragment = self.att_only_fragments[0]
        ttrna = self._get_ttrna_of_fragment(right_fragment)
        right_fragment, _ = self._orient_fragments_according_ttrna(ttrna, right_fragment)

        return self._create_pipolin(middle=middle_fragment, right=right_fragment)

    def _pipolb_plus_two_atts(self):
        # 2) ---att---...---pol---...---att(t)---        tRNA is required
        middle_fragment = self._orient_pipolb_only_fragment()
        left_fragment, right_fragment, ttrna = self._define_left_right_fragments_and_ttrna()
        right_fragment, (left_fragment,) = self._orient_fragments_according_ttrna(ttrna, right_fragment, left_fragment)

        return self._create_pipolin(left=left_fragment, middle=middle_fragment, right=right_fragment)

    # -----------------------------------
    # auxiliary functions
    # -----------------------------------

    def inflate_single_fragment_pipolin(self):
        atts_and_pipolbs = self._get_fragment_atts_and_pipolbs(self.pipolin.fragments[0])

        left = _BORDER_INFLATE if atts_and_pipolbs[0][1] == FeatureType.ATT else _NO_BORDER_INFLATE
        right = _BORDER_INFLATE if atts_and_pipolbs[-1][1] == FeatureType.ATT else _NO_BORDER_INFLATE
        fragment = self._inflate_fragment(self.pipolin.fragments[0], left, right)
        return Pipolin.from_fragments(fragment)

    def inflate_unscaffolded_pipolin(self):
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
        return PipolinFragment(new.location, new.contig_id, new.genome, new_features, new.orientation)

    def _create_pipolin(self, left=None, middle=None, right=None) -> Pipolin:
        fragments = []
        if left is not None:
            fragments.append(self._inflate_fragment(left, _BORDER_INFLATE, _NO_BORDER_INFLATE))
        if middle is not None:
            fragments.append(self._inflate_fragment(middle, _NO_BORDER_INFLATE, _NO_BORDER_INFLATE))
        if right is not None:
            fragments.append(self._inflate_fragment(right, _NO_BORDER_INFLATE, _BORDER_INFLATE))
        return Pipolin.from_fragments(*fragments)

    def _orient_pipolb_only_fragment(self) -> PipolinFragment:
        pipolb_only_fragment = self.pipolb_only_fragments[0]
        if pipolb_only_fragment.features[0][1] == FeatureType.PIPOLB:
            if pipolb_only_fragment.features[0][0].strand == Strand.REVERSE:
                pipolb_only_fragment = self._make_fragment_reverse(pipolb_only_fragment)
        else:
            raise AssertionError(f'The first feature of the pipolb_only_fragment is '
                                 f'{pipolb_only_fragment.features[0][1]} != FeatureType.PIPOLB')
        return pipolb_only_fragment

    def _define_left_right_fragments_and_ttrna(
            self) -> Tuple[PipolinFragment, PipolinFragment, Feature]:
        try:
            ttrna = self._get_ttrna_of_fragment(self.att_only_fragments[0])
            right_fragment = self.att_only_fragments[0]
            left_fragment = self.att_only_fragments[1]
        except CannotScaffoldError:
            ttrna = self._get_ttrna_of_fragment(self.att_only_fragments[1])
            right_fragment = self.att_only_fragments[1]
            left_fragment = self.att_only_fragments[0]
        return left_fragment, right_fragment, ttrna

    def _orient_fragment_according_main(
            self, main_fragment: PipolinFragment, dependent_fragment: PipolinFragment) -> PipolinFragment:
        att_ids = set(f[0].att_id for f in main_fragment.features if f[1] == FeatureType.ATT)
        if len(att_ids) != 1:
            raise CannotScaffoldError
        att_id = att_ids.pop()

        main_fragment_strand = self._get_fragment_atts_strand(main_fragment, att_id)
        dependent_fragment_strand = self._get_fragment_atts_strand(dependent_fragment, att_id)

        if main_fragment_strand != dependent_fragment_strand:
            dependent_fragment = self._make_fragment_reverse(dependent_fragment)
        return dependent_fragment

    def _orient_fragments_according_ttrna(self, ttrna: Feature, ttrna_containing_fragment: PipolinFragment,
                                          *other_fragment: PipolinFragment,
                                          ) -> Tuple[PipolinFragment, Iterable[PipolinFragment]]:
        new_other_fragments = []

        att = self._get_att_overlapping_ttrna(ttrna)
        other_fragment_strands = [self._get_fragment_atts_strand(f, att.att_id) for f in other_fragment]

        if ttrna.strand == Strand.FORWARD and ttrna.start < att.start:  # 3'-overlap
            ttrna_containing_fragment = self._make_fragment_reverse(ttrna_containing_fragment)
            for f, s in zip(other_fragment, other_fragment_strands):
                if att.strand == s:
                    new_other_fragments.append(self._make_fragment_reverse(f))
                else:
                    new_other_fragments.append(f)

        elif ttrna.strand == Strand.REVERSE and ttrna.end > att.end:  # 3'-overlap
            for f, s in zip(other_fragment, other_fragment_strands):
                if att.strand != s:
                    new_other_fragments.append(self._make_fragment_reverse(f))
                else:
                    new_other_fragments.append(f)

        else:
            raise CannotScaffoldError

        return ttrna_containing_fragment, new_other_fragments

    @staticmethod
    def _get_ttrna_of_fragment(ttrna_containing_fragment: PipolinFragment) -> Feature:
        ttrnas = [f[0] for f in ttrna_containing_fragment.features if f[1] == FeatureType.TARGET_TRNA]
        if len(ttrnas) != 1:
            raise CannotScaffoldError
        return ttrnas[0]

    def _get_att_overlapping_ttrna(self, ttrna):
        att = self.genome.features.get_features(FeatureType.ATT).get_overlapping(ttrna)
        if not att:
            raise AssertionError
        return att

    @staticmethod
    def _get_fragment_atts_strand(fragment: PipolinFragment, att_id):
        atts = set(f[0].strand for f in fragment.features if f[1] == FeatureType.ATT and f[0].att_id == att_id)
        if len(atts) != 1:
            raise AssertionError
        return atts.pop()

    @staticmethod
    def _make_fragment_reverse(fragment: PipolinFragment):
        return PipolinFragment(fragment.location, fragment.contig_id, fragment.genome, fragment.features,
                               orientation=Strand.REVERSE)

    def _classify_fragments(self, fragments: Sequence[PipolinFragment]):
        for fragment in fragments:
            fragment_atts_and_pipolbs = self._get_fragment_atts_and_pipolbs(fragment)

            if all([i[1] == FeatureType.PIPOLB for i in fragment_atts_and_pipolbs]):
                self.pipolb_only_fragments.append(fragment)
            elif all([i[1] == FeatureType.ATT for i in fragment_atts_and_pipolbs]):
                self.att_only_fragments.append(fragment)

            elif fragment_atts_and_pipolbs[0][1] == FeatureType.ATT and \
                    fragment_atts_and_pipolbs[-1][1] == FeatureType.ATT:
                self.att_pipolb_att_fragments.append(fragment)
            else:
                self.att_pipolb_fragments.append(fragment)

    @staticmethod
    def _get_fragment_atts_and_pipolbs(fragment: PipolinFragment):
        return [f for f in fragment.features if f[1] != FeatureType.TARGET_TRNA]


def _filter_overlapping(pipolins_to_filter: Sequence[Pipolin]) -> Sequence[Pipolin]:
    filtered = set()
    for comb in combinations(pipolins_to_filter, 2):
        p1_features = [chain.from_iterable(f.features) for f in comb[0].fragments]
        p2_features = [chain.from_iterable(f.features) for f in comb[1].fragments]
        if set(p1_features) in set(p2_features):
            filtered.add(comb[1])
        elif set(p2_features) in set(p1_features):
            filtered.add(comb[0])
        else:
            filtered.update(comb)
    return list(filtered)