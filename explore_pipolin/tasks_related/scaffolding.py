from enum import Enum, auto
from typing import Sequence, MutableSequence, Optional, Set

from prefect import task

from explore_pipolin.common import PipolinFragment, FeatureType, Contig, Feature, Pipolin, Genome, Range
from explore_pipolin.tasks_related.misc import get_ranges_around_pipolbs
from explore_pipolin.utilities.logging import genome_specific_logging


@task()
@genome_specific_logging
def scaffold_pipolins(genome: Genome, pipolins: Sequence[Pipolin]) -> Sequence[Pipolin]:
    # Useful link to check feature's qualifiers: https://www.ebi.ac.uk/ena/WebFeat/
    # https://github.com/biopython/biopython/issues/1755
    return scaffold(genome)


class Direction(Enum):
    LEFT = auto()
    RIGHT = auto()


def scaffold(genome: Genome) -> Sequence[Pipolin]:
    if genome.features.is_on_the_same_contig(FeatureType.PIPOLB, FeatureType.ATT, FeatureType.TARGET_TRNA):
        print('>>> Scaffolding is not required!')
        return [_create_pipolin_fragments_single_contig(genome)]
    else:
        print('>>> Scaffolding is required!')
        unchangeable_contigs = _get_unchangeable_contigs(genome)

        if len(unchangeable_contigs) == 1:
            unchangeable_contig = unchangeable_contigs[0]
            print(f'The unchangeable contig is {unchangeable_contig.id}!')

            att_only_contigs = _order_att_only_contigs(genome)
            if len(att_only_contigs) == 1:
                direction = _get_direction_of_unchangeable(genome, unchangeable_contig.id)

                orphan_atts = genome.features.atts_dict()[att_only_contigs[0].id]

                if direction is None:
                    contig_length = unchangeable_contig.length
                    pipolin_range = Range(
                        start=genome.features.atts_dict()[unchangeable_contig.id][0].start,
                        end=genome.features.atts_dict()[unchangeable_contig.id][-1].end
                    )
                    pipolin_range = pipolin_range.inflate(50, _max=contig_length)

                    pipolin = PipolinFragment(contig_id=unchangeable_contig.id, location=pipolin_range)
                    return [Pipolin.from_fragments(pipolin)]
                    # TODO: if att_only_contig has a target_trna, it could be added on the right
                if direction.RIGHT:
                    left_fragment = _create_unchangeable_fragment(genome, unchangeable_contig, direction)
                    right_fragment = _create_att_contig_fragment(orphan_atts, direction)
                    return [Pipolin.from_fragments(left_fragment, right_fragment)]
                else:
                    left_fragment = _create_att_contig_fragment(orphan_atts, direction)
                    right_fragment = _create_unchangeable_fragment(genome, unchangeable_contig, direction)
                    return [Pipolin.from_fragments(left_fragment, right_fragment)]
            elif len(att_only_contigs) == 2:
                # TODO: the order can be also [middle_fragment, left_fragment, right_fragment]
                middle_fragment = PipolinFragment(contig_id=unchangeable_contig.id,
                                                  location=Range(0, unchangeable_contig.length))

                left_atts = genome.features.atts_dict()[att_only_contigs[0].id]
                left_fragment = _create_att_contig_fragment(left_atts, Direction.LEFT)

                right_atts = genome.features.atts_dict()[att_only_contigs[1].id]
                right_fragment = _create_att_contig_fragment(right_atts, Direction.RIGHT)

                return [Pipolin.from_fragments(left_fragment, middle_fragment, right_fragment)]
            else:
                raise NotImplementedError
        elif len(unchangeable_contigs) == 2:
            target_trnas_contig0 = genome.features.target_trnas_dict()[unchangeable_contigs[0].id]
            target_trnas_contig1 = genome.features.target_trnas_dict()[unchangeable_contigs[1].id]

            # TODO: FIX THIS!!!!!! both contigs can have tRNAs => both != 0
            if len(target_trnas_contig0) != 0:
                right_contig = unchangeable_contigs[0]
                left_contig = unchangeable_contigs[1]
            elif len(target_trnas_contig1) != 0:
                right_contig = unchangeable_contigs[1]
                left_contig = unchangeable_contigs[0]
            else:
                raise NotImplementedError

            left_direction = _get_direction_of_unchangeable(genome, left_contig.id)
            left_fragment = _create_unchangeable_fragment(genome, left_contig, left_direction)

            right_direction = _get_direction_of_unchangeable(genome, right_contig.id)
            right_fragment = _create_unchangeable_fragment(genome, right_contig, right_direction)

            return [Pipolin.from_fragments(left_fragment, right_fragment)]
        elif len(unchangeable_contigs) > 2:
            raise NotImplementedError
        else:
            return [_try_finish_separate(genome)]


def _create_unchangeable_fragment(genome: Genome, contig: Contig, direction: Direction) -> PipolinFragment:
    if direction is Direction.RIGHT:
        fragment_range = Range(genome.features.atts_dict()[contig.id][0].start, contig.length)
        fragment_range = fragment_range.inflate(50, _max=contig.length)
        fragment = PipolinFragment(contig_id=contig.id, location=fragment_range)
    else:
        fragment_range = Range(0, genome.features.atts_dict()[contig.id][-1].end)
        fragment_range = fragment_range.inflate(50, _max=contig.length)
        fragment = PipolinFragment(contig_id=contig.id, location=fragment_range)

    return fragment


def _try_finish_separate(genome: Genome) -> Pipolin:
    polbs_fragments = _create_polbs_fragments(genome)

    if len(genome.features.get_features(FeatureType.TARGET_TRNA)) != 1:
        raise NotImplementedError('At the moment I can only handle one target tRNA')

    the_target_trna: Feature = next(iter(genome.features.get_features(FeatureType.TARGET_TRNA)))

    right_fragment = _create_att_contig_fragment(
        genome.features.atts_dict()[the_target_trna.contig_id], Direction.RIGHT)

    atts_contigs = set([i.contig for i in genome.features.get_features(FeatureType.ATT)])
    if len(atts_contigs) == 2:
        print('The single record can be created!!!\n')

        atts_contigs = list(atts_contigs)
        if atts_contigs[0].id == the_target_trna.contig_id:
            left_contig = atts_contigs[1]
        else:
            left_contig = atts_contigs[0]

        left_fragment = _create_att_contig_fragment(
            genome.features.atts_dict()[left_contig.id], Direction.LEFT)
    else:
        raise NotImplementedError

    return Pipolin.from_fragments(left_fragment, *polbs_fragments, right_fragment)


def _create_polbs_fragments(genome: Genome) -> Sequence[PipolinFragment]:
    polbs_contigs = set([i.contig for i in genome.features.get_features(FeatureType.PIPOLB)])
    if len(polbs_contigs) == 2:
        print('Two "polb only contigs" were found!')
        polbs_contigs = list(polbs_contigs)
        if polbs_contigs[0].length + polbs_contigs[1].length < 15000:
            polb0_fragment = PipolinFragment(contig_id=polbs_contigs[0].id,
                                             location=Range(0, polbs_contigs[0].length))
            polb1_fragment = PipolinFragment(contig_id=polbs_contigs[1].id,
                                             location=Range(0, polbs_contigs[1].length))

            contig_features = genome.features.get_features(FeatureType.PIPOLB).get_list_of_contig_sorted(
                polbs_contigs[0].id)
            polb0_length = sum((i.end - i.start) for i in contig_features)
            contig_features = genome.features.get_features(FeatureType.PIPOLB).get_list_of_contig_sorted(
                polbs_contigs[1].id)
            polb1_length = sum((i.end - i.start) for i in contig_features)
            # TODO: comparing just by length is an unreliable way! REWRITE if possible!
            if polb0_length < polb1_length:
                return [polb0_fragment, polb1_fragment]
            else:
                return [polb1_fragment, polb0_fragment]
        else:
            raise NotImplementedError

    elif len(polbs_contigs) == 1:
        the_polb_contig, = polbs_contigs
        return [PipolinFragment(contig_id=the_polb_contig.id, location=Range(0, the_polb_contig.length))]

    raise NotImplementedError


def _get_unchangeable_contigs(genome: Genome) -> MutableSequence[Contig]:
    polbs_contigs = set([i.contig for i in genome.features.get_features(FeatureType.PIPOLB)])

    contigs_to_return = []
    for contig in polbs_contigs:
        atts_next_polbs = genome.features.get_features(FeatureType.ATT).get_list_of_contig_sorted(contig.id)
        if len(atts_next_polbs) != 0:
            contigs_to_return.append(contig)

    return contigs_to_return


def _create_att_contig_fragment(contig_atts: Sequence[Feature], direction: Direction)\
        -> PipolinFragment:
    contig = contig_atts[0].contig
    contig_length = contig.length
    if direction is Direction.RIGHT:
        fragment_range = Range(0, contig_atts[-1].end + 50)
    else:
        fragment_range = Range(max(0, contig_atts[0].start - 50), contig_length)

    fragment = PipolinFragment(contig_id=contig.id, location=fragment_range.clamp(0, contig_length))
    return fragment


def _order_att_only_contigs(genome: Genome) -> MutableSequence[Contig]:
    att_only_contigs = list(_get_att_only_contigs(genome))
    if len(att_only_contigs) == 1:
        print('The single record can be created!!!\n')
        return att_only_contigs
    elif len(att_only_contigs) == 2:
        if _contig_has_feature_type(genome, att_only_contigs[0].id, FeatureType.TARGET_TRNA):
            # TODO: here sometimes fails for LREC241: KeyError: 'NODE_38'
            print('The single record can be created!!!\n')
            return [att_only_contigs[1], att_only_contigs[0]]
        elif _contig_has_feature_type(genome, att_only_contigs[1].id, FeatureType.TARGET_TRNA):
            print('The single record can be created!!!\n')
            return [att_only_contigs[0], att_only_contigs[1]]
        else:
            raise NotImplementedError
    else:
        raise NotImplementedError


def _get_att_only_contigs(genome: Genome) -> Set[Contig]:
    att_only_contigs = set()
    for att in genome.features.get_features(FeatureType.ATT):
        polbs_next_att = genome.features.get_features(FeatureType.PIPOLB).get_list_of_contig_sorted(
            att.contig_id)
        if len(polbs_next_att) == 0:
            att_only_contigs.add(att.contig)

    return att_only_contigs


def _get_direction_of_unchangeable(genome: Genome, contig_id: str) -> Optional[Direction]:
    polbs_sorted = genome.features.pipolbs_dict()[contig_id]
    atts_sorted = genome.features.atts_dict()[contig_id]
    # TODO: FIX THIS!!!!!!! for the case '---att---pol---att---pol---...---att(t)---'
    if polbs_sorted[0].start > atts_sorted[-1].end:
        return Direction.RIGHT
    elif polbs_sorted[-1].end < atts_sorted[0].start:
        return Direction.LEFT
    else:
        return None


def _contig_has_feature_type(genome: Genome, contig_id: str, feature_type: FeatureType) -> bool:
    feature_dict = genome.features.get_features(feature_type).get_dict_by_contig_sorted()
    return len(feature_dict[contig_id]) != 0


def _create_pipolin_fragments_single_contig(genome: Genome) -> Pipolin:
    if len(genome.features.atts_dict()) != 0:
        pipolin_range = _get_pipolin_bounds(genome)
        pipolin = PipolinFragment(contig_id=genome.features.get_features(FeatureType.PIPOLB).first.contig_id,
                                  location=pipolin_range)

        return Pipolin.from_fragments(pipolin)
    else:
        windows = get_ranges_around_pipolbs(genome)
        if len(windows) != 1:
            raise NotImplementedError('this method is for a single contig => should be only one Window!')
        pipolin = PipolinFragment(contig_id=genome.features.get_features(FeatureType.PIPOLB).first.contig_id,
                                  location=Range(windows[0].left_range.start, windows[0].right_range.end))
        return Pipolin.from_fragments(pipolin)


def _get_pipolin_bounds(genome: Genome) -> Range:
    pipolbs = sorted((i for i in genome.features.get_features(FeatureType.PIPOLB)), key=lambda p: p.start)
    atts = sorted((i for i in genome.features.get_features(FeatureType.ATT)), key=lambda p: p.start)

    contig_length = pipolbs[0].contig.length
    left_edge = max(atts[0].start - 50, 0) if atts[0].start < pipolbs[0].start else 0
    right_edge = min(atts[-1].end + 50, contig_length) if atts[-1].end > pipolbs[-1].end else contig_length
    return Range(left_edge, right_edge)
