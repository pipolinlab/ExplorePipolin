from __future__ import annotations

from explore_pipolin.common import Orientation, Contig, Genome, Feature, FeatureType, FeaturesContainer


def is_single_target_trna_per_contig(features_container: FeaturesContainer):
    # TODO: don't like this
    # there was one case with two target trnas per genome, although usually only one
    target_trnas_dict = features_container.get_features_dict_by_contig_normalized(feature_type=FeatureType.TARGET_TRNA)
    target_trnas = features_container.get_features(feature_type=FeatureType.TARGET_TRNA)
    if len(target_trnas) != len(target_trnas_dict):
        raise AssertionError("We are expecting a single tRNA to overlap with a single att per contig!")


def get_contig_orientation(contig: Contig, features_container: FeaturesContainer) -> Orientation:
    target_trnas = features_container.get_features_of_contig_normalized(contig_id=contig.contig_id,
                                                                        feature_type=FeatureType.TARGET_TRNA)
    atts = features_container.get_features_of_contig_normalized(contig_id=contig.contig_id,
                                                                feature_type=FeatureType.ATT)
    atts_strands = [att.strand for att in atts]
    polbs = features_container.get_features_of_contig_normalized(contig_id=contig.contig_id,
                                                                 feature_type=FeatureType.PIPOLB)
    polbs_strands = [polb.strand for polb in polbs]

    if len(target_trnas) != 0:
        if len(set(atts_strands)) != 1:
            raise AssertionError('ATTs are expected to be in the same frame, as they are direct repeats!')
        if set(atts_strands).pop() == target_trnas[0].strand:
            raise AssertionError('ATT and tRNA are expected to be on the different strands!')
        return - target_trnas[0].strand

    elif len(atts) != 0:
        if len(set(atts_strands)) != 1:
            raise AssertionError('ATTs are expected to be in the same frame, as they are direct repeats!')
        return atts[0].strand

    if len(polbs) != 0:
        if len(set(polbs_strands)) != 1:  # an ambiguous case
            return contig.contig_orientation
        return polbs[0].strand


def join_it(iterable, delimiter):
    it = iter(iterable)
    yield next(it)
    for x in it:
        yield delimiter
        yield x


def create_fragment_record(fragment, genome_dict):
    fragment_record = genome_dict[fragment.contig.contig_id][fragment.start:fragment.end]
    if fragment.contig.contig_orientation == Orientation.REVERSE:
        fragment_record = fragment_record.reverse_complement()
    return fragment_record


def feature_from_blasthit(hit, contig_id: str, genome: Genome) -> Feature:
    return Feature(start=hit.hit_start, end=hit.hit_end,
                   strand=Orientation.orientation_from_blast(hit.hit_strand),
                   contig_id=contig_id, genome=genome)
