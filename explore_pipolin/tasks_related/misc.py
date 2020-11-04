from explore_pipolin.common import Orientation, Contig, Genome, Feature, FeatureType


def is_single_target_trna_per_contig(genome: Genome):
    # TODO: don't like this
    # there was one case with two target trnas per genome, although usually only one
    target_trnas_dict = genome.features.get_features_dict_by_contig_normalized(feature_type=FeatureType.TARGET_TRNA)
    target_trnas = genome.features.get_features(feature_type=FeatureType.TARGET_TRNA)
    if len(target_trnas) != len(target_trnas_dict):
        raise AssertionError("We are expecting a single tRNA to overlap with a single att per contig!")


def get_contig_orientation(contig: Contig, genome: Genome) -> Orientation:
    target_trnas = genome.features.get_features_of_contig_normalized(contig_id=contig.contig_id,
                                                                     feature_type=FeatureType.TARGET_TRNA)
    atts = genome.features.get_features_of_contig_normalized(contig_id=contig.contig_id,
                                                             feature_type=FeatureType.ATT)
    atts_strands = [att.strand for att in atts]
    polbs = genome.features.get_features_of_contig_normalized(contig_id=contig.contig_id,
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


def get_left_right_windows(genome: Genome, feature_type) -> [Feature, Feature]:
    features = genome.features.get_features(feature_type=feature_type)
    features = sorted(features, key=lambda x: x.start)

    # if features[-1].start - features[0].start > 10000:   # TODO: This should be changed!
    #     raise AssertionError(f'You have several piPolBs per genome and they are too far from each other: '
    #                          f'within the region ({features[0].start}, {features[-1].end}). It might be, '
    #                          f'that you have two or more pipolins per genome, but we are expecting only one.')

    contig_id = genome.features.get_features(feature_type=feature_type)[0].contig_id
    contig_length = genome.get_contig_by_id(contig_id=contig_id).contig_length

    left_edge = features[0].start - 100000
    left_window = Feature(start=left_edge if left_edge >= 0 else 0, end=features[0].start,
                          strand=Orientation.FORWARD,
                          contig_id=contig_id, genome=genome)
    right_edge = features[-1].end + 100000
    right_window = Feature(start=features[-1].end, end=right_edge if right_edge <= contig_length else contig_length,
                           strand=Orientation.FORWARD,
                           contig_id=contig_id, genome=genome)

    return left_window, right_window


def add_features_from_blast_entries(entries, feature_type: FeatureType, genome: Genome,):
    for entry in entries:
        for hit in entry:
            feature = feature_from_blasthit(hit=hit, contig_id=entry.id, genome=genome)
            genome.features.add_feature(feature=feature, feature_type=feature_type)
