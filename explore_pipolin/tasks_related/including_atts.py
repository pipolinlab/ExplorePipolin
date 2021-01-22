from random import randrange
from typing import MutableSequence, Tuple

from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

from explore_pipolin.common import Strand, Pipolin, FeatureType, ContigID
from explore_pipolin.utilities.io import SeqIORecords


def include_atts_into_gb(gb_records: SeqIORecords, pipolin: Pipolin):
    att_seq_features = _generate_att_seq_features(record_format='gb', pipolin=pipolin)
    for att in att_seq_features:
        _add_att_seq_feature(att_seq_feature=att[0], seq_record=gb_records[att[1]])


def include_atts_into_gff(gff_records: SeqIORecords, pipolin: Pipolin):
    att_seq_features = _generate_att_seq_features(record_format='gff', pipolin=pipolin)
    for att in att_seq_features:
        _add_att_seq_feature(att_seq_feature=att[0], seq_record=gff_records[att[1]])


def _generate_att_seq_features(record_format: str, pipolin: Pipolin):
    att_seq_features: MutableSequence[Tuple[SeqFeature, ContigID]] = []

    for fragment in pipolin.fragments:
        fragment_shift = fragment.start

        for att in [f[0] for f in fragment.features if f[1] == FeatureType.ATT]:
            att_start, att_end = (att.start - fragment_shift), (att.end - fragment_shift)
            if record_format == 'gb':
                att_feature = _create_gb_att_seq_feature(start=att_start, end=att_end,
                                                         strand=att.strand, contig_id=fragment.contig_id)
            elif record_format == 'gff':
                att_feature = _create_gff_att_seq_feature(start=att_start, end=att_end,
                                                          strand=att.strand, contig_id=fragment.contig_id)
            else:
                raise AssertionError

            att_seq_features.append((att_feature, fragment.contig_id))

    return att_seq_features


def _create_gb_att_seq_feature(start: int, end: int, strand: Strand, contig_id: str) -> SeqFeature:
    random_number = randrange(10000, 99999)
    gb_qualifiers = {'inference': ['HMM:custom'], 'locus_tag': [f'{contig_id}_{random_number}'],
                     'rpt_family': ['Att'], 'rpt_type': ['direct']}
    att_seq_feature = SeqFeature(type='repeat_region',
                                 location=FeatureLocation(start=start, end=end, strand=strand.to_pm_one_encoding()),
                                 qualifiers=gb_qualifiers)
    return att_seq_feature


def _create_gff_att_seq_feature(start: int, end: int, strand: Strand, contig_id: str) -> SeqFeature:
    random_number = randrange(10000, 99999)
    gff_qualifiers = {'phase': ['.'], 'source': ['HMM:custom'],
                      'ID': [f'{contig_id}_{random_number}'], 'inference': ['HMM:custom'],
                      'locus_tag': [f'{contig_id}_{random_number}'],
                      'rpt_family': ['Att'], 'rpt_type': ['direct']}
    att_seq_feature = SeqFeature(type='repeat_region',
                                 location=FeatureLocation(start=start, end=end, strand=strand.to_pm_one_encoding()),
                                 qualifiers=gff_qualifiers)
    return att_seq_feature


def _add_att_seq_feature(att_seq_feature: SeqFeature, seq_record: SeqRecord):
    seq_record.features.append(att_seq_feature)
    seq_record.features.sort(key=lambda x: x.location.start)
