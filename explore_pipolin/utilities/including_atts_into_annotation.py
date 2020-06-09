from random import randrange
from typing import MutableSequence

from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

from explore_pipolin.common import Orientation
from explore_pipolin.utilities.io import SeqIORecords
from explore_pipolin.utilities.misc import GQuery


def include_atts_into_gb(gb_records: SeqIORecords, gquery: GQuery):
    att_seq_features = _generate_att_seq_features(record_format='gb', gquery=gquery)
    for att in att_seq_features:
        _add_att_seq_feature(att_seq_feature=att, seq_record=gb_records[gquery.genome.genome_id])


def include_atts_into_gff(gff_records: SeqIORecords, gquery: GQuery):
    att_seq_features = _generate_att_seq_features(record_format='gff', gquery=gquery)
    for att in att_seq_features:
        _add_att_seq_feature(att_seq_feature=att, seq_record=gff_records[gquery.genome.genome_id])


def _generate_att_seq_features(record_format: str, gquery: GQuery) -> MutableSequence[SeqFeature]:
    att_seq_features = []
    in_start = 0
    for fragment in gquery.pipolin_fragments:
        fragment_shift = fragment.start if fragment.contig.contig_orientation == Orientation.FORWARD else fragment.end
        for att in fragment.atts:
            att_start, att_end = sorted([abs(att.start - fragment_shift), abs(att.end - fragment_shift)])
            if record_format == 'gb':
                att_feature = _create_gb_att_seq_feature(start=att_start + in_start, end=att_end + in_start,
                                                         strand=att.strand, genome_id=gquery.genome.genome_id,)
            elif record_format == 'gff':
                att_feature = _create_gff_att_seq_feature(start=att_start + in_start, end=att_end + in_start,
                                                          strand=att.strand, genome_id=gquery.genome.genome_id)
            else:
                raise AssertionError
            att_seq_features.append(att_feature)
        in_start += (fragment.end - fragment.start) + 100

    return att_seq_features


def _create_gb_att_seq_feature(start: int, end: int, strand: Orientation, genome_id: str) -> SeqFeature:
    random_number = randrange(10000, 99999)
    gb_qualifiers = {'inference': ['HMM:custom'], 'locus_tag': [f'{genome_id}_{random_number}'],
                     'rpt_family': ['Att'], 'rpt_type': ['direct']}
    att_seq_feature = SeqFeature(type='repeat_region',
                                 location=FeatureLocation(start=start, end=end, strand=strand.to_pm_one_encoding()),
                                 qualifiers=gb_qualifiers)
    return att_seq_feature


def _create_gff_att_seq_feature(start: int, end: int, strand: Orientation, genome_id: str) -> SeqFeature:
    random_number = randrange(10000, 99999)
    gff_qualifiers = {'phase': ['.'], 'source': ['HMM:custom'],
                      'ID': [f'{genome_id}_{random_number}'], 'inference': ['HMM:custom'],
                      'locus_tag': [f'{genome_id}_{random_number}'],
                      'rpt_family': ['Att'], 'rpt_type': ['direct']}
    att_seq_feature = SeqFeature(type='repeat_region',
                                 location=FeatureLocation(start=start, end=end, strand=strand.to_pm_one_encoding()),
                                 qualifiers=gff_qualifiers)
    return att_seq_feature


def _add_att_seq_feature(att_seq_feature: SeqFeature, seq_record: SeqRecord):
    seq_record.features.append(att_seq_feature)
    seq_record.features.sort(key=lambda x: x.location.start)
