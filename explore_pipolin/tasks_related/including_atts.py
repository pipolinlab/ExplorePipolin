from random import randrange
from typing import MutableSequence

from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

from explore_pipolin.common import Strand, Pipolin, Genome
from explore_pipolin.utilities.io import SeqIORecords


def include_atts_into_gb(gb_records: SeqIORecords, genome: Genome, pipolin: Pipolin):
    att_seq_features = _generate_att_seq_features(record_format='gb', genome=genome, pipolin=pipolin)
    for att in att_seq_features:
        _add_att_seq_feature(att_seq_feature=att, seq_record=gb_records[genome.id])


def include_atts_into_gff(gff_records: SeqIORecords, genome: Genome, pipolin: Pipolin):
    att_seq_features = _generate_att_seq_features(record_format='gff', genome=genome, pipolin=pipolin)
    for att in att_seq_features:
        _add_att_seq_feature(att_seq_feature=att, seq_record=gff_records[genome.id])


def _generate_att_seq_features(record_format: str, genome: Genome, pipolin: Pipolin) -> MutableSequence[SeqFeature]:
    att_seq_features = []
    in_start = 0
    for fragment in pipolin.fragments:
        fragment_shift = fragment.start
        for att in fragment.atts:
            att_start, att_end = sorted([abs(att.start - fragment_shift), abs(att.end - fragment_shift)])
            if record_format == 'gb':
                att_feature = _create_gb_att_seq_feature(start=att_start + in_start, end=att_end + in_start,
                                                         strand=att.strand, genome_id=genome.id, )
            elif record_format == 'gff':
                att_feature = _create_gff_att_seq_feature(start=att_start + in_start, end=att_end + in_start,
                                                          strand=att.strand, genome_id=genome.id)
            else:
                raise AssertionError
            att_seq_features.append(att_feature)
        in_start += (fragment.end - fragment.start) + 100

    return att_seq_features


def _create_gb_att_seq_feature(start: int, end: int, strand: Strand, genome_id: str) -> SeqFeature:
    random_number = randrange(10000, 99999)
    gb_qualifiers = {'inference': ['HMM:custom'], 'locus_tag': [f'{genome_id}_{random_number}'],
                     'rpt_family': ['Att'], 'rpt_type': ['direct']}
    att_seq_feature = SeqFeature(type='repeat_region',
                                 location=FeatureLocation(start=start, end=end, strand=strand.to_pm_one_encoding()),
                                 qualifiers=gb_qualifiers)
    return att_seq_feature


def _create_gff_att_seq_feature(start: int, end: int, strand: Strand, genome_id: str) -> SeqFeature:
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
