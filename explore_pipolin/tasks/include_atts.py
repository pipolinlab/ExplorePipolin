import os
from random import randrange
from typing import MutableSequence, Tuple, Sequence

from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from prefect import task

from explore_pipolin.common import Strand, Pipolin, FeatureType, ContigID, Genome
from explore_pipolin.utilities.io import SeqIORecords, create_seqio_records_dict, write_genbank_records, \
    read_gff_records, write_gff_records
from explore_pipolin.utilities.logging import genome_specific_logging


@task()
@genome_specific_logging
def include_atts(genome: Genome, prokka_dir, pipolins: Sequence[Pipolin]):
    results_dir = os.path.join(genome.work_dir, 'pipolins')
    os.makedirs(results_dir, exist_ok=True)

    for prokka_file in os.listdir(prokka_dir):

        if prokka_file.startswith(genome.id):
            pipolin_index = int(os.path.splitext(prokka_file)[0].split(sep='_')[-1])
            cur_pipolin = pipolins[pipolin_index]

            if prokka_file.endswith('.gbk'):
                gb_records = create_seqio_records_dict(file=os.path.join(prokka_dir, prokka_file),
                                                       file_format='genbank')

                include_atts_into_gb(gb_records=gb_records, pipolin=cur_pipolin)

                output_file = os.path.join(results_dir, prokka_file)
                write_genbank_records(gb_records=gb_records, output_file=output_file)

            if prokka_file.startswith(genome.id) and prokka_file.endswith('.gff'):
                gff_records = read_gff_records(file=os.path.join(prokka_dir, prokka_file))
                include_atts_into_gff(gff_records=gff_records, pipolin=cur_pipolin)

                output_file = os.path.join(results_dir, prokka_file)
                write_gff_records(gff_records=gff_records, output_file=output_file)

    return results_dir


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
                                                         strand=att.strand, contig_id=att.contig_id)
            elif record_format == 'gff':
                att_feature = _create_gff_att_seq_feature(start=att_start, end=att_end,
                                                          strand=att.strand, contig_id=att.contig_id)
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
