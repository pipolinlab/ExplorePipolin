import os
from typing import MutableSequence, Tuple, Sequence

from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from prefect import task

from explore_pipolin.common import Strand, Pipolin, FeatureType, ContigID, Genome, PipolinVariants
from explore_pipolin.tasks.easyfig_coloring import easyfig_add_colours
from explore_pipolin.utilities.io import SeqIORecords, create_seqio_records_dict, write_seqio_records, \
    read_gff_records, write_gff_records
from explore_pipolin.utilities.logging import genome_specific_logging
import explore_pipolin.settings as settings


@task()
@genome_specific_logging
def generate_results(genome: Genome, prokka_dir, pipolins: Sequence[PipolinVariants]):
    results_dir = os.path.join(os.path.dirname(prokka_dir), 'pipolins')
    os.makedirs(results_dir, exist_ok=True)

    for prokka_file in os.listdir(prokka_dir):

        if prokka_file.startswith(genome.id) and (prokka_file.endswith('.gbk') or prokka_file.endswith('.gff')):
            # genome_N_vN.type.ext
            variant_index = int(os.path.splitext(os.path.splitext(prokka_file)[0])[0].split(sep='_')[-1][1:])
            pipolin_index = int(os.path.splitext(os.path.splitext(prokka_file)[0])[0].split(sep='_')[-2])

            cur_pipolin = pipolins[pipolin_index].variants[variant_index]

            if prokka_file.endswith('.gbk'):
                gb_records = create_seqio_records_dict(file=os.path.join(prokka_dir, prokka_file),
                                                       file_format='genbank')

                include_atts_into_gb(gb_records=gb_records, pipolin=cur_pipolin)
                if settings.get_instance().skip_colours is False:
                    easyfig_add_colours(gb_records=gb_records, pipolin=cur_pipolin)

                single_record = create_single_gb_record(gb_records=gb_records, pipolin=cur_pipolin)
                output_file_single_record = os.path.join(
                    results_dir, os.path.splitext(prokka_file)[0] + '.single_record.gbk'
                )
                write_seqio_records(single_record, output_file_single_record, 'genbank')

                output_file = os.path.join(results_dir, prokka_file)
                write_seqio_records(gb_records, output_file, 'genbank')

            if prokka_file.endswith('.gff'):
                gff_records = read_gff_records(file=os.path.join(prokka_dir, prokka_file))
                include_atts_into_gff(gff_records=gff_records, pipolin=cur_pipolin)

                output_file = os.path.join(results_dir, prokka_file)
                write_gff_records(gff_records=gff_records, output_file=output_file)

    return results_dir


def create_single_gb_record(gb_records: SeqIORecords, pipolin: Pipolin) -> SeqIORecords:
    record = revcompl_if_reverse(gb_records[pipolin.fragments[0].contig_id],
                                 pipolin.fragments[0].orientation)
    if len(pipolin.fragments) > 1:
        for fragment in pipolin.fragments[1:]:
            record += SeqRecord(seq='N' * 100)
            del gb_records[fragment.contig_id].features[0]   # delete source feature
            record += revcompl_if_reverse(gb_records[fragment.contig_id], fragment.orientation)

    old_source = record.features[0]
    new_source = SeqFeature(FeatureLocation(0, len(record), 1), type=old_source.type,
                            location_operator=old_source.location_operator, strand=old_source.strand,
                            id=old_source.id, qualifiers=old_source.qualifiers,
                            ref=old_source.ref, ref_db=old_source.ref_db)
    del record.features[0]
    record.features.insert(0, new_source)

    genome_id = pipolin.fragments[0].genome.id
    return {genome_id: record}


def revcompl_if_reverse(gb_record: SeqRecord, orientation: Strand) -> SeqRecord:
    if orientation == Strand.REVERSE:
        return gb_record.reverse_complement(id=True, name=True, description=True, annotations=True)
    return gb_record


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

        for att in [f for f in fragment.features if f.ftype == FeatureType.ATT]:
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
    gb_qualifiers = {'inference': ['HMM:custom'], 'locus_tag': [f'{contig_id}_00000'],
                     'rpt_family': ['Att'], 'rpt_type': ['direct']}
    att_seq_feature = SeqFeature(type='repeat_region',
                                 location=FeatureLocation(start=start, end=end, strand=strand.to_pm_one_encoding()),
                                 qualifiers=gb_qualifiers)
    return att_seq_feature


def _create_gff_att_seq_feature(start: int, end: int, strand: Strand, contig_id: str) -> SeqFeature:
    gff_qualifiers = {'phase': ['.'], 'source': ['HMM:custom'],
                      'ID': [f'{contig_id}_00000'], 'inference': ['HMM:custom'],
                      'locus_tag': [f'{contig_id}_00000'],
                      'rpt_family': ['Att'], 'rpt_type': ['direct']}
    att_seq_feature = SeqFeature(type='repeat_region',
                                 location=FeatureLocation(start=start, end=end, strand=strand.to_pm_one_encoding()),
                                 qualifiers=gff_qualifiers)
    return att_seq_feature


def _add_att_seq_feature(att_seq_feature: SeqFeature, seq_record: SeqRecord):
    seq_record.features.append(att_seq_feature)
    seq_record.features.sort(key=lambda x: x.location.start)
