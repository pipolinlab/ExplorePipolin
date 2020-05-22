import copy
import os
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from prefect import task
from typing import MutableSequence

from explore_pipolin.utilities import GQuery, Orientation
from explore_pipolin.utilities import read_gff_records
from explore_pipolin.utilities import read_seqio_records
from explore_pipolin.utilities import write_genbank_records
from explore_pipolin.utilities import write_gff_records


def add_new_gb_feature(new_feature: SeqFeature, record: SeqRecord):
    record.features.append(new_feature)
    record.features.sort(key=lambda x: x.location.start)


def create_att_seqfeatures(record_format: str, gquery: GQuery) -> MutableSequence[SeqFeature]:
    att_seqfeatures = []
    in_start = 0
    for fragment in gquery.pipolin_fragments:
        fragment_shift = fragment.start if fragment.contig.contig_orientation == Orientation.FORWARD else fragment.end
        for att in fragment.atts:
            att_start, att_end = sorted([abs(att.start - fragment_shift), abs(att.end - fragment_shift)])
            att_feature = gquery.create_att_feature(start=att_start + in_start, end=att_end + in_start,
                                                    frame=att.frame, records_format=record_format)
            att_seqfeatures.append(att_feature)
        in_start += (fragment.end - fragment.start) + 100

    return att_seqfeatures


@task
def include_atts_into_annotation(gquery, prokka_dir, root_dir):
    gb_records = read_seqio_records(file=os.path.join(prokka_dir, gquery.gquery_id + '.gbk'), file_format='genbank')
    gff_records = read_gff_records(file=os.path.join(prokka_dir, gquery.gquery_id + '.gff'))

    att_seqfeatures = create_att_seqfeatures(record_format='gb', gquery=gquery)
    for att in att_seqfeatures:
        add_new_gb_feature(new_feature=att, record=gb_records[gquery.gquery_id])
    att_seqfeatures = create_att_seqfeatures(record_format='gff', gquery=gquery)
    for att in att_seqfeatures:
        add_new_gb_feature(new_feature=att, record=gff_records[gquery.gquery_id])

    prokka_atts_dir = os.path.join(root_dir, 'prokka_atts')
    os.makedirs(prokka_atts_dir, exist_ok=True)

    write_genbank_records(gb_records=gb_records, out_dir=prokka_atts_dir, gquery=gquery)
    write_gff_records(gff_records=gff_records, out_dir=prokka_atts_dir, gquery=gquery)

    return prokka_atts_dir


@task
def set_correct_positions(gquery: GQuery, prokka_atts_dir, root_dir):
    gb_records = read_seqio_records(file=os.path.join(prokka_atts_dir, gquery.gquery_id + '.gbk'),
                                    file_format='genbank')
    # gff_records = read_gff_records(file=os.path.join(prokka_atts_dir, gquery.gquery_id + '.gff'))

    new_gb_records = {gquery.gquery_id: create_new_gb_record(gquery=gquery, gb_record=gb_records[gquery.gquery_id])}

    prokka_atts_positions_dir = os.path.join(root_dir, 'prokka_atts_positions')
    os.makedirs(prokka_atts_positions_dir, exist_ok=True)

    write_genbank_records(gb_records=new_gb_records, out_dir=prokka_atts_positions_dir, gquery=gquery)
    # write_gff_records(gff_records=gff_records, out_dir=prokka_atts_positions_dir, gquery=gquery)

    return prokka_atts_positions_dir


def create_new_gb_record(gquery: GQuery, gb_record: SeqRecord) -> SeqRecord:
    new_record = copy.deepcopy(gb_record)
    new_source_features = []

    in_start = 0
    next_start = 0
    for fragment in gquery.pipolin_fragments:
        fragment_shift = fragment.start
        next_start += (fragment.end - fragment.start) + 100

        for i_f, feature in enumerate(gb_record.features):
            if feature.location.start >= in_start and feature.location.end <= next_start:
                old_start, old_end, old_strand = feature.location.start, feature.location.end, feature.location.strand
                new_record.features[i_f].location = FeatureLocation(start=old_start - in_start + fragment_shift,
                                                                    end=old_end - in_start + fragment_shift,
                                                                    strand=old_strand)

        in_start += next_start

        source_location = FeatureLocation(start=fragment.start, end=fragment.end + 100)
        source_feature = SeqFeature(type='source', location=source_location,
                                    qualifiers=copy.deepcopy(gb_record.features[0].qualifiers))
        source_feature.qualifiers.update({'note': [fragment.contig.contig_id,
                                                   fragment.contig.contig_orientation.to_string()]})
        new_source_features.append(source_feature)

    del new_record.features[0]
    for feature in new_source_features:
        add_new_gb_feature(new_feature=feature, record=new_record)

    return new_record


# def create_assembly_gap_record(record):
#     source_feature = SeqFeature(type='source', location=FeatureLocation(1, 100, strand=+1),
#                                 qualifiers={'mol_type': record.features[0].qualifiers['mol_type'],
#                                             'organism': record.features[0].qualifiers['organism'],
#                                             'strain': record.features[0].qualifiers['strain']})
#     assembly_gap_seq = Seq('N' * 100, alphabet=IUPACAmbiguousDNA())
#     assembly_gap_qualifiers = {'estimated_length': ['unknown'],
#                                'gap_type': ['within_scaffolds'],
#                                'linkage_evidence': ['pipolin_structure']}
#     assembly_gap_feature = SeqFeature(type='assembly_gap',
#                                       location=FeatureLocation(1, 100, strand=+1),
#                                       qualifiers=assembly_gap_qualifiers)
#     assembly_gap_record = SeqRecord(seq=assembly_gap_seq, id=record.id, name=record.name,
#                                     description=record.description, features=[source_feature, assembly_gap_feature],
#                                     annotations=record.annotations)
#
#     return assembly_gap_record
