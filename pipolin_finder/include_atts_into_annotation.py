import os
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from prefect import task
from typing import MutableSequence

from pipolin_finder.utilities import GQuery, Orientation
from pipolin_finder.utilities import read_gff_records
from pipolin_finder.utilities import read_seqio_records
from pipolin_finder.utilities import write_genbank_records
from pipolin_finder.utilities import write_gff_records


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
