import os
from prefect import task
from explore_pipolin.utilities import create_att_seqfeatures
from explore_pipolin.utilities import read_gff_records
from explore_pipolin.utilities import read_seqio_records
from explore_pipolin.utilities import write_genbank_records
from explore_pipolin.utilities import write_gff_records
from explore_pipolin.utilities import add_new_gb_feature


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
