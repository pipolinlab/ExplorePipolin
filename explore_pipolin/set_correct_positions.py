import os
from prefect import task
from explore_pipolin.utilities import GQuery
from explore_pipolin.utilities import read_seqio_records
from explore_pipolin.utilities import write_genbank_records
from explore_pipolin.utilities import create_new_gb_record


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
