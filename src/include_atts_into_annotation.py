#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os

from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from prefect import task
import click
from typing import MutableSequence

from utilities import CONTEXT_SETTINGS
from utilities import read_gff_records
from utilities import read_seqio_records
from utilities import write_genbank_records
from utilities import write_gff_records
from utilities import GQuery


def add_new_gb_feature(new_feature: SeqFeature, record: SeqRecord):
    record.features.append(new_feature)
    record.features.sort(key=lambda x: x.location.start)


def create_att_seqfeatures(record_format: str, gquery: GQuery) -> MutableSequence[SeqFeature]:
    att_seqfeatures = []
    in_start = 0
    for fragment in gquery.pipolin_fragments:
        for att in fragment.atts:
            att_start = att.start - fragment.start + in_start
            att_end = att.end - fragment.start + in_start
            att_feature = gquery.create_att_feature(start=att_start, end=att_end, frame=att.frame,
                                                    records_format=record_format)
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


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('shelve-file', type=click.Path(exists=True))
@click.option('--object-name', required=True)
@click.argument('orig-annot-dir', type=click.Path(exists=True))
@click.argument('new-annot-dir', type=click.Path())
def main(shelve_file, object_name, orig_annot_dir, new_annot_dir):
    """
    Adds att regions to the annotations (*.gbk and *.gff files)
    TODO: add also *.faa and *.ffn files!
    """
    include_atts_into_annotation(shelve_file, object_name, orig_annot_dir)


if __name__ == '__main__':
    main()
