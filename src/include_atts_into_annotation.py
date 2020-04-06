#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
from prefect import task
from typing import Mapping
import click
from random import randrange
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from BCBio import GFF
from utilities import CONTEXT_SETTINGS
from utilities import read_from_shelve
from utilities import read_seqio_records, write_genbank_records
from utilities import write_gff_records
from utilities import GQuery, Feature


def read_gff_records(file):
    gff_records = {}
    with open(file) as inf:
        for entry in GFF.parse(inf):
            gff_records[entry.id] = entry

    return gff_records


def create_att_feature(att: Feature, gquery: GQuery, records_format: str):
    f_start = att.start - gquery.pipolin_fragment.start
    f_end = att.end - gquery.pipolin_fragment.start
    random_number = randrange(10000, 99999)
    gb_qualifiers = {'inference': ['HMM:custom'], 'locus_tag': [f'{gquery.gquery_id}_{random_number}'],
                     'rpt_family': ['Att'], 'rpt_type': ['direct']}
    gff_qualifiers = {'phase': ['.'], 'source': ['HMM:custom'],
                      'ID': [f'{gquery.gquery_id}_{random_number}'], 'inference': ['HMM:custom'],
                      'locus_tag': [f'{gquery.gquery_id}_{random_number}'],
                      'rpt_family': ['Att'], 'rpt_type': ['direct']}
    att_feature = SeqFeature(type='repeat_region',
                             location=FeatureLocation(start=f_start, end=f_end, strand=att.frame.to_pm_one_encoding()),
                             qualifiers=gb_qualifiers if records_format == 'gb' else gff_qualifiers)
    return att_feature


def add_new_gb_feature(new_feature, record):
    new_feature_index = None
    for i_f, feature in enumerate(record.features):
        if new_feature.location.start < feature.location.start:
            new_feature_index = i_f
            break

    if new_feature_index is None:
        new_feature_index = len(record.features)
    record.features.insert(new_feature_index, new_feature)


def add_atts(records, records_format, gquery: GQuery):
    for att in gquery.pipolin_fragment.atts:
        att_feature = create_att_feature(att=att, gquery=gquery, records_format=records_format)
        add_new_gb_feature(att_feature, records[gquery.pipolin_fragment.contig.contig_id])


@task
def include_atts_into_annotation(gquery, prokka_dir, root_dir):
    gb_records = read_seqio_records(file=os.path.join(prokka_dir, gquery.gquery_id + '.gbk'), file_format='genbank')
    gff_records = read_gff_records(file=os.path.join(prokka_dir, gquery.gquery_id + '.gff'))

    add_atts(gb_records, 'gb', gquery)
    add_atts(gff_records, 'gff', gquery)

    prokka_atts_dir = os.path.join(root_dir, 'prokka_atts')
    os.makedirs(prokka_atts_dir, exist_ok=True)

    write_genbank_records(gb_records=gb_records, out_dir=prokka_atts_dir, gquery=gquery)
    write_gff_records(in_records=gff_records, out_dir=prokka_atts_dir, gquery=gquery)

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
