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


def read_gff_records(annot_dir):
    gff_records: Mapping[str, Mapping[str, SeqIO.SeqRecord]] = {}
    for file in os.listdir(annot_dir):
        if file.endswith('.gff'):
            with open(os.path.join(annot_dir, file)) as inf:
                records = {}
                for entry in GFF.parse(inf):
                    records[entry.id] = entry
                gff_records[os.path.splitext(file)[0]] = records
    return gff_records


def get_bounds_with_strand(att_feature):
    if att_feature.start < att_feature.end:
        return (att_feature.start, att_feature.end), +1
    else:
        return (att_feature.end, att_feature.start), -1


def create_att_feature(att, pipolin, records_format):
    bounds, strand = get_bounds_with_strand(att)
    random_number = randrange(10000, 99999)
    gb_qualifiers = {'inference': ['HMM:custom'], 'locus_tag': [f'{pipolin.strain_id}_{random_number}'],
                     'rpt_family': ['Att'], 'rpt_type': ['direct']}
    gff_qualifiers = {'phase': ['.'], 'source': ['HMM:custom'],
                      'ID': [f'{pipolin.strain_id}_{random_number}'], 'inference': ['HMM:custom'],
                      'locus_tag': [f'{pipolin.strain_id}_{random_number}'],
                      'rpt_family': ['Att'], 'rpt_type': ['direct']}
    att_feature = SeqFeature(type='repeat_region',
                             location=FeatureLocation(bounds[0], bounds[1], strand=strand),
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


def add_atts(records, records_format, pipolins):
    for pipolin in pipolins:
        if pipolin.is_complete_genome():
            record = records[pipolin.strain_id][pipolin.strain_id]
            for att in pipolin.atts:
                att_feature = create_att_feature(att, pipolin, records_format)
                add_new_gb_feature(att_feature, record)
        else:
            for att in pipolin.atts:
                record = records[pipolin.strain_id][att.node]
                att_feature = create_att_feature(att, pipolin, records_format)
                add_new_gb_feature(att_feature, record)


@task
def include_atts_into_annotation(shelve_in_dir, object_name, orig_annot_dir):
    pipolins = read_from_shelve(os.path.join(shelve_in_dir, 'shelve.db'), object_name)
    gb_records = read_seqio_records(orig_annot_dir)
    gff_records = read_gff_records(orig_annot_dir)
    add_atts(gb_records, 'gb', pipolins)
    add_atts(gff_records, 'gff', pipolins)
    new_annot_dir = os.path.join(shelve_in_dir, 'prokka_atts')
    os.makedirs(new_annot_dir, exist_ok=True)
    write_genbank_records(gb_records, new_annot_dir)
    write_gff_records(gff_records, new_annot_dir)

    return new_annot_dir


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
