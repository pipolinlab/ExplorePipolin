#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
from typing import Mapping
import click
from random import randrange
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from BCBio import GFF
from utilities import CONTEXT_SETTINGS
from utilities import read_from_shelve
from utilities import read_genbank_records, write_genbank_records
from utilities import check_dir


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
                     'note': ['att repeat'], 'product': ['att repeat']}
    gff_qualifiers = {'ID': [f'{pipolin.strain_id}_{random_number}'],
                      'inference': ['HMM:custom'], 'locus_tag': [f'{pipolin.strain_id}_{random_number}'],
                      'product': ['att repeat'], 'source': ['HMM:custom'], 'phase': ['0']}
    att_feature = SeqFeature(type='repeat_region',
                             location=FeatureLocation(bounds[0], bounds[1], strand=strand),
                             qualifiers=gb_qualifiers if records_format == 'genbank' else gff_qualifiers)
    return att_feature


def add_atts(records, records_format, pipolins):
    for pipolin in pipolins:
        if pipolin.is_complete_genome():
            record = records[pipolin.strain_id][pipolin.strain_id]
            for att in pipolin.atts:
                att_feature = create_att_feature(att, pipolin, records_format)
                record.features.append(att_feature)
        else:
            for att in pipolin.atts:
                record = records[pipolin.strain_id][att.node]
                att_feature = create_att_feature(att, pipolin, records_format)
                record.features.append(att_feature)


def write_gff_records(gff_records, new_annot_dir):
    for key, value in gff_records.items():
        records = [record for record in value.values()]
        with open(os.path.join(new_annot_dir, f'{key}.gff'), 'w') as ouf:
            GFF.write(records, ouf)
            print('##FASTA', file=ouf)
            SeqIO.write(records, ouf, format='fasta')


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('shelve-file', type=click.Path(exists=True))
@click.option('--object-name', required=True)
@click.argument('orig-annot-dir', type=click.Path(exists=True))
@click.argument('new-annot-dir', type=click.Path())
def include_atts_into_annotation(shelve_file, object_name, orig_annot_dir, new_annot_dir):
    """
    Adds att regions to the annotations (*.gbk and *.gff files)
    """
    pipolins = read_from_shelve(shelve_file, object_name)

    gb_records = read_genbank_records(orig_annot_dir)
    gff_records = read_gff_records(orig_annot_dir)

    add_atts(gb_records, 'genbank', pipolins)
    add_atts(gff_records, 'gff', pipolins)

    check_dir(new_annot_dir)
    write_genbank_records(gb_records, new_annot_dir)
    write_gff_records(gff_records, new_annot_dir)


if __name__ == '__main__':
    include_atts_into_annotation()
