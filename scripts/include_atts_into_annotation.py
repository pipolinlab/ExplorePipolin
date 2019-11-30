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
from utilities import Feature, Pipolin
from utilities import read_from_shelve
from utilities import check_dir


def get_bounds_with_strand(att_feature):
    if att_feature.start < att_feature.end:
        return (att_feature.start, att_feature.end), +1
    else:
        return (att_feature.end, att_feature.start), -1


def create_att_feature(att, pipolin, format):
    bounds, strand = get_bounds_with_strand(att)
    random_number = randrange(10000, 99999)
    gb_qualifiers = {'inference': ['HMM:custom'], 'locus_tag': [f'{pipolin.strain_id}_{random_number}'],
                     'note': ['att repeat'], 'product': ['att repeat']}
    gff_qualifiers = {'ID': [f'{pipolin.strain_id}_{random_number}'],
                      'inference': ['HMM:custom'], 'locus_tag': [f'{pipolin.strain_id}_{random_number}'],
                      'product': ['att repeat'], 'source': ['HMM:custom'], 'phase': ['0']}
    att_feature = SeqFeature(type='repeat_region',
                             location=FeatureLocation(bounds[0], bounds[1], strand=strand),
                             qualifiers=gb_qualifiers if format == 'gb' else gff_qualifiers)
    return att_feature


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('shelve-file', type=click.Path(exists=True))
@click.argument('annot-dir', type=click.Path(exists=True))
@click.argument('new-annot-dir')
def include_atts_into_annotation(shelve_file, annot_dir, new_annot_dir):
    """
    TODO
    """
    pipolins = read_from_shelve(shelve_file, 'pipolins')

    genbank_records: Mapping[str, Mapping[str, SeqIO.SeqRecord]] = {}
    for file in os.listdir(annot_dir):
        if file.endswith('.gbk'):
            # TODO: fixed LOCUS name (see extract_pipolin_regions.py)
            record = SeqIO.to_dict(SeqIO.parse(os.path.join(annot_dir, file), format='genbank'))
            genbank_records[os.path.splitext(file)[0]] = record

    gff_records: Mapping[str, Mapping[str, SeqIO.SeqRecord]] = {}
    for file in os.listdir(annot_dir):
        if file.endswith('.gff'):
            with open(os.path.join(annot_dir, file)) as inf:
                records = {}
                for entry in GFF.parse(inf):
                    records[entry.id] = entry
                gff_records[os.path.splitext(file)[0]] = records

    for pipolin in pipolins:
        if pipolin.is_complete_genome():
            gb_record = genbank_records[pipolin.strain_id][pipolin.strain_id]
            gff_record = gff_records[pipolin.strain_id][pipolin.strain_id]
            for att in pipolin.atts:
                att_gb_feature = create_att_feature(att, pipolin, 'gb')
                gb_record.features.append(att_gb_feature)
                att_gff_feature = create_att_feature(att, pipolin, 'gff')
                gff_record.features.append(att_gff_feature)
        else:
            for att in pipolin.atts:
                node = f'{att.node.split("_")[0]}_{att.node.split("_")[1]}'
                gb_record = genbank_records[pipolin.strain_id][node]
                gff_record = gff_records[pipolin.strain_id][node]
                att_gb_feature = create_att_feature(att, pipolin, 'gb')
                gb_record.features.append(att_gb_feature)
                att_gff_feature = create_att_feature(att, pipolin, 'gff')
                gff_record.features.append(att_gff_feature)

    check_dir(new_annot_dir)
    for key, value in genbank_records.items():
        records = [record for record in value.values()]
        with open(os.path.join(new_annot_dir, f'{key}.gbk'), 'w') as ouf:
            SeqIO.write(records, ouf, 'genbank')

    for key, value in gff_records.items():
        records = [record for record in value.values()]
        with open(os.path.join(new_annot_dir, f'{key}.gff'), 'w') as ouf:
            GFF.write(records, ouf)
            print('##FASTA', file=ouf)
            SeqIO.write(records, ouf, format='fasta')


if __name__ == '__main__':
    include_atts_into_annotation()