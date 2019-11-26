#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
import click
from Bio import SeqIO
from utilities import CONTEXT_SETTINGS
from utilities import Feature, Pipolin
from utilities import read_pipolins_from_shelve
from utilities import check_dir


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('shelve-file', type=click.Path(exists=True))
@click.argument('annot-dir', type=click.Path(exists=True))
@click.argument('new-annot-dir')
def include_atts_into_annotation(shelve_file, annot_dir, new_annot_dir):
    """
    TODO
    """
    pipolins = read_pipolins_from_shelve(shelve_file)

    genbank_records = {}
    for file in os.listdir(annot_dir):
        if file.endswith('.gbk'):
            # TODO: fixed LOCUS name (see extract_pipolin_regions.py)
            record = SeqIO.to_dict(SeqIO.parse(os.path.join(annot_dir, file), format='genbank'))
            genbank_records[os.path.splitext(file)[0]] = record

    check_dir(new_annot_dir)
    for key, value in genbank_records.items():
        records = [record for record in value.values()]
        with open(os.path.join(new_annot_dir, f'{key}.gbk'), 'w') as ouf:
            SeqIO.write(records, ouf, 'genbank')


if __name__ == '__main__':
    include_atts_into_annotation()
