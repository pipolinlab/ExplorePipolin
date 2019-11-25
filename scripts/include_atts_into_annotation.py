#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
import click
from Bio import SeqIO
from utilities import CONTEXT_SETTINGS


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('annot-dir', type=click.Path(exists=True))
def include_atts_into_annotation(annot_dir):
    """
    TODO
    """
    annotations = {}
    for file in os.listdir(annot_dir):
        strain_id = os.path.splitext(file)[0]
        if file.endswith('.gbk'):
            annotations[strain_id] = SeqIO.parse(os.path.join(annot_dir, file), format='genbank')

    print(next(iter(annotations['NZ_JNMI01000006.1'])))


if __name__ == '__main__':
    include_atts_into_annotation()
