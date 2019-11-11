#!/usr/bin/env python3
# -*- encoding:utf-8 -*-

import os
import click
from Bio import SeqIO
from utilities import CONTEXT_SETTINGS


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('annot-dir', type=click.Path(exists=True))
def get_all_pipolB_seqs(annot_dir):
    """
    TODO
    """
    os.chdir(os.path.dirname(annot_dir))
    faa_files = [file for file in os.listdir(annot_dir) if file.endswith('.faa')]
    all_proteins = {}
    for file in faa_files:
        all_proteins.update(SeqIO.to_dict(SeqIO.parse(os.path.join(annot_dir, file), 'fasta')))

    pipolBs = []
    for _, protein in all_proteins.items():
        if protein.description == f'{protein.id} Primer-independent DNA polymerase PolB':
            pipolBs.append(protein)

    with open('pipolBs.faa', 'w') as ouf:
        SeqIO.write(pipolBs, ouf, 'fasta')


if __name__ == '__main__':
    get_all_pipolB_seqs()