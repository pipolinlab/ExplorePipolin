#!/usr/bin/env python3
# -*- encoding:utf-8 -*-

import os
import click
from Bio import SeqIO
from utilities import CONTEXT_SETTINGS


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('annot-dir', type=click.Path(exists=True))
@click.argument('out-dir', type=click.Path(exists=True))
@click.option('--ext', type=click.Choice(['faa', 'ffn']), required=True,
              help='Put "faa" if you need amino acid sequences and "ffn" if you need nucleotide sequences')
def get_all_pipolB_seqs(annot_dir, out_dir, ext):
    """
    TODO
    """
    files = [file for file in os.listdir(annot_dir) if file.endswith(ext)]
    all_seqs = {}
    for file in files:
        all_seqs.update(SeqIO.to_dict(SeqIO.parse(os.path.join(annot_dir, file), 'fasta')))

    pipolBs = []
    for _, seq in all_seqs.items():
        if seq.description == f'{seq.id} Primer-independent DNA polymerase PolB':
            pipolBs.append(seq)

    with open(os.path.join(out_dir, f'pipolBs.{ext}'), 'w') as ouf:
        SeqIO.write(pipolBs, ouf, 'fasta')


if __name__ == '__main__':
    get_all_pipolB_seqs()
