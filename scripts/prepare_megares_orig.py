#!/usr/bin/evn python3
# -*- encoding: utf-8 -*-

import click
from Bio import SeqIO
from utilities import CONTEXT_SETTINGS


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('megares-in', type=click.Path(exists=True))
@click.argument('megares-out', type=click.Path())
def prepare_megares_orig(megares_in, megares_out):
    """
    The script takes MEGARES_IN fasta file, downloaded from MEGARes database,
    and creates a MEGARES_OUT fasta file, from which the sequences with the flag
    "RequiresSNPConfirmation" are removed.
    """
    records_in = SeqIO.parse(megares_in, format='fasta')
    records_out = []
    for record in records_in:
        if 'RequiresSNPConfirmation' not in record.id:
            records_out.append(record)
    with open(megares_out, 'w') as ouf:
        SeqIO.write(records_out, ouf, format='fasta')


if __name__ == '__main__':
    prepare_megares_orig()
