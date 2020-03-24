#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import click
from Bio import SeqIO
from utilities import CONTEXT_SETTINGS, read_from_prokka_dir


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('ids_file', type=click.Path(exists=True))
@click.argument('prokka-dir', type=click.Path(exists=True))
@click.argument('out-file', type=click.Path())
@click.option('--ext', type=click.Choice(['faa', 'ffn']), required=True,
              help='Put "faa" if you need amino acid sequences and "ffn" if you need nucleotide sequences')
def extract_sequences_from_prokka(ids_file, prokka_dir, out_file, ext):
    """
    TODO
    """
    seq_ids = []

    with open(ids_file) as inf:
        for line in inf:
            seq_ids.append(line.strip())

    records = read_from_prokka_dir(prokka_dir, ext)

    seqs = [records[s_id] for s_id in seq_ids if s_id in records]
    with open(out_file, 'w') as ouf:
        SeqIO.write(seqs, ouf, 'fasta')


if __name__ == '__main__':
    extract_sequences_from_prokka()
