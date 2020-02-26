#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
import click
from Bio import SeqIO
from utilities import CONTEXT_SETTINGS, read_from_prokka_dir
from utilities import check_dir


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('roary-dir', type=click.Path(exists=True))
@click.argument('prokka-dir', type=click.Path(exists=True))
@click.argument('out-dir', type=click.Path())
@click.option('--ext', type=click.Choice(['faa', 'ffn']), required=True,
              help='Put "faa" if you need amino acid sequences and "ffn" if you need nucleotide sequences')
def extract_roary_groups(roary_dir, prokka_dir, out_dir, ext):
    """
    TODO
    """
    roary_groups = {}
    with open(os.path.join(roary_dir, 'gene_presence_absence.csv')) as inf:
        _ = inf.readline()  # skip a header line
        for line in inf:
            group = line.strip().replace('\t', ',').split(sep=',')
            group = [i.replace('"', '') for i in group]  # remove quotation marks from items
            roary_groups[group[0]] = [group[i] for i in range(14, len(group)) if len(group[i]) > 0]

    records = read_from_prokka_dir(prokka_dir, ext)

    check_dir(out_dir)
    for group_id, genes in roary_groups.items():
        seqs = [records[gene] for gene in genes if gene in records]
        with open(os.path.join(out_dir, f'{group_id}.{ext}'), 'w') as ouf:
            SeqIO.write(seqs, ouf, 'fasta')


if __name__ == '__main__':
    extract_roary_groups()
