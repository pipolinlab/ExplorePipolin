#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import click
from rename_files import create_dict_w_strainnames
from Bio import AlignIO

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('in_alignment', type=click.Path(exists=True))
@click.argument('out_alignment', type=click.Path())
@click.argument('names_file', type=click.Path(exists=True))
def rename_fasta_ids(in_alignment, out_alignment, names_file):
    """
    TODO
    """
    accs_to_strainnames = create_dict_w_strainnames(names_file)

    alignment = AlignIO.read(in_alignment, format="fasta")
    for record in alignment:
        if record.id[:4] != 'chr_':
            record.id = accs_to_strainnames[record.id]
        else:
            record.id = record.id[4:]

    AlignIO.write(alignment, out_alignment, format="fasta")


if __name__ == '__main__':
    rename_fasta_ids()
