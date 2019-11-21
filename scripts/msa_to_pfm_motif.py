#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import click
from Bio import AlignIO
from Bio import motifs
from Bio.Alphabet import IUPAC
from Bio.Alphabet import Gapped
from utilities import CONTEXT_SETTINGS


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('msa', type=click.Path(exists=True))
@click.argument('out-file')
def att_motif_scan(msa, out_file):
    """
    TODO
    """
    alignment = AlignIO.read(msa, format='fasta', alphabet=Gapped(IUPAC.unambiguous_dna, gap_char='-'))
    motif = motifs.create([i.seq for i in alignment])
    with open(out_file, 'w') as ouf:
        ouf.write(motifs.write([motif], 'pfm'))


if __name__ == '__main__':
    att_motif_scan()
