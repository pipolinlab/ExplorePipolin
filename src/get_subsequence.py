#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
import click
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from utilities import CONTEXT_SETTINGS


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('in-fasta', type=click.Path(exists=True))
@click.argument('start', type=click.INT)
@click.argument('end', type=click.INT)
@click.argument('out-fasta', type=click.Path())
@click.option('--contig-name', help='If IN_FASTA is not a complete genome, but rather '
                                    'several contigs, specify the contig name.')
def get_subsequence(in_fasta, start, end, contig_name, out_fasta):
    """
    This script creates an OUT_FASTA file with the subsequence
    that you wish to extract from IN_FASTA file.
    If the START position is greater than the END position,
    the subsequence is considered to come from the minus strand
    and the reverse-complement will be returned.
    """
    in_seqs = SeqIO.to_dict(SeqIO.parse(in_fasta, 'fasta'))
    in_seq = in_seqs[contig_name] if contig_name else next(iter(in_seqs.values()))

    subseq = in_seq.seq[end:start].reverse_complement() if start > end else in_seq.seq[start:end]
    header = os.path.splitext(os.path.basename(out_fasta))[0]
    record = SeqRecord(seq=subseq, id=header, description='')

    with open(out_fasta, 'w') as ouf:
        SeqIO.write(record, ouf, 'fasta')


if __name__ == '__main__':
    get_subsequence()
