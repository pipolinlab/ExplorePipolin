#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
import click
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from utilities import CONTEXT_SETTINGS


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('reference', type=click.Path(exists=True))
@click.argument('out-dir', type=click.Path(exists=True))
def extract_ref_att_seq(reference, out_dir):
    """
    This is a stupid script to generate ref_att.fa file, which will be used to roughly
    estimate pipolin region boundaries. In future, it will be replaced by att profile.
    REFERENCE is a NZ_JNMI01000006.1.fa file.
    """
    ref_seq = SeqIO.read(reference, 'fasta')
    att_record = SeqRecord(seq=ref_seq.seq[64241: 64373],
                           id='att-site', description='')
    with open(os.path.join(out_dir, 'ref_att.fa'), 'w') as ouf:
        SeqIO.write(att_record, ouf,'fasta')


if __name__ == '__main__':
    extract_ref_att_seq()
