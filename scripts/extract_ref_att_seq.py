#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
import click
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from utilities import CONTEXT_SETTINGS


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('reference', type=click.Path(exists=True))
def extract_ref_att_seq(reference):
    """
    This is a stupid script to generate ref_att.fa file, which will be used to roughly
    estimate pipolin region boundaries. In future, it will be replaced by att profile.
    REFERENCE is a NZ_JNMI01000006.1.fa file.
    """
    ref_seq = SeqIO.read(reference, 'fasta')
    att_record = SeqRecord(seq=ref_seq.seq[64241: 64373],
                           id='att-site', description='')
    os.chdir(os.path.dirname(reference))
    with open('ref_att.fa', 'w') as ouf:
        SeqIO.write(att_record, ouf,'fasta')


if __name__ == '__main__':
    extract_ref_att_seq()