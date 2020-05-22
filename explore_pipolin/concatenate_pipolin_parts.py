#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import click
from explore_pipolin.utilities import CONTEXT_SETTINGS
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import collections


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('fasta-file',type=click.Path(exists=True))
@click.argument('new_fasta', type=click.Path())
def concatenate(fasta_file, new_fasta):
    """
    TODO
    """
    records = list(SeqIO.parse(fasta_file, format='fasta'))

    for record in records:
        record.id = '_'.join(record.id.split(sep='_')[:-1])
    ids = [record.id for record in records]
    duplicates = [item for item, count in collections.Counter(ids).items() if count > 1]

    new_records = []
    for duplicate in duplicates:
        records_to_concat = []
        for record in records:
            if record.id == duplicate:
                records_to_concat.append(record)
        if len(records_to_concat) <= 1:
            raise AssertionError('Something wrong!!!')
        records_to_concat.sort(key=lambda x: x.name)
        new_seq = ''.join(str(record.seq) for record in records_to_concat)
        new_records.append(SeqRecord(seq=Seq(new_seq), id=duplicate,
                                     description='Primer-independent DNA polymerase PolB'))

    new_records.extend(record for record in records if record.id not in duplicates)
    with open(new_fasta, 'w') as ouf:
        SeqIO.write(new_records, ouf, format='fasta')


if __name__ == '__main__':
    concatenate()
