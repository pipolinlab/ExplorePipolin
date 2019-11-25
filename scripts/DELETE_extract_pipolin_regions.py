#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
import click
import pandas
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from utilities import CONTEXT_SETTINGS


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('csv-file', type=click.Path(exists=True))
@click.argument('genomes-dir', type=click.Path(exists=True))
@click.argument('out-dir', type=click.Path(exists=True))
def extract_pipolin_regions(csv_file, genomes_dir, out_dir):
    """
    This script reads the information about atts from the CSV_FILE
    and creates fasta files with pipolin regions
    """
    pipolins = pandas.read_csv(csv_file, na_values='None')

    genome_seqs = {}
    for file in os.listdir(genomes_dir):
        genome_seqs.update(SeqIO.to_dict(SeqIO.parse(os.path.join(genomes_dir, file), 'fasta')))

    for _, row in pipolins.iterrows():
        if not pandas.isna(row['att3_e']):
            start, end = [int(i) for i in [row['att1_s'], row['att3_e']]]
        else:
            start, end = [int(i) for i in [row['att1_s'], row['att2_e']]]

        with open(os.path.join(out_dir, f'{row["id"]}-pipolin.fa'), 'w') as ouf:
            record = SeqRecord(seq=genome_seqs[row['id']].seq[start - 50:end + 50],
                               id=row['id'], description='')
            SeqIO.write(record, ouf, 'fasta')


if __name__ == '__main__':
    extract_pipolin_regions()
