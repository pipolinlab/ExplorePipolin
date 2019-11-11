#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
import click
import pandas
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from utilities import CONTEXT_SETTINGS   # TODO fix this!


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('csv-file', type=click.Path(exists=True))
@click.argument('contigs-dir', type=click.Path(exists=True))
def extract_pipolin_regions(csv_file, contigs_dir):
    """
    This script reads the information about atts from the CSV_FILE
    and creates fasta files with pipolin regions
    """
    pipolins = pandas.read_csv(csv_file, na_values='None')

    contig_seqs = {}
    for file in os.listdir(contigs_dir):
        # TODO: all nodes are unique, but it'd be better to have a test for this
        contig_seqs.update(SeqIO.to_dict(SeqIO.parse(os.path.join(contigs_dir, file), 'fasta')))

    os.chdir(os.path.dirname(csv_file))
    if not os.path.exists('pipolin_regions'):
        os.mkdir('pipolin_regions')

    for _, row in pipolins.iterrows():
        # TODO: the assumption will be: it cannot be less than 50000 bp between pi-polB and atts
        # For example, make a combinations of two to compare all the nodes and identify those that are identical
        records = []
        if row['polB1_node'] != row['att1_node']:
            records.append(SeqRecord(seq=contig_seqs[row['polB1_node']].seq,
                              id=row['polB1_node'], description=''))
            records.append(SeqRecord(seq=contig_seqs[row['att1_node']].seq,
                                     id=row['att1_node'], description=''))
            if not pandas.isna(row['polB2_node']):
                records.append(SeqRecord(seq=contig_seqs[row['polB2_node']].seq,
                                         id=row['polB2_node'], description=''))
            if not pandas.isna(row['att2_node']):
                records.append(SeqRecord(seq=contig_seqs[row['att2_node']].seq,
                                         id=row['att2_node'], description=''))
            if not pandas.isna(row['att3_node']):
                records.append(SeqRecord(seq=contig_seqs[row['att3_node']].seq,
                                         id=row['att3_node'], description=''))
        elif row['polB1_node'] != row['att2_node']:
            records.append(SeqRecord(seq=contig_seqs[row['polB1_node']].seq,
                                     id=row['polB1_node'], description=''))
            if not pandas.isna(row['polB2_node']):
                records.append(SeqRecord(seq=contig_seqs[row['polB2_node']].seq,
                                         id=row['polB2_node'], description=''))
            if not pandas.isna(row['att2_node']):
                records.append(SeqRecord(seq=contig_seqs[row['att2_node']].seq,
                                         id=row['att2_node'], description=''))
            if not pandas.isna(row['att3_node']):
                records.append(SeqRecord(seq=contig_seqs[row['att3_node']].seq,
                                         id=row['att3_node'], description=''))
        else:
            records.append(SeqRecord(seq=contig_seqs[row['polB1_node']].seq,
                                     id=row['polB1_node'], description=''))
            if not pandas.isna(row['polB2_node']):
                records.append(SeqRecord(seq=contig_seqs[row['polB2_node']].seq,
                                         id=row['polB2_node'], description=''))
            if not pandas.isna(row['att3_node']):
                records.append(SeqRecord(seq=contig_seqs[row['att3_node']].seq,
                                         id=row['att3_node'], description=''))

        with open(f'pipolin_regions/{row["id"]}-pipolin.fa', 'w') as ouf:
            SeqIO.write(records, ouf, 'fasta')


if __name__ == '__main__':
    extract_pipolin_regions()