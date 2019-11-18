#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
import click
import pandas
from collections import namedtuple
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from utilities import CONTEXT_SETTINGS
from utilities import save_as_csv


FIELDS = ('id', 'node', 'att_s', 'att_e')
AttRegions = namedtuple('AttRegions', FIELDS, defaults=(None, ) * len(FIELDS))


def get_att_regions(pipolins, long=False):
    att_regions = []
    for _, row in pipolins.iterrows():
        if long:
            if row['polB1_node'] != row['att2_node']:
                continue

        if pandas.isna(row['att3_s']) or (row['polB1_s'] > row['att1_e'] and row['polB1_e'] < row['att2_s']):
            att_regions.append(AttRegions(id=row['id'], att_s=row['att1_s'], att_e=row['att1_e']))
            att_regions.append(AttRegions(id=row['id'], att_s=row['att2_s'], att_e=row['att2_e']))
        else:
            # pi-polB within att2 and att3
            att_regions.append(AttRegions(id=row['id'], att_s=row['att2_s'], att_e=row['att2_e']))
            att_regions.append(AttRegions(id=row['id'], att_s=row['att3_s'], att_e=row['att3_e']))

        if long:
            att_regions[-2] = att_regions[-2]._replace(node=row['att1_node'])
            att_regions[-1] = att_regions[-1]._replace(node=row['att1_node'])

    return att_regions


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('ncbi-csv', type=click.Path(exists=True))
@click.argument('genomes-dir', type=click.Path(exists=True))
@click.argument('new-csv', type=click.Path(exists=True))
@click.argument('contigs-dir', type=click.Path(exists=True))
@click.argument('out-dir', type=click.Path(exists=True))
def prepare_atts_for_msa(ncbi_csv, genomes_dir, new_csv, contigs_dir, out_dir):
    """
    TODO
    """
    ncbi_pipolins = pandas.read_csv(ncbi_csv, na_values='None')
    ncbi_att_regions = get_att_regions(ncbi_pipolins)
    print(f'Number of atts from NCBI_CSV file is {len(ncbi_att_regions)}')

    new_pipolins = pandas.read_csv(new_csv, na_values='None')
    new_att_regions = get_att_regions(new_pipolins, long=True)
    print(f'Number of atts from NEW_CSV file is {len(new_att_regions)}')

    all_att_regions = ncbi_att_regions + new_att_regions
    lengths = [att.att_e - att.att_s for att in all_att_regions]
    print(f'> Maximum att length is {max(lengths)}')
    print(f'> Minimum att length is {min(lengths)}')

    att_records = []
    for att in all_att_regions:
        if att.node:
            contig_seq = SeqIO.to_dict(SeqIO.parse(os.path.join(contigs_dir, f'{att.id}.fa'), 'fasta'))[att.node]
            record = SeqRecord(seq=contig_seq.seq[int(att.att_s) - 10:int(att.att_e) + 10],
                               id=f'{att.id} {att.node}', description='')   # TODO: add node to the name!
        else:
            genome_seq = SeqIO.read(os.path.join(genomes_dir, f'{att.id}.fa'), 'fasta')
            record = SeqRecord(seq=genome_seq.seq[int(att.att_s) - 10:int(att.att_e) + 10],
                               id=att.id, description='')
        att_records.append(record)

    save_as_csv(AttRegions, all_att_regions, os.path.join(out_dir, 'atts.csv'))
    with open(os.path.join(out_dir, 'att_sequences.fa'), 'w') as ouf:
        SeqIO.write(att_records, ouf, 'fasta')

    # TODO: after extracting att sequences they need to be sorted somehow as attL and attR



if __name__ == '__main__':
    prepare_atts_for_msa()
