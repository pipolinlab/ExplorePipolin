#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import click
import pandas
from utilities import CONTEXT_SETTINGS


def get_att_regions(pipolins, long=False):
    att_regions = []
    for _, row in pipolins.iterrows():
        if long:
            if row['polB1_node'] != row['att2_node']:
                continue

        if pandas.isna(row['att3_s']) or (row['polB1_s'] > row['att1_e'] and row['polB1_e'] < row['att2_s']):
            att_regions.append((row['att1_s'], row['att1_e']))
            att_regions.append((row['att2_s'], row['att2_e']))
        else:
            # pi-polB within att2 and att3
            att_regions.append((row['att2_s'], row['att2_e']))
            att_regions.append((row['att3_s'], row['att3_e']))

    return att_regions


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('ncbi-csv', type=click.Path(exists=True))
@click.argument('new-csv', type=click.Path(exists=True))
def prepare_atts_for_msa(ncbi_csv, new_csv):
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
    lengths = [att[1] - att[0] for att in all_att_regions]
    print(f'> Maximum att length is {max(lengths)}')
    print(f'> Minimum att length is {min(lengths)}')

    # TODO: to extract att sequences I also need sequence ids!

    # TODO: after extracting att sequences they need to be sorted somehow as attL and attR


if __name__ == '__main__':
    prepare_atts_for_msa()