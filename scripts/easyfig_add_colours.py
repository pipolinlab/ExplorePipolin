#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import click
from utilities import CONTEXT_SETTINGS, GenBankRecords
from utilities import read_genbank_records, write_genbank_records

red = '255 0 0'   # Primer-independent DNA polymerase PolB
brick_red = '139 58 58'   # Tyrosine recombinase XerC
yellow = '255 255 0'   # Type I site-specific deoxyribonuclease (hsdR)
# Type I restriction modification enzyme
# Type I restriction modification system methyltransferase (hsdM)
magenta = '255 0 255'   # metallohydrolase
purple = '178 58 238'   # excisionase
cyan = '0 255 255'   # Uracil-DNA glycosylase
green = '0 255 0'   # tRNA-Leu
blue = '0 0 255'   # att repeat
light_grey = '170 170 170'   # others

products_to_colour = {'Primer-independent DNA polymerase PolB': red,
                      'Tyrosine recombinase XerC': brick_red,
                      'Type I site-specific deoxyribonuclease (hsdR)': yellow,
                      'Type I restriction modification enzyme': yellow,
                      'Type I restriction modification system methyltransferase (hsdM)': yellow,
                      'metallohydrolase': magenta, 'excisionase': purple,
                      'Uracil-DNA glycosylase': cyan, 'tRNA-Leu': green, 'att repeat': blue}


def colour_feature(qualifiers):
    if 'product' in qualifiers:
        for product in qualifiers['product']:
            if product in products_to_colour:
                qualifiers['colour'] = [products_to_colour[product]]
            else:
                qualifiers['colour'] = [light_grey]


def add_colours(gb_records: GenBankRecords):
    for record_set in gb_records.values():
        for record in record_set.values():
            for feature in record.features:
                colour_feature(feature.qualifiers)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('in-dir', type=click.Path(exists=True))
def easyfig_add_colours(in_dir):
    """
    IN_DIR contains *.gbk file to modify.
    """
    gb_records = read_genbank_records(in_dir)
    add_colours(gb_records)
    write_genbank_records(gb_records, in_dir)


if __name__ == '__main__':
    easyfig_add_colours()
