#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
import click
from utilities import CONTEXT_SETTINGS, GenBankRecords
from utilities import read_genbank_records, write_genbank_records
from Bio.SeqIO import SeqRecord

red = '255 0 0'   # Primer-independent DNA polymerase PolB
brick_red = '139 58 58'   # Tyrosine recombinase XerC
brown = '200 150 100'   # Prophage integrase IntS
yellow = '255 255 0'   # Type I site-specific deoxyribonuclease (hsdR)
# Type I restriction modification enzyme
# Type I restriction modification system methyltransferase (hsdM)
magenta = '255 0 255'   # metallohydrolase
purple = '178 58 238'   # excisionase
cyan = '0 255 255'   # Uracil-DNA glycosylase
green = '0 255 0'   # tRNA-Leu
blue = '0 0 255'   # repeat_region
floral_white = '255 250 240'   # others
black = '0 0 0'   # pipolin_structure
pink = '255 200 200'   # paired-ends

orange = '255 165 0'   # virulence gene
light_steel_blue = '176 196 222'   # AMR gene

products_to_colours = {'Primer-independent DNA polymerase PolB': red,
                       'Tyrosine recombinase XerC': brick_red,
                       'Type I site-specific deoxyribonuclease (hsdR)': yellow,
                       'Type I restriction modification enzyme': yellow,
                       'Type I restriction modification system methyltransferase (hsdM)': yellow,
                       'metallohydrolase': magenta, 'excisionase': purple,
                       'Uracil-DNA glycosylase': cyan, 'tRNA-Leu': green,
                       'repeat_region': blue, 'pipolin_structure': black,
                       'paired-ends': pink, 'Prophage integrase IntS': brown,
                       'other': floral_white}

VIRULENCE_DBs = ['vfdb', 'ecoli_vf', 'ecoli_virfinder']
AMR_DBs = ['megares_noSNPs', 'argannot', 'ncbi', 'card', 'resfinder']


def colour_feature(qualifiers):
    if 'product' in qualifiers:
        for product in qualifiers['product']:
            if product in products_to_colours:
                qualifiers['colour'] = [products_to_colours[product]]
            else:
                qualifiers['colour'] = [products_to_colours['other']]
    if 'linkage_evidence' in qualifiers:
        for evidence in qualifiers['linkage_evidence']:
            if evidence in products_to_colours:
                qualifiers['colour'] = [products_to_colours[evidence]]
            else:
                qualifiers['colour'] = [products_to_colours['other']]


def add_colours(gb_records: GenBankRecords):
    for record_set in gb_records.values():
        for record in record_set.values():
            for feature in record.features:
                if feature.type in products_to_colours:
                    feature.qualifiers['colour'] = products_to_colours[feature.type]
                else:
                    colour_feature(feature.qualifiers)


def read_record_ranges(record: SeqRecord):
    ranges = []
    for feature in record.features:
        ranges.append((feature.location.start, feature.location.end, feature.location.strand))
    return ranges


def find_feature_position(s_ranges, q_range):
    all_ranges = s_ranges + [q_range]
    all_ranges.sort(key=lambda x: (x[0], x[1]))
    q_position = all_ranges.index(q_range)
    q = set(range(q_range[0], q_range[1] + 1))
    if all_ranges[q_position][0] < all_ranges[q_position - 1][1]:
        if all_ranges[q_position][2] == all_ranges[q_position - 1][2]:
            s = set(range(all_ranges[q_position - 1][0], all_ranges[q_position - 1][1]))
            qs = q.intersection(s)
            # if len(qs) / len(s) >= 0.9:
            s_position = s_ranges.index(all_ranges[q_position - 1])
            return s_position
    if all_ranges[q_position][1] > all_ranges[q_position + 1][0]:
        if all_ranges[q_position][2] == all_ranges[q_position + 1][2]:
            s = set(range(all_ranges[q_position + 1][0], all_ranges[q_position + 1][1]))
            qs = q.intersection(s)
            # if len(qs) / len(s) >= 0.9:
            s_position = s_ranges.index(all_ranges[q_position + 1])
            return s_position
    return None


def color_feature(record: SeqRecord, feature_position, db_type):
    if record.features[feature_position].qualifiers['colour'] == ['255 250 240']:
        if db_type == 'vir':
            record.features[feature_position].qualifiers['colour'] = orange
        elif db_type == 'amr':
            record.features[feature_position].qualifiers['colour'] = light_steel_blue
        else:
            raise AssertionError('Something is wrong here!')


def add_info_from_summary(gb_records: GenBankRecords, summary_path, db_type):
    with open(summary_path) as inf:
        _ = inf.readline()   # skip header
        for line in inf:
            entry_line = line.strip().split(sep='\t')
            q_id = entry_line[1]
            q_location = (int(entry_line[2]), int(entry_line[3]), 1 if entry_line[4] == '+' else -1)
            q_cov = float(entry_line[10])
            if q_cov >= 10.0:
                record = gb_records[os.path.basename(summary_path)[:-4]][q_id]
                record_ranges = read_record_ranges(record)
                feature_position = find_feature_position(record_ranges, q_location)
                if feature_position is not None:
                    color_feature(record, feature_position, db_type)
                gb_records[os.path.basename(summary_path)[:-4]][q_id] = record


def find_and_color_amr_and_virulence(gb_records: GenBankRecords, abricate_dir):
    contents = os.listdir(abricate_dir)
    for content in contents:
        if content in VIRULENCE_DBs or content in AMR_DBs:
            content_path = os.path.join(abricate_dir, content)
            db_type = 'vir' if content in VIRULENCE_DBs else 'amr'
            summaries = os.listdir(content_path)
            for summary in summaries:
                add_info_from_summary(gb_records, os.path.join(content_path, summary), db_type)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.argument('in-dir', type=click.Path(exists=True))
@click.argument('abricate_dir', type=click.Path(exists=True), required=False)
def easyfig_add_colours(in_dir, abricate_dir):
    """
    IN_DIR contains *.gbk files to modify.
    """
    gb_records = read_genbank_records(in_dir)
    add_colours(gb_records)
    if abricate_dir is not None:
        find_and_color_amr_and_virulence(gb_records, abricate_dir)
    write_genbank_records(gb_records, in_dir)


if __name__ == '__main__':
    easyfig_add_colours()
