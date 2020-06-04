import os
from explore_pipolin.utilities.misc import GQuery
from explore_pipolin.utilities.io import SeqIORecords
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

orange = '255 165 0'   # AMR gene
light_steel_blue = '176 196 222'   # virulence gene

products_to_colours = {'Primer-independent DNA polymerase PolB': red,
                       'Tyrosine recombinase XerC': brick_red,
                       'Type I site-specific deoxyribonuclease (hsdR)': yellow,
                       'Type I restriction modification enzyme': yellow,
                       'Type I restriction modification system methyltransferase (hsdM)': yellow,
                       'metallohydrolase': magenta, 'excisionase': purple,
                       'Uracil-DNA glycosylase': cyan, 'tRNA-Leu': green,  # 'tRNA-Arg': green,
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
    elif 'linkage_evidence' in qualifiers:
        if qualifiers['estimated_length'] == ['100']:
            qualifiers['linkage_evidence'] = ['pipolin_structure']
            qualifiers['estimated_length'] = ['unknown']
            qualifiers['colour'] = [products_to_colours['pipolin_structure']]
        else:
            qualifiers['color'] = [products_to_colours[qualifiers['linkage_evidence'][0]]]
    else:
        qualifiers['colour'] = [products_to_colours['other']]


def add_colours(record: SeqRecord):
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

    if all_ranges[q_position][0] < all_ranges[q_position - 1][1]:
        if all_ranges[q_position][2] == all_ranges[q_position - 1][2]:
            s_position = s_ranges.index(all_ranges[q_position - 1])
            return s_position

    if all_ranges[q_position][1] > all_ranges[q_position + 1][0]:
        if all_ranges[q_position][2] == all_ranges[q_position + 1][2]:
            s_position = s_ranges.index(all_ranges[q_position + 1])
            return s_position

    return None


def color_feature(record: SeqRecord, feature_position, db_type):
    if record.features[feature_position].qualifiers['colour'] == ['255 250 240']:
        if db_type == 'vir':
            record.features[feature_position].qualifiers['colour'] = light_steel_blue
        elif db_type == 'amr':
            record.features[feature_position].qualifiers['colour'] = orange
        else:
            raise AssertionError('Something is wrong here!')


def add_info_from_summary(gb_records: SeqIORecords, summary_path, db_type):
    with open(summary_path) as inf:
        _ = inf.readline()   # skip header
        for line in inf:
            entry_line = line.strip().split(sep='\t')
            q_id = entry_line[1]
            q_location = (int(entry_line[2]), int(entry_line[3]), 1 if entry_line[4] == '+' else -1)
            q_cov = float(entry_line[9])

            if q_cov >= 10.0:
                record = gb_records[q_id]
                record_ranges = read_record_ranges(record)
                feature_position = find_feature_position(record_ranges, q_location)
                if feature_position is not None:
                    color_feature(record, feature_position, db_type)
                gb_records[q_id] = record


def find_and_color_amr_and_virulence(gquery: GQuery, gb_records: SeqIORecords, abricate_dir):
    contents = os.listdir(abricate_dir)
    for content in contents:
        if content in VIRULENCE_DBs or content in AMR_DBs:
            content_path = os.path.join(abricate_dir, content)
            db_type = 'vir' if content in VIRULENCE_DBs else 'amr'
            summary_path = os.path.join(content_path, gquery.genome.id + '.tab')
            add_info_from_summary(gb_records, summary_path, db_type)
