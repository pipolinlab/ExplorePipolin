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
