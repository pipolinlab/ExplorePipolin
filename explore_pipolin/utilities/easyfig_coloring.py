from enum import Enum
from Bio.SeqIO import SeqRecord


class EasyfigColour(Enum):
    RED = '255 0 0'   # Primer-independent DNA polymerase PolB
    BRICK_RED = '139 58 58'   # Tyrosine recombinase XerC
    BROWN = '200 150 100'   # Prophage integrase IntS
    YELLOW = '255 255 0'   # Type I site-specific deoxyribonuclease (hsdR)
    # Type I restriction modification enzyme
    # Type I restriction modification system methyltransferase (hsdM)
    MAGENTA = '255 0 255'   # metallohydrolase
    PURPLE = '178 58 238'   # excisionase
    CYAN = '0 255 255'   # Uracil-DNA glycosylase
    GREEN = '0 255 0'   # tRNA-Leu
    BLUE = '0 0 255'   # repeat_region
    FLORAL_WHITE = '255 250 240'   # others
    BLACK = '0 0 0'   # pipolin_structure
    PINK = '255 200 200'   # paired-ends


products_to_colours = {'Primer-independent DNA polymerase PolB': EasyfigColour.RED,
                       'Tyrosine recombinase XerC': EasyfigColour.BRICK_RED,
                       'Type I site-specific deoxyribonuclease (hsdR)': EasyfigColour.YELLOW,
                       'Type I restriction modification enzyme': EasyfigColour.YELLOW,
                       'Type I restriction modification system methyltransferase (hsdM)': EasyfigColour.YELLOW,
                       'metallohydrolase': EasyfigColour.MAGENTA, 'excisionase': EasyfigColour.PURPLE,
                       'Uracil-DNA glycosylase': EasyfigColour.CYAN, 'tRNA-Leu': EasyfigColour.GREEN,
                       'repeat_region': EasyfigColour.BLUE, 'pipolin_structure': EasyfigColour.BLACK,
                       'paired-ends': EasyfigColour.PINK, 'Prophage integrase IntS': EasyfigColour.BROWN,
                       'other': EasyfigColour.FLORAL_WHITE}


def colour_feature(qualifiers):
    if 'product' in qualifiers:
        for product in qualifiers['product']:
            if product in products_to_colours:
                qualifiers['colour'] = [products_to_colours[product].value]
            else:
                qualifiers['colour'] = [products_to_colours['other'].value]
    elif 'linkage_evidence' in qualifiers:
        if qualifiers['estimated_length'] == ['100']:
            qualifiers['linkage_evidence'] = ['pipolin_structure']
            qualifiers['estimated_length'] = ['unknown']
            qualifiers['colour'] = [products_to_colours['pipolin_structure'].value]
        else:
            qualifiers['color'] = [products_to_colours[qualifiers['linkage_evidence'][0]].value]
    else:
        qualifiers['colour'] = [products_to_colours['other'].value]


def add_colours(record: SeqRecord):
    for feature in record.features:
        if feature.type in products_to_colours:
            feature.qualifiers['colour'] = products_to_colours[feature.type].value
        else:
            colour_feature(feature.qualifiers)
