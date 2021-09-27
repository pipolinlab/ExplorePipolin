from enum import Enum
from Bio.SeqIO import SeqRecord

from explore_pipolin.common import Pipolin, FeatureType, Range
from explore_pipolin.utilities.io import SeqIORecords


def easyfig_add_colours(gb_records: SeqIORecords, pipolin: Pipolin):
    for record in gb_records.values():
        add_colours(record)

    for fragment in pipolin.fragments:
        fragment_shift = fragment.start

        for att in [f for f in fragment.features if f.ftype == FeatureType.ATT]:
            att_start, att_end = (att.start - fragment_shift), (att.end - fragment_shift)

            for feature in gb_records[fragment.contig_id].features:
                if feature.location.start == att_start and feature.location.end == att_end:
                    feature.qualifiers['colour'] = _products_to_colours['FeatureType.ATT'].value

        for ttrna in [f for f in fragment.features if f.ftype == FeatureType.TARGET_TRNA]:
            ttrna_start, ttrna_end = (ttrna.start - fragment_shift), (ttrna.end - fragment_shift)

            for feature in gb_records[fragment.contig_id].features:
                feature_range = Range(start=feature.location.start, end=feature.location.end)
                if feature.type == 'tRNA' and feature_range.is_overlapping(Range(start=ttrna_start, end=ttrna_end)):
                    feature.qualifiers['colour'] = _products_to_colours['FeatureType.TARGET_TRNA'].value

        for pipolb in [f for f in fragment.features if f.ftype == FeatureType.PIPOLB]:
            pipolb_start, pipolb_end = (pipolb.start - fragment_shift), (pipolb.end - fragment_shift)

            for feature in gb_records[fragment.contig_id].features:
                feature_range = Range(start=feature.location.start, end=feature.location.end)
                if feature.type == 'CDS' and feature_range.is_overlapping(Range(start=pipolb_start, end=pipolb_end)):
                    feature.qualifiers['colour'] = _products_to_colours['FeatureType.PIPOLB'].value


class EasyfigColour(Enum):
    RED = '255 0 0'   # FeatureType.PIPOLB
    BLUE = '0 0 255'  # FeatureType.ATT
    GREEN = '0 255 0'  # FeatureType.TARGET_TRNA

    # gaps
    BLACK = '0 0 0'   # pipolin_structure (gap from reconstruction)
    PINK = '255 200 200'   # paired-ends (scaffolding gap)

    # other useful features
    BRICK_RED = '139 58 58'   # Tyrosine recombinase XerC
    BROWN = '200 150 100'   # Prophage integrase IntS
    YELLOW = '255 255 0'   # Type I site-specific deoxyribonuclease (hsdR)
                           # Type I restriction modification enzyme
                           # Type I restriction modification system methyltransferase (hsdM)
    MAGENTA = '255 0 255'   # metallohydrolase
    PURPLE = '178 58 238'   # excisionase
    CYAN = '0 255 255'   # Uracil-DNA glycosylase

    FLORAL_WHITE = '255 250 240'   # other


_products_to_colours = {'FeatureType.PIPOLB': EasyfigColour.RED,
                        'FeatureType.ATT': EasyfigColour.BLUE,
                        'FeatureType.TARGET_TRNA': EasyfigColour.GREEN,
                        # gaps
                        'pipolin_structure': EasyfigColour.BLACK, 'paired-ends': EasyfigColour.PINK,
                        # other useful features
                        'Tyrosine recombinase XerC': EasyfigColour.BRICK_RED,
                        'Type I site-specific deoxyribonuclease (hsdR)': EasyfigColour.YELLOW,
                        'Type I restriction modification enzyme': EasyfigColour.YELLOW,
                        'Type I restriction modification system methyltransferase (hsdM)': EasyfigColour.YELLOW,
                        'metallohydrolase': EasyfigColour.MAGENTA,
                        'excisionase': EasyfigColour.PURPLE,
                        'Uracil-DNA glycosylase': EasyfigColour.CYAN,
                        'Prophage integrase IntS': EasyfigColour.BROWN,
                        'other': EasyfigColour.FLORAL_WHITE}


def add_colours(record: SeqRecord):
    for feature in record.features:
        _colour_feature(feature.qualifiers)


def _colour_feature(qualifiers):
    if 'product' in qualifiers:
        for product in qualifiers['product']:
            if product in _products_to_colours:
                qualifiers['colour'] = [_products_to_colours[product].value]
            else:
                qualifiers['colour'] = [_products_to_colours['other'].value]
    elif 'linkage_evidence' in qualifiers:
        if qualifiers['linkage_evidence'] == ['paired-ends']:
            qualifiers['colour'] = [_products_to_colours['paired-ends'].value]
        elif qualifiers['estimated_length'] == ['100']:
            qualifiers['linkage_evidence'] = ['pipolin_structure']
            qualifiers['estimated_length'] = ['unknown']
            qualifiers['colour'] = [_products_to_colours['pipolin_structure'].value]
        else:
            qualifiers['colour'] = [_products_to_colours['other'].value]
    else:
        qualifiers['colour'] = [_products_to_colours['other'].value]
